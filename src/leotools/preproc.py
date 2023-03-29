"""This module contains tools used for preprocessing satellite images."""

from pathlib import Path
from shutil import copy
import fnmatch
from tempfile import TemporaryDirectory
import tarfile
import zipfile
import xml.etree.ElementTree as ET
import re
import operator

import numpy as np
import pandas as pd
from osgeo import gdal
import rasterio
from rasterio.enums import Resampling
from rasterio.merge import merge
import geopandas as gpd
from shapely.geometry import Polygon, box

from leotools.constants import EOV, GTIFF_UINT16, GTIFF_UINT16_COMPATIBLE
from leotools.basetools import ProcessTimer, check_path, load_files
from leotools.gistools import round_extent, extract_bands, image_to_array, make_aux, make_ovr, get_tags

### Constants
L47_BANDS = { ### number: name
    'B1': 'Blue',
    'B2': 'Green',
    'B3': 'Red',
    'B4': 'NIR',
    'B5': 'SWIR',
    'B6': 'TIR', ### Normally unused
    'B7': 'SWIR',
    'B8': 'Panchromatic', ### Only in Landsat 7, Normally unused
}

L89_BANDS = { ### number: name
    'B1': 'Coastal aerosol',
    'B2': 'Blue',
    'B3': 'Green',
    'B4': 'Red',
    'B5': 'NIR',
    'B6': 'SWIR 1',
    'B7': 'SWIR 2',
    'B8': 'Panchromatic', ### Normally unused
    'B9': 'Cirrus', ### Normally unused
    'B10': 'TIRS 1', ### Normally unused
    'B11': 'TIRS 2', ### Normally unused
}

S2_BANDS = { ### number: (resolution, name)
    'B01': (60, 'Coastal aerosol'),
    'B02': (10, 'Blue'),
    'B03': (10, 'Green'),
    'B04': (10, 'Red'),
    'B05': (20, 'Red edge 1'),
    'B06': (20, 'Red edge 2'),
    'B07': (20, 'Red edge 3'),
    'B08': (10, 'NIR'),
    'B8A': (20, 'NIR narrow'),
    'B09': (60, 'Water vapor'), ### Normally unused
    'B10': (60, 'Cirrus'), ### Normally unused
    'B11': (20, 'SWIR 1'),
    'B12': (20, 'SWIR 2'),
}

def calc_transform(src_bounding_box, dst_resolution, src_crs=None, dst_round=300):
    """Calculates source and destination transforms from provided parameters.

    Args:
        src_bounding_box (list): Coordinates of the source bounding box.
        dst_resolution (int): Destination resolution in meters.
        src_crs (str, optional): Projection of the source.
        dst_round (number, optional): The resolution to round to.

    Returns:
        dst_transform (dict): Destination transform.
        dst_width (int): Width of the destination in pixels.
        dst_height (int): Height of the destination in pixels.
    """

    if src_crs:
        extent_utm = gpd.GeoSeries(src_bounding_box, crs=src_crs)
        extent_eov = round_extent(extent_utm.to_crs(EOV).total_bounds, dst_round)
    else:
        extent_eov =  src_bounding_box.bounds

    ### Reprojection
    dst_width = int((extent_eov[2]-extent_eov[0])/dst_resolution)
    dst_height = int((extent_eov[3]-extent_eov[1])/dst_resolution)

    dst_transform = rasterio.transform.from_bounds(
        *extent_eov,
        width=dst_width,
        height=dst_height,
    )

    return dst_transform, dst_width, dst_height

def ls_tile(tar_file, image_dir, meta_dir=None, used_bands=None, mode=0):
    """Exports a GeoTIFF and a metadata from a Landsat 4-9 archive.

    Args:
        tar_file (path): A Landsat 4-9 archive. Has to be a .tar or .tar.gz file.
        image_dir (path): Directory the GeoTIFF is extracted to.
        meta_dir (path, optional): Directory the metadata .txt is extracted to
            in a subdirectory with the same name as the image.
        used_bands (list, optional): Bands included in the final image.
        mode (int, optional): Mode used for checking the output paths. Based on
            the `check_path` function.

    Returns:
        None
    """
    
    dtype = 'uint16'
    resolution = 30
    src_nodata = 0
    dst_nodata = 0

    ### Checking paths
    check_path(image_dir, mode)
    if meta_dir:
        check_path(meta_dir, mode)

    tar_file = Path(tar_file)

    ### Printing and timing
    print(f"Processing {tar_file.name}")
    timer = ProcessTimer()
    timer.start()

    if tar_file.match("*.tar.gz"):
        mode = 'r:gz'
    elif tar_file.match("*.tar"):
        mode = 'r:'
    else:
        raise ValueError(f"{tar_file} is not a .tar or .tar.gz file.")

    ### Handling the tar file
    with tarfile.open(tar_file, mode) as t:

        ### Loading in metadata
        meta_name = next((i for i in t.getnames() if i.endswith("_MTL.xml")))
        meta = ET.parse(t.extractfile(meta_name))
        root = meta.getroot()

        ### Collecting name parameters
        date = root.find("IMAGE_ATTRIBUTES/DATE_ACQUIRED").text.replace('-', '')
        sat = f'L{root.find("IMAGE_ATTRIBUTES/SPACECRAFT_ID").text.split("_", 1)[1]}'
        tier = root.find("PRODUCT_CONTENTS/COLLECTION_CATEGORY").text
        path = root.find("IMAGE_ATTRIBUTES/WRS_PATH").text
        row = root.find("IMAGE_ATTRIBUTES/WRS_ROW").text
        refl = 'boa' if root.find("PRODUCT_CONTENTS/PROCESSING_LEVEL").text == 'L2SP' else 'toa'

        output_name = f"{date}_{sat}_{tier}_{path}_{row}_{refl}_eov"
        image_path = Path(image_dir, output_name + ".tif")

        ### Getting the ID of the source image
        source = root.find("PRODUCT_CONTENTS/LANDSAT_PRODUCT_ID").text

        ### Loading bands into memory and stacking them into an array
        if sat in ['L8', 'L9']:
            if not used_bands:
                used_bands = ['B1', 'B2', 'B3', 'B4', 'B5', 'B6', 'B7']
            bands = {i: L89_BANDS[i] for i in used_bands}
        
        elif sat in ['L4', 'L5', 'L7']:
            if not used_bands:
                used_bands = ['B1', 'B2', 'B3', 'B4', 'B5', 'B7']
            bands = {i: L47_BANDS[i] for i in used_bands}

        image_list = [fnmatch.filter(t.getnames(), f"*_{i}.TIF")[0] for i in bands]
        band_list = [image_to_array(t.extractfile(i), 1) for i in image_list]
        src_array = np.stack(band_list, axis=0)
        src_shape = src_array.shape ### bands, height, width

        ### Extracting metadata
        if meta_dir:
            t.extract(meta_name, Path(meta_dir, output_name))

    bounding_box = Polygon([
        (
            float(root.find(f"PROJECTION_ATTRIBUTES/CORNER_{i}_PROJECTION_X_PRODUCT").text),
            float(root.find(f"PROJECTION_ATTRIBUTES/CORNER_{i}_PROJECTION_Y_PRODUCT").text)
        )
        for i in ['UL', 'UR', 'LR', 'LL']
    ])

    utm = f"EPSG:{32600 + int(root.find('PROJECTION_ATTRIBUTES/UTM_ZONE').text)}"

    src_transform = rasterio.transform.from_bounds(*(bounding_box.bounds), width=src_shape[2], height=src_shape[1])

    dst_transform, dst_width, dst_height = calc_transform(bounding_box, resolution, utm)

    dst_array, dst_transform = rasterio.warp.reproject(
        source=src_array,
        src_transform=src_transform,
        src_crs=utm,
        src_nodata=src_nodata,

        destination=np.zeros((src_shape[0], dst_height, dst_width), dtype),
        dst_transform=dst_transform,
        dst_crs=EOV,
        dst_nodata=dst_nodata,
        dst_resolution=(resolution, resolution),

        resampling=Resampling.nearest,
        num_threads=4,
        warp_mem_limit=0,
    )

    profile = GTIFF_UINT16(
        width=dst_width,
        height=dst_height,
        count=src_shape[0],
        transform=dst_transform,
        nodata=dst_nodata,
    )

    ### Writing output image
    with rasterio.open(image_path, 'w', **profile) as dst:
        dst.write(dst_array)
        dst.descriptions = [f"{k} {v}" for k, v in bands.items()]
        dst.update_tags(
            date=date,
            sat=sat,
            path=path,
            refl=refl,
            source=source,
            comp_pred=profile['predictor'],
            comp_level=profile['zstd_level']
        )

    timer.stop(f"Created {output_name}.tif")

def s2_tile(zip_file, image_dir, meta_dir=None, incl_gt=False, used_bands=None, mode=0):
    """Exports a GeoTIFF and a metadata from a Sentinel-2 archive.

    Args:
        zip_file (path): A Sentinel 2 archive. Has to be a .zip file.
        image_dir (path): Directory the GeoTIFF is extracted to.
        meta_dir (path, optional): Directory the metadata .xml is extracted to
            in a subdirectory with the same name as the image.
        incl_gt (bool, optional): If true, the image's name will include
            generation time.
        used_bands (list, optional): Bands included in the final image.
        mode (int, optional): Mode used for checking the output paths. Based on
            the `check_path` function.

    Returns:
        None
    """

    if not used_bands:
        used_bands = ['B01', 'B02', 'B03', 'B04', 'B05', 'B06', 'B07', 'B08', 'B8A', 'B11', 'B12']

    bands = {i: S2_BANDS[i] for i in used_bands}
    band_idxs = [list(S2_BANDS.keys()).index(i) for i in used_bands] ### Indexing from 0

    dtype = 'uint16'
    resolution = 10
    src_nodata = 0
    dst_nodata = 0

    ### Checking paths
    check_path(image_dir, mode)
    if meta_dir:
        check_path(meta_dir, mode)

    zip_file = Path(zip_file)

    ### Printing and timing
    print(f"Processing {zip_file.name}")
    timer = ProcessTimer()
    timer.start()

    if not zip_file.match("*.zip"):
        raise ValueError(f"{zip_file} is not a .zip file.")

    ### Handling the zip file
    with zipfile.ZipFile(zip_file, 'r') as z:
        meta_name = fnmatch.filter(z.namelist(), "*MTD_MSI???.xml")[0]
        meta = ET.parse(z.open(meta_name))
        root = meta.getroot()
        namespace = '{' + root.attrib['{http://www.w3.org/2001/XMLSchema-instance}schemaLocation'] + '}'
        xml_path = f"{namespace}General_Info/Product_Info"

        ### Collecting name parameters
        date = root.find(f"{xml_path}/PRODUCT_START_TIME").text[:10].replace('-', '')
        sat = root.find(f"{xml_path}/PRODUCT_URI").text[:3].lower()
        orbit = root.find(f"{xml_path}/PRODUCT_URI").text[33:37].lower()
        tile = root.find(f"{xml_path}/PRODUCT_URI").text[41:44].lower()
        refl = 'boa' if root.find(f"{xml_path}/PROCESSING_LEVEL").text == 'Level-2A' else 'toa'
        gen_time = root.find(f"{xml_path}/PRODUCT_URI").text[44:-5] if incl_gt else '' ### This includes underscore

        output_name = f"{date}_{sat}_{orbit}_{tile}_{refl}_eov{gen_time}"
        image_path = Path(image_dir, output_name + ".tif")

        ### Getting the ID of the source image
        source = root.find(f"{xml_path}/PRODUCT_URI").text[:-5]

        ### Getting the processing baseline and offset values for recalculating to the old baseline
        proc_baseline = float(root.find(f"{xml_path}/PROCESSING_BASELINE").text)
        
        if proc_baseline >= 4:
            xml_path2 = f"{namespace}General_Info/Product_Image_Characteristics" ### Changing xml path

            ### Path differences based on processing level
            if refl == 'boa': ### L2A
                path_end = "/BOA_ADD_OFFSET_VALUES_LIST/BOA_ADD_OFFSET"
            else: ### L1C
                path_end = "/Radiometric_Offset_List/RADIO_ADD_OFFSET"

            offset_path = f"{xml_path2}{path_end}"
            offsets = [float(root.find(f"{offset_path}[@band_id='{i}']").text) for i in band_idxs]

        ### Getting the bands
        if refl == 'boa': ### L2A
            band_list = [fnmatch.filter(z.namelist(), f"*IMG_DATA/*{k}_{v[0]}m.jp2")[0] for k, v in bands.items()]
        else: ### L1C
            band_list = [fnmatch.filter(z.namelist(), f"*IMG_DATA/*{k}.jp2")[0] for k in bands.keys()]

        print(band_list)

        ### Extracting metadata
        if meta_dir:
            zip_info = z.getinfo(meta_name)
            zip_info.filename = root.find(f"{xml_path}/PRODUCT_URI").text[:-5] + ".xml"
            z.extract(zip_info, Path(meta_dir, output_name))

        with rasterio.open(z.open(band_list[0])) as src:
            bounding_box = box(*(src.bounds))
            utm = src.crs

        dst_transform, dst_width, dst_height = calc_transform(bounding_box, resolution, utm)

        profile = GTIFF_UINT16(
            width=dst_width,
            height=dst_height,
            count=len(band_list),
            transform=dst_transform,
            nodata=dst_nodata,
            BIGTIFF='YES',
        )

        ### Writing output image
        with rasterio.open(image_path, 'w', **profile) as dst:
            for i in range(len(band_list)):
                with rasterio.open(z.open(band_list[i])) as src:
                        dst_array, dst_transform = rasterio.warp.reproject(
                            source=src.read(1),
                            src_transform=src.transform,
                            src_crs=utm,
                            src_nodata=src_nodata,

                            destination=np.zeros((dst_height, dst_width), dtype),
                            dst_transform=dst_transform,
                            dst_crs=EOV,
                            dst_nodata=dst_nodata,
                            dst_resolution=(resolution, resolution),

                            resampling=Resampling.nearest,
                            num_threads=4,
                            warp_mem_limit=0,
                        )

                        ### Recalculating baseline 4.0 images
                        if proc_baseline >= 4:
                            dst_array = np.add(dst_array, offsets[i])
                            dst_array[dst_array < 0] = 0

                        dst.write(dst_array.astype(dtype), i+1)

            ### Adding band names and tags
            dst.descriptions = [f"{k} {v[1]}" for k, v in bands.items()]
            dst.update_tags(
                date=date,
                sat=sat,
                path=orbit,
                refl=refl,
                source=source,
                comp_pred=profile['predictor'],
                comp_level=profile['zstd_level']
            )

    timer.stop(f"Created {output_name}.tif")

def reproj_tile(input_paths, image_dir, meta_dir=None, ls_kwargs=None, s2_kwargs=None, recursive=False, mode=0):
    """Processes image archives. Handles Sentinel-2 and Landsat 4-9 images.

    Args:
        input_paths (path or list): Archive or a list or direcotry of archives.
        image_dir (path): Output directory for GeoTIFF images.
        meta_dir (path, optional): Output directory for metadata files.
        ls_kwargs (dict, optional): Keyword args used for Landsat processing.
        s2_kwargs (dict, optional): Keyword args used for Sentinel-2 processing.
        recursive (bool, optional): Whether to scan the input paths recursively.
        mode (int, optional): Mode used for checking the output paths. Based on
            the `check_path` function.

    Returns:
        None
    """

    file_list = load_files(input_paths, '*.*', recursive)

    tar_list = fnmatch.filter(file_list, "*.tar.gz") + fnmatch.filter(file_list, "*.tar")
    zip_list = fnmatch.filter(file_list, "*.zip")

    ### Initializing process counting and timing
    print("\nProcessing tiles:")
    timer = ProcessTimer()
    timer.start("\nTiles processed")
    len_full = len(tar_list + zip_list)
    count = 0

    for i in tar_list:
        count += 1
        print(f"\n{count}/{len_full}")
        ls_tile(i, image_dir, meta_dir, mode=mode, **(ls_kwargs or {}))

    for i in zip_list:
        count += 1
        print(f"\n{count}/{len_full}")
        s2_tile(i, image_dir, meta_dir, mode=mode, **(s2_kwargs or {}))

    timer.stop()

def make_extras(input_paths, ovr=True, aux=True, recursive=False):
    """Creates external overview and statistics files to images.

    Args:
        input_paths (path or list): Images to work on.
        ovr (bool, optional): Whether to create external overviews.
        aux (bool, optional): Whether to create external statistics files.
        recursive (bool, optional): Whether to scan the input paths recursively.

    Returns:
        None
    """

    image_list = load_files(input_paths, '*.tif', recursive)

    print("\nMaking external files:")
    timer = ProcessTimer()
    subtimer = ProcessTimer()
    timer.start("\nExternal files made")

    if ovr or aux:
        for i in image_list:
            print(f"\nProcessing {i.name}")
            subtimer.start(f"Extras made")
            if ovr:
                make_ovr(i)
            if aux:
                make_aux(i)

    subtimer.stop()
    timer.stop()


def merge_datatake(input_paths, image_output_dir, meta_dir=None, meta_output_dir=None, ovr=True, aux=True, recursive=False, mode=0):
    """Merges individual images by satellite, date and datatake.

    Args:
        input_paths (path or list): Images to work on.
        image_output_dir (path): Output directory for merged images.
        meta_dir (path, optional): Directory containing metadata files in
            subdirectories matching image names.
        meta_output_dir (path, optional): Output directory for collected
            metadata files.
        ovr (bool, optional): Whether to create external overviews.
        aux (bool, optional): Whether to create external statistics files.
        recursive (bool, optional): Whether to scan the input paths recursively.
        mode (int, optional): Mode used for checking the output paths. Based on
            the `check_path` function.

    Returns:
        None
    """

    ### Checking paths
    check_path(image_output_dir, mode)
    image_list = load_files(input_paths, '*.tif', recursive)

    pattern = re.compile(r".*(\d{8}_(?:L[45789]|s2[ab])(?:_T\d)?_(?:\d{3}|r\d{2,3}))_(?:\d{2}|\w{3})_(boa_eov)(?:_\d{8}T\d{6})?.tif")
    matches = [pattern.fullmatch(str(i)) for i in image_list]
    affixes = {i.group(1, 2) for i in matches if i} ### Extracting groups 1 and 2 from the match object

    ### Initializing process counting and timing
    print("\nMerging Datatakes:")
    timer = ProcessTimer()
    subtimer = ProcessTimer()
    timer.start("\nDatatakes merged")
    len_affixes = len(affixes)
    count = 0

    ### Merging the individual images into a single datatake
    for a in affixes:

        ### Creating path name
        filter_ex = r".*" + a[0] + r"_(\d{2}|\w{3})_" + a[1] + r"(_\d{8}T\d{6})?.tif"
        image_pattern = re.compile(filter_ex)
        images = [i for i in image_list if image_pattern.fullmatch(str(i))]
        dst_name = f"{a[0]}_{a[1]}.tif"
        dst_path = Path(image_output_dir, dst_name)

        ### Process counting and timing
        count += 1
        subtimer.start(dst_name)
        print(f"\n{count}/{len_affixes}")
        print(f"Merging {dst_name}")

        ### Collecting metadata files
        if meta_dir and meta_output_dir:
            meta_dir = Path(meta_dir)
            check_path(meta_output_dir, mode)

            meta_pattern = re.compile(filter_ex[:-4])
            meta_subdirs = [i for i in meta_dir.iterdir() if meta_pattern.fullmatch(str(i))] ### The filter expression without the '.tif'
            meta_output_path = Path(meta_output_dir, dst_name[:-4])
            check_path(meta_output_path, 1)

            for m in meta_subdirs:
                meta = next(m.iterdir()) ### Only copies the first file it finds
                copy(meta, meta_output_path / meta.name)

        ### Merging the images
        opened_images = [rasterio.open(i) for i in images]

        try:
            sources = [i.tags()['source'] for i in opened_images] ### Storing the source image IDs
        except:
            sources = []

        try:
            dst_desc = opened_images[0].descriptions
        except:
            dst_desc = []
        
        try:
            dst_tags = {i: opened_images[0].tags()[i] for i in ['date', 'sat', 'path', 'refl']}
        except:
            dst_tags = {}

        bounds = list(zip(*[i.bounds for i in opened_images]))
        bounding_box = box(min(bounds[0]), min(bounds[1]), max(bounds[2]), max(bounds[3]))

        dst_transform, dst_width, dst_height = calc_transform(bounding_box, opened_images[0].res[0])

        profile = opened_images[0].profile
        profile.update(
            GTIFF_UINT16(
                transform=dst_transform,
                width=dst_width,
                height=dst_height,
                BIGTIFF='YES',
            )
        )

        dst_tags.update(
            comp_pred=profile['predictor'],
            comp_level=profile['zstd_level']
        )

        dst_count = opened_images[0].count
        with rasterio.open(dst_path, 'w', **profile) as dst:

            ### Merging images
            for i in range(dst_count):
                dst_array, dst_transform = merge(opened_images, indexes=[i+1])
                ### https://rasterio.readthedocs.io/en/latest/api/rasterio.merge.html
                dst_array = dst_array[0, :, :]
                dst.write(dst_array, i+1)
                del dst_array
                print(f"Bands wirtten: {i+1}/{dst_count}", end='\r')
            print()

            ### Closing opened images
            for i in opened_images:
                i.close()

            ### Attaching band names and metadata tags to the datatake
            dst.descriptions = dst_desc
            dst.update_tags(**dst_tags)
            if len(sources) == len(images):
                dst.update_tags(sources=', '.join(sources))

        ### Making additional files
        if ovr:
            make_ovr(dst_path)
        if aux:
            make_aux(dst_path)

        subtimer.stop()

    timer.stop()

def preproc(input_paths, image_output_dir, meta_output_dir=None, temp_dir=None, recursive=False, mode=0):
    """Handles the full preprocessing of Landsat 4-9 and Sentinel 2 archives.

    It first creates tiled GeoTIFF images and metadata files in a temp
    directory and then merges and collects them respectively and places them
    into the output directories.

    Args:
        input_paths (path or list): List or directory of archives.
        image_output_dir (path): Output directory for merged images.
        meta_output_dir (path, optional): Output directory for collected
            metadata files.
        temp_dir (path, optional): The specific directory for temporary files.
        recursive (bool, optional): Whether to scan the input paths recursively.
        mode (int, optional): Mode used for checking the temp and output paths.
            Based on the `check_path` function.

    Returns:
        None
    """

    timer = ProcessTimer()
    timer.start("\nFull preprocessing")

    with TemporaryDirectory(dir=temp_dir) as tmp:
        reproj_tile(input_paths, tmp, tmp, s2_kwargs={'incl_gt': True}, recursive=recursive, mode=mode)
        merge_datatake(tmp, image_output_dir, tmp, meta_output_dir, mode=mode)

    timer.stop()

def reformat(input_paths, output_dir, compatible=False, suffix='', recursive=False, mode=0):
    """Reformats images based on the current format standards.

    To avoid any unintended consequences, make sure that the input files and
    the initial contents of the output directory don't have matching file names.
    Consider using the `suffix` arg.

    Args:
        input_paths (path or list): Images to reformat.
        output_dir (path): Directory to put the reformated files. These will
            have the same names as the originals.
        compatible (bool, otpional): Whether to use a compatible format or the
            most recent one.
        suffix (str, optional): A suffix to append to filenames.
        recursive (bool, optional): Whether to scan the input paths recursively.
        mode (int, optional): Mode used for checking the temp and output paths.
            Based on the `check_path` function.

    Returns:
        None
    """

    check_path(output_dir, mode)
    image_list = load_files(input_paths, '*.tif', recursive)

    print("\nReformating images:")
    timer = ProcessTimer()
    subtimer = ProcessTimer()
    timer.start("\nImages reformated")
    len_image_list = len(image_list)
    count = 0

    versions = [GTIFF_UINT16, GTIFF_UINT16_COMPATIBLE]
    levels = ['zstd_level', 'zlevel']

    for i in image_list:
        dst_path = Path(output_dir, f"{i.stem}_{suffix}.tif")
        count += 1
        print(f"{count}/{len_image_list}")
        print(i.name)
        subtimer.start("Writing time")

        with rasterio.open(i) as src:
            profile = src.profile
            profile.update(**versions[compatible](), BIGTIFF='YES')
            dst_count = src.count
            with rasterio.open(dst_path, 'w', **profile) as dst:
                for j in range(1, dst_count+1):
                    array = src.read(j)
                    dst.write(array, j)

                dst.update_tags(**src.tags())
                dst.update_tags(
                    comp_pred=profile['predictor'],
                    comp_level=profile[levels[compatible]]
                )

        subtimer.stop()

    timer.stop()

### Data pool management

def op_cc(a, b):
    """Contains operator function."""
    return a.str.contains(b, regex=False)

def op_nc(a, b):
    """Excludes operator function."""
    return ~op_cc(a, b)

def filter_images(file_list, expressions):
    """Filters a list of images based on tags.

    The expressoin should follow the `"<tag> <operator> <value>"` format eg.
    `"date >= 20220101"`. Beside the usual `==`, `!=`, `>`, `<`, `>=` and `<=`
    operators, the `cc` and `!c` operators are also available, standing in for
    "contains" and "does not contain" operations.

    Args:
        file_list (list): List of image paths.
        expressions (list): List of str filter expressions.

    Returns:
        filtered_list (list): Lits of images that match the filter conditions.
    """

    ops = {
        '==': operator.eq,
        '!=': operator.ne,
        '>': operator.gt,
        '<': operator.lt,
        '>=': operator.ge,
        '<=': operator.le,
        'cc': op_cc,
        '!c': op_nc,
    }

    tags = pd.DataFrame([get_tags(i) for i in file_list]).astype(pd.StringDtype())

    for i in expressions:
        i_split = i.split(maxsplit=2)
        tags = tags[ops[i_split[1]](tags[i_split[0]], i_split[2])]

    return np.array(file_list)[[i for i in tags.index]]

def list_unique(file_list, tag):
    """Lists unique values for a tag in the images in the list.

    Args:
        file_list (list): List of image paths.
        tag (str): Tag to list values for.

    Returns:
        unique_values (list): List of unique values.
    """

    tags = pd.DataFrame([get_tags(i) for i in file_list])
    return tags[tag].unique()