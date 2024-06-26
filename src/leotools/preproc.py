"""This module contains tools used for preprocessing satellite images."""

from pathlib import Path
from shutil import copy
import fnmatch
from tempfile import TemporaryDirectory
import xml.etree.ElementTree as ET
import re
import operator

import numpy as np
import pandas as pd
import rasterio
from rasterio.enums import Resampling
from shapely.geometry import Polygon, box

from .constants import EOV, GTIFF_UINT16, GTIFF_UINT16_COMPATIBLE
from .basetools import ProcessTimer, check_path, load_files, FileContainer
from .gistools import calc_transform, image_to_array, make_aux, make_ovr, get_tags, get_desc, merge_images

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
    'QA_PIXEL': 'Pixel QA', ### Normally unused
    'QA_RADSAT': 'Radiometric saturation QA', ### Normally unused
    'CLOUD_QA': 'Surface reflectance cloud QA', ### Only in BOA, Normally unused
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
    'QA_PIXEL': 'Pixel QA', ### Normally unused
    'QA_RADSAT': 'Radiometric saturation QA', ### Normally unused
    'QA_AEROSOL': 'Surface reflectance QA', ### Normally unused
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
    'AOT': (10, 'Aerosol Optical Thickness'), ### Normally unused
    'WVP': (10, 'Water Vapor Pressure'), ### Normally unused
    'SCL': (20, 'Scene Classification'), ### Normally unused
    ### TCI band not implemented
}

def ls_tile(input_file, image_dir, meta_dir=None, used_bands=None, dtype='uint16', mode=0):
    """Exports a GeoTIFF and a metadata from a Landsat 4-9 archive.

    Args:
        input_file (path): A Landsat 4-9 archive as a .tar file or directory.
        image_dir (path): Directory the GeoTIFF is extracted to.
        meta_dir (path, optional): Directory the metadata .txt is extracted to
            in a subdirectory with the same name as the image.
        used_bands (list, optional): Bands included in the final image. Bands
            are not duplicated, only the first mention is used.
        dtype (str): Data type of the output image.
        mode (int, optional): Mode used for checking the output paths. Based on
            the `check_path` function.

    Returns:
        None
    """
    
    resolution = 30
    src_nodata = 0
    dst_nodata = 0

    check_path(image_dir, mode)
    if meta_dir:
        check_path(meta_dir, mode)

    input_file = Path(input_file)

    ### Printing and timing
    print(f"Processing {input_file.name}")
    timer = ProcessTimer()
    timer.start()

    ### Handling the tar file
    with FileContainer(input_file) as c:

        ### Loading in metadata
        meta_name = c.filter("*_MTL.xml")[0]
        meta = ET.parse(c.get(meta_name))
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

        image_list = [c.filter(f"*_{i}.TIF")[0] for i in bands]
        band_list = [image_to_array(c.get(i), 1) for i in image_list]
        src_array = np.stack(band_list, axis=0)
        src_shape = src_array.shape ### bands, height, width

        ### Extracting metadata
        if meta_dir:
            c.extract(meta_name, Path(meta_dir, output_name + '.xml'))

    bounding_box = Polygon([
        (
            float(root.find(f"PROJECTION_ATTRIBUTES/CORNER_{i}_PROJECTION_X_PRODUCT").text),
            float(root.find(f"PROJECTION_ATTRIBUTES/CORNER_{i}_PROJECTION_Y_PRODUCT").text)
        )
        for i in ['UL', 'UR', 'LR', 'LL']
    ])

    utm = f"EPSG:{32600 + int(root.find('PROJECTION_ATTRIBUTES/UTM_ZONE').text)}"

    src_transform = rasterio.transform.from_bounds(*(bounding_box.bounds), width=src_shape[2], height=src_shape[1])

    dst_transform, dst_width, dst_height = calc_transform(bounding_box, resolution, utm, EOV, 300)

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
        dtype=dtype,
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

def s2_tile(input_file, image_dir, meta_dir=None, used_bands=None, dtype='uint16', mode=0, incl_gt=False, ignore_pb=False):
    """Exports a GeoTIFF and a metadata from a Sentinel-2 archive.

    Args:
        input_file (path): A Sentinel 2 archive as a .zip file or directory.
        image_dir (path): Directory the GeoTIFF is extracted to.
        meta_dir (path, optional): Directory the metadata .xml is extracted to
            in a subdirectory with the same name as the image.
        used_bands (list, optional): Bands included in the final image. Bands
            are not duplicated, only the first mention is used.
        dtype (str): Data type of the output image.
        mode (int, optional): Mode used for checking the output paths. Based on
            the `check_path` function.
        incl_gt (bool, optional): If true, the image's name will include
            generation time.
        ignore_pb (bool): Ignore the processing baseline of the image and the
            associated calculations.

    Returns:
        None
    """

    if not used_bands:
        used_bands = ['B01', 'B02', 'B03', 'B04', 'B05', 'B06', 'B07', 'B08', 'B8A', 'B11', 'B12']

    bands = {i: S2_BANDS[i] for i in used_bands}
    band_idxs = [list(S2_BANDS.keys()).index(i) for i in used_bands] ### Indexing from 0

    resolution = 10
    src_nodata = 0
    dst_nodata = 0

    check_path(image_dir, mode)
    if meta_dir:
        check_path(meta_dir, mode)

    input_file = Path(input_file)

    ### Printing and timing
    print(f"Processing {input_file.name}")
    timer = ProcessTimer()
    timer.start()

    ### Handling the zip file
    with FileContainer(input_file) as c:
        meta_name = c.filter("*MTD_MSI???.xml")[0]
        root = ET.parse(c.get(meta_name)).getroot()
        
        namespace = '{' + root.attrib['{http://www.w3.org/2001/XMLSchema-instance}schemaLocation'] + '}General_Info'
        general_info = ET.parse(c.get(meta_name)).getroot().find(namespace)
        product_info = [i for i in general_info if i.tag.endswith('Product_Info')][0]
        product_uri = [i for i in product_info if i.tag.startswith('PRODUCT_URI')][0].text

        ### Collecting name parameters
        date = product_info.find('PRODUCT_START_TIME').text[:10].replace('-', '')
        sat = product_uri[:3].lower()
        orbit = product_uri[33:37].lower()
        tile = product_uri[41:44].lower()
        refl = 'boa' if product_info.find('PROCESSING_LEVEL').text.startswith('Level-2A') else 'toa'
        gen_time = product_uri[44:-5] if incl_gt else '' ### This includes underscore

        output_name = f"{date}_{sat}_{orbit}_{tile}_{refl}_eov{gen_time}"
        image_path = Path(image_dir, output_name + ".tif")

        ### Getting the ID of the source image
        source = product_uri[:-5]

        offset_paths = [
            "Product_Image_Characteristics/BOA_ADD_OFFSET_VALUES_LIST/BOA_ADD_OFFSET", ### L2A
            "Product_Image_Characteristics/Radiometric_Offset_List/RADIO_ADD_OFFSET", ### L1C
        ]

        offset_values = None

        for i in offset_paths:
            if general_info.find(i) is not None:
                offset_elements = [general_info.find(f"{i}[@band_id='{j}']") for j in band_idxs]
                offset_values = [float(j.text) if j is not None else 0.0 for j in offset_elements]
                break

        ### Getting the bands
        if refl == 'boa': ### L2A
            band_list = [c.filter(f"*IMG_DATA/R{v[0]}m/*{k}_{v[0]}m.jp2")[0] for k, v in bands.items()]
        else: ### L1C
            band_list = [c.filter(f"*IMG_DATA/*{k}.jp2")[0] for k in bands.keys()]
        
        ### Extracting metadata
        if meta_dir:
            c.extract(meta_name, Path(meta_dir, output_name + '.xml'))

        with rasterio.open(c.get(band_list[0])) as src:
            bounding_box = box(*(src.bounds))
            utm = src.crs

        dst_transform, dst_width, dst_height = calc_transform(bounding_box, resolution, utm, EOV, 300)

        profile = GTIFF_UINT16(
            dtype=dtype,
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
                with rasterio.open(c.get(band_list[i])) as src:
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
                        if offset_values and not ignore_pb:
                            dst_array = np.add(dst_array, offset_values[i])
                            dst_array[dst_array<0] = 0

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

def reproj_tiles(input_paths, image_dir, meta_dir=None, ls_kwargs=None, s2_kwargs=None, recursive=False, mode=0):
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
    tar_list = fnmatch.filter(file_list, '*.tar.gz') + fnmatch.filter(file_list, '*.tar')
    zip_list = fnmatch.filter(file_list, '*.zip')

    ### Initializing process counting and timing
    print("\nProcessing tiles:")
    timer = ProcessTimer()
    timer.start("\nTiles processed")
    len_full = len(tar_list + zip_list)
    count = 0

    for i in tar_list:
        count += 1
        print(f"\n[{count}/{len_full}]")
        ls_tile(i, image_dir, meta_dir, mode=mode, **(ls_kwargs or {}))

    for i in zip_list:
        count += 1
        print(f"\n[{count}/{len_full}]")
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

def merge_datatakes(input_paths, image_output_dir, meta_dir=None, meta_output_dir=None, ovr=True, aux=True, recursive=False, mode=0):
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

    check_path(image_output_dir, mode)
    if meta_dir and meta_output_dir:
        check_path(meta_output_dir, mode)
    
    ### Finding related images
    image_list = load_files(input_paths, '*.tif', recursive)
    pattern = re.compile(r".*(\d{8}_(?:L[45789]|s2[ab])(?:_T\d)?_(?:\d{3}|r\d{2,3}))_(?:\d{2}|\w{3})_([bt]oa_eov)(?:_\d{8}T\d{6})?.tif")
    matches = [pattern.fullmatch(str(i)) for i in image_list]
    affixes = {i.group(1, 2) for i in matches if i} ### Extracting groups 1 and 2 from the match object

    ### Initializing process counting and timing
    print("\nMerging Datatakes:")
    timer = ProcessTimer()
    subtimer = ProcessTimer()
    timer.start("\nDatatakes merged")
    len_affixes = len(affixes)
    count = 0

    ### Merging related images into a single datatake
    for a in affixes:

        ### Assembling path name
        filter_expression = r".*" + a[0] + r"_(\d{2}|\w{3})_" + a[1] + r"(_\d{8}T\d{6})?.tif"
        image_pattern = re.compile(filter_expression)
        images = [i for i in image_list if image_pattern.fullmatch(str(i))]
        dst_name = f"{a[0]}_{a[1]}.tif"
        dst_path = Path(image_output_dir, dst_name)

        ### Process counting and timing
        count += 1
        subtimer.start(dst_name)
        print(f"\n[{count}/{len_affixes}]")
        print(f"Merging {dst_name}")

        ### Collecting metadata files
        if meta_dir and meta_output_dir:
            meta_dir = Path(meta_dir)

            meta_pattern = re.compile(filter_expression[:-4]+'.xml') ### The filter expression without the '.tif'
            meta_files = [i for i in meta_dir.iterdir() if meta_pattern.fullmatch(str(i))]
            meta_output_path = Path(meta_output_dir, dst_name[:-4])
            check_path(meta_output_path, 1)

            for i in meta_files:
                copy(i, meta_output_path / i.name)

        profile = GTIFF_UINT16(BIGTIFF='YES')

        merge_images(images, dst_path, profile, 300, verbose=True)

        ### Attaching band names and metadata tags to the datatake
        with rasterio.open(dst_path, 'r+') as dst:
            tags = [get_tags(i) for i in images]
            incomplete = False

            try:
                dst.descriptions = get_desc(images[0])
            except:
                incomplete = True

            try:
                dst.update_tags(
                    **{i: tags[0][i] for i in ['date', 'sat', 'path', 'refl']},
                    comp_pred=profile['predictor'],
                    comp_level=profile['zstd_level'],
                )
            except:
                incomplete = True

            try:
                dst.update_tags(sources=', '.join([i['source'] for i in tags]))
            except:
                incomplete = True

            if incomplete:
                print("WARNING: Incomplete metadata.")

        ### Creating external files
        if ovr:
            print("Creating external overview")
            make_ovr(dst_path)
        if aux:
            print("Creating statistics file")
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
        reproj_tiles(input_paths, tmp, tmp, s2_kwargs={'incl_gt': True}, recursive=recursive, mode=mode)
        merge_datatakes(tmp, image_output_dir, tmp, meta_output_dir, mode=mode)

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

    ### Process counting and timing
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
        print(f"[{count}/{len_image_list}]")
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