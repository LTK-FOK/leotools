"""This module contains tools directly handling geoinformatic processes."""

from pathlib import Path
from math import ceil

import numpy as np
from osgeo import gdal
import rasterio
from rasterio.mask import mask
from rasterio.merge import merge
from rasterio.enums import Resampling
from rasterio.warp import aligned_target, reproject
import geopandas as gpd
from shapely.geometry import Polygon, box

from .basetools import check_path, load_files, ProcessTimer

### Metadata manipulation

def get_profile(image):
    """Returns the rasterio profile of an image as a dictionary."""
    with rasterio.open(image) as src:
        return src.profile

def get_tags(image):
    """Returns the GeoTiff tags of an image."""
    with rasterio.open(image) as src:
        return src.tags()

def get_desc(image):
    """Returns the decriptions of the image bands."""
    with rasterio.open(image) as src:
        return src.descriptions

def strip_tags(image):
    """Removes the GeoTIFF tags from an image."""
    ds = gdal.Open(image, gdal.GA_Update)
    ds.SetMetadata({})
    del ds

def get_file_size(image):
    """Returns file size in gigabytes."""
    size = Path(image).stat().st_size #>> 20
    return (size/1024**3)

def get_read_time(image):
    """Returns the read time of an image."""
    timer = ProcessTimer(verbose=False)
    timer.start()
    with rasterio.open(image) as src:
        src.read()
    return timer.stop()

def print_profile(image, read_time=True):
    """Prints the rasterio profile and additional statistics of an image.
    
    Args:
        image (path): The image to print profile and tags for.
        read_time (bool, optional): If True (default), prints reading time.

    Returns:
        None
    """
    
    image = Path(image)

    print(f"\n{image.name}\n")

    with rasterio.open(image) as src:

        ### Printing profile
        for k, v in src.profile.items():
            print(f"{k}:  {v}")

        print()

        ### Printing tags
        for k, v in src.tags().items():
            print(f"{k}:  {v}")

    ### Printing file size
    print(f"\nFile size: {get_file_size(image):.3f} GB")
    
    ### Measuring reading time
    if read_time:
        print(f"\nReading time: {get_read_time(image)}")

def strip_proj(image):
    """Removes the projection from the metadata of an image."""
    ds = gdal.Open(image, gdal.GA_Update)
    ds.SetProjection('')
    del ds

def set_nodata(image, nodata):
    """Sets the nodata value of an image to the specified number."""
    with rasterio.open(image, 'r+') as src:
        src.nodata = nodata

### Importing and exporting arrays

def image_to_array(image, bands=None):
    """Loads selected bands of an image into a numpy array.

    The dimensions of the output array are `count, rows(y), cols(x)`, or
        `rows(y), cols(x)` for a single-band image.

    Args:
        image (path): The image to load.
        bands (int or list, optional): Selected bands. Bands can be omitted,
            duplicated and rearranged.

    Returns
        array (ndarray): The loaded array.
    """

    with rasterio.open(image) as src:
        return src.read(bands)

def array_to_image(array, output, sample_image=None, profile=None, mode=0):
    """Exports a numpy array to an image with the properties of a sample file.

    The band count of the array is automatically detected. Dimensions should be
    `count, rows(y), cols(x)`, or `rows(y), cols(x)` for a single-band image.

    Args:
        array (ndarray): The array to export.
        output (path): Path to the output image.
        sample_image (path, optional): Path to the sample image.
        profile (dict, optional): A dict of additional metadata parameters.
        mode (int, optional): Mode used for checking the output paths. Based on
            the `check_path` function.

    Returns:
        None
    """

    ### https://rasterio.readthedocs.io/en/latest/topics/image_processing.html#imageorder
    
    check_path(Path(output).parent, mode)
    shape = array.shape
    len_shape = len(shape)

    if not sample_image and not profile:
        raise ValueError("Neither sample_tif nor profile were given.")
    
    if sample_image:
        sample_profile = get_profile(sample_image)
    else:
        sample_profile = {}
    
    if not profile:
        profile = {}

    ### Multi-band image shape = count, rows(y), cols(x)
    if len_shape == 3:
        shape_profile = {
            'count': shape[0],
            'width': shape[2],
            'height': shape[1]
        }

    ### Single-band image shape = rows(y), cols(x)
    elif len_shape == 2:
        shape_profile = {
            'count': 1,
            'width': shape[1],
            'height': shape[0]
        }
    
    ### Last resort
    else:
        raise ValueError("The array should be 2 or 3 dimensional.")
    
    dst_profile = {**sample_profile, **profile, **shape_profile}
    
    with rasterio.open(output, 'w', **dst_profile) as dst:
        if len_shape == 3:
            dst.write(array)
        elif len_shape == 2:
            dst.write(array, 1)
        else:
            raise ValueError("The array should be 2 or 3 dimensional.")

### Array operations

def remap_array(array, map_dict, default_value=None):
    """Replaces values in an array with a dictionary.

    Args:
        array (array): The array to remap.
        map_dict (dict): The dictionary that contains the remapping instructions.
        default_value (number, optional): Replace values that are not assigned
            in the dict. Keeps the original value if not given.

    Returns:
        map_array (ndarray): The remapped array.
    """

    if default_value is None:
        new_array = np.copy(array)
    else:
        new_array = np.full(array.shape, default_value)

    for k, v in map_dict.items():
        new_array[array == k] = v

    return new_array

def mask_array(array, mask, nodata=0, invert=False):
    """Replace values in the array with the nodata value where the mask is 0.
    
    Args:
        array (ndarray): The array to mask.
        mask (ndarray): The mask that decides which values to mask out.
        nodata (number, optional): The value which all masked values get.
        invert (bool, optonal): If True, True values in the mask will be
            replaced with nodata, If False (default), the opposite.

    Returns:
        masked_array (ndarray): The masked array.
    """
    
    masked_array = array.copy()

    if invert:
        masked_array[mask.astype('bool')] = nodata
    else:
        masked_array[~mask.astype('bool')] = nodata

    return masked_array

def stack_arrays(arrays):
    """Stacks multiple 2D or 3D arrays into a single 3D array.
    
    Args:
        arrays (list): The arrays to stack.

    Returns:
        stacked_array (array): The stacked array.
    """

    # arrays = [i for i in arrays if len(i.shape) == 2 else np.expand_dims(i, 0)]
    new_arrays = []
    for i in arrays:
        len_shape = len(i.shape)

        ### Multi-band image shape = count, rows(y), cols(x)
        if len_shape == 3:
            new_arrays += [i]

        ### Single-band image shape = rows(y), cols(x)
        elif len_shape == 2:
            new_arrays += [np.expand_dims(i, 0)]
        
        ### Last resort
        else:
            raise ValueError("The arrays sould be 2 or 3 dimensional.")

    return np.concatenate(new_arrays)

### Band operations and NDXI calculation

class BandSelector:
    """Wrapper class for rasterio's DatasetReader that allows indexing of bands.

    With `b[1]`, you'd get the first band of `b` as an array. `b` is an opened
    rasterio image. This is equivalent to `b.read(1)`.
    """
    def __init__(self, src):
        self.src = src

    def __getitem__(self, index):
        return self.src.read(index) ### Indexing from 1

### Raster operations

def ndiff(e1, e2):
    """Normalized difference function. Returns x, where x = (e1-e2) / (e1+e2)."""
    return (e1-e2) / (e1+e2)

def extract_bands(image, output_dir, bands=None, dtype=None, mode=0):
    """Extract bands or band combinations of an image.

    Args:
        image (path): The image that contains the bands.
        output_dir (path): The directory to put the output images into.
        bands (list or dict, optional): The bands and band combinations to
            output.
        dtype (str, optional): The datatype of the output images. Defaults to
            the datatype of the input.
        mode (int, optional): Mode used for checking the output paths. Based on
            the `check_path` function.

    Returns:
        None
    """

    check_path(output_dir, mode)
    image = Path(image)
    output_dir = Path(output_dir)

    with rasterio.open(image) as src:
        profile = src.profile
        b = BandSelector(src)

        if dtype:
            profile['dtype'] = dtype

        ### If bands arg is empty get every band of the image
        if not bands:
            bands = [i+1 for i in range(src.count)]

        ### If bands are given as a list
        if isinstance(bands, list):
            for i in bands:
                out_path = output_dir / f"{image.stem}_b{i}.tif"
                array = b[i].astype(profile['dtype'])
                array_to_image(array, out_path, profile=profile)

        ### If bands are given as a dictionary
        elif isinstance(bands, dict):
            for k, v in bands.items():

                ### Evaluating the expression
                with np.errstate(divide='ignore', invalid='ignore'):
                    array = eval(v, {'__builtins__': None, 'ndiff': ndiff}, {'b': b}).astype(profile['dtype'])
                    
                    ### Filling in nodata values
                    array[np.isinf(array)] = profile['nodata'] or 0
                    array[np.isnan(array)] = profile['nodata'] or 0
            
                ### Exporting the array
                out_path = output_dir / f"{image.stem}_{k}.tif"
                array_to_image(array, out_path, profile=profile)

        else:
            raise TypeError("Bands should be list, dict or None.")

def stack_images(input_paths, output, band=0, bands_csv=None, recursive=False, mode=0, verbose=False):
    """Stacks multiple images into a single image.

    Args:
        input_paths (path or list): A list of the source images.
        output (path): Output image.
        band (int, optional): Which band of the source images to use. Uses all
            bands of each image if 0.
        bands_csv (path, optional): A .csv file listing the source (image and
            band) for each band of the output image.
        recursive (bool, optional): Whether to scan the input directory
            recursively.
        mode (int, optional): Mode used for checking the output paths. Based on
            the `check_path` function.
        verbose (bool, optional): Whether to print current process status.

    Returns:
        None
    """

    check_path(Path(output).parent, mode)
    image_list = load_files(input_paths, '*.tif', recursive)
    dst_count = sum([get_profile(i)['count'] for i in image_list]) if not band else len(image_list)
    profile = get_profile(image_list[0])
    profile.update(count=dst_count)
    first_band = 1

    if verbose:
        print(f"Bands wirtten: 0/{dst_count}", end='\r')
    
    with rasterio.open(output, 'w', **profile) as dst:
        for i in image_list:
            with rasterio.open(i) as src:
                src_count = 1 if band else src.count
                dst_bands = first_band if band else list(range(first_band, first_band+src_count))
                dst.write(src.read(band or None), dst_bands)

                ### Creating csv
                if bands_csv:
                    with open(bands_csv, 'a') as f:
                        lines = [f"{first_band+j},{i},{band or j+1}\n" for j in range(src_count)]
                        f.writelines(lines)

                first_band += src_count

            if verbose:
                print(f"Bands wirtten: {first_band-1}/{dst_count}", end='\r')

        if verbose:
            print()

def merge_images(input_paths, output, profile=None, round=0, recursive=False, mode=0, verbose=False):
    """Merges multiple images into one.

    Args:
        input_paths (path or list): Images to work on.
        output (path): Output image.
        profile (dict, optional): A dict of additional image parameters.
        round (number, optional): The grid to round the extent to. 
        recursive (bool, optional): Whether to scan the input paths recursively.
        mode (int, optional): Mode used for checking the output paths. Based on
            the `check_path` function.
        verbose (bool, optional): Whether to print current process status.

    Returns:
        None
    """

    check_path(Path(output).parent, mode)
    image_list = load_files(input_paths, '*.tif', recursive)
    opened_images = [rasterio.open(i) for i in image_list]

    all_bounds = list(zip(*[i.bounds for i in opened_images]))
    final_bounds = box(min(all_bounds[0]), min(all_bounds[1]), max(all_bounds[2]), max(all_bounds[3]))

    dst_transform, dst_width, dst_height = calc_transform(
        final_bounds,
        opened_images[0].res[0],
        opened_images[0].crs,
        opened_images[0].crs,
        round,
    )

    if not profile:
        profile = {}

    dst_profile = opened_images[0].profile

    dst_profile.update(
        **profile,
        transform=dst_transform,
        width=dst_width,
        height=dst_height,
    )

    count = opened_images[0].count
    with rasterio.open(output, 'w', **dst_profile) as dst:

        ### Merging images
        if verbose:
            print(f"Bands wirtten: 0/{count}", end='\r')

        for i in range(count):
            ### https://rasterio.readthedocs.io/en/latest/api/rasterio.merge.html
            dst_array, dst_transform = merge(opened_images, indexes=[i+1])
            dst_array = dst_array[0, :, :]
            dst.write(dst_array, i+1)
            del dst_array

            if verbose:
                print(f"Bands wirtten: {i+1}/{count}", end='\r')
        
        if verbose:
            print()

        ### Closing opened images
        for i in opened_images:
            i.close()

### Making external files

def make_ovr(image, levels=(2, 8, 32, 128), compress='zstd', compress_level=None, predictor=2):
    """Creates pyramid layers for an image in an external .ovr file.
    
    Args:
        image (path): The image to make an overview for.
        levels (tuple, optional): Overview levels.
        compress (str, optonal): Compression method used.
        compress_level (int, optional): Compression level used. Uses default for
            the given compression if not provided.
        predictor (int): Predictor used.

    Returns:
        None
    """

    with rasterio.Env(
        TIFF_USE_OVR=True,
        COMPRESS_OVERVIEW=compress,
        ZLEVEL_OVERVIEW=compress_level or 6,
        ZSTD_LEVEL_OVERVIEW=compress_level or 9,
        PREDICTOR_OVERVIEW=predictor,
        GDAL_NUM_THREADS=4,
    ):
        with rasterio.open(image, 'r+') as src:
            src.build_overviews(levels, Resampling.average)

def make_aux(image):
    """Calculates image statistics and puts them into an .aux.xml file."""

    Path(image).with_suffix('.tif.aux.xml').unlink(missing_ok=True)
    src = gdal.Open(str(image))

    for i in range(src.RasterCount):
        src.GetRasterBand(i + 1).ComputeStatistics(0)

    src = None

### Vector operations

def load_polygons(polygons):
    """Returns provided vector data as a list of polygons."""
    
    ### If given a path
    if isinstance(polygons, (str, Path)):
        polygons = gpd.read_file(polygons) ### Loads it as GeoDataFrame

    ### If given a geopandas GeoDataFrame
    if isinstance(polygons, gpd.GeoDataFrame):
        polygons = polygons['geometry'] ### Loads it as GeoSeries

    ### If given a geopandas GeoSeries
    if isinstance(polygons, gpd.geoseries.GeoSeries):
        return list(polygons)

    ### if given a shapely Polygon
    elif isinstance(polygons, Polygon):
        return [polygons]
    
    else:
        raise TypeError("Polygons should be paths or geometries.")

def multi_mask(image, output, nodata=0, area_masks=None, bound_masks=None, all_touched=True, invert=False, mode=0):
    """Pixels of the raster that touch the mask shapes are filled with nodata.

    A mask can be:
        - Path to a file containing vector data.
        - A geopandas `GeoSeries` object.
        - A shapely `Polygon` object.
    
    Args:
        image (path): The image to mask.
        output (path): The masked output image.
        nodata (number, optional): Masked areas will get this nodata value.
        area_masks (list, optional): List of masked areas.
        bound_masks (list, optional): List of masked boundaries.
        all_touched (bool, optional): If True (default), every pixel touching
            the geometries is masked out.
        invert (bool, optional): If True, pixels under the geometries are
            discarded, if False (default), they are kept.
        mode (int, optional): Mode used for checking the output paths. Based on
            the `check_path` function.

    Returns:
        None
    """

    check_path(Path(output).parent, mode)

    ### Combining vector mask
    sum_mask = []
    if area_masks:
        for i in area_masks:
            sum_mask += load_polygons(i)

    if bound_masks:
        for i in bound_masks:
            sum_mask += [j.boundary for j in load_polygons(i)]

    if not sum_mask:
        raise ValueError("None of the provided masks contain polygons!")

    ### Applying vector mask
    ### https://rasterio.readthedocs.io/en/latest/api/rasterio.mask.html
    with rasterio.open(image) as src:
        profile = src.profile
        out_img, transform = mask(src, sum_mask, all_touched=all_touched, invert=invert, nodata=nodata)

    profile.update(transform=transform, nodata=nodata)

    with rasterio.open(output, 'w', **profile) as dst:
        dst.write(out_img)

### Coordinate operations

def round_to(x, y):
    """Rounds x to the closest multiple of y.
    
    A positive/negative y rounds up/down respectively.
    """

    return (ceil(x / y)) * y

def round_extent(bounds, grid):
    """Rounds extent bounds to the grid.
    
    Args:
        bounds (list): Extent bounds in this order: xmin, ymin, xmax, ymax.
        grid (int): The grid to round to. A positive/negative number
            expands/shrinks the extent respectively.

    Returns:
        new_bounds (list): The new, rounded extent bounds.
    """
    
    if grid:
        new_bounds = [
            round_to(bounds[0], -grid), # xmin
            round_to(bounds[1], -grid), # ymin
            round_to(bounds[2], +grid), # xmax
            round_to(bounds[3], +grid) # ymax
        ]

        return new_bounds
    
    else: ### grid == 0
        return bounds

def calc_transform(bounds, dst_resolution, src_crs, dst_crs, dst_round=0, tap=True):
    """Calculates source and destination transforms from provided parameters.

    Args:
        bounds (list, tuple or polygon): A series of extent bounds (xmin, ymin,
            xmax, ymax) or a polygon to calculate the extent from.
        dst_resolution (int): Destination resolution in meters.
        src_crs (str): Projection of the source.
        dst_src (str): Projection of the destination.
        dst_round (number, optional): The grid to round the extent to.
        tap (bool, optional): Whether to align pixels to target resolution.

    Returns:
        dst_transform (dict): Destination transform.
        dst_width (int): Width of the destination in pixels.
        dst_height (int): Height of the destination in pixels.
    """

    if isinstance(bounds, (list, tuple)):
        bounds = box(*bounds) ### Convert into polygon

    if isinstance(bounds, Polygon):
        if not src_crs:
            src_crs = dst_crs

        src_extent = gpd.GeoSeries(bounds, crs=src_crs)
        dst_extent = round_extent(src_extent.to_crs(dst_crs).total_bounds, dst_round)

    else:
        raise TypeError("The bounding box needs to be list, tuple or polygon.")

    ### Reprojection
    dst_width = int((dst_extent[2]-dst_extent[0])/dst_resolution)
    dst_height = int((dst_extent[3]-dst_extent[1])/dst_resolution)

    dst_transform = rasterio.transform.from_bounds(
        *dst_extent,
        width=dst_width,
        height=dst_height,
    )

    if tap:
        return aligned_target(dst_transform, dst_width, dst_height, dst_resolution)
    else:
        return dst_transform, dst_width, dst_height

def reproj_image(image, output, resolution, crs=None, round=0, tap=True, mode=0, verbose=False):
    """Reprojects an image.

    Args:
        image (path): Input image. Its projection information has to be correct.
        output (path): Output image.
        resolution (int): Resolution of output image in meters.
        crs (str, optional): Projection of the output image. If not given, uses
            input projection.
        round (number, optional): The grid to round the extent to.
        tap (bool, optional): Whether to align pixels to target resolution.
        mode (int, optional): Mode used for checking the output paths. Based on
            the `check_path` function.
        verbose (bool, optional): Whether to print current process status.

    Returns:
        None
    """

    check_path(Path(output).parent, mode)

    with rasterio.open(image) as src:
        profile = src.profile

        dst_transform, dst_width, dst_height = calc_transform(
            src.bounds,
            resolution,
            src.crs,
            crs or src.crs,
            round,
            tap,
        )
        
        count = src.count

        profile.update(
            transform=dst_transform,
            width=dst_width,
            height=dst_height,
            crs=crs
        )

        with rasterio.open(output, 'w', **profile) as dst:
            if verbose:
                print(f"Bands wirtten: 0/{count}", end='\r')

            for i in range(count):
                dst_array, dst_transform = reproject(
                    source=src.read(i+1),
                    src_transform=src.transform,
                    src_crs=src.crs,
                    src_nodata=profile['nodata'],

                    destination=np.zeros((dst_height, dst_width), profile['dtype']),
                    dst_transform=dst_transform,
                    dst_crs=crs,
                    dst_nodata=profile['nodata'],
                    dst_resolution=(resolution, resolution),

                    resampling=Resampling.nearest,
                    num_threads=4,
                    warp_mem_limit=0,
                )

                dst.write(dst_array, i+1)

                if verbose:
                    print(f"Bands wirtten: {i+1}/{count}", end='\r')
                
            if verbose:
                print()