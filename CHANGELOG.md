# v0.9.0

## Added

- `basetools.format_time` function.
- `gistools.validate_image` function.
- `gistools.cut_image` function.
- `crs` and `suffix` args to core preprocessing functions.
- `sort` arg to `gistools.extract_bands` function.
- `ovr` and `aux` arg to `preproc.reformat` function.
- - Image procession now validates the bands of the output image.

## Changed

- `basetools.timestamp` now returns its value instead of printing it.

## Removed

- `crs` key from `constants.GTIFF_UINT16` and `constants.GTIFF_UINT16_COMPATIBLE`.
- `preproc.make_extras` function.

## Fixed

- Sentinel-2 processing can now deal with S2C images.
- `gistools.round_to` now accepts 0 as `y` arg.

# v0.8.0

## Added

- `basetools.timestamp` function.
- `gistools.get_desc` function.
- `gistools.merge_images` function. It has also been added as a command line function.
- `gistools.reproj_image` function. It has also been added as a command line function.
- `DirInterface`, `ZIPInterface`, `TARInterface` and `FileContainer` classes to `basetools`.
- Path checking for outputs and the `mode` arg to `array_to_image`, `extract_bands`, `stack_images`, `merge_images` and `multi_mask` functions in the `gistools` submodule.
- `dst_crs` arg to `gistools.calc_transform` function.
- `verbose` arg to `gistools.stack_images` function.
- `dtype` arg to `preproc.ls_tile`. `dtype` and `ignore_pb` args to `preproc.s2_tile`.
- QA_PIXEL, QA_RADSAT, QA_CLOUD_QA bands to the list of extractable Landsat 4-7 bands. QA_PIXEL, QA_RADSAT, QA_AEROSOL bands to the list of extractable Landsat 8-9 bands. AOT, WVP, SCL bands to the list of extractable Sentinel-2 bands.

## Changed

- `calc_transform` function was relocated from the `preproc` to the `gistools` module.
- `reproj_tile` and `merge_datatake` functions in the `preproc` submodule have been renamed to `reproj_tiles` and `merge_datatakes` respectively.
- `src_bounding box` arg of the `gistools.calc_transform` function has been renamed to `bounding_box`. Now it also accepts a list of extent bounds as an input.
- `tar_file` arg of the `preproc.ls_tile` function has been renamed to `input_file`.
- `zip_file` arg of the `preproc.l2_tile` function has been renamed to `input_file`.
- `bounding_box` arg of `gistools.calc_transform` has been renamed to `bounds`.
- `num` and `resolution` args of the `gistools.round_to` function have been renamed to `x` and `y`.
- `resolution` arg of the `gistools.round_extent` function has been renamed to `grid`.
- `ls_tile` and `s2_tile` functions don't put metadata xml files into separate folders anymore.
- `src_crs` and `dst_crs` args of the `gistools.calc_transform` function are now required. The default of the `dst_round` arg is now 0 instead of 300.

## Fixed

- `basetools.load_files` can now handle a list of strings as an input.
- `gistools.calc_transform` now does the rounding if a `src_crs` is not given.
- `gistools.BandSelector` now retains the datatype of the original image`

# v0.7.2

## Fixed

- Sentinel-2 processing no longer prints the bands of the images.
- The merging of S2 L1C tiles is now functional.

# v0.7.1

## Added

- `preproc.s2_tile` now handles the processing baseline 4.0 change of Level-1C images.
- The MIT License notice is now included in `LICENSE.txt`.
- `__license__` module level variable.
- More detailed install instructions to the readme.

## Changed

- Dependencies now include all directly imported packages.
- `leotools` directory has been moved under `src`.
- Build metadata has been moved from `setup.py` to `pyproject.toml`.

## Fixed

- Package building and installation using pip now works.
- Functionality of `preproc.op_nc`.
- `ovr` and `aux` args are now properly implemented for `preproc.make_extras`.
- The `preproc` functions handling Landsat images have had their descriptions updated and now correctly say "Landsat 4-9". 
- The functionality of Sentinel-2 L1C processing.

# v0.7.0

## Added

- An external documentation in Hungarian in the README.md. English section has been removed.
- `levels`, `compress`, `compress_level`, `predictor` args to `gistools.make_ovr` function.
- `gistools.load_polygons` function.
- `all_touched` arg to `gistools.multi_mask` function.
- `verbose` init arg and `pause` method to `basetools.ProcessTimer` class.
- `get_tags`, `get_file_size` and `get_read_time` functions to the `gistools` module.
- `gistools.stack_arrays` function.
- `recursive` and `mode` args to `preproc.reformat`
- `default_value` arg to `remap_array`.
- `gistools.stack_images` to `__main__` as a command line function.
- `gistools.strip_tags` function.
- `op_cc`, `op_nc`, `filter_images` and `list_unique` functions to `preproc` module.
- `used_bands` arg to `preproc.ls_tile` and `preproc.s2_tile` functions.
- `ls_kwargs` and `s2_kwargs` args to `preproc.reproj_tile` function.

## Changed

- `basetools.load_inputs` has been renamed to `basetools.load_files`. It has also had its functionality overhauled. Now it only returns paths that exist.
- `gistools.recode_array` has been renamed to `gistools.remap_array`. Its `rec_dict` arg has been renamed to `map_dict`.
- The `all_bands` arg of `gistools.stack_images` has been replaced by the `band` arg, which decides which band to stack.
- The `make_ovr` and `make_aux` functons were relocated from the `preproc` to the `gistools` module.
- The `gistools.multi_mask` function can now handle shapely Polygons and geopandas GeoSeries but does not use the `load_files` function.
- The default value of the `invert` arg of `gistools.multi_mask` function is now `False`.
- The `invert` tag of `gistools.mask_array` works in reverse now.
- Sentinel-2 preprocessing now recalculates image values from the 4.0 processing baseline to the former 03.01 baseline.

## Removed

- `basetools.load_params` and `basetools.interpret_params` functions. Functions using these have the functionality integrated into them. If you need to load in a dictionary from a plaintext file, use the `yaml` module.
- `preproc.extract_ndxi` function and the associated commandline option.
- Redundant `gistools.band_calc` function.
- `gistools.recalc_dtype` function. Change the values and datatpye of the array in code. Also removed `multiply` and `add` args from `gistools.extract_bands`.

## Fixed

- `basetools.ProcessTimer` is now imported into the `gistools` module.
- `gistools.BandSelector` now casts read band to float, making its behavior more consistent.

# v0.6.1

## Fixed

- `image_to_array` error in Landsat processing.

# v0.6.0

## Added

- `__version__`, `__author__`, `__email__`, `__description__` and `__license__` module level variables.
- Support for the preprocessing of Landsat 4, 5 and 7 images.
- The commandline help function now raises a ValueError if asked to help with a function that doesn't exist in the catalog.
- Package building now includes testing.
- `dtype`, `multiply` and `add` args to `gistools.extract_bands`.
- `suffix` arg to `preproc.reformat`.
- `info` to the commandline function catalog.
- `bands_csv`and `recursive` args to `gistools.stack_images`.
- The original band designations are now added to preprocessed images as band descriptions.

## Changed

- `preproc.merge_datatake` now uses `basetools.load_inputs` for inputs.
- `gistools.stack_images` now uses `basetools.load_inputs`  for inputs.
- The `sample_tif` arg of `gistools.array_to_image` is now `sample_image`.
- Certain functions got updated docstrings.

## Fixed

- The `multi_mask` function now works if area an bound masks are not given.
- `gistools.extract_bands` now evaluates expressions safely.

# v0.5.2

## Fixed

- Landsat 9 is now supported by the `preproc.merge_datatake` and `preproc.extract_ndxi`.

# v0.5.1

## Added

- Support for the preprocessing of Landsat 9 images.

# v0.5.0

## Changed

- leotools is now a proper Python package.
- The `leotools` command is now an entry point of the package, it can be used from anywhere as long as the installation environment is active.
- Function kwargs in the commandline can now be specified with one or two hyphens.

## Removed

- .conf files used for ndxi calculation. Now this data resides in the `constants` submodule.

# v0.4.4

## Changed

- inside the `__main__` module, the `main` functions handles the functions instead of the old `run_func` function.

## Fixed

- `preproc.reformat` now runs without an error, correctly copies metadata for images.

# v0.4.3

## Added

- `preproc.reformat` now has the `compatible` arg, offering an alternate compression.

# v0.4.2

## Fixed

- `preproc.preproc` creates ovr files correctly.

# v0.4.1

## Fixed

- `preproc.reproj_tiles` input.

# v0.4.0

## Added

- `gistools.print_profile` function. It can be accessed from the commandline as `profile`.
- `basetools.load_inputs` function. This extends the ways input files can be provided to certain functions.
- Additional functionality to `gistools.stack_images` function.
- `preproc.make_ovr`, `preproc.make_aux` and `preproc.make_extras` functions.
- `preproc.reproj_tile` now has an `incl_tg` arg, coming from `s2_tile`. `preproc.preproc` uses this set to `True`.

## Changed

- The standard profile of preprocessed images now follows a new standard.
- `constants.LEOTIFF` has been renamed to `GTIFF_UINT16`.
- Processed tiles now contain their own metadata tags.
- The names of several leotools commandline function have been shortened. Use `leotools list` to see them all.
- The commandline interpreter now accepts `None`, `True` and `False` arguments and converts them into None-type and Bool-type objects respectively.
- The name of `gistools.stack_bands` funtion to `gistools.stack_images`.
- Some functions now accept inputs in more diverse formats. These are `preproc.reproj_tile`, `gistools.stack_images`, `gistools.multi_mask`.
- `basetools.ProcessTimer` now measures days and centiseconds.

## Removed

- `basetools.path_matches` and `basetools.list_dir` functions.

# v0.3.2

## Added

- `preproc.s2_tile` now has an `incl_gt` arg.

## Fixed

- Two-part S2 tiles are now handled correctly.

# v0.3.1

## Added

- The `preproc.prepcorc` function now has a `temp_dir` arg.

## Changed

- Datatake merging is now done band-by-band. This should deal with memory problems.

# v0.3.0

## Added

- leotools is now usable from the commandline.
- `constants` submodule.
- `preproc.preproc` and `preproc.extract_ndxi` functions.
- The package now contains `ndxi_l8.conf` and `ndxi_s2.conf` ndxi configuration files.
- Optinal args are indicated in function docstrings.
- Functions in the `preproc` submodule now have a `mode` arg that decides their behavior when encountering a non-existant directory.

## Changed

- Functions in the `preproc` submodule are now verbose.
- `preproc.mass_proc` was renamed to `preproc.reproj_tile`. It now runs on individual files as well as directories.
- `preproc.l8_tile` now uses and collects the .xml metadata file instead of .txt.
- Nodata values used in the `preproc` submodule have been set to 0 across the board.
- Preprocessed tiles now contain their own `source` tags. Merged datatakes receive these sources instead of composing them from metadata filenames.

## Removed

- Unnecessary args in `preproc.l8_tile` and `preproc.s2_tile` functions.

## Fixed

- `preproc.merge_datatake` now accepts the same path as both an input and output directory. When starting, the directory should only contain unmerged tiles.

# v0.2.0

## Added

- New satellite image preprocessing tools.

## Changed

- Instead of a single module, leotools is now a package containing three submodules.

# v0.1.0

## Added

- `leo_def` is the first iteration of a unified code base