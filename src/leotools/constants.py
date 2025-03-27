"""This module contains often used constants."""

from rasterio.profiles import Profile

### Constants
EOV = 'EPSG:23700'

class GTIFF_UINT16(Profile):
    """Standard GTiff profile for int images."""

    defaults = {
        'driver': 'GTiff',
        'interleave': 'band',
        'tiled': True,
        'blockxsize': 256,
        'blockysize': 256,
        'compress': 'zstd',
        'zstd_level': 1,
        'predictor': 2,
        'dtype': 'uint16',
        'nodata': 0,
    }

class GTIFF_UINT16_COMPATIBLE(Profile):
    """Standard GTiff profile for int images."""

    defaults = {
        'driver': 'GTiff',
        'interleave': 'band',
        'tiled': True,
        'blockxsize': 256,
        'blockysize': 256,
        'compress': 'deflate',
        'zlevel': 1,
        'predictor': 2,
        'dtype': 'uint16',
        'nodata': 0,
    }

NDXI_L8 = {
    'ndvi': r"ndiff(b[5], b[4])",
    'ndwi': r"ndiff(b[3], b[5])",
    'ndmi': r"ndiff(b[5], b[6])",
    'ndwig': r"ndiff(b[6], b[4])",
    'bsi11': r"ndiff(b[6]+b[4], b[5]+b[2])",
    'bsi12': r"ndiff(b[7]+b[4], b[5]+b[2])",
    'ndsig': r"ndiff(b[6], b[5])",
    'nmdi': r"ndiff(b[5], b[6]-b[7])",
    # 'grvi': r"b[8] / b[3]",
    'tct1': r"b[1]*0.3029 + b[2]*0.2786 + b[3]*0.4733 + b[4]*0.5599 + b[5]*0.5080 + b[6]*0.1872",
    'tct2': r"b[1]*-0.2941 + b[2]*-0.2430 + b[3]*-0.5424 + b[4]*0.7276 + b[5]*0.0713 + b[6]*-0.1608",
    'tct3': r"b[1]*0.1511 + b[2]*0.1973 + b[3]*0.3283 + b[4]*0.3407 + b[5]*-0.7117 + b[6]*-0.4559",
    'tct4': r"b[1]*-0.8239 + b[2]*0.0849 + b[3]*0.4396 + b[4]*-0.0580 + b[5]*0.2013 + b[6]*-0.2773",
    'tct5': r"b[1]*-0.3294 + b[2]*0.0557 + b[3]*0.1056 + b[4]*0.1855 + b[5]*-0.4349 + b[6]*0.8085",
    'tct6': r"b[1]*0.1079 + b[2]*-0.9023 + b[3]*0.4119 + b[4]*0.0575 + b[5]*-0.0259 + b[6]*0.0252",
}

NDXI_S2 = {
    'ndvi': r"ndiff(b[8], b[4])",
    'ndwi': r"ndiff(b[3], b[8])",
    'ndmi': r"ndiff(b[8], b[10])",
    'ndmid': r"ndiff(b[9], b[10])",
    'ndwig': r"ndiff(b[10], b[4])",
    'bsi11': r"ndiff(b[10]+b[4], b[8]+b[2])",
    'bsi12': r"ndiff(b[11]+b[4], b[8]+b[2])",
    'rlai': r"ndiff(b[5], b[4])",
    'psri': r"(b[4]-b[2]) / b[6]",
    'ndsig': r"ndiff(b[10], b[8])",
    'nmdi': r"ndiff(b[8], b[10]-b[11])",
    'rendvi': r"ndiff(b[6], b[5])",
    'grvi': r"b[5] / b[3]",
    'nvia': r"ndiff(b[7], b[10])",
}