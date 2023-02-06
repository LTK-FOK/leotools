from setuptools import setup #, ssl_support
from leotools import __version__, __author__, __email__, __description__, __license__

# ssl_support.cert_paths = [os.environ.get('SSL_CERT_FILE')]

name = "leotools"

setup(
    name=name,
    description=__description__,
    version=__version__,
    author=__author__,
    author_email=__email__,
    keywords=["python", "GIS"],

    packages=[name],
    install_requires=[
        "rasterio",
        "geopandas",
    ],
    entry_points={
        "console_scripts": [f"{name} = {name}.__main__:main"]
    }
)