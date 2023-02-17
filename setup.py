"""linz-stationcoordmodel setup module
"""

# Comments from samply python packaging script

# Always prefer setuptools over distutils
from setuptools import setup, find_packages

# To use a consistent encoding
from codecs import open
from os import path

here = path.abspath(path.dirname(__file__))

# Get the long description from the relevant file
with open(path.join(here, "DESCRIPTION.rst"), encoding="utf-8") as f:
    long_description = f.read()

setup(
    name="linz-stationcoordmodel",
    # Versions should comply with PEP440.  For a discussion on single-sourcing
    # the version across setup.py and the project code, see
    # https://packaging.python.org/en/latest/single_source_version.html
    version="1.1.1",
    description="LINZ.station_coord_model module - model time series of station coordinates",
    long_description=long_description,
    # The project's main homepage.
    url="https://github.com/linz/python-linz-stationcoordmodel",
    # Author details
    author="Chris Crook",
    author_email="ccrook@linz.govt.nz",
    # Choose your license
    license="GPL",
    # See https://pypi.python.org/pypi?%3Aaction=list_classifiers
    classifiers=[
        # How mature is this project? Common values are
        #   3 - Alpha
        #   4 - Beta
        #   5 - Production/Stable
        "Development Status :: 5 - Stable",
        # Indicate who your project is intended for
        "Intended Audience :: LINZ coders",
        "Topic :: Geodetic",
        # Pick your license as you wish (should match "license" above)
        "License :: OSI Approved :: GPL License",
        # Specify the Python versions you support here. In particular, ensure
        # that you indicate whether you support Python 2, Python 3 or both.
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.5",
    ],
    # What does your project relate to?
    keywords="geodesy",
    # You can just specify the packages manually here if your project is
    # simple. Or you can use find_packages().
    packages=find_packages(),
    # List run-time dependencies here.  These will be installed by pip when
    # your project is installed. For an analysis of "install_requires" vs pip's
    # requirements files see:
    # https://packaging.python.org/en/latest/requirements.html
    # Note: matplotlib dependency not included as causes probrams not to run,
    # forces dependency on tornado (python web server)
    # install_requires=['linz-geodetic','numpy','scipy','matplotlib','pandas'],
    install_requires=["linz-geodetic", "numpy", "scipy", "pandas", "requests"],
    # List additional groups of dependencies here (e.g. development
    # dependencies). You can install these using the following syntax,
    # for example:
    # $ pip install -e .[dev,test]
    extras_require={},
    # If there are data files included in your packages that need to be
    # installed, specify them here.  If using Python 2.6 or less, then these
    # have to be included in MANIFEST.in as well.
    package_data={
        "LINZ": [
            "CORS_Analyst.yaml",
            "report_datum_status.cfg",
        ]
    },
    # Namespace package - other modules may include into these packages
    namespace_packages=["LINZ"],
    # To provide executable scripts, use entry points in preference to the
    # "scripts" keyword. Entry points provide cross-platform support and allow
    # pip to create the appropriate form of executable for the target platform.
    entry_points={
        "console_scripts": [
            "calc_refstation_coords=LINZ.calc_refstation_coords:main",
            "station_coord_model=LINZ.StationCoordinateModel:main",
            "spm_update=LINZ.spm_update:main",
            "spm_clean=LINZ.StationCoordinateModel:clean",
            "analyse_cors=LINZ.CORS_Analyst:main",
            "report_datum_status=LINZ.report_datum_status:main",
            "cors_station_coordinates=LINZ.cors_station_coordinates:main",
            "positionzpp_spm_upload=LINZ.positionzpp_spm_upload:main",
        ]
    },
)
