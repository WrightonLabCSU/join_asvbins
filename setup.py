"""Setup file for package"""
from setuptools import setup, find_packages
# from some_python_init_file import __version__ as version
from os import path

__author__ = 'rmflynn'
__version__ = '0.0.1'

here = path.abspath(path.dirname(__file__))

with open(path.join(here, 'README.md'), encoding='utf-8') as f:
        long_description = f.read()

setup(
    name="asvbin_joiner_rmflynn",
    version=__version__,
    scripts=['scripts/join_asvbins.py'],
    package_dir={"": "src/"},
    packages=find_packages(where='src/'),
    description="Asv to bin joining tool",
    long_description=long_description,
    long_description_content_type='text/markdown',  # Optional (see note above)
    python_requires='>=3',
    install_requires=['scikit-bio', 'pandas', 'numpy', 'snakemake', 'graphviz'],
    author="Rory Flynn",
    author_email='Rory.Flynn@colostate.edu',
    url="https://github.com/rmFlynn/16s_to_bins_project",  # this will change
    download_url="None for %s" % __version__)

