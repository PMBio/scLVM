
"""scLVM - single cell latent variable model
See:
https://github.com/PMBio/scLVM
"""

# Always prefer setuptools over distutils
from setuptools import setup, find_packages
# To use a consistent encoding
from codecs import open
from os import path

here = path.abspath(path.dirname(__file__))

# Get the long description from the relevant file
with open(path.join(here, 'README.md'), encoding='utf-8') as f:
    long_description = f.read()
    
reqs = ['numpy', 'scipy', 'pygp >=1.1.07', 'matplotlib >=1.2','h5py','limix >=0.6.6']

setup(
    name='scLVM',

    # Versions should comply with PEP440.  For a discussion on single-sourcing
    # the version across setup.py and the project code, see
    # https://packaging.python.org/en/latest/single_source_version.html
    version='0.1.5',

    description='scLVM',
    long_description=long_description,

    # The project's main homepage.
    url='https://github.com/PMBio/scLVM',

    # Author details
    author='Florian Buettner, F Paolo Casale, Oliver Stegle',
    author_email='scLVM-dev@ebi.ac.uk',

    # Choose your license
    license='Apache 2.0',


    # What does your project relate to?
    keywords=['single-cell RNA-seq, latent variable model',' Variance compoennt modelling'],

    # You can just specify the packages manually here if your project is
    # simple. Or you can use find_packages().
    packages=['scLVM', 'scLVM.utils'],#find_packages(),

    # List run-time dependencies here.  These will be installed by pip when
    # your project is installed. For an analysis of "install_requires" vs pip's
    # requirements files see:
    # https://packaging.python.org/en/latest/requirements.html
    
    install_requires=reqs,


    py_modules = ['scLVM.core','scLVM.gp_clvm']

)
