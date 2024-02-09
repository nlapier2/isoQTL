from setuptools import setup
#from distutils.core import setup, Extension
from Cython.Build import cythonize
import numpy

setup(
    name='IsoQTL',
    ext_modules=cythonize("cis_pass_opposites.pyx"), 
    include_dirs=[numpy.get_include()],
    zip_safe=False,
)
