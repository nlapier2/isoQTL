from setuptools import setup
#from distutils.core import setup, Extension
from Cython.Build import cythonize
import numpy

setup(
    name='IsoQTL',
    #ext_modules=cythonize("hashimoto_ohtani_cython.pyx"),
    ext_modules=cythonize("isoqtl_cis_pass.pyx"), 
    include_dirs=[numpy.get_include()],
    zip_safe=False,
)
