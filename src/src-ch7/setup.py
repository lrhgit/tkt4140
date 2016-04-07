from distutils.core import setup
from Cython.Build import cythonize
import numpy

setup(
    ext_modules = cythonize("gauss_seidel_sor.pyx"),
    include_dirs=[numpy.get_include()]
)


