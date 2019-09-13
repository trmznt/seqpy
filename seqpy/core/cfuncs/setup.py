from distutils.core import setup
from Cython.Build import cythonize
import numpy

setup(name="fastdx"
	, ext_modules=cythonize('fastdx.pyx', language_level=3)
	, include_dirs=[numpy.get_include()]
	,)
setup(name="genoutils"
	, ext_modules=cythonize('genoutils.pyx', language_level=3)
	, include_dirs=[numpy.get_include()]
	,)

