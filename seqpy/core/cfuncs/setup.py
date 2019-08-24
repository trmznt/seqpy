from distutils.core import setup
from Cython.Build import cythonize
import numpy

setup(name="fastdx"
	, ext_modules=cythonize('fastdx.pyx')
	, include_dirs=[numpy.get_include()]
	,)
setup(name="genoutils"
	, ext_modules=cythonize('genoutils.pyx')
	, include_dirs=[numpy.get_include()]
	,)

