from distutils.core import setup
from Cython.Build import cythonize

setup(name="fastdx", ext_modules=cythonize('fastdx.pyx'),)
setup(name="genoutils", ext_modules=cythonize('genoutils.pyx'),)

