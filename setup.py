import os, sys

from distutils.core import setup, Extension
from distutils import sysconfig

cpp_args = ['-std=c++11', '-stdlib=libc++']

sfc_module = Extension(
    'crp', sources=['source.cpp', "Geophone.cpp", "box.cpp", "ImageP.cpp","IpToGeoGeometry.cpp"],
    include_dirs=['c:/pybind11/include','eigen'],
    language='c++',
    extra_compile_args=cpp_args,
    )

setup(
    name='crp',
    version='1.0',
    description='Python package with superfastcode C++ extension (PyBind11)',
    ext_modules=[sfc_module],
)