import os
from distutils.core import setup
from distutils.extension import Extension
import pathlib
import warnings

import numpy
from Cython.Distutils import build_ext

import settings

###########################################
# Get the third-party include directories #
###########################################
include_dirs = [numpy.get_include(), str(settings.CPP_SOURCE_DIRECTORY), settings.EIGEN_DIR, str(settings.CONDA_ENVIRONMENT_INCLUDE)]

############################
# Get the static libraries #
############################
# Find current project static library
project_static_library = settings.BUILD_DIRECTORY / settings.CPP_SOURCE_SUBDIRECTORY / f"lib{settings.PROJECT_NAME}.a"
static_libraries = [str(project_static_library.resolve())]

# Get all of the upstream static libraries
# TODO: make the static library path more robust. Does it need to be more robust or is this logic consistent with Conda
# practices?
if settings.CONDA_ENVIRONMENT_LIB64.exists() and settings.CONDA_ENVIRONMENT_LIB64.is_dir():
    libdir="lib64"
else:
    libdir="lib"

for upstream_project in settings.UPSTREAM_PROJECTS:
    upstream_installed = settings.CONDA_ENVIRONMENT / f"{libdir}/lib{upstream_project}.a"
    upstream_insource = settings.BUILD_DIRECTORY / f"_deps/{upstream_project}-build" / settings.CPP_SOURCE_SUBDIRECTORY / f"lib{upstream_project}.a"
    if upstream_installed.exists() and upstream_installed.is_file():
        static_libraries.append(str(upstream_installed.resolve()))
    elif upstream_insource.exists() and upstream_insource.is_file():
        static_libraries.append(str(upstream_insource.resolve()))
    else:
        warnings.warn(f"Could not find upstream static library from '{upstream_project}'", RuntimeWarning)

################################################
# Get the upstream project include directories #
################################################

# Get all of the possible in-source build include locations
for upstream_project in settings.UPSTREAM_PROJECTS:
    upstream_insource = settings.BUILD_DIRECTORY / f"_deps/{upstream_project}-src" / settings.CPP_SOURCE_SUBDIRECTORY 
    if upstream_insource.exists() and upstream_insource.is_dir():
        include_dirs.append(upstream_insource.resolve())
    upstream_insource = settings.BUILD_DIRECTORY / f"_deps/{upstream_project}-src" / settings.PYTHON_SOURCE_SUBDIRECTORY
    if upstream_insource.exists() and upstream_insource.is_dir():
        include_dirs.append(upstream_insource.resolve())

# Define the build configuration
ext_modules = [Extension(settings.PROJECT_NAME,
                     sources=["main.pyx"],
                     language='c++',
                     extra_objects=static_libraries,
                     include_dirs=include_dirs,
                     extra_compile_args=[f"-std=c++{settings.CXX_STANDARD}"],
                     extra_link_args=[f"-std=c++{settings.CXX_STANDARD}"]
                     )]

setup(
  name = settings.PROJECT_NAME,
  cmdclass = {'build_ext': build_ext},
  ext_modules = ext_modules
)
