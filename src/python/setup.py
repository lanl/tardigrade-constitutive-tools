import os
import re
from distutils.core import setup
from distutils.extension import Extension
import pathlib

import numpy
from Cython.Distutils import build_ext

import settings


def return_group_or_error(regex, contents):
    """
    Return a regex group or raise a ValueError

    :param str regex: the Python 3 re module regex
    :param str contents: the string to search for the regex group
    :returns: string of regex group
    """
    search_results = re.search(regex, contents)
    if search_results:
        return search_results.group(0).strip()
    else:
        raise ValueError(f"'{regex}' pattern not found in CMake string contents")

# Search operations on the cmake cache file
# Open the project cmake cache
with open(settings.PROJECT_CMAKE_CACHE, 'r') as cmake_cache_file:
    cmake_cache_contents = cmake_cache_file.read()

# Get the sub-project source directories if the fetch type is local
local_libraries = [settings.CPP_BUILD_DIRECTORY]
library_search_string = "**/*-src*/"

###########################################
# Get the third-party include directories #
###########################################
include_dirs = [numpy.get_include(), str(settings.CPP_SOURCE_DIRECTORY), settings.EIGEN_DIR, str(settings.CONDA_ENVIRONMENT_INCLUDE)]

############################
# Get the static libraries #
############################
# Find current project static library
project_static_library = settings.CPP_BUILD_DIRECTORY / settings.CPP_SOURCE_SUBDIRECTORY / f"lib{settings.PROJECT_NAME}.a"
static_libraries = [str(project_static_library.resolve())]

# Get all of the upstream static libraries
for upstream_project in settings.STATIC_LIBRARY_LINKING_ORDER[1:]:
    upstream_installed = pathlib.Path(settings.CONDA_ENVIRONMENT) / f"lib/lib{upstream_project}.a"
    upstream_insource = settings.CPP_BUILD_DIRECTORY / f"_deps/{upstream_project}-build" / settings.CPP_SOURCE_SUBDIRECTORY / f"lib{upstream_project}.a"
    if upstream_installed.exists() and upstream_installed.is_file():
        static_libraries.append(str(upstream_installed.resolve()))
    elif upstream_insource.exists() and upstream_insource.is_file():
        static_libraries.append(str(upstream_insource.resolve()))
    else:
        raise RuntimeError(f"Could not find upstream static library from '{upstream_project}'")

################################################
# Get the upstream project include directories #
################################################

# Get all of the possible in-source build include locations
for upstream_project in settings.UPSTREAM_PROJECTS:
    upstream_insource = settings.CPP_BUILD_DIRECTORY / f"_deps/{upstream_project}-src" / settings.CPP_SOURCE_SUBDIRECTORY 
    if upstream_insource.exists() and upstream_insource.is_dir():
        include_dirs.append(upstream_insource.resolve())
    upstream_insource = settings.CPP_BUILD_DIRECTORY / f"_deps/{upstream_project}-src" / settings.PYTHON_SOURCE_SUBDIRECTORY
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
