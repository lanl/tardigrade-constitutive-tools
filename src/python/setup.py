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

# Search operations on the cmake lists file
# Open the project root CMake configuration file
with open(settings.PROJECT_CMAKE_FILE, 'r') as cmake_lists_file:
    cmake_lists_contents = cmake_lists_file.read()

# Get the CXX standard
project_cxx_std_regex = '(?<=CMAKE_CXX_STANDARD)(.*)(?=\))'
cxx_standard = return_group_or_error(project_cxx_std_regex, cmake_lists_contents)

# Search operations on the cmake cache file
# Open the project cmake cache
with open(settings.PROJECT_CMAKE_CACHE, 'r') as cmake_cache_file:
    cmake_cache_contents = cmake_cache_file.read()

# Get the project name
project_name_regex = '(?<=CMAKE_PROJECT_NAME:STATIC=).*'
project_name = return_group_or_error(project_name_regex, cmake_cache_contents)

# Get the sub-project source directories if the fetch type is local
local_libraries = [settings.CPP_BUILD_DIRECTORY]
library_search_string = "**/*-src*/"

###############################
# Get the include directories #
###############################
# FIXME: VIP-648 - use the installed upstream packages for the "include_dirs" whenever possible

include_dirs = [numpy.get_include(), settings.CPP_SOURCE_DIRECTORY]

# Get the Eigen library
eigen_regex = '(?<=EIGEN_DIR:PATH=).*'
include_dirs.append(return_group_or_error(eigen_regex, cmake_cache_contents))

############################
# Get the static libraries #
############################
# FIXME: VIP-648 - use the installed upstream packages for the "ordered_static_libraries" whenever possible
# Find current project static library
project_static_library = pathlib.Path(settings.CPP_BUILD_DIRECTORY) / settings.CPP_SOURCE_SUBDIRECTORY / f"lib{project_name}.a"
static_libraries = [str(project_static_library.resolve())]

# Get all of the upstream static libraries
for upstream_project in settings.STATIC_LIBRARY_LINKING_ORDER[1:]:
    upstream_installed = pathlib.Path(settings.CONDA_ENVIRONMENT) / f"lib/lib{upstream_project}.a"
    upstream_insource = pathlib.Path(settings.CPP_BUILD_DIRECTORY) / f"_deps/{upstream_project}-build" / settings.CPP_SOURCE_DIRECTORY / f"lib{upstream_project}.a"
    if upstream_installed.exists() and upstream_installed.is_file():
        static_libraries.append(str(upstream_installed.resolve()))
    elif upstream_insource.exists() and upstream_insource.is_file():
        static_libraries.append(str(upstream_in_source.resolve()))
    else:
        raise RuntimeError(f"Could not find upstream static library from '{upstream_project}'")

import pdb
pdb.set_trace()

###################################
# Get all of the pyx source files #
###################################

# Get all of the include locations
for source_subpath in (settings.PYTHON_SOURCE_SUBDIRECTORY, settings.CPP_SOURCE_SUBDIRECTORY):

    for local_library in local_libraries:

        for dir in pathlib.Path(local_library).glob(library_search_string + source_subpath):
            if not dir.is_dir():
                continue

            include_dirs.append(str(dir))

# Define the build configuration
ext_modules = [Extension(project_name,
                     sources=["main.pyx"],
                     language='c++',
                     extra_objects=ordered_static_libraries,
                     include_dirs=include_dirs,
                     extra_compile_args=[f"-std=c++{cxx_standard}"],
                     extra_link_args=[f"-std=c++{cxx_standard}"]
                     )]

setup(
  name = project_name,
  cmdclass = {'build_ext': build_ext},
  ext_modules = ext_modules
)
