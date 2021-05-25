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

# Get the project fetch source type
project_fetch_source_regex = '(?<=FETCH_SOURCE:STRING=).*'
fetch_source = return_group_or_error(project_fetch_source_regex, cmake_cache_contents)

# Get the sub-project source directories if the fetch type is local
if fetch_source == "REPO":
    local_libraries = [settings.CPP_BUILD_DIRECTORY]
    library_search_string = "**/*-src*/"

elif fetch_source == "LOCAL":
    local_libraries = []
    library_search_string = "**/"
    for source_variable_name in settings.LIBRARY_SOURCE_VARIABLE_NAMES:
        regex = f'(?<={source_variable_name}:PATH=).*'
        local_libraries.append(return_group_or_error(regex, cmake_cache_contents))
else:
    raise ValueError(f"FETCH_SOURCE {fetch_source} not recognized")

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

static_libraries = []

# Get all of the static libraries
static_libraries = [str(lib.resolve()) for lib in pathlib.Path(settings.CPP_BUILD_DIRECTORY).glob("**/*.a")]

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

# Re-order the static libraries so they are in the correct include order
if len(settings.STATIC_LIBRARY_LINKING_ORDER) != len(static_libraries):
    raise ValueError("The expected static libraries and the detected static libraries have different sizes")

ordered_static_libraries = [None for _ in static_libraries]
for index, library in enumerate(settings.STATIC_LIBRARY_LINKING_ORDER):

    ordered_static_libraries[index] = [library_path for library_path in static_libraries if ('lib' + library) in library_path][0]

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
