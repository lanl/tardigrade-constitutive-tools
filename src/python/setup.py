import os
import re
from distutils.core import setup
from distutils.extension import Extension
import pathlib

import numpy
from Cython.Distutils import build_ext

import settings


# Open the project cmake cache
with open(settings.PROJECT_CMAKE_CACHE, 'r') as cmake_file:
    cmake_contents = cmake_file.read()

# Get the project name
project_name_regex = '(?<=CMAKE_PROJECT_NAME:STATIC=).*'
project_name_search = re.search(project_name_regex, cmake_contents)

if project_name_search.group:
    project_name = project_name_search.group(0)
else:
    raise ValueError("CMAKE_PROJECT_NAME was not found in CMakeCache.txt")

# Get the project fetch source type
project_fetch_source_regex = '(?<=FETCH_SOURCE:STRING=).*'
project_fetch_source_search = re.search(project_fetch_source_regex, cmake_contents)

if project_fetch_source_search.group:
    fetch_source = project_fetch_source_search.group(0)
else:
    raise ValueError("FETCH_SOURCE was not found in CMakeCache.txt")

# Get the sub-project source directories if the fetch type is local
if fetch_source == "REPO":
    local_libraries = [settings.CPP_BUILD_DIRECTORY]
    library_search_string = "**/*-src*/"
    
elif fetch_source == "LOCAL":
    local_libraries = []
    library_search_string = "**/"
    for source_variable_name in settings.LIBRARY_SOURCE_VARIABLE_NAMES:
        regex = f'(?<={source_variable_name}:PATH=).*'
        search = re.search(regex, cmake_contents)

        if search.group:
            local_libraries.append(search.group(0))
        else:
            raise ValueError(f"{source_variable_name} was not found in CMakeCache.txt")
else:
    raise ValueError(f"CMAKE_FETCH_SOURCE {fetch_source} not recognized")

###############################
# Get the include directories #
###############################

include_dirs = [numpy.get_include(), settings.CPP_SOURCE_DIRECTORY]

# Get the Eigen library
eigen_regex = '(?<=CMAKE_EIGEN_DIR:PATH=).*'
eigen_search = re.search(eigen_regex, cmake_contents)

if eigen_search.group:
    include_dirs.append(eigen_search.group(0))
else:
    raise ValueError("CMAKE_EIGEN_DIR was not found in CMakeCache.txt")

############################
# Get the static libraries #
############################

static_libraries = []

# Get all of the static libraries
static_libraries = [str(lib.resolve()) for lib in pathlib.Path(settings.CPP_BUILD_DIRECTORY).glob("**/*.a")]

# Ignore all of the libraries except for the one associated with this project
static_libaries = [sl for sl in static_libraries if ("lib" + project_name) in sl]

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
ext_modules = [Extension("constitutive_tools",
                     sources=["main.pyx"],
                     language='c++',
                     extra_objects=static_libraries,
                     include_dirs=include_dirs,
                     extra_compile_args=[f"-std=c++{settings.CMAKE_CXX_STANDARD}"],
                     extra_link_args=[f"-std=c++{settings.CMAKE_CXX_STANDARD}"]
                     )]

setup(
  name = 'constitutive_tools',
  cmdclass = {'build_ext': build_ext},
  ext_modules = ext_modules
)
