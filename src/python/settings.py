import os
import re
import pathlib


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


CONDA_ENVIRONMENT = pathlib.Path(os.environ['CONDA_DEFAULT_ENV'])
CONDA_ENVIRONMENT_INCLUDE = CONDA_ENVIRONMENT / "include"
ROOT_DIRECTORY = pathlib.Path(".").resolve().parent.parent
CPP_BUILD_DIRECTORY = ROOT_DIRECTORY / "build"
PROJECT_CMAKE_CACHE = CPP_BUILD_DIRECTORY / "CMakeCache.txt"
PROJECT_CMAKE_FILE = ROOT_DIRECTORY / "CMakeLists.txt"
PYTHON_SOURCE_SUBDIRECTORY = pathlib.Path("src/python")
CPP_SOURCE_SUBDIRECTORY = pathlib.Path("src/cpp")
CPP_SOURCE_DIRECTORY = ROOT_DIRECTORY / CPP_SOURCE_SUBDIRECTORY
PYTHON_SOURCE_DIRECTORY = ROOT_DIRECTORY / PYTHON_SOURCE_SUBDIRECTORY
UPSTREAM_PROJECTS = ['vector_tools', 'error_tools']
LIBRARY_SOURCE_VARIABLE_NAMES = [f"{project.upper()}_PATH" for project in UPSTREAM_PROJECTS]
STATIC_LIBRARY_LINKING_ORDER = ['constitutive_tools', 'error_tools']

# Open the project root CMake configuration file
with open(PROJECT_CMAKE_FILE, 'r') as cmake_lists_file:
    cmake_lists_contents = cmake_lists_file.read()

# Open the project cmake cache
with open(PROJECT_CMAKE_CACHE, 'r') as cmake_cache_file:
    cmake_cache_contents = cmake_cache_file.read()

# Get the project name
project_name_regex = '(?<=CMAKE_PROJECT_NAME:STATIC=).*'
PROJECT_NAME = return_group_or_error(project_name_regex, cmake_cache_contents)

# Get the CXX standard
project_cxx_std_regex = '(?<=CMAKE_CXX_STANDARD)(.*)(?=\))'
CXX_STANDARD = return_group_or_error(project_cxx_std_regex, cmake_lists_contents)

# Get the Eigen header files directory
eigen_regex = '(?<=EIGEN_DIR:PATH=).*'
EIGEN_DIR = return_group_or_error(eigen_regex, cmake_cache_contents)
