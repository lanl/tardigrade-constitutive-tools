import os
import pathlib

CONDA_ENVIRONMENT = pathlib.Path(os.environ['CONDA_DEFAULT_ENV'])
CONDA_ENVIRONMENT_INCLUDE = CONDA_ENVIRONMENT / "include"
ROOT_DIRECTORY = os.path.abspath(os.path.join("..", ".."))
PROJECT_CMAKE_FILE = os.path.join(ROOT_DIRECTORY, "CMakeLists.txt")
PYTHON_SOURCE_SUBDIRECTORY = os.path.join("src", "python")
CPP_SOURCE_SUBDIRECTORY = os.path.join("src", "cpp")
CPP_SOURCE_DIRECTORY = os.path.join(ROOT_DIRECTORY, CPP_SOURCE_SUBDIRECTORY)
PYTHON_SOURCE_DIRECTORY = os.path.join(ROOT_DIRECTORY, PYTHON_SOURCE_SUBDIRECTORY)
CPP_BUILD_DIRECTORY = os.path.join(ROOT_DIRECTORY, "build")
PROJECT_CMAKE_CACHE = os.path.join(CPP_BUILD_DIRECTORY, "CMakeCache.txt")
UPSTREAM_PROJECTS = ['vector_tools', 'error_tools']
LIBRARY_SOURCE_VARIABLE_NAMES = [f"{project.upper()}_PATH" for project in UPSTREAM_PROJECTS]
STATIC_LIBRARY_LINKING_ORDER = ['constitutive_tools', 'error_tools']
