import os


CMAKE_CXX_STANDARD = 11
ROOT_DIRECTORY = os.path.abspath(os.path.join("..", ".."))
PROJECT_CMAKE_FILE = os.path.join(ROOT_DIRECTORY, "CMakeLists.txt")
PYTHON_SOURCE_SUBDIRECTORY = os.path.join("src", "python")
CPP_SOURCE_SUBDIRECTORY = os.path.join("src", "cpp")
CPP_SOURCE_DIRECTORY = os.path.join(ROOT_DIRECTORY, CPP_SOURCE_SUBDIRECTORY)
PYTHON_SOURCE_DIRECTORY = os.path.join(ROOT_DIRECTORY, PYTHON_SOURCE_SUBDIRECTORY)
CPP_BUILD_DIRECTORY = os.path.join(ROOT_DIRECTORY, "build")
PROJECT_CMAKE_CACHE = os.path.join(CPP_BUILD_DIRECTORY, "CMakeCache.txt")
PYTHON_LIBS = ['error_tools']
