cmake_minimum_required (VERSION 3.14)
project (constitutive_tools VERSION 0.0.1)
enable_testing()

# Configure version
configure_file("src/cpp/config.h.in" "${PROJECT_BINARY_DIR}/config.h")

# Find eigen
set(CMAKE_PREFIX_PATH ${CMAKE_PREFIX_PATH} $ENV{HOME}/.local $ENV{HOME}/local/ $ENV{HOME})
find_package (Eigen3 3.3 REQUIRED NO_MODULE)
include_directories(${EIGEN3_INCLUDE_DIR})

# Add the cmake folder to locate the FindSphinx module
set(CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake" ${CMAKE_MODULE_PATH})

# Find related, but external, projects
# TODO: Decide how to handle master/dev build differences. bash env variables for GIT TAG?
include(FetchContent)
FetchContent_Declare(
    error_tools
    GIT_REPOSITORY ssh://git@xcp-stash.lanl.gov:7999/mm/error_tools.git
    GIT_TAG        origin/dev
)
FetchContent_Declare(
    vector_tools
    GIT_REPOSITORY ssh://git@xcp-stash.lanl.gov:7999/mm/vector_tools.git
    GIT_TAG        origin/dev
)
FetchContent_MakeAvailable(error_tools vector_tools)

# Add project source directories
add_subdirectory ("src/cpp")
if(CMAKE_PROJECT_NAME STREQUAL PROJECT_NAME)
    # Find Boost. Required for tests and the python interface
    find_package(Boost 1.53.0 REQUIRED COMPONENTS unit_test_framework python)
    # Add the current conda environment to the CMAKE_PREFIX_PATH
    set(Python_ROOT_DIR $ENV{$CONDA_PREFIX})
    set(Python_FIND_STRATEGY "LOCATION")
    find_package(Python REQUIRED COMPONENTS Interpreter Development)
    message("PythonINCLUDE_DIRS: ${Python_INCLUDE_DIRS}")
    message("Python_ROOT_DIR: ${Python_ROOT_DIR}")
    # Add tests and docs
    add_subdirectory ("src/cpp/tests")
    add_subdirectory ("docs")
    # Link the libraries for the python interface
    target_link_libraries(constitutive_tools ${Boost_LIBRARIES} ${PYTHON_LIBRARIES})
    target_include_directories(constitutive_tools PRIVATE ${Python_INCLUDE_DIRS})
    message("PYTHON LIBRARIES: ${Python_LIBRARIES}")
    message("PYTHON INCLUDE: ${Python_INCLUDE_DIRS}")
endif()

# Set the c++ standard
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED True)
