cmake_minimum_required(VERSION 3.14)
project(constitutive_tools VERSION 0.1.0)
enable_testing()

# Make the code position independent
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

# Configure version
configure_file("src/cpp/config.h.in" "${PROJECT_BINARY_DIR}/config.h")

# Find eigen
find_package(Eigen3 3.3 REQUIRED NO_MODULE)
include_directories(${EIGEN3_INCLUDE_DIR})
if(EIGEN3_FOUND)
    message(STATUS "Found Eigen3: ${EIGEN3_INCLUDE_DIR}")
endif()

# Save the eigen directory for use in the python interface
set(EIGEN_DIR ${EIGEN3_INCLUDE_DIR} CACHE PATH "The path to the eigen include directory")

# Add the cmake folder to locate the FindSphinx module
set(CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake" ${CMAKE_MODULE_PATH})

# Set build type checks 
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE "None")
endif()
string(TOLOWER "${CMAKE_BUILD_TYPE}" cmake_build_type_lower)
set(upstream_required "")
if(cmake_build_type_lower STREQUAL "release")
    set(upstream_required "REQUIRED")
endif()
# Find related, but external, projects in installed environment
set(upstream_packages "error_tools" "vector_tools")
foreach(package ${upstream_packages})
    find_package(${package} ${upstream_required} CONFIG)
    if(${package}_FOUND)
        message(STATUS "Found ${package}: ${${package}_DIR}")
    else()
        message(WARNING "Did not find ${package} package")
    endif()
endforeach(package)

# Find related, but external, projects using FetchContent and building locally
if(NOT error_tools_FOUND AND NOT vector_tools_FOUND)
    message(STATUS "Upstream package(s) not found. Attempting local dependencies build(s).")
    # Set a variable indicating if the other projects should come from the git repository or a local source
    set(FETCH_SOURCE "REPO" CACHE STRING "The source of the external projects. Options are REPO and LOCAL if LOCAL then the paths must be defined")
    set(ERROR_TOOLS_PATH "" CACHE PATH "The path to the local version of error_tools")
    set(VECTOR_TOOLS_PATH "" CACHE PATH "The path to the local version of vector_tools")
    # TODO: Decide how to handle master/dev build differences. bash env variables for GIT TAG?
    include(FetchContent)
    if(FETCH_SOURCE STREQUAL "REPO")
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
    elseif(FETCH_SOURCE STREQUAL "LOCAL")
        FetchContent_Declare(
            error_tools
            SOURCE_DIR ${ERROR_TOOLS_PATH}
        )
        FetchContent_Declare(
            vector_tools
            SOURCE_DIR ${VECTOR_TOOLS_PATH}
        )
    else()
        message(FATAL_ERROR "Fetch source ${FETCH_SOURCE} is not recognized. Options are REPO and LOCAL")
    endif()
    FetchContent_MakeAvailable(error_tools vector_tools)
endif()

# Add project source directories
add_subdirectory("src/cpp")
if(CMAKE_PROJECT_NAME STREQUAL PROJECT_NAME)
    # Find Boost. Required for tests
    find_package(Boost 1.53.0 REQUIRED COMPONENTS unit_test_framework)
    # Add tests and docs
    add_subdirectory("src/cpp/tests")
    add_subdirectory("docs")
endif()

# Set the c++ standard
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED True)
