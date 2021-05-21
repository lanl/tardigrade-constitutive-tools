cmake_minimum_required(VERSION 3.14)
project(constitutive_tools VERSION 0.1.0)
enable_testing()

# Make the code position independent
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

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
string(TOLOWER "${CMAKE_BUILD_TYPE}" cmake_build_type_lower)
set(upstream_required "")
if(cmake_build_type_lower STREQUAL "release")
    set(upstream_required "REQUIRED")
endif()

# Set a variable indicating if upstream projects should come from the git repository or a local source
include(FetchContent)
set(FETCH_SOURCE "REPO" CACHE STRING "The source of the external projects. Options are REPO and LOCAL if LOCAL then the paths must be defined")
set(ERROR_TOOLS_PATH "" CACHE PATH "The path to the local version of error_tools")
set(VECTOR_TOOLS_PATH "" CACHE PATH "The path to the local version of vector_tools")

# Find related, but external, projects in installed environment
set(upstream_packages "error_tools" "vector_tools")
foreach(package ${upstream_packages})
    # FIXME: Always search for installed upstream packages after fixing Cython build include/library paths
    # Returning to conditional search for upstream deletes the following lines
    #   if(upstream_required STREQUAL "REQUIRED")'
    #   endif()
    # And removes the indents from lines
    #   find_package(${package} ${upstream_required} CONFIG)
    #   if(${package}_FOUND)
    #       message(STATUS "Found ${package}: ${${package}_DIR}")
    if(upstream_required STREQUAL "REQUIRED")
        find_package(${package} ${upstream_required} CONFIG)
        if(${package}_FOUND)
            message(STATUS "Found ${package}: ${${package}_DIR}")
        endif()
    else()
        # Find related, but external, projects using FetchContent and building locally
        message(WARNING "Did not find an installed ${package} package. Attempting local build with FetchContent.")
        # TODO: Decide how to handle master/dev build differences. bash env variables for GIT TAG?
        if(FETCH_SOURCE STREQUAL "REPO")
            FetchContent_Declare(
                ${package}
                GIT_REPOSITORY ssh://git@xcp-stash.lanl.gov:7999/mm/${package}.git
                GIT_TAG        origin/dev
            )
        elseif(FETCH_SOURCE STREQUAL "LOCAL")
            FetchContent_Declare(
                ${package}
                SOURCE_DIR ${${error_tools_PATH}
            )
        else()
            message(FATAL_ERROR "Fetch source ${FETCH_SOURCE} is not recognized. Options are REPO and LOCAL")
        endif()
        FetchContent_MakeAvailable(${package})
    endif()
endforeach(package)

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
