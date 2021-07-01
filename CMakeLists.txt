#===================================================================================================== PROJECT SETUP ===
cmake_minimum_required(VERSION 3.14)
project(constitutive_tools VERSION 0.2.0)

# Set common project paths relative to project root directory
set(CPP_SRC_PATH "src/cpp")
set(PYTHON_SRC_PATH "src/python")
set(CMAKE_SRC_PATH "src/cmake")

# Add the cmake folder to locate project CMake module(s)
set(CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/${CMAKE_SRC_PATH}" ${CMAKE_MODULE_PATH})

# Get version number from Git
set(VERSION_UPDATE_FROM_GIT True)
# FIXME: Figure out why include() doesn't pick up CMAKE_MODULE_PATH correctly
include(${CMAKE_SRC_PATH}/GetVersionFromGitTag.cmake)
project(${PROJECT_NAME} VERSION ${${PROJECT_NAME}_VERSION})

# Add installation directory variables
include(GNUInstallDirs)

# Make the code position independent
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

# Set the c++ standard
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED True)

# Enable CTest
enable_testing()

#================================================================================================= FIND DEPENDENCIES ===
# Find eigen
find_package(Eigen3 3.3 REQUIRED NO_MODULE)
include_directories(${EIGEN3_INCLUDE_DIR})
if(EIGEN3_FOUND)
    message(STATUS "Found Eigen3: ${EIGEN3_INCLUDE_DIR}")
endif()

# Save the eigen directory for use in the python interface
set(EIGEN_DIR ${EIGEN3_INCLUDE_DIR} CACHE PATH "The path to the eigen include directory")

# Add the cmake folder to locate the FindSphinx module
set(CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/${CMAKE_SRC_PATH}" ${CMAKE_MODULE_PATH})

# Set build type checks
string(TOLOWER "${CMAKE_BUILD_TYPE}" cmake_build_type_lower)
set(upstream_required "")
if(cmake_build_type_lower STREQUAL "release")
    set(upstream_required "REQUIRED")
endif()

# Set a variable indicating if upstream projects should come from the git repository or a local source
set(ERROR_TOOLS_PATH "" CACHE PATH "The path to the local version of error_tools")
set(ERROR_TOOLS_GITTAG "" CACHE PATH "The version to build of error_tools")
set(VECTOR_TOOLS_PATH "" CACHE PATH "The path to the local version of vector_tools")
set(VECTOR_TOOLS_GITTAG "" CACHE PATH "The version to build of of vector_tools")

# Find related, but external, projects in installed environment
include(FetchContent)
set(upstream_packages "error_tools" "vector_tools")
foreach(package ${upstream_packages})
    string(TOUPPER "${package}" package_upper)
    find_package(${package} ${upstream_required} CONFIG)
    if(${package}_FOUND)
        message(STATUS "Found ${package}: ${${package}_DIR}")
    else()
        # Find related, but external, projects using FetchContent and building locally
        message(WARNING "Did not find an installed ${package} package. Attempting local build with FetchContent.")
        # Pick repository path to build
        if(${package_upper}_PATH)
            set(PACKAGE_REPOSITORY ${${package_upper}_PATH})
        else()
            set(PACKAGE_REPOSITORY ssh://git@xcp-stash.lanl.gov:7999/mm/${package}.git)
        endif()
        message(STATUS "Building ${package} from ${PACKAGE_REPOSITORY}")
        # Pick repository Git tag to build
        if(${package_upper}_GITTAG)
            set(PACKAGE_GITTAG ${${package_upper}_GITTAG})
        else()
            set(PACKAGE_GITTAG origin/dev)
        endif()
        message(STATUS "Building ${package} version ${PACKAGE_GITTAG}")
        FetchContent_Declare(
            ${package}
            GIT_REPOSITORY ${PACKAGE_REPOSITORY}
            GIT_TAG        ${PACKAGE_GITTAG}
        )
        FetchContent_MakeAvailable(${package})
    endif()
endforeach(package)

#=============================================================================================== ADD PROJECT TARGETS ===
# MUST COME AFTER DEPENDENCY LOCATING
# Add project source directories
include_directories("${CPP_SRC_PATH}")
add_subdirectory(${CPP_SRC_PATH})
add_subdirectory(${PYTHON_SRC_PATH})

# Only add tests and documentation for current project builds. Protects downstream project builds.
if(CMAKE_PROJECT_NAME STREQUAL PROJECT_NAME)
    # Find Boost. Required for tests
    find_package(Boost 1.53.0 REQUIRED COMPONENTS unit_test_framework)
    # Add tests and docs
    add_subdirectory("${CPP_SRC_PATH}/tests")
    add_subdirectory("docs")
endif()

#==================================================================================== SETUP INSTALLATION CMAKE FILES ===
include(CMakePackageConfigHelpers)
write_basic_package_version_file("${PROJECT_NAME}ConfigVersion.cmake"
                                 VERSION ${PROJECT_VERSION}
                                 COMPATIBILITY SameMajorVersion)
configure_package_config_file(
  "${PROJECT_SOURCE_DIR}/${CMAKE_SRC_PATH}/Config.cmake.in"
  "${PROJECT_BINARY_DIR}/${PROJECT_NAME}Config.cmake"
  INSTALL_DESTINATION
  ${CMAKE_INSTALL_DATAROOTDIR}/${PROJECT_NAME}/cmake)

# CMake won't build the targets for local builds of upstream projects
if(cmake_build_type_lower STREQUAL release)
    install(EXPORT ${PROJECT_NAME}_Targets
            FILE ${PROJECT_NAME}Targets.cmake
            DESTINATION ${CMAKE_INSTALL_DATAROOTDIR}/${PROJECT_NAME}/cmake)
endif()

install(FILES "${PROJECT_BINARY_DIR}/${PROJECT_NAME}Config.cmake"
              "${PROJECT_BINARY_DIR}/${PROJECT_NAME}ConfigVersion.cmake"
        DESTINATION ${CMAKE_INSTALL_DATAROOTDIR}/${PROJECT_NAME}/cmake)
