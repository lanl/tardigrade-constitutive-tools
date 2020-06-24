project (constitutive_tools)

cmake_minimum_required (VERSION 3.14)

# Find eigen
set(CMAKE_PREFIX_PATH $ENV{HOME} $ENV{HOME}/local/ $ENV{HOME}/.local ${CMAKE_PREFIX_PATH})
set(ENV{EIGEN_ROOT} "$ENV{HOME}/include")
set(ENV{EIGEN_ROOT_DIR} "$ENV{HOME}/include")
find_package (Eigen3 3.3 REQUIRED NO_MODULE)
include_directories(${EIGEN3_INCLUDE_DIR})

# Add the cmake folder to locate the FindSphinx module
set(CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake" ${CMAKE_MODULE_PATH})

# Find related, but external, projects
#find_package (VECTOR_TOOLS REQUIRED)
#find_package (ERROR_TOOLS REQUIRED)

# Add project source directories
add_subdirectory ("src/cpp")
add_subdirectory ("docs")

# Set the c++ standard
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED True)
