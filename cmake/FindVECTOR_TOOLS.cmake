# - Try to find vector_tools
# Once done this will define
#  VECTOR_TOOLS_FOUND - System has vector_tools
#  VECTOR_TOOLS_INCLUDE_DIRS - The vector_tools include directories
#  VECTOR_TOOLS_LIBRARIES - The libraries needed to use vector_tools

# Set hard coded relative path
set(VECTOR_TOOLS_RELATIVE_PATH "${PROJECT_SOURCE_DIR}/../vector_tools/src/cpp")
set(VECTOR_TOOLS_ABS_PATH
    /Users/projects/kbrindley/e13repos/constitutive_models/vector_tools/src/cpp)

# Include dir
find_path(VECTOR_TOOLS_INCLUDE_DIR
          NAMES vector_tools.h
          PATHS ${VECTOR_TOOLS_RELATIVE_PATH} ${VECTOR_TOOLS_ABS_PATH})

find_library(VECTOR_TOOLS_LIBRARY
             NAMES vector_tools
             PATHS ${VECTOR_TOOLS_RELATIVE_PATH} ${VECTOR_TOOLS_ABS_PATH})

include(FindPackageHandleStandardArgs)
# handle the QUIETLY and REQUIRED arguments and set VECTOR_TOOLS_FOUND to TRUE
# if all listed variables are TRUE
find_package_handle_standard_args(vector_tools DEFAULT_MSG
                                  VECTOR_TOOLS_LIBRARY VECTOR_TOOLS_INCLUDE_DIR)

mark_as_advanced(VECTOR_TOOLS_INCLUDE_DIR VECTOR_TOOLS_LIBRARY )

set(VECTOR_TOOLS_LIBRARIES ${VECTOR_TOOLS_LIBRARY} )
set(VECTOR_TOOLS_INCLUDE_DIRS ${VECTOR_TOOLS_INCLUDE_DIR} )
