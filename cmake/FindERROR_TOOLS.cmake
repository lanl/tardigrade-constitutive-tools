# - Try to find error_tools
# Once done this will define
#  ERROR_TOOLS_FOUND - System has error_tools
#  ERROR_TOOLS_INCLUDE_DIRS - The error_tools include directories
#  ERROR_TOOLS_LIBRARIES - The libraries needed to use error_tools

# Set hard coded relative path
set(ERROR_TOOLS_RELATIVE_PATH "${PROJECT_SOURCE_DIR}/../error_tools/src/cpp")
set(ERROR_TOOLS_ABS_PATH
    /Users/projects/kbrindley/e13repos/constitutive_models/error_tools/src/cpp)

# Include dir
find_path(ERROR_TOOLS_INCLUDE_DIR
          NAMES error_tools.h
          PATHS ${ERROR_TOOLS_RELATIVE_PATH} ${ERROR_TOOLS_ABS_PATH})

find_library(ERROR_TOOLS_LIBRARY
             NAMES error_tools
             PATHS ${ERROR_TOOLS_RELATIVE_PATH} ${ERROR_TOOLS_ABS_PATH})

include(FindPackageHandleStandardArgs)
# handle the QUIETLY and REQUIRED arguments and set ERROR_TOOLS_FOUND to TRUE
# if all listed variables are TRUE
find_package_handle_standard_args(error_tools DEFAULT_MSG
                                  ERROR_TOOLS_LIBRARY ERROR_TOOLS_INCLUDE_DIR)

mark_as_advanced(ERROR_TOOLS_INCLUDE_DIR ERROR_TOOLS_LIBRARY )

set(ERROR_TOOLS_LIBRARIES ${ERROR_TOOLS_LIBRARY} )
set(ERROR_TOOLS_INCLUDE_DIRS ${ERROR_TOOLS_INCLUDE_DIR} )
