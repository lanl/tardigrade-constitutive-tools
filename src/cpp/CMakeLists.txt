add_library (constitutive_tools "constitutive_tools.cpp" "constitutive_tools.h")

target_include_directories(constitutive_tools PUBLIC .)
set_target_properties(constitutive_tools PROPERTIES CXX_STANDARD 11)
target_compile_options(constitutive_tools PUBLIC
                       -Wall -ansi -pedantic -O3 -fmax-errors=5 -ggdb)

target_link_libraries(constitutive_tools Eigen3::Eigen)
target_include_directories(constitutive_tools PUBLIC "${error_tools_SOURCE_DIR}/src/cpp" "${vector_tools_SOURCE_DIR}/src/cpp")
