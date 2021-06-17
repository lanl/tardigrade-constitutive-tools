set(CYTHON_SOURCE_FILES
    main.pyx
    constitutive_tools_python.pxd
    constitutive_tools_python.pyx
    setup.py
    conftest.py
)
add_custom_target(CYTHON ALL
                  DEPENDS ${PROJECT_NAME})
add_custom_command(TARGET CYTHON PRE_BUILD
                   COMMAND ${CMAKE_COMMAND} -E copy_directory
                   "${PROJECT_SOURCE_DIR}/${PYTHON_SRC_PATH}" "${PROJECT_BINARY_DIR}/${PYTHON_SRC_PATH}"
                   COMMAND python setup.py build_ext --inplace
                   DEPENDS ${CYTHON_SOURCE_FILES}
                  )
# Add pytests as a ctest function for automated testing under unified CMake/CTest tools
add_test(NAME pytest
         COMMAND pytest --cov-report=xml:$(pwd)/cobertura.xml --cov=$(pwd)/. --cov-report=html:$(pwd)/htmlcov
         WORKING_DIRECTORY ${PROJECT_BINARY_DIR}/${PYTHON_SRC_PATH}
        )
