set(CYTHON_SOURCE_FILES
    conftest.py
    constitutive_tools_python.pxd
    constitutive_tools_python.pyx
    main.pyx
    settings.py.in
    setup.py
)
configure_file(settings.py.in settings.py)
add_custom_target(${PROJECT_NAME}_CYTHON ALL
                  DEPENDS ${PROJECT_NAME})
add_custom_command(TARGET ${PROJECT_NAME}_CYTHON PRE_BUILD
                   COMMAND ${CMAKE_COMMAND} -E copy_directory
                   "${PROJECT_SOURCE_DIR}/${PYTHON_SRC_PATH}" "${PROJECT_BINARY_DIR}/${PYTHON_SRC_PATH}"
                   COMMAND python setup.py build_ext --inplace
                   DEPENDS ${CYTHON_SOURCE_FILES}
                  )
install(FILES 
            "${PROJECT_SOURCE_DIR}/${PYTHON_SRC_PATH}/${PROJECT_NAME}_python.pyx"
            "${PROJECT_SOURCE_DIR}/${PYTHON_SRC_PATH}/${PROJECT_NAME}_python.pxd"
        DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})
# Add pytests as a ctest function for automated testing under unified CMake/CTest tools
add_test(NAME pytest
         COMMAND
             ${CMAKE_COMMAND} -E env LD_LIBRARY_PATH=${CMAKE_BINARY_DIR}/${CPP_SRC_PATH}:$ENV{LD_LIBRARY_PATH}
             pytest --cov-report=xml:$(pwd)/cobertura.xml --cov=$(pwd)/. --cov-report=html:$(pwd)/htmlcov
         WORKING_DIRECTORY ${PROJECT_BINARY_DIR}/${PYTHON_SRC_PATH}
        )
