set(CYTHON_SOURCE_FILES
    main.pyx
    constitutive_tools_python.pxd
    constitutive_tools_python.pyx
    setup.py
    conftest.py
)
add_custom_target(CYTHON ALL
                  DEPENDS ${CYTHON_SOURCE_FILES})
add_custom_command(TARGET CYTHON PRE_BUILD
                   COMMAND ${CMAKE_COMMAND} -E copy_directory
                   "${PROJECT_SOURCE_DIR}/${PYTHON_SRC_PATH}" "${PROJECT_BINARY_DIR}/${PYTHON_SRC_PATH}"
                   COMMAND python setup.py build_ext --inplace
                   DEPENDS ${CYTHON_SOURCE_FILES}
                  )
