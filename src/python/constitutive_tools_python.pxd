from libcpp.vector cimport vector

import numpy as np
cimport numpy as np

cimport error_tools_python


# Internal functions
cdef map_1D_array_to_vector(np.ndarray array)

cdef map_vector_to_1D_array(vector[double] cvector)

cdef map_2D_array_to_matrix(np.ndarray array)

cdef map_matrix_to_2D_array(vector[vector[double]] matrix)

cdef extern from "constitutive_tools.h" namespace "constitutiveTools":

    error_tools_python.Node* decomposeGreenLagrangeStrain(const vector[double] &, vector[double] &, double &)
