from libcpp.vector cimport vector

import numpy as np
cimport numpy as np

cimport error_tools_python


cdef extern from "constitutive_tools.h" namespace "constitutiveTools":

    error_tools_python.Node* decomposeGreenLagrangeStrain(const vector[double] &, vector[double] &, double &)
