from libcpp.vector cimport vector

import numpy as np
cimport numpy as np

cimport error_tools_python
cimport constitutive_tools_python


cdef map_1D_array_to_vector(np.ndarray array):
    """
    Map a 1D numpy array to a vector

    TODO: This is not a very efficient way to do this

    TODO: This should *probably* go into a wrapper for vector tools

    :param np.ndarray array: The array to map to the standard vector
    """

    cdef vector[double] cvector;
    cdef int i

    cvector.reserve(array.size)

    for i in range(array.size):
        cvector.push_back(array[i])

    return cvector


cdef map_vector_to_1D_array(vector[double] cvector):
    """
    Map a c vector to a numpy array

    TODO: This is not a very efficient way to do this

    TODO: This should *probably* go into a wrapper for vector tools

    :param vector[double] cvector: The c vector to map to a numpy array
    """

    array = np.zeros(cvector.size())

    for i in range(array.size):
        array[i] = cvector[i]

    return array


cdef map_2D_array_to_matrix(np.ndarray array):
    """
    Map a 2D numpy array to a vector

    TODO: This is not a very efficient way to do this

    :param np.ndarray array: The array to map to the matrix
    """

    if ((array.shape[0] == 0 ) or (array.shape[1] == 0)):
        raise ValueError("The array is empty")

    cdef vector[vector[double]] cmatrix;
    cdef int i

    cmatrix.resize(array.shape[0])

    for i in range(array.shape[0]):
        cmatrix[i] = map_1D_array_to_vector(array[i])

    return cmatrix


cdef map_matrix_to_2D_array(vector[vector[double]] matrix):
    """
    Map a C++ matrix to an array

    TODO: This is not a very efficient way to do this

    :param np.ndarray array: The matrix to map to the array
    """

    cdef rows = matrix.size()
    if (rows == 0):
        return np.empty((0, 0))

    cols = matrix[0].size()

    array = np.zeros((rows, cols))

    for i in range(rows):
        array[i, :] = map_vector_to_1D_array(matrix[i])

    return array


def py_decomposeGreenLagrangeStrain(np.ndarray greenLagrangeStrain):
    """
    Wrapper for the C++ function constitutiveTools::decomposeGreenLagrangeStrain
    that breaks the strain into isochoric and volumetric parts.

    :param np.ndarray greenLagrangeStrain: The Green-Lagrange strain in vector
        notation [E11, E12, E13, E21, E22, E23, E31, E32, E33]
    """

    cdef vector[double] c_greenLagrangeStrain
    cdef vector[double] c_isochoricGreenLagrangeStrain
    cdef double c_volumetricGreenLagrangeStrain = 0
    cdef error_tools_python.Node *error

    cdef np.ndarray isochoricGreenLagrangeStrain

    c_greenLagrangeStrain = map_1D_array_to_vector(greenLagrangeStrain)

    error = constitutive_tools_python.decomposeGreenLagrangeStrain(c_greenLagrangeStrain, c_isochoricGreenLagrangeStrain, c_volumetricGreenLagrangeStrain)

    if error:
        error.c_print(True)
        raise ValueError("Error in decompose Green-Lagrange strain")

    isochoricGreenLagrangeStrain = map_vector_to_1D_array(c_isochoricGreenLagrangeStrain)

    return isochoricGreenLagrangeStrain, c_volumetricGreenLagrangeStrain

def py_midpointEvolution(Dt, np.ndarray Ap, np.ndarray DApDt, np.ndarray DADt, np.ndarray alpha,
                         compute_jacobians=False):
    """
    Wrapper for the c++ function constitutiveTools::midpointEvolution that computes the integration of a
    function using an implicit midpoint rule

    :param float Dt: The change in time
    :param np.ndarray Ap: The previous value of the vector
    :param np.ndarray DApDt: The time rate of change of the vector at the previous timestep
    :param np.ndarray DADt: The time rate of change of the vector at the current timestep
    :param np.ndarray alpha: The vector of integration parameters
    :param bool compute_jacobians: A flag indicating if the Jacobians should be computed

    :returns: The updated value of A and, if compute_jacobians is True, the Jacobian w.r.t. the current
        time rate of change of A

    :rtype: np.ndarray and, if compute_jacobians is True, np.ndarray
    """

    cdef double c_Dt = Dt
    cdef vector[double] c_Ap = map_1D_array_to_vector(Ap)
    cdef vector[double] c_DApDt = map_1D_array_to_vector(DApDt)
    cdef vector[double] c_DADt = map_1D_array_to_vector(DADt)
    cdef vector[double] c_alpha = map_1D_array_to_vector(alpha)
    cdef vector[double] c_A
    cdef error_tools_python.Node *error

    cdef vector[vector[double]] c_DADADT

    cdef np.ndarray A
    cdef np.ndarray DADADT

    if compute_jacobians:

        error = constitutive_tools_python.midpointEvolution(c_DT, c_Ap, c_DApDt, c_DADt, c_A, c_DADADT, alpha)

        if error:
            error.c_print(True)
            raise ValueError("Error in the midpoint evolution function")

        A = map_vector_to_1D_array(c_A)

        DADADT = map_matrix_to_2D_array(c_DADADT)

        return A, DADADT

    else:

        error = constitutive_tools_python.midpointEvolution(c_DT, c_Ap, c_DApDt, c_DADt, c_A, alpha)

        if error:
            error.c_print(True)
            raise ValueError("Error in the midpoint evolution function")

        A = map_vector_to_1D_array(c_A)

        return A
