import pytest
import numpy as np

import constitutive_tools


# Implement a numeric calculation of the Jacobian of a vector function

def approx_fprime(f, x0, relative_epsilon=1e-6, absolute_epsilon=1e-6, args=[], kwargs={}):
    """
    Approximate the derivative of a vector-valued function

    :param function f: The vector valued function to evaluate
    :param np.ndarray x0: The point at which to evaluate the Jacobian
    :param float relative_epsilon: The relative delta value
    :param float absolute_epsilon: The absolute minimum delta value
    :param list args: The arguments for the function
    :param dict kwargs: The keyword arguments for the function
    """

    # Check that the function can be evaluated at the x0 point
    f0 = f(x0, *args, **kwargs)

    # Initialize the Jacobian
    jacobian = np.zeros((f0.size, x0.size))

    for i, xi in enumerate(x0):

        # Perturb x0
        delta = np.zeros(x0.size)
        delta[i] = relative_epsilon * np.abs(xi) + absolute_epsilon

        fm = f(x0 - delta, *args, **kwargs)
        fp = f(x0 + delta, *args, **kwargs)

        jacobian[:, i] = (fp - fm) / (2 * np.linalg.norm(delta))

    return jacobian


# Test the decomposition of the Green-Lagrange strain

deformationGradient = np.array([0.69646919, 0.28613933, 0.22685145,\
                                0.55131477, 0.71946897, 0.42310646,\
                                0.98076420, 0.68482974, 0.4809319]).reshape((3,3))

J = np.linalg.det(deformationGradient)
Fbar = deformationGradient / (J**(1./3))
Ebar = 0.5 * (np.dot(Fbar.T, Fbar) - np.eye(3)).flatten()

greenLagrangeStrain = 0.5 * (np.dot(deformationGradient.T, deformationGradient) - np.eye(3)).flatten()

data = [(greenLagrangeStrain, {'isochoric':Ebar, 'volumetric':J})]

@pytest.mark.parametrize("greenLagrangeStrain,answers",data)
def test_decomposeGreenLagrangeStrain(greenLagrangeStrain, answers):
    """
    Test the decomposition of the Green-Lagrange strain to parts arising
    from the isochoric and volumetric parts of the deformation gradient.

    :param np.ndarray greenLagrangeStrain: The Green-Lagrange strain in
        vector form.
    :param dict answers: The answer dictionary. Expected keys are
        'isochoric' and 'volumetric'
    """

    isochoric_result, volumetric_result = constitutive_tools.py_decomposeGreenLagrangeStrain(greenLagrangeStrain)

    assert np.allclose(isochoric_result, answers['isochoric'])

    assert np.isclose(volumetric_result, answers['volumetric'])


# Test the computation of the midpoint evolution


Dt = 2.5

Ap = np.array([9, 10, 11, 12])

DApDt = np.array([1, 2, 3, 4])

DADt = np.array([5, 6, 7, 8])

alphaVec = np.array([0.1, 0.2, 0.3, 0.4])

answers = {'A':np.array([20.5, 23. , 25.5, 28.])}

data = [(Dt, Ap, DApDt, DADt, alphaVec, answers)]

@pytest.mark.parametrize("Dt, Ap, DApDt, DADt, alphaVec, answers", data)
def test_midpointEvolution(Dt, Ap, DApDt, DADt, alphaVec, answers):
    """
    Test the midpoint evolution function

    :param float Dt: The timestep
    :param np.ndarray Ap: The previous value of A
    :param np.ndarray DApDt: The previous rate of change of A
    :param np.ndarray DADT: The current rate of change of A
    :param np.ndarray alphaVec: The integration parameter vector
    :param dict answers: The answer dictionary
    """

    A = constitutive_tools.py_midpointEvolution(Dt, Ap, DApDt, DADt, alphaVec, False)

    assert np.allclose(A, answers['A'])

    A, DADADT = constitutive_tools.py_midpointEvolution(Dt, Ap, DApDt, DADt, alphaVec, True)

    def function_wrapper(x0, *args, **kwargs):
        """
        Wrapper for the evaluation of the midpoint evolution function

        :param np.ndarray x0: The point about which to compute the Jacobian
        :param list args: A list of arguments
        :param dict kwargs: A dictionary of keyword arguments
        """

        DADt = x0
        
        return constitutive_tools.py_midpointEvolution(Dt, Ap, DApDt, DADt, alphaVec, False)

    assert np.allclose(A, answers['A'])

    x0 = np.copy(DADt)

    jacobian = approx_fprime(function_wrapper, x0)

    assert np.allclose(jacobian, DADADT)
