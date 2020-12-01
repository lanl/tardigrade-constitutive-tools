import pytest
import numpy as np

import constitutive_tools


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
