/*
===============================================================================
|                            constitutive_tools.cpp                           |
===============================================================================
| A collection of tools useful for constitutive models. These tools are       |
| intended to be generalized expressions which perform operations commonly    |
| encountered in the development of constitutive models. This will enable     |
| users to develop new models quickly and (in principle) with less errors     |
| resulting in delays.                                                        |
|                                                                             |
| Developers should use vector_tools to perform vector multiplications and    |
| matrix solves since this library can be independently checked. Also using   |
| templates and typedef for the data is strongly encouraged so that single    |
| and double precision values (in addition to other types) can be used        |
| simply without significant reworking of the code. Passing by reference is   |
| also encouraged if appropriate to allow users of FORTRAN to interface with  |
| the code without extensive modifications.                                   |
===============================================================================
*/

#include<constitutive_tools.h>

namespace constitutiveTools{

    floatType deltaDirac(const unsigned int i, const unsigned int j){
        /*!
         * The delta dirac function
         * 
         * if i==j return 1
         * if i!=j return 0
         * 
         * :params const unsigned int i: The first index
         * :params const unsigned int j: The second index
         */
    
        if (i==j){
            return 1.;
        }
        return 0;
    }

    errorOut rotateMatrix(const floatVector &A, const floatVector &Q, floatVector &rotatedA){
        /*!
         * Rotate a matrix A using the orthogonal matrix Q with the form
         * A'_ij = Q_{Ii} A_{IJ} Q_{Jj}
         * 
         * TODO: Generalize to non square matrices
         * 
         * :param const floatVector &A: The matrix to be rotated
         * :param const floatVector &Q: The rotation matrix
         * :param floatVector &rotatedA: The rotated matrix
         */

        //Check the size of A
        if (A.size() != Q.size()){
            return new errorNode("rotateMatrix", "A and Q must have the same number of values");
        }

        //Set the dimension to be the square-root of the size of A
        unsigned int dim = std::sqrt(A.size());
        if (A.size() % dim != 0){
            return new errorNode("rotateMatrix", "A must be square");
        }

        //Resize rotated A
        rotatedA.resize(A.size());

        for (unsigned int i=0; i<dim; i++){
            for (unsigned int j=0; j<dim; j++){
                for (unsigned int I=0; I<dim; I++){
                    for (unsigned int J=0; J<dim; J++){
                        rotatedA[i*dim + j] += Q[I*dim + i] * A[I*dim + J] * Q[J*dim + j];
                    }
                }
            }
        }

        return NULL;
    }
    
    errorOut computeGreenLagrangeStrain(const std::vector< floatType > &F,
                                        std::vector< floatType > &E){
        /*!
         * Compute the Green-Lagrange strain from the deformation gradient. The operation is:
         * E = 0.5 (F_{iI} F_{iJ} - delta_{IJ})
         * 
         * Where F is the deformation gradient and delta is the kronecker delta. 
         * 
         * :params std::vector< floatType > &F: A reference to the deformation gradient.
         * :params std::vector< floatType > &E: The resulting Green-Lagrange strain.
         * 
         * The deformation gradient is organized as  F11, F12, F13, F21, F22, F23, F31, F32, F33
         * The Green-Lagrange strain is organized as E11, E12, E13, E21, E22, E23, E31, E32, E33
         */
   
        if (F.size() != 9){
            return new errorNode("computeGreenLagrangeStrain", "The deformation gradient must be 3D.");
        }
    
        const unsigned int dim=3;
        E.resize(dim*dim);
        
        for (unsigned int I=0; I<dim; I++){
            for (unsigned int J=0; J<dim; J++){
                E[dim*I + J] = -deltaDirac(I, J);
                for (unsigned int i=0; i<dim; i++){
                    E[dim*I + J] += F[dim*i + I]*F[dim*i + J];
                }
                E[dim*I + J] *= 0.5;
            }
        }
        return NULL;
    }

    errorOut decomposeGreenLagrangeStrain(const floatVector &E, floatVector &Ebar, floatType &J){
        /*!
         * Decompose the Green-Lagrange strain tensor into isochoric and volumetric parts.
         * where J    = det(F) = sqrt(det(2*E + I))
         *       Ebar_IJ = 0.5*((1/(J**(2/3))) F_iI F_iJ - I_IJ) = (1/(J**(2/3)))*E_IJ
         * 
         * :param const floatVector &E: The Green-Lagrange strain tensor
         * :param floatVector &Ebar: The isochoric Green-Lagrange strain tensor.
         *     format = E11, E12, E13, E21, E22, E23, E31, E32, E33
         * :param floatType &J: The Jacobian of deformation (det(F))
         */

        if (E.size() != 9){
            return new errorNode("decomposeGreenLagrangeStrain", "the Green-Lagrange strain must be 3D");
        }

        //Construct the identity tensor
        floatVector eye = {1, 0, 0, 0, 1, 0, 0, 0, 1};

        J = sqrt(vectorTools::determinant(2*E + eye, 3, 3));
        Ebar = E/(pow(J, 2./3)) + 0.5*(1/pow(J, 2./3) - 1)*eye;
        return NULL;
    }

    errorOut mapPK2toCauchy(const floatVector &PK2Stress, const floatVector &deformationGradient, floatVector &cauchyStress){
        /*!
         * Map the PK2 stress to the current configuration resulting in the Cauchy stress.
         * cauchy_ij = (1/det(F)) F_{iI} PK2_{IJ} F_{jJ}
         * where F is the deformation gradient
         * 
         * :param const floatVector &PK2Stress: The Second Piola-Kirchoff stress
         * :param const floatVector &deformationGradient: The total deformation gradient.
         * :param floatVector &cauchyStress: The Cauchy stress.
         */

        if (PK2Stress.size() != 9){
            return new errorNode("mapPK2toCauchy", "The cauchy stress must have nine components (3D)");
        }

        if (deformationGradient.size() != PK2Stress.size()){
            return new errorNode("mapPK2toCauchy", "The deformation gradient and the PK2 stress don't have the same size");
        }

        //Compute the determinant of the deformation gradient
        floatType detF = vectorTools::determinant(deformationGradient, 3, 3);

        //Initialize the Cauchy stress
        cauchyStress = floatVector(PK2Stress.size(), 0);

        for (unsigned int i=0; i<3; i++){
            for (unsigned int j=0; j<3; j++){
                for (unsigned int I=0; I<3; I++){
                    for (unsigned int J=0; J<3; J++){
                        cauchyStress[3*i + j] += deformationGradient[3*i + I]*PK2Stress[3*I + J]*deformationGradient[3*j + J]/detF;
                    }
                }
            }
        }
        return NULL;
    }

}
