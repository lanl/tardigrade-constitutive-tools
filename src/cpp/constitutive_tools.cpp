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
    
    errorOut computeGreenLagrangeStrain(const floatVector &F,
                                        floatVector &E){
        /*!
         * Compute the Green-Lagrange strain from the deformation gradient. The operation is:
         * E = 0.5 (F_{iI} F_{iJ} - delta_{IJ})
         * 
         * Where F is the deformation gradient and delta is the kronecker delta. 
         * 
         * :params floatVector &F: A reference to the deformation gradient.
         * :params floatVector &E: The resulting Green-Lagrange strain.
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

    errorOut computeDGreenLagrangeStrainDF(const floatVector &F,
                                        floatMatrix &dEdF){
        /*!
         * Compute the derivative of the Green-Lagrange strain w.r.t. the deformation gradient.
         * dE_IJdFkK = 0.5 ( delta_{IK} F_{kJ} + F_{kI} delta_{JK})
         * 
         * Where F is the deformation gradient and delta is the kronecker delta. 
         * 
         * :params floatVector &F: A reference to the deformation gradient.
         * :params floatVector &dEdF: The resulting gradient.
         * 
         * The deformation gradient is organized as  F11, F12, F13, F21, F22, F23, F31, F32, F33
         */

        if (F.size() != 9){
            return new errorNode("decomposeGreenLagrangeStrain", "the Green-Lagrange strain must be 3D");
        }

        dEdF = floatMatrix(F.size(), floatVector(F.size(), 0));
        for (unsigned int I=0; I<3; I++){
            for (unsigned int J=0; J<3; J++){
                for (unsigned int k=0; k<3; k++){
                    for (unsigned int K=0; K<3; K++){
                        dEdF[3*I + J][3*k + K] = 0.5*(deltaDirac(I, K)*F[3*k + J] + F[3*k + I]*deltaDirac(J, K));
                    }
                }
            }
        }

        return NULL;
    }

    errorOut decomposeGreenLagrangeStrain(const floatVector &E, floatVector &Ebar, floatType &J){
        /*!
         * Decompose the Green-Lagrange strain tensor into isochoric and volumetric parts.
         * where J    = det(F) = sqrt(det(2*E + I))
         *       Ebar_IJ = 0.5*((1/(J**(2/3))) F_iI F_iJ - I_IJ) = (1/(J**(2/3)))*E_IJ + 0.5(1/(J**(2/3)) - 1)*I_{IJ}
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

        floatType Jsq = vectorTools::determinant(2*E + eye, 3, 3);
        if (Jsq<=0){
            return new errorNode("decomposeGreenLagrangeStrain", "the determinant of the Green-Lagrange strain is negative");
        }

        J = sqrt(Jsq);
        Ebar = E/(pow(J, 2./3)) + 0.5*(1/pow(J, 2./3) - 1)*eye;
        return NULL;
    }

    errorOut decomposeGreenLagrangeStrain(const floatVector &E, floatVector &Ebar, floatType &J,
                                          floatMatrix &dEbardE, floatVector &dJdE){
        /*!
         * Decompute the Green-Lagrange strain tensor into isochoric and volumetric parts.
         * where J    = det(F) = sqrt(det(2*E + I))
         *       Ebar_IJ = 0.5*((1/(J**(2/3))) F_iI F_iJ - I_IJ) = (1/(J**(2/3)))*E_IJ + 0.5(1/(J**(2/3)) - 1)*I_{IJ}
         * 
         * :param const floatVector &E: The Green-Lagrange strain tensor
         * :param floatVector &Ebar: The isochoric Green-Lagrange strain tensor.
         *     format = E11, E12, E13, E21, E22, E23, E31, E32, E33
         * :param floatType &J: The Jacobian of deformation (det(F))
         * :param floatMatrix &dEbardE: The derivative of the isochoric Green-Lagrange strain 
         *     tensor w.r.t. the total strain tensor.
         * :param floatMatrix &dJdE: The derivative of the jacobian of deformation w.r.t. the 
         *     Green-Lagrange strain tensor.
         */

        errorOut error = decomposeGreenLagrangeStrain(E, Ebar, J);
        if (error){
            errorOut result = new errorNode("decomposeGreenLagrangeStrain", "Error in computation of the isochoric volumetric split");
            result->addNext(error);
            return result;
        }

        //Compute the derivative of the jacobian of deformation w.r.t. the Green-Lagrange strain
        floatVector eye = {1, 0, 0, 0, 1, 0, 0, 0, 1};
        dJdE = J * vectorTools::inverse(2*E + eye, 3, 3);

        //Compute the derivative of the isochoric part of the Green-Lagrange strain w.r.t. the Green-Lagrange strain
        floatMatrix EYE = vectorTools::eye<floatType>(9);

        floatType invJ23 = 1./pow(J, 2./3);
        floatType invJ53 = 1./pow(J, 5./3);

        dEbardE = invJ23*EYE - (2./3)*invJ53*vectorTools::dyadic(E, dJdE) - (1./3)*invJ53*vectorTools::dyadic(eye, dJdE);

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

    errorOut WLF(const floatType &temperature, const floatVector &WLFParameters, floatType &factor){
        /*!
         * An implementation of the Williams-Landel-Ferry equation.
         * 
         * factor = 10**((-C1*(T - Tr))/(C2 + T - Tr))
         * 
         * where T is the temperature, Tr is the reference temperature, and C1 and C2 are parameters
         * 
         * :param const floatType &temperature: The temperature
         * :param const floatVector &WLFParameters: The parameters for the function [Tr, C1, C2]
         * :param floatType &factor: The shift factor
         */

        if (WLFParameters.size() != 3){
            return new errorNode("WLF", "The parameters have the wrong number of terms");
        }

        floatType Tr = WLFParameters[0];
        floatType C1 = WLFParameters[1];
        floatType C2 = WLFParameters[2];

        if (vectorTools::fuzzyEquals(C2 + (temperature - Tr), 0.)){
            return new errorNode("WLF", "Zero in the denominator");
        }

        factor = pow(10., -C1*(temperature - Tr)/(C2 + (temperature - Tr)));

        return NULL;
    }

    errorOut WLF(const floatType &temperature, const floatVector &WLFParameters, floatType &factor, floatType &dfactordT){
        /*!
         * An implementation of the Williams-Landel-Ferry equation that also returns the gradient w.r.t. T
         *
         * :param const floatType &temperature: The temperature
         * :param const floatVector &WLFParameters: The parameters for the function [Tr, C1, C2]
         * :param floatType &factor: The shift factor
         * :param floatType &dfactordT: The derivative of the shift factor w.r.t. the temperature.
         */

        errorOut error = WLF(temperature, WLFParameters, factor);
        if (error){
            errorOut result = new errorNode("WLF", "error in computation of WLF factor");
            result->addNext(error);
            return result;
        }

        floatType Tr = WLFParameters[0];
        floatType C1 = WLFParameters[1];
        floatType C2 = WLFParameters[2];

        dfactordT = std::log(10)*factor*(-C1/(C2 + temperature - Tr) + (C1*(temperature - Tr)/pow(C2 + temperature - Tr, 2)));

        return NULL;
    }

    errorOut midpointEvolution(const floatType &Dt, const floatVector &Ap, const floatVector &DApDt, const floatVector &DADt,
                               floatVector &A, const floatType alpha){
        /*!
         * Perform midpoint rule based evolution of a vector. Defaults to the trapezoidal rule.
         * alpha=0 (implicit)
         * alpha=1 (explicit)
         * 
         * :param const floatType &Dt: The change in time.
         * :param const floatVector &Ap: The previous value of the vector
         * :param const floatVector &DApDt: The previous time rate of change of the vector.
         * :param const floatVector *DADt: The current time rate of change of the vector.
         * :param floatVector &A: The current value of the vector.
         * :param const floatType alpha: The integration parameter.
         */

        if ((Ap.size() != DApDt.size()) || (Ap.size() != DADt.size())){
            return new errorNode("midpointEvolution", "The size of the previous value of the vector and the two rates are not equal");
        }

        if ((alpha > 1) || (alpha < 0)){
            return new errorNode("midpointEvolution", "Alpha must be between 0 and 1");
        }

        A = Ap + Dt*(alpha*DApDt + (1 - alpha)*DADt);
        return NULL;
    }
}
