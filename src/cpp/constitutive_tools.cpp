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

    errorOut computeGreenLagrangeStrain(const floatVector &F, floatVector &E, floatMatrix &dEdF){
        /*!
         * Compute the Green-Lagrange strain from the deformation gradient and it's jacobian.
         * 
         * :param floatVector &F: A reference to the deformation gradient.
         * :param floatVector &E: The resulting Green-Lagrange strain.
         * :param floatMatrix &dEdF: The jacobian of the Green-Lagrange strain w.r.t. the 
         *     deformation gradient.
         * 
         * The deformation gradient is organized as  F11, F12, F13, F21, F22, F23, F31, F32, F33
         * The Green-Lagrange strain is organized as E11, E12, E13, E21, E22, E23, E31, E32, E33
         */

        errorOut error = computeGreenLagrangeStrain(F, E);

        if (error){
            errorOut result = new errorNode("computeGreenLagrangeStrain (jacobian)", "Error in computation of Green-Lagrange strain");
            result->addNext(error);
            return result;
        }

        error = computeDGreenLagrangeStrainDF(F, dEdF);

        if (error){
            errorOut result = new errorNode("computeGreenLagrangeStrain (jacobian)", "Error in computation of Green-Lagrange strain jacobian");
            result->addNext(error);
            return result;
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

    errorOut computeDFDt(const floatVector &velocityGradient, const floatVector &deformationGradient, floatVector &DFDt){
        /*!
         * Compute the total time derivative of the deformation gradient.
         * 
         * \dot{F}_{iI} = L_{ij} F_{jI}
         * 
         * :param const floatVector &velocityGradient: The velocity gradient L_{ij}
         * :param const floatVector &deformationGradient: The deformation gradient F_{iI}
         * :param floatVector &DFDt: The total time derivative of the deformation gradient
         */

        //Assume 3D
        unsigned int dim = 3;

        if (velocityGradient.size() != deformationGradient.size()){
            return new errorNode("computeDFDt", "The velocity gradient and deformation gradient must have the same size");
        }

        if (velocityGradient.size() != dim*dim){
            return new errorNode("computeDFDt", "The velocity gradient doesn't have enough entries");
        }

        DFDt = floatVector(velocityGradient.size(), 0);

        for (unsigned int i=0; i<dim; i++){
            for (unsigned int I=0; I<dim; I++){
                for (unsigned int j=0; j<dim; j++){
                    DFDt[dim*i + I] += velocityGradient[dim*i + j] * deformationGradient[dim*j + I];
                }
            }
        }
        return NULL;
    }

    errorOut computeDFDt(const floatVector &velocityGradient, const floatVector &deformationGradient, floatVector &DFDt,
                         floatMatrix &dDFDtdL, floatMatrix &dDFDtdF){
        /*!
         * Compute the total time derivative of the deformation gradient
         * and return the partial derivatives w.r.t. L and F.
         * 
         * \dot{F}_{iI} = L_{ij} F_{jI}
         * \frac{\partial \dot{F}_{iI}}{\partial L_{kl}} = \delta_{ik} F{lI}
         * \frac{\partial \dot{F}_{iI}}{\partial F_{kK}} = L_{ik} \delta{IK}
         * 
         * :param const floatVector &velocityGradient: The velocity gradient L_{ij}
         * :param const floatVector &deformationGradient: The deformation gradient F_{iI}
         * :param floatVector &DFDt: The total time derivative of the deformation gradient
         */

        //Assume 3D
        unsigned int dim = 3;

        errorOut error = computeDFDt(velocityGradient, deformationGradient, DFDt);

        if (error){
            errorOut result = new errorNode("computeDFDt (jacobian)", "error in computation of DFDt");
            result->addNext(error);
            return result;
        }

        //Form the identity tensor
        floatVector eye(dim*dim, 0);
        vectorTools::eye(eye);

        //Form the partial w.r.t. L and F
        dDFDtdL = floatMatrix(dim*dim, floatVector(dim*dim, 0));
        dDFDtdF = dDFDtdL;

        for (unsigned int i=0; i<dim; i++){
            for (unsigned int I=0; I<dim; I++){
                for (unsigned int k=0; k<dim; k++){
                    for (unsigned int l=0; l<dim; l++){ //Note that we will use l = K for efficiency
                        dDFDtdL[dim*i + I][dim*k + l] = eye[dim*i + k] * deformationGradient[dim*l + I];
                        dDFDtdF[dim*i + I][dim*k + l] = velocityGradient[dim*i + k] * eye[dim*I + l];
                    }
                }
            }
        }

        return NULL;
    }

    errorOut midpointEvolution(const floatType &Dt, const floatVector &Ap, const floatVector &DApDt, const floatVector &DADt,
                               floatVector &A, const floatVector &alpha){
        /*!
         * Perform midpoint rule based evolution of a vector.
         * alpha=0 (implicit)
         * alpha=1 (explicit)
         * 
         * :param const floatType &Dt: The change in time.
         * :param const floatVector &Ap: The previous value of the vector
         * :param const floatVector &DApDt: The previous time rate of change of the vector.
         * :param const floatVector &DADt: The current time rate of change of the vector.
         * :param floatVector &A: The current value of the vector.
         * :param const floatVector &alpha: The integration parameter.
         */

        if ((Ap.size() != DApDt.size()) || (Ap.size() != DADt.size())){
            return new errorNode("midpointEvolution", "The size of the previous value of the vector and the two rates are not equal");
        }
        if (Ap.size() != alpha.size()){
            return new errorNode("midpointEvolution", "The size of the alpha vector is not the same size as the previous vector value");
        }

        A = floatVector(Ap.size(), 0);
        unsigned int i = 0;
        for (auto ai = alpha.begin(); ai != alpha.end(); ai++, i++){
            if (((*ai) < 0) || ((*ai) > 1)){
                return new errorNode("midpointEvolution", "Alpha must be between 0 and 1");
            }
            A[i] = Ap[i] + Dt * (*ai * DApDt[i] + (1 - *ai) * DADt[i]);
        }
        
        return NULL;
    }

    errorOut midpointEvolution(const floatType &Dt, const floatVector &Ap, const floatVector &DApDt, const floatVector &DADt,
                               floatVector &A, floatMatrix &DADADt, const floatVector &alpha){
        /*!
         * Perform midpoint rule based evolution of a vector and return the jacobian.
         * alpha=0 (implicit)
         * alpha=1 (explicit)
         * 
         * :param const floatType &Dt: The change in time.
         * :param const floatVector &Ap: The previous value of the vector
         * :param const floatVector &DApDt: The previous time rate of change of the vector.
         * :param const floatVector &DADt: The current time rate of change of the vector.
         * :param floatVector &A: The current value of the vector.
         * :param floatMatrix &DADADt: The gradient of A w.r.t. the current rate of change.
         * :param const floatVector &alpha: The integration parameter.
         */

        errorOut error = midpointEvolution(Dt, Ap, DApDt, DADt, A, alpha);

        if (error){
            errorOut result = new errorNode("midpointEvolution (jacobian)", "Error in computation of the integrated term");
            result->addNext(error);
            return result;
        }

        DADADt = floatMatrix(A.size(), floatVector(A.size(), 0));
        unsigned int i=0;
        for (auto ai = alpha.begin(); ai != alpha.end(); ai++, i++){
            DADADt[i][i] = Dt * (1 - *ai);
        }
        
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
         * :param const floatVector &DADt: The current time rate of change of the vector.
         * :param floatVector &A: The current value of the vector.
         * :param const floatType alpha: The integration parameter.
         */

        return midpointEvolution(Dt, Ap, DApDt, DADt, A, alpha*floatVector(Ap.size(), 1));
    }

    errorOut midpointEvolution(const floatType &Dt, const floatVector &Ap, const floatVector &DApDt, const floatVector &DADt,
                               floatVector &A, floatMatrix &DADADt, const floatType alpha){
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
         * :param floatMatrix &DADADt: The derivative of the vector w.r.t. the rate of change of the vector.
         * :param const floatType alpha: The integration parameter.
         */

        return midpointEvolution(Dt, Ap, DApDt, DADt, A, DADADt, alpha*floatVector(Ap.size(), 1));
    }

    errorOut evolveF(const floatType &Dt, const floatVector &Fp, const floatVector &Lp, const floatVector &L,
                     floatVector &F, const floatType alpha, const unsigned int mode){
        /*!
         * Evolve F using the midpoint integration method.
         * 
         * mode 1:
         * F_{iI}^{t + 1} = \left[\delta_{ij} - \Delta t \left(1 - \alpha \right) L_{ij}^{t+1} \right]^{-1} \left[F_{iI}^{t} + \Delta t \alpha \dot{F}_{iI}^{t} \right]
         * 
         * mode 2:
         * F_{iI}^{t + 1} = \left[F_{iJ}^{t} + \Delta t \alpha \dot{F}_{iJ}^{t} \right] \left[\delta_{IJ} - \Delta T \left( 1- \alpha \right) L_{IJ}^{t+1} \right]^{-1}
         *
         * :param const floatType &Dt: The change in time.
         * :param const floatVector &Fp: The previous value of the deformation gradient
         * :param const floatVector &Lp: The previous velocity gradient in the current configuration (mode 1) or 
         *     reference configuration (mode 2).
         * :param const floatVector &L: The current velocity gradient in the current configuration (mode 1) or 
         *     reference configuration (mode 2).
         * :param floatVector &F: The computed current deformation gradient.
         * :param const floatType alpha: The integration parameter.
         * :param const unsigned int mode: The mode of the ODE. Whether the velocity gradient is known in the 
         *     current (mode 1) or reference (mode 2) configuration. 
         */

        //Assumes 3D
        const unsigned int dim = 3;
        if (Fp.size() != dim*dim){
            return new errorNode("evolveF", "The deformation gradient doesn't have enough terms (require 9 for 3D)");
        }

        if (Lp.size() != Fp.size()){
            return new errorNode("evolveF", "The previous velocity gradient and deformation gradient aren't the same size");
        }

        if (Fp.size() != L.size()){
            return new errorNode("evolveF", "The previous deformation gradient and the current velocity gradient aren't the same size");
        }

        if ((mode != 1) && (mode != 2)){
            return new errorNode("evolveF", "The mode of evolution is not recognized");
        }

        //Compute the time-rate of change of the previous deformation gradient from the velocity gradient.
        floatVector DFpDt;
        if (mode == 1){
            DFpDt = vectorTools::matrixMultiply(Lp, Fp, dim, dim, dim, dim);
        }
        if (mode == 2){
            DFpDt = vectorTools::matrixMultiply(Fp, Lp, dim, dim, dim, dim);
        }

        //Compute the left-hand side
        floatVector eye(dim*dim);
        vectorTools::eye(eye);
        floatVector LHS = eye - Dt*(1 - alpha)*L;
        
        //Compute the inverse of the left-hand side
        floatVector invLHS = vectorTools::inverse(LHS, dim, dim);

        //Compute the right-hand side
        floatVector RHS = Fp + Dt*alpha*DFpDt;
        F = floatVector(dim*dim, 0);

        //Compute the new value of F
        if (mode == 1){
            F = vectorTools::matrixMultiply(invLHS, RHS, dim, dim, dim, dim);
        }
        if (mode == 2){
            F = vectorTools::matrixMultiply(RHS, invLHS, dim, dim, dim, dim);
        }
        return NULL;
    }

    errorOut evolveF(const floatType &Dt, const floatVector &Fp, const floatVector &Lp, const floatVector &L,
                     floatVector &F, floatMatrix &dFdL, const floatType alpha, const unsigned int mode){
        /*!
         * Evolve F using the midpoint integration method and return the jacobian w.r.t. L.
         * 
         * mode 1:
         * F_{iI}^{t + 1} = \left[\delta_{ij} - \Delta t \left(1 - \alpha \right) L_{ij}^{t+1} \right]^{-1} \left[F_{iI}^{t} + \Delta t \alpha \dot{F}_{iI}^{t} \right]
         * \frac{\partial F_{jI}^{t + 1}}{\partial L_{kl}^{t+1}} &= \left[\delta_{kj} - \Delta t \left(1 - \alpha\right) L_{kj}\right]^{-1} \Delta t \left(1 - \alpha\right) F_{lI]^{t + 1}
         * 
         * mode 2:
         * F_{iI}^{t + 1} = \left[F_{iJ}^{t} + \Delta t \alpha \dot{F}_{iJ}^{t} \right] \left[\delta_{IJ} - \Delta T \left( 1- \alpha \right) L_{IJ}^{t+1} \right]^{-1}
         * \frac{\partial F_{iJ}^{t + 1}{\partial L_{KL}} = \Delta t (1 - \alpha) F_{iK}^{t + 1} \left[\delta_{JL} - 
         * 
         * :param const floatType &Dt: The change in time.
         * :param const floatVector &Fp: The previous value of the deformation gradient
         * :param const floatVector &Lp: The previous velocity gradient.
         * :param const floatVector &L: The current velocity gradient.
         * :param floatVector &F: The computed current deformation gradient.
         * :param const floatType alpha: The integration parameter.
         */

        //Assumes 3D
        const unsigned int dim = 3;
        if (Fp.size() != dim*dim){
            return new errorNode("evolveF", "The deformation gradient doesn't have enough terms (require 9 for 3D)");
        }

        if (Lp.size() != Fp.size()){
            return new errorNode("evolveF", "The previous velocity gradient and deformation gradient aren't the same size");
        }

        if (Fp.size() != L.size()){
            return new errorNode("evolveF", "The previous deformation gradient and the current velocity gradient aren't the same size");
        }

        //Compute the time-rate of change of the previous deformation gradient from the velocity gradient.
        floatVector DFpDt;
        if (mode == 1){
            DFpDt = vectorTools::matrixMultiply(Lp, Fp, dim, dim, dim, dim);
        }
        if (mode == 2){
            DFpDt = vectorTools::matrixMultiply(Fp, Lp, dim, dim, dim, dim);
        }

        //Compute the left-hand side
        floatVector eye(dim*dim);
        vectorTools::eye(eye);
        floatVector LHS = eye - Dt*(1 - alpha)*L;
        
        //Compute the inverse of the left-hand side
        floatVector invLHS = vectorTools::inverse(LHS, dim, dim);

        //Compute the right-hand side
        floatVector RHS = Fp + Dt*alpha*DFpDt;

        //Evolve the deformation gradient
        if (mode == 1){
            F = vectorTools::matrixMultiply(invLHS, RHS, dim, dim, dim, dim);
        }
        if (mode == 2){
            F = vectorTools::matrixMultiply(RHS, invLHS, dim, dim, dim, dim);
        }

        //Compute the jacobian
        dFdL = floatMatrix(F.size(), floatVector(L.size(), 0));
        if (mode == 1){
            for (unsigned int j=0; j<dim; j++){
                for (unsigned int I=0; I<dim; I++){
                    for (unsigned int k=0; k<dim; k++){
                        for (unsigned int l=0; l<dim; l++){
                            dFdL[dim*j + I][dim*k + l] = invLHS[dim*j + k] * Dt * (1 - alpha)*F[dim*l + I];
                        }
                    }
                }
            }
        }
        if (mode == 2){
            for (unsigned int j=0; j<dim; j++){
                for (unsigned int I=0; I<dim; I++){
                    for (unsigned int K=0; K<dim; K++){
                        for (unsigned int L=0; L<dim; L++){
                            dFdL[dim*j + I][dim*K + L] = invLHS[dim*L + I] * Dt * (1 - alpha)*F[dim*j + K];
                        }
                    }
                }
            }
        }
        return NULL;
    }

    floatType mac(const floatType &x){
        /*!
         * Compute the Macaulay brackets of a scalar x
         * 
         * returns x if x>0, 0 otherwise
         * 
         * :param const floatType &x: The incoming scalar.
         */

        return 0.5*(fabs(x) + x);
    }

    floatType mac(const floatType &x, floatType &dmacdx){
        /*!
         * Compute the Macaulay brackets of the scalar x and 
         * return the jacobian as well.
         */

        if (x<0){ dmacdx = 0; }
        else {dmacdx = 1;}
        return mac(x);
    }

    errorOut computeUnitNormal(const floatVector &A, floatVector &Anorm){
        /*!
         * Compute the unit normal of a second order tensor (or strictly speaking 
         * any tensor).
         * 
         * :param const floatVector &A: The second order tensor
         * :param const floatVector &Anorm: The unit normal in the direction of A
         */

        floatType norm = sqrt(vectorTools::inner(A, A));
        
        Anorm = A/norm;

        return NULL;
    }

    errorOut computeUnitNormal(const floatVector &A, floatVector &Anorm, floatMatrix &dAnormdA){
        /*!
         * Compute the unit normal of a second order tensor (or strictly speaking any 
         * tensor) and the gradient of that unit normal w.r.t. the tensor.
         * 
         * :param const floatVector &A: The second order tensor
         * :param const floatVector &Anorm: The unit normal in the direction of A
         * :param const floatMatrix &dAnormdA: The gradient of the unit normal w.r.t. A
         */

        floatType norm = sqrt(vectorTools::inner(A, A));
        
        Anorm = A/norm;

        floatMatrix eye = vectorTools::eye<floatType>(A.size());

        dAnormdA = (eye - vectorTools::dyadic(Anorm, Anorm))/norm;

        return NULL;
    }

    errorOut pullBackVelocityGradient(const floatVector &velocityGradient, const floatVector &deformationGradient,
                                      floatVector &pullBackVelocityGradient){
        /*!
         * Pull back the velocity gradient to the configuration indicated by deformationGradient.
         * i.e. $totalDeformationGradient_{iI} = deformationGradient_{i \bar{I}} remainingDeformationGradient_{\bar{I} I}$
         * This is done via
         * $L_{\bar{I} \bar{J}} = deformationGradient_{\bar{I} i}^{-1} velocityGradient_{ij} deformationGradient_{j \bar{J}}$
         * 
         * :param const floatVector &velocityGradient: The velocity gradient in the current configuration.
         * :param const floatVector &deformationGradient: The deformation gradient between the desired configuration 
         *     and the current configuration.
         * :param floatVector &pullBackVelocityGradient: The pulled back velocity gradient.
         */

        //Assume 3D
        unsigned int dim = 3;

        //Invert the deformation gradient
        floatVector inverseDeformationGradient = vectorTools::inverse(deformationGradient, dim, dim);

        //Pull back the velocity gradient
        pullBackVelocityGradient = vectorTools::matrixMultiply(inverseDeformationGradient, velocityGradient, dim, dim, dim, dim);
        pullBackVelocityGradient = vectorTools::matrixMultiply(pullBackVelocityGradient, deformationGradient, dim, dim, dim, dim);

        return NULL;
    }

    errorOut pullBackVelocityGradient(const floatVector &velocityGradient, const floatVector &deformationGradient,
                                      floatVector &pullBackVelocityGradient, floatMatrix &dPullBackLdL, 
                                      floatMatrix &dPullBackLdF){
        /*!
         * Pull back the velocity gradient to the configuration indicated by deformationGradient.
         * i.e. $totalDeformationGradient_{iI} = deformationGradient_{i \bar{I}} remainingDeformationGradient_{\bar{I} I}$
         * This is done via
         * $L_{\bar{I} \bar{J}} = deformationGradient_{\bar{I} i}^{-1} velocityGradient_{ij} deformationGradient_{j \bar{J}}$
         * 
         * :param const floatVector &velocityGradient: The velocity gradient in the current configuration.
         * :param const floatVector &deformationGradient: The deformation gradient between the desired configuration 
         *     and the current configuration.
         * :param floatVector &pullBackVelocityGradient: The pulled back velocity gradient.
         * :param floatMatrix &dPullBackLdL: The gradient of the pulled back velocity gradient 
         *     w.r.t. the velocity gradient.
         * :param floatMatrix &dPullBackLdF: The gradient of the pulled back velocity gradient 
         *     w.r.t. the deformation gradient.
         */

        //Assume 3D
        unsigned int dim = 3;

        //Invert the deformation gradient
        floatVector inverseDeformationGradient = vectorTools::inverse(deformationGradient, dim, dim);

        //Pull back the velocity gradient
        pullBackVelocityGradient = vectorTools::matrixMultiply(inverseDeformationGradient, velocityGradient, dim, dim, dim, dim);
        pullBackVelocityGradient = vectorTools::matrixMultiply(pullBackVelocityGradient, deformationGradient, dim, dim, dim, dim);

        //Construct the gradients
        dPullBackLdL = floatMatrix(pullBackVelocityGradient.size(), floatVector( velocityGradient.size(), 0));
        dPullBackLdF = floatMatrix(pullBackVelocityGradient.size(), floatVector( deformationGradient.size(), 0));

        floatVector term1 = vectorTools::matrixMultiply(inverseDeformationGradient, velocityGradient, dim, dim, dim, dim);
        term1 = vectorTools::matrixMultiply(term1, deformationGradient, dim, dim, dim, dim);

        floatVector term2 = vectorTools::matrixMultiply(inverseDeformationGradient, velocityGradient, dim, dim, dim, dim);

        for (unsigned int I=0; I<dim; I++){
            for (unsigned int J=0; J<dim; J++){
                for (unsigned int k=0; k<dim; k++){
                    for (unsigned int l=0; l<dim; l++){
                        dPullBackLdL[dim*I + J][dim*k + l] = inverseDeformationGradient[dim*I + k] * deformationGradient[dim*l + J];
                    }

                    for (unsigned int K=0; K<dim; K++){
                        dPullBackLdF[dim*I + J][dim*k + K] += -inverseDeformationGradient[dim*I + k] * term1[dim*K + J]
                                                              + term2[dim*I + k] * deltaDirac(J, K);
                    }
                }
            }
        }

        return NULL;
    }

    errorOut quadraticThermalExpansion(const floatType &temperature, const floatType &referenceTemperature, 
                                       const floatVector &linearParameters, const floatVector &quadraticParameters, 
                                       floatVector &thermalExpansion){
        /*! 
         * Define a quadratic equation for the thermal expansion. This could be the 
         * thermal strain or the value of the stretch tensor.
         * 
         * :param const floatType &temperature: The temperature
         * :param const floatType &referenceTemperature: The reference temperature
         * :param const floatVector &linearParameters: The linear thermal expansion parameters.
         * :param const floatVector &quadraticParameters: The quadratic thermal expansion parameters.
         * :param floatVector &thermalExpansion: The resulting thermal expansion.
         */

        if (linearParameters.size() != quadraticParameters.size()){
            return new errorNode("quadraticThermalExpansion", "The linear and quadratic parameters must have the same length");
        }

        floatType relativeTemperature  = temperature - referenceTemperature;
        thermalExpansion = linearParameters*relativeTemperature + quadraticParameters * relativeTemperature * relativeTemperature;

        return NULL;
    }

    errorOut quadraticThermalExpansion(const floatType &temperature, const floatType &referenceTemperature, 
                                       const floatVector &linearParameters, const floatVector &quadraticParameters, 
                                       floatVector &thermalExpansion, floatVector &thermalExpansionJacobian){
        /*! 
         * Define a quadratic equation for the thermal expansion. This could be the 
         * thermal strain or the value of the stretch tensor.
         * 
         * :param const floatType &temperature: The temperature
         * :param const floatType &referenceTemperature: The reference temperature
         * :param const floatVector &linearParameters: The linear thermal expansion parameters.
         * :param const floatVector &quadraticParameters: The quadratic thermal expansion parameters.
         * :param floatVector &thermalExpansion: The resulting thermal expansion.
         * :param floatVector &thermalExpansionJacobian: The gradient of the thermal expansion w.r.t. 
         *     the temperature.
         */

        errorOut error = quadraticThermalExpansion(temperature, referenceTemperature, linearParameters, quadraticParameters, 
                                                   thermalExpansion);

        if (error){
            errorOut result = new errorNode("quadraticThermalExpansion (jacobian)", "Error in computation of thermal expansion");
            result->addNext(error);
            return result;
        }

        thermalExpansionJacobian = linearParameters + 2 * quadraticParameters * (temperature - referenceTemperature);

        return NULL;
    }

    errorOut pushForwardGreenLagrangeStrain(const floatVector &greenLagrangeStrain, const floatVector &deformationGradient,
                                            floatVector &almansiStrain){
        /*!
         * Push forward the Green-Lagrange strain to the current configuration.
         * 
         * $e_{ij} = F_{Ii}^{-1} E_{IJ} F_{Jj}^{-1}$
         *
         * where $e_{ij}$ is the Almansi strain (the strain in the current configuration, $F_{iI}^{-1}$ is the 
         * inverse of the deformation gradient, and $E_{IJ}$ is the Green-Lagrange strain.
         * 
         * :param const floatVector &greenLagrangeStrain: The Green-Lagrange strain.
         * :param const floatVector &deformationGradient: The deformation gradient mapping between configurations.
         * :param floatVector &almansiStrain: The strain in the current configuration indicated by the deformation gradient.
         */

        //Assume 3D
        unsigned int dim = 3;

        //Compute the inverse deformation gradient
        floatVector inverseDeformationGradient = vectorTools::inverse(deformationGradient, dim, dim);

        //Map the Green-Lagrange strain to the current configuration
        almansiStrain = vectorTools::matrixMultiply(greenLagrangeStrain, inverseDeformationGradient,
                                                    dim, dim, dim, dim, 0, 0);
        almansiStrain = vectorTools::matrixMultiply(inverseDeformationGradient, almansiStrain,
                                                    dim, dim, dim, dim, 1, 0);
        return NULL;
    }

    errorOut pushForwardGreenLagrangeStrain(const floatVector &greenLagrangeStrain, const floatVector &deformationGradient,
                                            floatVector &almansiStrain, floatMatrix &dalmansiStraindE, floatMatrix &dalmansiStraindF){
        /*!
         * Push forward the Green-Lagrange strain to the current configuration 
         * and return the jacobians.
         * 
         * $e_{ij} = F_{Ii}^{-1} E_{IJ} F_{Jj}^{-1}$
         * $\frac{\partial e_{ij}}{\partial E_{KL}} = F_{Ki}^{-1} F_{Kj}^{-1}$
         * $\frac{\partial e_{ij}}{\partial F_{kK}} = -F_{Ik}^{-1} F_{Ki}^{-1} E_{IJ} F_{J j}^{-1} - F_{Ii}^{-1} E_{IJ} F_{Jk}^{-1} F_{Kj}^{-1}$
         *
         * where $e_{ij}$ is the Almansi strain (the strain in the current configuration, $F_{iI}^{-1}$ is the 
         * inverse of the deformation gradient, and $E_{IJ}$ is the Green-Lagrange strain.
         * 
         * :param const floatVector &greenLagrangeStrain: The Green-Lagrange strain.
         * :param const floatVector &deformationGradient: The deformation gradient mapping between configurations.
         * :param floatVector &almansiStrain: The strain in the current configuration indicated by the deformation gradient.
         * :param floatMatrix &dalmansiStraindE: Compute the derivative of the almansi strain w.r.t. the Green-Lagrange strain.
         * :param floatMatrix &dalmansiStraindF: Compute the derivative of the almansi strain w.r.t. the deformation gradient.
         */

        //Assume 3D
        unsigned int dim = 3;

        //Compute the inverse deformation gradient
        floatVector inverseDeformationGradient = vectorTools::inverse(deformationGradient, dim, dim);

        //Compute the jacobian of the inverse deformation gradient
        floatMatrix dFinvdF(dim*dim, floatVector(dim*dim, 0));
        for (unsigned int I=0; I<dim; I++){
            for (unsigned int l=0; l<dim; l++){
                for (unsigned int k=0; k<dim; k++){
                    for (unsigned int K=0; K<dim; K++){
                        dFinvdF[dim*I + l][dim*k + K] = -inverseDeformationGradient[dim*I + k] *
                                                         inverseDeformationGradient[dim*K + l];
                    }
                }
            }
        }

        //Map the Green-Lagrange strain to the current configuration
        almansiStrain = vectorTools::matrixMultiply(greenLagrangeStrain, inverseDeformationGradient,
                                                    dim, dim, dim, dim, 0, 0);
        almansiStrain = vectorTools::matrixMultiply(inverseDeformationGradient, almansiStrain,
                                                    dim, dim, dim, dim, 1, 0);

        //Compute the jacobians
        dalmansiStraindE = floatMatrix(dim*dim, floatVector(dim*dim, 0));
        for (unsigned int i=0; i<dim; i++){
            for (unsigned int j=0; j<dim; j++){
                for (unsigned int K=0; K<dim; K++){
                    for (unsigned int L=0; L<dim; L++){
                        dalmansiStraindE[dim*i + j][dim*K + L] = inverseDeformationGradient[dim*K + i] * 
                                                                 inverseDeformationGradient[dim*L + j];
                    }
                }
            }
        }

        dalmansiStraindF = floatMatrix(dim*dim, floatVector(dim*dim, 0));
        floatVector term1 = vectorTools::matrixMultiply(greenLagrangeStrain, inverseDeformationGradient,
                                                        dim, dim, dim, dim, 0, 0);
        floatVector term2 = vectorTools::matrixMultiply(inverseDeformationGradient, greenLagrangeStrain,
                                                        dim, dim, dim, dim, 1, 0);

        for (unsigned int i=0; i<dim; i++){
            for (unsigned int j=0; j<dim; j++){
                for (unsigned int k=0; k<dim; k++){
                    for (unsigned int K=0; K<dim; K++){
                        for (unsigned int I=0; I<dim; I++){
                            dalmansiStraindF[dim*i + j][dim*k + K] += dFinvdF[dim*I + i][dim*k + K] * term1[dim*I + j]
                                                                    + term2[dim*i + I] * dFinvdF[dim*I + j][dim*k + K];
                        }
                    }
                }
            }
        }

        return NULL;
    }

    errorOut pullBackAlmansiStrain( const floatVector &almansiStrain, const floatVector &deformationGradient,
                                    floatVector &greenLagrangeStrain ){
        /*!
         * Pull back the almansi strain to the configuration indicated by the deformation gradient.
         * 
         * :param const floatVector &almansiStrain: The strain in the deformation gradient's current configuration.
         * :param const floatVector &deformationGradient: The deformation gradient between configurations.
         * :param floatVector &greenLagrangeStrain: The Green-Lagrange strain which corresponds to the reference 
         *     configuration of the deformation gradient.
         */

        //Assume 3d
        unsigned int dim = 3;

        greenLagrangeStrain = vectorTools::matrixMultiply( deformationGradient, almansiStrain, dim, dim, dim, dim, 1, 0 );
        greenLagrangeStrain = vectorTools::matrixMultiply( greenLagrangeStrain, deformationGradient, dim, dim, dim, dim, 0, 0 );

        return NULL;
    }

    errorOut pullBackAlmansiStrain( const floatVector &almansiStrain, const floatVector &deformationGradient,
                                    floatVector &greenLagrangeStrain, floatMatrix &dEde, floatMatrix &dEdF ){
        /*!
         * Pull back the almansi strain to the configuration indicated by the deformation gradient.
         * 
         * Also return the Jacobians.
         * 
         * :param const floatVector &almansiStrain: The strain in the deformation gradient's current configuration.
         * :param const floatVector &deformationGradient: The deformation gradient between configurations.
         * :param floatVector &greenLagrangeStrain: The Green-Lagrange strain which corresponds to the reference 
         *     configuration of the deformation gradient.
         */

        //Assume 3d
        unsigned int dim = 3;

        errorOut error = pullBackAlmansiStrain( almansiStrain, deformationGradient, greenLagrangeStrain );

        if ( error ){
            errorOut result = new errorNode( "pullBackAlmansiStrain (jacobian)",
                                             "Error in computation of Green-Lagrange strain" );
            result->addNext( error );
            return result;
        }

        floatVector eye( dim * dim );
        vectorTools::eye( eye );

        dEde = floatMatrix( dim * dim, floatVector( dim * dim, 0 ) );
        dEdF = floatMatrix( dim * dim, floatVector( dim * dim, 0 ) );

        for ( unsigned int I = 0; I < dim; I++ ){
            for ( unsigned int J = 0; J < dim; J++ ){
                for ( unsigned int K = 0; K < dim; K++ ){
                    for ( unsigned int L = 0; L < dim; L++ ){
                        dEde[ dim * I + J ][ dim * K + L ] = deformationGradient[ dim * K + I ] * deformationGradient[ dim * L + J ];
                        for ( unsigned int j = 0; j < dim; j++ ){
                            dEdF[ dim * I + J ][ dim * K + L ] += eye[ dim * I + L ] * almansiStrain[ dim * K + j ] * deformationGradient[ dim * j + J ]
                                                                + deformationGradient[ dim * j + I ] * almansiStrain[ dim * j + K ] * eye[ dim * J + L ];
                        }
                    }
                }
            }
        }

        return NULL;
    }
}
