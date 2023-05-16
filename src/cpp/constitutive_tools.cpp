/**
  *****************************************************************************
  * \file constitutive_tools.cpp
  *****************************************************************************
  * A collection of tools useful for constitutive models. These tools are
  * intended to be generalized expressions which perform operations commonly
  * encountered in the development of constitutive models. This will enable
  * users to develop new models quickly and (in principle) with less errors
  * resulting in delays.
  *
  * Developers should use vector_tools to perform vector multiplications and
  * matrix solves since this library can be independently checked. Also using
  * templates and typedef for the data is strongly encouraged so that single
  * and double precision values (in addition to other types) can be used
  * simply without significant reworking of the code. Passing by reference is
  * also encouraged if appropriate to allow users of FORTRAN to interface with
  * the code without extensive modifications.
  *****************************************************************************
  */

#include<constitutive_tools.h>

namespace constitutiveTools{

    floatType deltaDirac(const unsigned int i, const unsigned int j){
        /*!
         * The delta dirac function \f$\delta\f$
         *
         * if i==j return 1
         * if i!=j return 0
         *
         * \param i: The first index
         * \param j: The second index
         */

        if (i==j){
            return 1.;
        }
        return 0;
    }

    errorOut rotateMatrix(const floatVector &A, const floatVector &Q, floatVector &rotatedA){
        /*!
         * Rotate a matrix \f$A\f$ using the orthogonal matrix \f$Q\f$ with the form
         * 
         * \f$A'_{ij} = Q_{Ii} A_{IJ} Q_{Jj}\f$
         *
         * TODO: Generalize to non square matrices
         *
         * \param &A: The matrix to be rotated ( \f$A\f$ )
         * \param &Q: The rotation matrix ( \f$Q\f$Q )
         * \param &rotatedA: The rotated matrix ( \f$A'\f$ )
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

    errorOut computeRightCauchyGreen( const floatVector &deformationGradient, floatVector &C ){
        /*!
         * Compute the Right Cauchy-Green deformation tensor ( \f$C\f$ )
         *
         * \f$C_{IJ} = F_{iI} F_{iJ}\f$
         *
         * \param &deformationGradient: A reference to the deformation gradient ( \f$F\f$ )
         * \param &C: The resulting Right Cauchy-Green deformation tensor ( \f$C\f$ )
         *
         * The deformation gradient is organized as F11, F12, F13, F21, F22, F23, F31, F32, F33
         *
         * The Right Cauchy-Green deformation tensor is organized as C11, C12, C13, C21, C22, C23, C31, C32, C33
         */

        //Assume 3D
        unsigned int dim = 3;

        if ( deformationGradient.size( ) != dim * dim ){
            return new errorNode( "computeRightCauchyGreen",
			          "The deformation gradient must be 3D" );
        }

        C = vectorTools::matrixMultiply( deformationGradient, deformationGradient, dim, dim, dim, dim, 1, 0 );

        return NULL;
    }

    errorOut computeRightCauchyGreen( const floatVector &deformationGradient, floatVector &C, floatMatrix &dCdF ){
        /*!
         * Compute the Right Cauchy-Green deformation tensor ( \f$C\f$ ) from the deformation gradient ( \f$F\f$ )
         * 
         * \f$C_{IJ} = F_{iI} F_{iJ}\f$
         * 
         * \param &deformationGradient: A reference to the deformation gradient ( \f$F\f$ )
         * \param &C: The resulting Right Cauchy-Green deformation tensor ( \f$C\f$ )
         * \param &dCdF: The Jacobian of the Right Cauchy-Green deformation tensor
         *     with regards to the deformation gradient ( \f$\frac{\partial C}{\partial F}\f$ ).
         *
         * The deformation gradient is organized as F11, F12, F13, F21, F22, F23, F31, F32, F33
         *
         * The Right Cauchy-Green deformation tensor is organized as C11, C12, C13, C21, C22, C23, C31, C32, C33
         */

        //Assume 3D
        unsigned int dim = 3;

        errorOut error = computeRightCauchyGreen( deformationGradient, C );

        if ( error ){
            errorOut result = new errorNode( "computeRightCauchyGreen (jacobian)",
                                             "Error in computation of Right Cauchy green deformation tensor" );
            result->addNext( error );
            return result;
        }
        
        //Assemble the Jacobian
        floatVector eye( dim * dim );
        vectorTools::eye( eye );
        
        dCdF = floatMatrix( C.size( ), floatVector( deformationGradient.size( ), 0 ) );
        
        for ( unsigned int I = 0; I < dim; I++ ){
            for ( unsigned int J = 0; J < dim; J++ ){
                for ( unsigned int k = 0; k < dim; k++ ){
                    for ( unsigned int K =0 ; K < dim; K++ ){
                        dCdF[ dim * I + J ][ dim * k + K ] = eye[ dim * I + K ] * deformationGradient[ dim * k + J ]
                                                           + deformationGradient[ dim * k + I ] * eye[ dim * J + K ];
                    }
                }
            }
        }

        return NULL;

    }

    errorOut computeGreenLagrangeStrain( const floatVector &deformationGradient,
                                         floatVector &E ){
        /*!
         * Compute the Green-Lagrange strain ( \f$E\f$ ) from the deformation gradient ( \f$F\f$ ). The operation is:
         *
         * \f$E = 0.5 (F_{iI} F_{iJ} - \delta_{IJ})\f$
         *
         * Where \f$F\f$ is the deformation gradient and \f$\delta\f$ is the kronecker delta.
         *
         * \param &deformationGradient: A reference to the deformation gradient ( \f$F\f$ ).
         * \param &E: The resulting Green-Lagrange strain ( \f$E\f$ ).
         *
         * The deformation gradient is organized as  F11, F12, F13, F21, F22, F23, F31, F32, F33
         *
         * The Green-Lagrange strain is organized as E11, E12, E13, E21, E22, E23, E31, E32, E33
         */

        if ( deformationGradient.size( ) != 9 ){
            return new errorNode( "computeGreenLagrangeStrain", "The deformation gradient must be 3D." );
        }

        const unsigned int dim=3;
        E.resize( dim * dim );

        for ( unsigned int I = 0; I < dim; I++ ){
            for ( unsigned int J = 0; J < dim; J++ ){
                E[ dim * I + J ] = -deltaDirac( I, J );
                for ( unsigned int i = 0; i < dim; i++ ){
                    E[ dim * I + J ] += deformationGradient[ dim * i + I ] * deformationGradient[ dim * i + J ];
                }
                E[ dim * I + J ] *= 0.5;
            }
        }
        return NULL;
    }

    errorOut computeGreenLagrangeStrain( const floatVector &deformationGradient, floatVector &E, floatMatrix &dEdF){
        /*!
         * Compute the Green-Lagrange strain ( \f$E\f$ ) from the deformation gradient ( \f$F\f$ ) and it's jacobian.
         *
         * \param &deformationGradient: A reference to the deformation gradient ( \f$F\f$ ).
         * \param &E: The resulting Green-Lagrange strain ( \f$E\f$ ).
         * \param &dEdF: The jacobian of the Green-Lagrange strain w.r.t. the
         *     deformation gradient ( \f$\frac{\partial E}{\partial F}\f$ ).
         *
         * The deformation gradient is organized as  F11, F12, F13, F21, F22, F23, F31, F32, F33
         *
         * The Green-Lagrange strain is organized as E11, E12, E13, E21, E22, E23, E31, E32, E33
         */

        errorOut error = computeGreenLagrangeStrain( deformationGradient, E );

        if ( error ){
            errorOut result = new errorNode( "computeGreenLagrangeStrain (jacobian)", "Error in computation of Green-Lagrange strain" );
            result->addNext( error );
            return result;
        }

        error = computeDGreenLagrangeStrainDF( deformationGradient, dEdF );

        if ( error ){
            errorOut result = new errorNode( "computeGreenLagrangeStrain (jacobian)", "Error in computation of Green-Lagrange strain jacobian" );
            result->addNext( error );
            return result;
        }
        return NULL;
    }

    errorOut computeDGreenLagrangeStrainDF(const floatVector &deformationGradient, floatMatrix &dEdF){
        /*!
         * Compute the derivative of the Green-Lagrange strain ( \f$E\f$ )w.r.t. the deformation gradient ( \f$F\f$ ).
         *
         * \f$\frac{\partial E_{IJ}}{\partial F_{kK}} = 0.5 ( \delta_{IK} F_{kJ} + F_{kI} \delta_{JK})\f$
         *
         * Where \f$F\f$ is the deformation gradient and \f$\delta\f$ is the kronecker delta.
         *
         * \param &deformationGradient: A reference to the deformation gradient ( \f$F\f$ ).
         * \param &dEdF: The resulting gradient ( \f$\frac{\partial E}{\partial F}\f$ ).
         *
         * The deformation gradient is organized as  F11, F12, F13, F21, F22, F23, F31, F32, F33
         */

        if ( deformationGradient.size( ) != 9 ){
            return new errorNode( "decomposeGreenLagrangeStrain", "the Green-Lagrange strain must be 3D" );
        }

        dEdF = floatMatrix( deformationGradient.size(), floatVector( deformationGradient.size( ), 0 ) );
        for ( unsigned int I = 0; I < 3; I++ ){
            for ( unsigned int J = 0; J < 3; J++ ){
                for ( unsigned int k = 0; k < 3; k++ ){
                    for (unsigned int K = 0; K < 3; K++ ){
                        dEdF[ 3 * I + J ][ 3 * k + K ] = 0.5 * ( deltaDirac( I, K ) * deformationGradient[ 3 * k + J ]
                                                       + deformationGradient[ 3 * k + I ] * deltaDirac( J, K ) );
                    }
                }
            }
        }

        return NULL;
    }

    errorOut decomposeGreenLagrangeStrain( const floatVector &E, floatVector &Ebar, floatType &J ){
        /*!
         * Decompose the Green-Lagrange strain tensor ( \f$E\f$ ) into isochoric ( \f$\bar{E}\f$ ) and volumetric ( \f$J\f$ ) parts where
         *
         * \f$J = det(F) = sqrt(det(2*E + I))\f$
         *
         * \f$\bar{E}_{IJ} = 0.5*((1/(J**(2/3))) F_{iI} F_{iJ} - I_{IJ}) = (1/(J**(2/3)))*E_{IJ} + 0.5(1/(J**(2/3)) - 1)*I_{IJ}\f$
         *
         * \param &E: The Green-Lagrange strain tensor ( \f$E\f$ )
         * \param &Ebar: The isochoric Green-Lagrange strain tensor ( \f$\bar{E}\f$ ).
         *     format = E11, E12, E13, E21, E22, E23, E31, E32, E33
         * \param &J: The Jacobian of deformation ( \f$J\f$ )
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
         * Decompose the Green-Lagrange strain tensor ( \f$E\f$ ) into isochoric ( \f$\bar{E}\f$ ) and volumetric ( \f$J\f$ ) parts where
         *
         * \f$J = det(F) = sqrt(det(2*E + I))\f$
         *
         * \f$\bar{E}_{IJ} = 0.5*((1/(J**(2/3))) F_{iI} F_{iJ} - I_{IJ}) = (1/(J**(2/3)))*E_{IJ} + 0.5(1/(J**(2/3)) - 1)*I_{IJ}\f$
         *
         * \param &E: The Green-Lagrange strain tensor ( \f$E\f$ )
         * \param &Ebar: The isochoric Green-Lagrange strain tensor ( \f$\bar{E}\f$ ).
         *     format = E11, E12, E13, E21, E22, E23, E31, E32, E33
         * \param &J: The Jacobian of deformation ( \f$J\f$ )
         * \param &dEbardE: The derivative of the isochoric Green-Lagrange strain
         *     tensor w.r.t. the total strain tensor ( \f$\frac{\partial \bar{E}}{\partial E}\f$ ).
         * \param &dJdE: The derivative of the jacobian of deformation w.r.t. the
         *     Green-Lagrange strain tensor ( \f$\frac{\partial J}{\partial E}\f$ ).
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
         * Map the PK2 stress ( \f$P^{II}\f$ ) to the current configuration resulting in the Cauchy stress ( \f$\sigma\f$ ).
         *
         * \f$\sigma_{ij} = (1/det(F)) F_{iI} P^{II}_{IJ} F_{jJ}\f$
         *
         * where \f$F\f$ is the deformation gradient
         *
         * \param &PK2Stress: The Second Piola-Kirchoff stress ( \f$P^{II}\f$ )
         * \param &deformationGradient: The total deformation gradient ( \f$F\f$ ).
         * \param &cauchyStress: The Cauchy stress (\f$\sigma\f$ ).
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
         * \f$factor = 10**((-C_1*(T - T_r))/(C_2 + T - T_r))\f$
         *
         * where \f$T\f$ is the temperature, \f$T_r\f$ is the reference temperature, and \f$C_1\f$ and \f$C_2\f$ are parameters
         *
         * \param &temperature: The temperature \f$T\f$
         * \param &WLFParameters: The parameters for the function [\f$T_r\f$, \f$C_1\f$, \f$C_2\f$]
         * \param &factor: The shift factor
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
         * An implementation of the Williams-Landel-Ferry equation that also returns the gradient w.r.t. \f$T\f$
         *
         * \param &temperature: The temperature ( \f$T\f$ )
         * \param &WLFParameters: The parameters for the function [\f$T_r\f$, \f$C_1\f$, \f$C_2\f$]
         * \param &factor: The shift factor
         * \param &dfactordT: The derivative of the shift factor w.r.t. the temperature ( \f$\frac{\partial factor}{\partial T}\f$ )
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
         * \f$\dot{F}_{iI} = L_{ij} F_{jI}\f$
         *
         * \param &velocityGradient: The velocity gradient \f$L_{ij}\f$
         * \param &deformationGradient: The deformation gradient \f$F_{iI}\f$
         * \param &DFDt: The total time derivative of the deformation gradient
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
         * \f$\dot{F}_{iI} = L_{ij} F_{jI}\f$
         * \f$\frac{\partial \dot{F}_{iI}}{\partial L_{kl}} = \delta_{ik} F{lI}\f$
         * \f$\frac{\partial \dot{F}_{iI}}{\partial F_{kK}} = L_{ik} \delta{IK}\f$
         *
         * \param &velocityGradient: The velocity gradient \f$L_{ij}\f$
         * \param &deformationGradient: The deformation gradient \f$F_{iI}\f$
         * \param &DFDt: The total time derivative of the deformation gradient
         * \param &dDFDtdL: The derivative of the total time derivative of the deformation gradient
         *     with respect to the velocity gradient.
         * \param &dDFDtdF: The derivative of the total time derivative of the deformation gradient
         *     with respect to the deformation gradient.
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

    errorOut midpointEvolution( const floatType &Dt, const floatVector &Ap, const floatVector &DApDt, const floatVector &DADt,
                                floatVector &dA, floatVector &A, const floatVector &alpha ){
        /*!
         * Perform midpoint rule based evolution of a vector.
         *
         * alpha=0 (implicit)
         *
         * alpha=1 (explicit)
         *
         * \param &Dt: The change in time.
         * \param &Ap: The previous value of the vector
         * \param &DApDt: The previous time rate of change of the vector.
         * \param &DADt: The current time rate of change of the vector.
         * \param &dA: The change in the vector
         * \param &A: The current value of the vector.
         * \param &alpha: The integration parameter.
         */

        if ( ( Ap.size( ) != DApDt.size( ) ) || ( Ap.size( ) != DADt.size( ) ) ){

            return new errorNode( __func__, "The size of the previous value of the vector and the two rates are not equal" );

        }

        if ( Ap.size( ) != alpha.size( ) ){

            return new errorNode( __func__, "The size of the alpha vector is not the same size as the previous vector value" );

        }

        dA = floatVector( Ap.size( ), 0 );

        A = floatVector( Ap.size( ), 0 );

        unsigned int i = 0;

        for ( auto ai = alpha.begin( ); ai != alpha.end( ); ai++, i++ ){

            if ( ( ( *ai ) < 0) || ( ( *ai ) > 1 ) ){

                return new errorNode( __func__, "Alpha must be between 0 and 1" );

            }

            dA[ i ] = Dt * ( *ai * DApDt[ i ] + ( 1 - *ai ) * DADt[ i ] );

            A[ i ]  = Ap[ i ] + dA[ i ];

        }

        return NULL;

    }

    errorOut midpointEvolution( const floatType &Dt, const floatVector &Ap, const floatVector &DApDt, const floatVector &DADt,
                                floatVector &dA, floatVector &A, floatMatrix &DADADt, const floatVector &alpha ){
        /*!
         * Perform midpoint rule based evolution of a vector and return the jacobian.
         *
         * alpha=0 (implicit)
         *
         * alpha=1 (explicit)
         *
         * \param &Dt: The change in time.
         * \param &Ap: The previous value of the vector
         * \param &DApDt: The previous time rate of change of the vector.
         * \param &DADt: The current time rate of change of the vector.
         * \param &A: The change in value of the vector.
         * \param &A: The current value of the vector.
         * \param &DADADt: The gradient of A w.r.t. the current rate of change.
         * \param &alpha: The integration parameter.
         */

        errorOut error = midpointEvolution( Dt, Ap, DApDt, DADt, dA, A, alpha );

        if ( error ){

            errorOut result = new errorNode( __func__, "Error in computation of the integrated term" );

            result->addNext( error );

            return result;

        }

        DADADt = floatMatrix( A.size( ), floatVector( A.size( ), 0 ) );

        unsigned int i = 0;

        for ( auto ai = alpha.begin( ); ai != alpha.end( ); ai++, i++ ){

            DADADt[ i ][ i ] = Dt * ( 1 - *ai );

        }

        return NULL;

    }

    errorOut midpointEvolution( const floatType &Dt, const floatVector &Ap, const floatVector &DApDt, const floatVector &DADt,
                                floatVector &dA, floatVector &A, floatMatrix &DADADt, floatMatrix &DADADtp,
                                const floatVector &alpha ){
        /*!
         * Perform midpoint rule based evolution of a vector and return the jacobian.
         *
         * alpha=0 (implicit)
         *
         * alpha=1 (explicit)
         *
         * Note that the gradient of A w.r.t. Ap is identity and the gradient of dA w.r.t. Ap is zero
         *
         * \param &Dt: The change in time.
         * \param &Ap: The previous value of the vector
         * \param &DApDt: The previous time rate of change of the vector.
         * \param &DADt: The current time rate of change of the vector.
         * \param &A: The change in value of the vector.
         * \param &A: The current value of the vector.
         * \param &DADADt: The gradient of A w.r.t. the current rate of change.
         * \param &DADADtp: The gradient of A w.r.t. the previous rate of change.
         * \param &alpha: The integration parameter.
         */

        errorOut error = midpointEvolution( Dt, Ap, DApDt, DADt, dA, A, DADADt, alpha );

        if ( error ){

            errorOut result = new errorNode( __func__, "Error in computation of the integrated term" );

            result->addNext( error );

            return result;

        }

        DADADtp = floatMatrix( A.size( ), floatVector( A.size( ), 0 ) );

        unsigned int i = 0;

        for ( auto ai = alpha.begin( ); ai != alpha.end( ); ai++, i++ ){

            DADADtp[ i ][ i ] = Dt * ( *ai );

        }

        return NULL;

    }

    errorOut midpointEvolution( const floatType &Dt, const floatVector &Ap, const floatVector &DApDt, const floatVector &DADt,
                                floatVector &dA, floatVector &A, const floatType alpha ){
        /*!
         * Perform midpoint rule based evolution of a vector. Defaults to the trapezoidal rule.
         *
         * alpha=0 (implicit)
         *
         * alpha=1 (explicit)
         *
         * \param &Dt: The change in time.
         * \param &Ap: The previous value of the vector
         * \param &DApDt: The previous time rate of change of the vector.
         * \param &DADt: The current time rate of change of the vector.
         * \param &dA: The change in the value of the vector.
         * \param &A: The current value of the vector.
         * \param alpha: The integration parameter.
         */

        return midpointEvolution( Dt, Ap, DApDt, DADt, dA, A, alpha * floatVector( Ap.size( ), 1 ) );

    }

    errorOut midpointEvolution( const floatType &Dt, const floatVector &Ap, const floatVector &DApDt, const floatVector &DADt,
                                floatVector &dA, floatVector &A, floatMatrix &DADADt, const floatType alpha ){
        /*!
         * Perform midpoint rule based evolution of a vector. Defaults to the trapezoidal rule.
         *
         * alpha=0 (implicit)
         *
         * alpha=1 (explicit)
         *
         * \param &Dt: The change in time.
         * \param &Ap: The previous value of the vector
         * \param &DApDt: The previous time rate of change of the vector.
         * \param *DADt: The current time rate of change of the vector.
         * \param &dA: The change in the vector
         * \param &A: The current value of the vector.
         * \param &DADADt: The derivative of the vector w.r.t. the rate of change of the vector.
         * \param alpha: The integration parameter.
         */

        return midpointEvolution( Dt, Ap, DApDt, DADt, dA, A, DADADt, alpha * floatVector( Ap.size( ), 1 ) );

    }

    errorOut midpointEvolution( const floatType &Dt, const floatVector &Ap, const floatVector &DApDt, const floatVector &DADt,
                                floatVector &dA, floatVector &A, floatMatrix &DADADt, floatMatrix &DADADtp,
                                const floatType alpha ){
        /*!
         * Perform midpoint rule based evolution of a vector. Defaults to the trapezoidal rule.
         *
         * alpha=0 (implicit)
         *
         * alpha=1 (explicit)
         * 
         * Note that the gradient of A w.r.t. Ap is identity and the gradient of dA w.r.t. Ap is zero
         *
         * \param &Dt: The change in time.
         * \param &Ap: The previous value of the vector
         * \param &DApDt: The previous time rate of change of the vector.
         * \param *DADt: The current time rate of change of the vector.
         * \param &dA: The change in the vector
         * \param &A: The current value of the vector.
         * \param &DADADt: The derivative of the vector w.r.t. the rate of change of the vector.
         * \param &DADADtp: The derivative of the vector w.r.t. the previous rate of change of the vector.
         * \param alpha: The integration parameter.
         */

        return midpointEvolution( Dt, Ap, DApDt, DADt, dA, A, DADADt, DADADtp, alpha * floatVector( Ap.size( ), 1 ) );

    }

    errorOut evolveF(const floatType &Dt, const floatVector &previousDeformationGradient, const floatVector &Lp, const floatVector &L,
                     floatVector &dF, floatVector &deformationGradient, const floatType alpha, const unsigned int mode){
        /*!
         * Evolve the deformation gradient ( F ) using the midpoint integration method.
         *
         * mode 1:
         * \f$F_{iI}^{t + 1} = \left[\delta_{ij} - \Delta t \left(1 - \alpha \right) L_{ij}^{t+1} \right]^{-1} \left[F_{iI}^{t} + \Delta t \alpha \dot{F}_{iI}^{t} \right]\f$
         *
         * mode 2:
         * \f$F_{iI}^{t + 1} = \left[F_{iJ}^{t} + \Delta t \alpha \dot{F}_{iJ}^{t} \right] \left[\delta_{IJ} - \Delta T \left( 1- \alpha \right) L_{IJ}^{t+1} \right]^{-1}\f$
         *
         * \param &Dt: The change in time.
         * \param &previousDeformationGradient: The previous value of the deformation gradient
         * \param &Lp: The previous velocity gradient in the current configuration (mode 1) or
         *     reference configuration (mode 2).
         * \param &L: The current velocity gradient in the current configuration (mode 1) or
         *     reference configuration (mode 2).
         * \param &dF: The change in the deformation gradient \f$\Delta \bf{F}\f$ such that \f$F_{iI}^{t+1} = F_{iI}^t + \Delta F_{iI}\f$
         * \param &deformationGradient: The computed current deformation gradient.
         * \param alpha: The integration parameter.
         * \param mode: The mode of the ODE. Whether the velocity gradient is known in the
         *     current (mode 1) or reference (mode 2) configuration.
         */

        //Assumes 3D
        const unsigned int dim = 3;
        if ( previousDeformationGradient.size( ) != dim * dim ){
            return new errorNode( "evolveF", "The deformation gradient doesn't have enough terms (require 9 for 3D)" );
        }

        if ( Lp.size( ) != previousDeformationGradient.size( ) ){
            return new errorNode( "evolveF", "The previous velocity gradient and deformation gradient aren't the same size" );
        }

        if ( previousDeformationGradient.size( ) != L.size( ) ){
            return new errorNode( "evolveF", "The previous deformation gradient and the current velocity gradient aren't the same size" );
        }

        if ( ( mode != 1 ) && ( mode != 2 ) ){
            return new errorNode( "evolveF", "The mode of evolution is not recognized" );
        }

        //Compute L^{t + \alpha}
        floatVector LtpAlpha = alpha * Lp + ( 1 - alpha ) * L;

        //Compute the right hand side

        floatVector RHS;
        if ( mode == 1 ){
            RHS = vectorTools::matrixMultiply( LtpAlpha, previousDeformationGradient, dim, dim, dim, dim );
        }
        if ( mode == 2 ){
            RHS = vectorTools::matrixMultiply( previousDeformationGradient, LtpAlpha, dim, dim, dim, dim );
        }

        RHS *= Dt;

        //Compute the left-hand side
        floatVector eye( dim * dim );

        vectorTools::eye( eye );

        floatVector LHS = eye - Dt * ( 1 - alpha ) * L;
        floatVector invLHS = vectorTools::inverse( LHS, dim, dim );

        if ( mode == 1 ){

            dF = vectorTools::matrixMultiply( invLHS, RHS, dim, dim, dim, dim );

        }
        if ( mode == 2 ){

            dF = vectorTools::matrixMultiply( RHS, invLHS, dim, dim, dim, dim );

        }

        deformationGradient = previousDeformationGradient + dF;

        return NULL;

    }

    errorOut evolveF(const floatType &Dt, const floatVector &previousDeformationGradient, const floatVector &Lp, const floatVector &L,
                     floatVector &deformationGradient, const floatType alpha, const unsigned int mode){
        /*!
         * Evolve the deformation gradient ( F ) using the midpoint integration method.
         *
         * mode 1:
         * \f$F_{iI}^{t + 1} = \left[\delta_{ij} - \Delta t \left(1 - \alpha \right) L_{ij}^{t+1} \right]^{-1} \left[F_{iI}^{t} + \Delta t \alpha \dot{F}_{iI}^{t} \right]\f$
         *
         * mode 2:
         * \f$F_{iI}^{t + 1} = \left[F_{iJ}^{t} + \Delta t \alpha \dot{F}_{iJ}^{t} \right] \left[\delta_{IJ} - \Delta T \left( 1- \alpha \right) L_{IJ}^{t+1} \right]^{-1}\f$
         *
         * \param &Dt: The change in time.
         * \param &previousDeformationGradient: The previous value of the deformation gradient
         * \param &Lp: The previous velocity gradient in the current configuration (mode 1) or
         *     reference configuration (mode 2).
         * \param &L: The current velocity gradient in the current configuration (mode 1) or
         *     reference configuration (mode 2).
         * \param &deformationGradient: The computed current deformation gradient.
         * \param alpha: The integration parameter.
         * \param mode: The mode of the ODE. Whether the velocity gradient is known in the
         *     current (mode 1) or reference (mode 2) configuration.
         */

        floatVector dF;

        return evolveF( Dt, previousDeformationGradient, Lp, L, dF, deformationGradient, alpha, mode );

    }

    errorOut evolveF( const floatType &Dt, const floatVector &previousDeformationGradient, const floatVector &Lp, const floatVector &L,
                      floatVector &dF, floatVector &deformationGradient, floatMatrix &dFdL, const floatType alpha, const unsigned int mode ){
        /*!
         * Evolve the deformation gradient ( F ) using the midpoint integration method and return the jacobian w.r.t. L.
         *
         * mode 1:
         * \f$F_{iI}^{t + 1} = \left[\delta_{ij} - \Delta t \left(1 - \alpha \right) L_{ij}^{t+1} \right]^{-1} \left[F_{iI}^{t} + \Delta t \alpha \dot{F}_{iI}^{t} \right]\f$
         * \f$\frac{\partial F_{jI}^{t + 1}}{\partial L_{kl}^{t+1}} = \left[\delta_{kj} - \Delta t \left(1 - \alpha\right) L_{kj}\right]^{-1} \Delta t \left(1 - \alpha\right) F_{lI}^{t + 1}\f$
         *
         * mode 2:
         * \f$F_{iI}^{t + 1} = \left[F_{iJ}^{t} + \Delta t \alpha \dot{F}_{iJ}^{t} \right] \left[\delta_{IJ} - \Delta T \left( 1- \alpha \right) L_{IJ}^{t+1} \right]^{-1}\f$
         * \f$\frac{\partial F_{iJ}^{t + 1}}{\partial L_{KL}} = \Delta t (1 - \alpha) F_{iK}^{t + 1} \left[\delta_{JL} - \right ]\f$
         *
         * \param &Dt: The change in time.
         * \param &previousDeformationGradient: The previous value of the deformation gradient
         * \param &Lp: The previous velocity gradient.
         * \param &L: The current velocity gradient.
         * \param &dF: The change in the deformation gradient \f$\Delta \bf{F}\f$ such that \f$F_{iI}^{t+1} = F_{iI}^t + \Delta F_{iI}\f$
         * \param &deformationGradient: The computed current deformation gradient.
         * \param &dFdL: The derivative of the deformation gradient w.r.t. the velocity gradient
         * \param alpha: The integration parameter ( 0 for implicit, 1 for explicit )
         * \param mode: The form of the ODE. See above for details.
         */

        //Assumes 3D
        const unsigned int dim = 3;
        errorOut error = evolveF( Dt, previousDeformationGradient, Lp, L, dF, deformationGradient, alpha, mode);

        if ( error ){

            errorOut result = new errorNode( __func__, "Error when computing the evolved deformation gradient" );
            result->addNext( error );
            return result;

        }

        //Compute the left hand side
        floatVector eye( dim * dim );

        vectorTools::eye( eye );

        floatVector LHS = eye - Dt * ( 1 - alpha ) * L;

        //Compute the inverse of the left-hand side
        floatVector invLHS = vectorTools::inverse( LHS, dim, dim );

        //Compute the jacobian
        dFdL = floatMatrix( deformationGradient.size( ), floatVector( L.size( ), 0 ) );
        if ( mode == 1 ){
            for ( unsigned int j = 0; j < dim; j++ ){
                for ( unsigned int I = 0; I < dim; I++ ){
                    for ( unsigned int k = 0; k < dim; k++ ){
                        for ( unsigned int l = 0; l < dim; l++ ){
                            dFdL[ dim * j + I ][ dim * k + l ] += Dt * ( 1 - alpha ) * invLHS[ dim * j + k ] * deformationGradient[ dim * l + I ];
                        }
                    }
                }
            }
        }
        if ( mode == 2 ){
            for ( unsigned int j = 0; j < dim; j++ ){
                for ( unsigned int I = 0; I < dim; I++ ){
                    for ( unsigned int K = 0; K < dim; K++ ){
                        for ( unsigned int _L = 0; _L < dim; _L++ ){
                            dFdL[ dim * j + I ][ dim * K + _L ] += Dt * ( 1 - alpha ) * invLHS[ dim * _L + I ] * deformationGradient[ dim * j + K ];
                        }
                    }
                }
            }
        }
        return NULL;
    }

    errorOut evolveF( const floatType &Dt, const floatVector &previousDeformationGradient, const floatVector &Lp, const floatVector &L,
                      floatVector &deformationGradient, floatMatrix &dFdL, const floatType alpha, const unsigned int mode ){
        /*!
         * Evolve the deformation gradient ( F ) using the midpoint integration method and return the jacobian w.r.t. L.
         *
         * mode 1:
         * \f$F_{iI}^{t + 1} = \left[\delta_{ij} - \Delta t \left(1 - \alpha \right) L_{ij}^{t+1} \right]^{-1} \left[F_{iI}^{t} + \Delta t \alpha \dot{F}_{iI}^{t} \right]\f$
         * \f$\frac{\partial F_{jI}^{t + 1}}{\partial L_{kl}^{t+1}} = \left[\delta_{kj} - \Delta t \left(1 - \alpha\right) L_{kj}\right]^{-1} \Delta t \left(1 - \alpha\right) F_{lI}^{t + 1}\f$
         *
         * mode 2:
         * \f$F_{iI}^{t + 1} = \left[F_{iJ}^{t} + \Delta t \alpha \dot{F}_{iJ}^{t} \right] \left[\delta_{IJ} - \Delta T \left( 1- \alpha \right) L_{IJ}^{t+1} \right]^{-1}\f$
         * \f$\frac{\partial F_{iJ}^{t + 1}}{\partial L_{KL}} = \Delta t (1 - \alpha) F_{iK}^{t + 1} \left[\delta_{JL} - \right ]\f$
         *
         * \param &Dt: The change in time.
         * \param &previousDeformationGradient: The previous value of the deformation gradient
         * \param &Lp: The previous velocity gradient.
         * \param &L: The current velocity gradient.
         * \param &deformationGradient: The computed current deformation gradient.
         * \param &dFdL: The derivative of the deformation gradient w.r.t. the velocity gradient
         * \param alpha: The integration parameter ( 0 for implicit, 1 for explicit )
         * \param mode: The form of the ODE. See above for details.
         */

        floatVector dF;

        return evolveF( Dt, previousDeformationGradient, Lp, L, dF, deformationGradient, dFdL, alpha, mode );

    }

    errorOut evolveF( const floatType &Dt, const floatVector &previousDeformationGradient, const floatVector &Lp, const floatVector &L,
                      floatVector &dF, floatVector &deformationGradient, floatMatrix &dFdL, floatMatrix &ddFdFp, floatMatrix &dFdFp, floatMatrix &dFdLp, const floatType alpha, const unsigned int mode ){
        /*!
         * Evolve the deformation gradient ( F ) using the midpoint integration method and return the jacobian w.r.t. L.
         *
         * mode 1:
         * \f$F_{iI}^{t + 1} = \left[\delta_{ij} - \Delta t \left(1 - \alpha \right) L_{ij}^{t+1} \right]^{-1} \left[F_{iI}^{t} + \Delta t \alpha \dot{F}_{iI}^{t} \right]\f$
         * \f$\frac{\partial F_{jI}^{t + 1}}{\partial L_{kl}^{t+1}} = \left[\delta_{kj} - \Delta t \left(1 - \alpha\right) L_{kj}\right]^{-1} \Delta t \left(1 - \alpha\right) F_{lI}^{t + 1}\f$
         *
         * mode 2:
         * \f$F_{iI}^{t + 1} = \left[F_{iJ}^{t} + \Delta t \alpha \dot{F}_{iJ}^{t} \right] \left[\delta_{IJ} - \Delta T \left( 1- \alpha \right) L_{IJ}^{t+1} \right]^{-1}\f$
         * \f$\frac{\partial F_{iJ}^{t + 1}}{\partial L_{KL}} = \Delta t (1 - \alpha) F_{iK}^{t + 1} \left[\delta_{JL} - \right ]\f$
         *
         * \param &Dt: The change in time.
         * \param &previousDeformationGradient: The previous value of the deformation gradient
         * \param &Lp: The previous velocity gradient.
         * \param &L: The current velocity gradient.
         * \param &dF: The change in the deformation gradient \f$\Delta \bf{F}\f$ such that \f$F_{iI}^{t+1} = F_{iI}^t + \Delta F_{iI}\f$
         * \param &deformationGradient: The computed current deformation gradient.
         * \param &dFdL: The derivative of the deformation gradient w.r.t. the velocity gradient
         * \param &ddFdFp: The derivative of the change in the deformation gradient w.r.t. the previous deformation gradient
         * \param &dFdFp: The derivative of the deformation gradient w.r.t. the previous deformation gradient
         * \param &dFdLp: The derivative of the deformation gradient w.r.t. the previous velocity gradient
         * \param alpha: The integration parameter ( 0 for implicit, 1 for explicit )
         * \param mode: The form of the ODE. See above for details.
         */

        //Assumes 3D
        const unsigned int dim = 3;
        errorOut error = evolveF( Dt, previousDeformationGradient, Lp, L, dF, deformationGradient, alpha, mode);

        if ( error ){

            errorOut result = new errorNode( __func__, "Error when computing the evolved deformation gradient" );
            result->addNext( error );
            return result;

        }

        //Compute L^{t + \alpha}
        floatVector LtpAlpha = alpha * Lp + ( 1 - alpha ) * L;

        //Compute the left hand side
        floatVector eye( dim * dim );

        vectorTools::eye( eye );

        floatVector LHS = eye - Dt * ( 1 - alpha ) * L;

        //Compute the inverse of the left-hand side
        floatVector invLHS = vectorTools::inverse( LHS, dim, dim );

        //Compute the jacobian
        dFdL   = floatMatrix( deformationGradient.size( ), floatVector( previousDeformationGradient.size( ), 0 ) ); 
        ddFdFp = floatMatrix( deformationGradient.size( ), floatVector( previousDeformationGradient.size( ), 0 ) );
        dFdLp  = floatMatrix( deformationGradient.size( ), floatVector( Lp.size( ), 0 ) );
        if ( mode == 1 ){
            for ( unsigned int j = 0; j < dim; j++ ){
                for ( unsigned int I = 0; I < dim; I++ ){
                    for ( unsigned int k = 0; k < dim; k++ ){
                        for ( unsigned int l = 0; l < dim; l++ ){
                            dFdL[ dim * j + I ][ dim * k + l ]  += Dt * ( 1 - alpha ) * invLHS[ dim * j + k ] * deformationGradient[ dim * l + I ];
                            dFdLp[ dim * j + I ][ dim * k + l ] += invLHS[ dim * j + k ] * Dt * alpha * previousDeformationGradient[ dim * l + I ];
                            for ( unsigned int a = 0; a < dim; a++ ){
                                ddFdFp[ dim * j + I ][ dim * k + l ] += Dt * invLHS[ dim * j + a ] * LtpAlpha[ dim * a + k ] * eye[ dim * I + l ];
                            }
                        }
                    }
                }
            }
        }
        if ( mode == 2 ){
            for ( unsigned int j = 0; j < dim; j++ ){
                for ( unsigned int I = 0; I < dim; I++ ){
                    for ( unsigned int K = 0; K < dim; K++ ){
                        for ( unsigned int _L = 0; _L < dim; _L++ ){
                            dFdL[ dim * j + I ][ dim * K + _L ]  += Dt * ( 1 - alpha ) * invLHS[ dim * _L + I ] * deformationGradient[ dim * j + K ];
                            dFdLp[ dim * j + I ][ dim * K + _L ] += Dt * alpha * invLHS[ dim * _L + I ] * previousDeformationGradient[ dim * j + K ];
                            for ( unsigned int A = 0; A < dim; A++ ){
                                ddFdFp[ dim * j + I ][ dim * K + _L ] += Dt * eye[ dim * j + K ] * LtpAlpha[ dim * _L + A ] * invLHS[ dim * A + I ];
                            }
                        }
                    }
                }
            }
        }

        dFdFp = ddFdFp + vectorTools::eye< floatType >( previousDeformationGradient.size( ) );

        return NULL;
    }

    errorOut evolveF( const floatType &Dt, const floatVector &previousDeformationGradient, const floatVector &Lp, const floatVector &L,
                      floatVector &deformationGradient, floatMatrix &dFdL, floatMatrix &dFdFp, floatMatrix &dFdLp, const floatType alpha, const unsigned int mode ){
        /*!
         * Evolve the deformation gradient ( F ) using the midpoint integration method and return the jacobian w.r.t. L.
         *
         * mode 1:
         * \f$F_{iI}^{t + 1} = \left[\delta_{ij} - \Delta t \left(1 - \alpha \right) L_{ij}^{t+1} \right]^{-1} \left[F_{iI}^{t} + \Delta t \alpha \dot{F}_{iI}^{t} \right]\f$
         * \f$\frac{\partial F_{jI}^{t + 1}}{\partial L_{kl}^{t+1}} = \left[\delta_{kj} - \Delta t \left(1 - \alpha\right) L_{kj}\right]^{-1} \Delta t \left(1 - \alpha\right) F_{lI}^{t + 1}\f$
         *
         * mode 2:
         * \f$F_{iI}^{t + 1} = \left[F_{iJ}^{t} + \Delta t \alpha \dot{F}_{iJ}^{t} \right] \left[\delta_{IJ} - \Delta T \left( 1- \alpha \right) L_{IJ}^{t+1} \right]^{-1}\f$
         * \f$\frac{\partial F_{iJ}^{t + 1}}{\partial L_{KL}} = \Delta t (1 - \alpha) F_{iK}^{t + 1} \left[\delta_{JL} - \right ]\f$
         *
         * \param &Dt: The change in time.
         * \param &previousDeformationGradient: The previous value of the deformation gradient
         * \param &Lp: The previous velocity gradient.
         * \param &L: The current velocity gradient.
         * \param &deformationGradient: The computed current deformation gradient.
         * \param &dFdL: The derivative of the deformation gradient w.r.t. the velocity gradient
         * \param &dFdL: The derivative of the deformation gradient w.r.t. the previous deformation gradient
         * \param &dFdL: The derivative of the deformation gradient w.r.t. the previous velocity gradient
         * \param alpha: The integration parameter ( 0 for implicit, 1 for explicit )
         * \param mode: The form of the ODE. See above for details.
         */

        floatVector dF;
        floatMatrix ddFdFp;

        return evolveF( Dt, previousDeformationGradient, Lp, L, dF, deformationGradient, dFdL, ddFdFp, dFdFp, dFdLp, alpha, mode );

    }

    floatType mac(const floatType &x){
        /*!
         * Compute the Macaulay brackets of a scalar x
         *
         * returns x if x>0, 0 otherwise
         *
         * \param &x: The incoming scalar.
         */

        return 0.5 * ( fabs( x ) + x );
    }

    floatType mac( const floatType &x, floatType &dmacdx ){
        /*!
         * Compute the Macaulay brackets of the scalar x and
         * return the jacobian as well.
         * 
         * returns x if x>0, 0 otherwise
         * 
         * The Jacobian is the Heaviside function
         * 
         * \param &x: The incoming scalar
         * \param &dmacdx: The returned jacobian
         */

        dmacdx = 0;
        if ( x >= 0 ){ dmacdx = 1; }
        return mac( x );
    }

    errorOut computeUnitNormal(const floatVector &A, floatVector &Anorm){
        /*!
         * Compute the unit normal of a second order tensor (or strictly speaking
         * any tensor).
         *
         * \param &A: The second order tensor
         * \param &Anorm: The unit normal in the direction of A
         */

        floatType norm = sqrt(vectorTools::inner(A, A));

        if ( vectorTools::fuzzyEquals( norm, 0. ) ){
            Anorm = floatVector( A.size(), 0 );
        }
        else {
            Anorm = A/norm;
        }

        return NULL;
    }

    errorOut computeUnitNormal(const floatVector &A, floatVector &Anorm, floatMatrix &dAnormdA){
        /*!
         * Compute the unit normal of a second order tensor (or strictly speaking any
         * tensor) and the gradient of that unit normal w.r.t. the tensor.
         *
         * \param &A: The second order tensor
         * \param &Anorm: The unit normal in the direction of A
         * \param &dAnormdA: The gradient of the unit normal w.r.t. A
         */

        floatType norm = sqrt(vectorTools::inner(A, A));

        if ( vectorTools::fuzzyEquals( norm, 0. ) ){
            Anorm = floatVector( A.size(), 0 );
        }
        else {
            Anorm = A/norm;
        }

        floatMatrix eye = vectorTools::eye<floatType>(A.size());

        dAnormdA = (eye - vectorTools::dyadic(Anorm, Anorm))/norm;

        return NULL;
    }

    errorOut pullBackVelocityGradient(const floatVector &velocityGradient, const floatVector &deformationGradient,
                                      floatVector &pullBackVelocityGradient){
        /*!
         * Pull back the velocity gradient to the configuration indicated by deformationGradient, i.e.
         *
         * \f$totalDeformationGradient_{iI} = deformationGradient_{i \bar{I}} remainingDeformationGradient_{\bar{I}I}\f$
         *
         * This is done via
         *
         * \f$L_{\bar{I} \bar{J}} = deformationGradient_{\bar{I} i}^{-1} velocityGradient_{ij} deformationGradient_{j\bar{J}}\f$
         *
         * \param &velocityGradient: The velocity gradient in the current configuration.
         * \param &deformationGradient: The deformation gradient between the desired configuration
         *     and the current configuration.
         * \param &pullBackVelocityGradient: The pulled back velocity gradient.
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
         * Pull back the velocity gradient to the configuration indicated by deformationGradient, i.e.
         *
         * \f$totalDeformationGradient_{iI} = deformationGradient_{i \bar{I}} remainingDeformationGradient_{\bar{I}I}\f$
         *
         * This is done via
         *
         * \f$L_{\bar{I} \bar{J}} = deformationGradient_{\bar{I} i}^{-1} velocityGradient_{ij} deformationGradient_{j\bar{J}}\f$
         *
         * \param &velocityGradient: The velocity gradient in the current configuration.
         * \param &deformationGradient: The deformation gradient between the desired configuration
         *     and the current configuration.
         * \param &pullBackVelocityGradient: The pulled back velocity gradient.
         * \param &dPullBackLdL: The gradient of the pulled back velocity gradient
         *     w.r.t. the velocity gradient.
         * \param &dPullBackLdF: The gradient of the pulled back velocity gradient
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
         * \f$ e^{\theta}_{ij} = a_{ij} \left(\theta - \theta_0\right) + b_{ij} \left(\theta^2 - \theta_0^2\right)\f$
         * 
         * Where \f$e^{\theta}_{ij}\f$ is the thermal expansion, \f$a_{ij}\f$ are the linear parameters,
         * \f$b_{ij}\f$ are the quadratic parameters, \f$\theta\f$ is the current temperature, and \f$\theta_0\f$
         * is the reference temperature.
         *
         * \param &temperature: The temperature
         * \param &referenceTemperature: The reference temperature
         * \param &linearParameters: The linear thermal expansion parameters.
         * \param &quadraticParameters: The quadratic thermal expansion parameters.
         * \param &thermalExpansion: The resulting thermal expansion.
         */

        if (linearParameters.size() != quadraticParameters.size()){
            return new errorNode("quadraticThermalExpansion", "The linear and quadratic parameters must have the same length");
        }

        thermalExpansion = linearParameters * temperature          + quadraticParameters * temperature * temperature
                         - linearParameters * referenceTemperature - quadraticParameters * referenceTemperature * referenceTemperature;

        return NULL;
    }

    errorOut quadraticThermalExpansion(const floatType &temperature, const floatType &referenceTemperature,
                                       const floatVector &linearParameters, const floatVector &quadraticParameters,
                                       floatVector &thermalExpansion, floatVector &thermalExpansionJacobian){
        /*!
         * Define a quadratic equation for the thermal expansion. This could be the
         * thermal strain or the value of the stretch tensor.
         * 
         * \f$ e^{\theta}_{ij} = a_{ij} \left(\theta - \theta_0\right) + b_{ij} \left(\theta^2 - \theta_0^2\right)\f$
         * 
         * Where \f$e^{\theta}_{ij}\f$ is the thermal expansion, \f$a_{ij}\f$ are the linear parameters,
         * \f$b_{ij}\f$ are the quadratic parameters, \f$\theta\f$ is the current temperature, and \f$\theta_0\f$
         * is the reference temperature.
         *
         * \param &temperature: The temperature
         * \param &referenceTemperature: The reference temperature
         * \param &linearParameters: The linear thermal expansion parameters.
         * \param &quadraticParameters: The quadratic thermal expansion parameters.
         * \param &thermalExpansion: The resulting thermal expansion.
         * \param &thermalExpansionJacobian: The gradient of the thermal expansion w.r.t.
         *     the temperature.
         */

        errorOut error = quadraticThermalExpansion(temperature, referenceTemperature, linearParameters, quadraticParameters,
                                                   thermalExpansion);

        if (error){
            errorOut result = new errorNode("quadraticThermalExpansion (jacobian)", "Error in computation of thermal expansion");
            result->addNext(error);
            return result;
        }

        thermalExpansionJacobian = linearParameters + 2 * quadraticParameters * temperature;

        return NULL;
    }

    errorOut pushForwardGreenLagrangeStrain(const floatVector &greenLagrangeStrain, const floatVector &deformationGradient,
                                            floatVector &almansiStrain){
        /*!
         * Push forward the Green-Lagrange strain to the current configuration.
         *
         * \f$e_{ij} = F_{Ii}^{-1} E_{IJ} F_{Jj}^{-1}\f$
         *
         * where \f$e_{ij}\f$ is the Almansi strain (the strain in the current configuration, \f$F_{iI}^{-1}\f$ is the
         * inverse of the deformation gradient, and \f$E_{IJ}\f$ is the Green-Lagrange strain.
         *
         * \param &greenLagrangeStrain: The Green-Lagrange strain.
         * \param &deformationGradient: The deformation gradient mapping between configurations.
         * \param &almansiStrain: The strain in the current configuration indicated by the deformation gradient.
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
         * \f$e_{ij} = F_{Ii}^{-1} E_{IJ} F_{Jj}^{-1}\f$
         *
         * \f$\frac{\partial e_{ij}}{\partial E_{KL}} = F_{Ki}^{-1} F_{Kj}^{-1}\f$
         *
         * \f$\frac{\partial e_{ij}}{\partial F_{kK}} = -F_{Ik}^{-1} F_{Ki}^{-1} E_{IJ} F_{J j}^{-1} - F_{Ii}^{-1} E_{IJ} F_{Jk}^{-1} F_{Kj}^{-1}\f$
         *
         * where \f$e_{ij}\f$ is the Almansi strain (the strain in the current configuration, \f$F_{iI}^{-1}\f$ is the
         * inverse of the deformation gradient, and \f$E_{IJ}\f$ is the Green-Lagrange strain.
         *
         * \param &greenLagrangeStrain: The Green-Lagrange strain.
         * \param &deformationGradient: The deformation gradient mapping between configurations.
         * \param &almansiStrain: The strain in the current configuration indicated by the deformation gradient.
         * \param &dalmansiStraindE: Compute the derivative of the almansi strain w.r.t. the Green-Lagrange strain.
         * \param &dalmansiStraindF: Compute the derivative of the almansi strain w.r.t. the deformation gradient.
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
         * \param &almansiStrain: The strain in the deformation gradient's current configuration.
         * \param &deformationGradient: The deformation gradient between configurations.
         * \param &greenLagrangeStrain: The Green-Lagrange strain which corresponds to the reference
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
         * \param &almansiStrain: The strain in the deformation gradient's current configuration.
         * \param &deformationGradient: The deformation gradient between configurations.
         * \param &greenLagrangeStrain: The Green-Lagrange strain which corresponds to the reference
         *     configuration of the deformation gradient.
         * \param &dEde: The derivative of the Green-Lagrange strain w.r.t. the Almansi strain.
         * \param &dEdF: The derivative of the Green-Lagrange strain w.r.t. the deformation gradient
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

    errorOut computeSymmetricPart( const floatVector &A, floatVector &symmA, unsigned int &dim ){
        /*!
         * Compute the symmetric part of a second order tensor ( \f$A\f$ ) and return it.
         *
         * \f$symm( A )_ij = \frac{1}{2}\left(A_{ij} + A_{ji}\right)\f$
         *
         * \param &A: A constant reference to the second order tensor to process ( \f$A\f$ )
         * \param &symmA: The symmetric part of A ( \f$A^{symm}\f$ )
         * \param &dim: The dimension of A. Note that this is an output used for help
         *     with computing the Jacobian. If you don't need dim as an output use the
         *     version of this function without it.
         */
        
        //Get the dimension of A
        dim = ( unsigned int )( std::sqrt( ( double )A.size( ) ) + 0.5 );
    
        if ( dim * dim != A.size( ) ){
            return new errorNode( "computeSymmetricPart", "A is not a square matrix" );
        }
 
        symmA = floatVector( A.size( ), 0 );

        for ( unsigned int i = 0; i < dim; i++ ){
            for ( unsigned int j = 0; j < dim; j++ ){
                symmA[ dim * i + j ] = 0.5 * ( A[ dim * i + j ] + A[ dim * j + i ] );
            }
        }
    
        return NULL;
    }

    errorOut computeSymmetricPart( const floatVector &A, floatVector &symmA ){
        /*!
         * Compute the symmetric part of a second order tensor ( \f$A\f$ ) and return it.
         *
         * \f$symm( A )_ij = \frac{1}{2}\left(A_{ij} + A_{ji}\right)\f$
         *
         * \param &A: A constant reference to the second order tensor to process ( \f$A\f$ )
         * \param &symmA: The symmetric part of A ( \f$A^{symm}\f$ )
         */
    
        unsigned int dim;
        return computeSymmetricPart( A, symmA, dim );
    }

    errorOut computeSymmetricPart( const floatVector &A, floatVector &symmA, floatMatrix &dSymmAdA ){
        /*!
         * Compute the symmetric part of a second order tensor ( \f$A\f$ ) and return it.
         *
         * \f$( A )^{symm}_{ij} = \frac{1}{2}\left(A_{ij} + A_{ji}\right)\f$
         *
         * Also computes the jacobian
         * 
         * \f$\frac{\partial A^{symm}_{ij}}{\partial A_{kl}} = \frac{1}{2}\left( \delta_{ik} \delta_{jl} + \delta_{jk}\delta_{il} \right)
         *
         * \param &A: A constant reference to the second order tensor to process ( \f$A\f$ )
         * \param &symmA: The symmetric part of A ( \f$A^{symm}\f$ )
         * \param &dSymmAdA: The Jacobian of the symmetric part of A w.r.t. A ( \f$\frac{\partial A^{symm}}{\partial A}\f$ )
         */
        
        unsigned int dim;
        errorOut error = computeSymmetricPart( A, symmA, dim );
        
        if ( error ){
            errorOut result = new errorNode( "computeSymmetricPart (jacobian)",
                                             "Error in computation of the symmetric part of A" );
            result->addNext( error );
            return result;
        }
        
        floatVector eye( A.size( ) );
        vectorTools::eye( eye );
        
        dSymmAdA = floatMatrix( symmA.size( ), floatVector( A.size( ), 0 ) );
        
        for ( unsigned int i = 0; i < dim; i++ ){
            for ( unsigned int j = 0; j < dim; j++ ){
                for ( unsigned int k = 0; k < dim; k++ ){
                    for ( unsigned int l = 0; l < dim; l++ ){
                        dSymmAdA[ dim * i + j ][ dim * k + l ] = 0.5 * ( eye[ dim * i + k ] * eye[ dim * j + l ] + eye[ dim * j + k ] * eye[ dim * i + l ] );
                    }
                }
            }
        }
        
        return NULL;
    }

    errorOut pushForwardPK2Stress( const floatVector &PK2, const floatVector &F, floatVector &cauchyStress ){
        /*!
         * Push the Second Piola-Kirchhoff stress forward to the current configuration resulting in the Cauchy stress
         * 
         * \f$ \sigma_{ij} = \frac{1}{J} F_{iI} S_{IJ} F_{jJ} \f$
         * 
         * \param &PK2: The Second Piola-Kirchhoff stress \f$ S_{IJ} \f$
         * \param &F: The deformation gradient \f$ F_{iI} \f$
         * \param &cauchyStress: The Cauchy stress \f$ \sigma_{ij} \f$
         */

        unsigned int dim = ( unsigned int )( std::sqrt( ( double )PK2.size( ) ) + 0.5 );
        
        if ( dim * dim != PK2.size( ) ){

            return new errorNode( "The PK2 stress must have a size of " + std::to_string( dim * dim ) + " and has a size of " + std::to_string( PK2.size( ) ) );

        }

        if ( PK2.size( ) != F.size( ) ){

            return new errorNode( "The deformation gradient must have a size of " + std::to_string( PK2.size( ) ) + " and has a size of " + std::to_string( F.size( ) ) );

        }

        cauchyStress = floatVector( dim * dim, 0 );

        floatType J = vectorTools::determinant( F, dim, dim );

        cauchyStress = vectorTools::matrixMultiply( F, PK2, dim, dim, dim, dim );

        cauchyStress = vectorTools::matrixMultiply( cauchyStress, F, dim, dim, dim, dim, false, true );

        cauchyStress /= J;

        return NULL;

    }

    errorOut pushForwardPK2Stress( const floatVector &PK2, const floatVector &F, floatVector &cauchyStress,
                                   floatMatrix &dCauchyStressdPK2, floatMatrix &dCauchyStressdF ){
        /*!
         * Push the Second Piola-Kirchhoff stress forward to the current configuration resulting in the Cauchy stress
         * 
         * \f$ \sigma_{ij} = \frac{1}{J} F_{iI} S_{IJ} F_{jJ} \f$
         * 
         * \param &PK2: The Second Piola-Kirchhoff stress \f$ S_{IJ} \f$
         * \param &F: The deformation gradient \f$ F_{iI} \f$
         * \param &cauchyStress: The Cauchy stress \f$ \sigma_{ij} \f$
         * \param &dCauchyStressdPK2: The gradient of the Cauchy stress w.r.t. the PK2 stress
         * \param &dCauchyStressdF: The gradient of the Cauchy stress w.r.t. the deformation gradient
         */

        unsigned int dim = ( unsigned int )( std::sqrt( ( double )PK2.size( ) ) + 0.5 );
        
        if ( dim * dim != PK2.size( ) ){

            return new errorNode( "The PK2 stress must have a size of " + std::to_string( dim * dim ) + " and has a size of " + std::to_string( PK2.size( ) ) );

        }

        if ( PK2.size( ) != F.size( ) ){

            return new errorNode( "The deformation gradient must have a size of " + std::to_string( PK2.size( ) ) + " and has a size of " + std::to_string( F.size( ) ) );

        }

        cauchyStress = floatVector( dim * dim, 0 );

        floatType J = vectorTools::determinant( F, dim, dim );

        floatVector dJdF = vectorTools::computeDDetADA( F, dim, dim );

        cauchyStress = vectorTools::matrixMultiply( F, PK2, dim, dim, dim, dim );

        cauchyStress = vectorTools::matrixMultiply( cauchyStress, F, dim, dim, dim, dim, false, true );

        dCauchyStressdF = vectorTools::dyadic( -cauchyStress / ( J * J ), dJdF );

        cauchyStress /= J;

        dCauchyStressdPK2 = floatMatrix( dim * dim, floatVector( dim * dim, 0 ) );

        floatVector eye( dim * dim, 0 );
        vectorTools::eye( eye );

        for ( unsigned int i = 0; i < dim; i++ ){

            for ( unsigned int j = 0; j < dim; j++ ){

                for ( unsigned int A = 0; A < dim; A++ ){

                    for ( unsigned int B = 0; B < dim; B++ ){

                        dCauchyStressdPK2[ dim * i + j ][ dim * A + B ] += J * F[ dim * i + A ] * F[ dim * j + B ];

                        for ( unsigned int I = 0; I < dim; I++ ){

                            dCauchyStressdF[ dim * i + j ][ dim * A + B ] += eye[ dim * i + A ] * PK2[ dim * B + I ] * F[ dim * j + I ] / J
                                                                           + F[ dim * i + I ] * PK2[ dim * I + B ] * eye[ dim * j + A ] / J;

                        }

                    }

                }

            }

        } 

        return NULL;

    }

}
