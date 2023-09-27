/**
  *****************************************************************************
  * \file tardigrade_constitutive_tools.h
  *****************************************************************************
  * A collection of tools useful for constitutive models. These tools are
  * intended to be generalized expressions which perform operations commonly
  * encountered in the development of constitutive models. This will enable
  * users to develop new models quickly and (in principle) with less errors
  * resulting in delays.
  *
  * Developers should use tardigrade_vector_tools to perform vector multiplications and
  * matrix solves since this library can be independently checked. Also using
  * templates and typedef for the data is strongly encouraged so that single
  * and double precision values (in addition to other types) can be used
  * simply without significant reworking of the code. Passing by reference is
  * also encouraged if appropriate to allow users of FORTRAN to interface with
  * the code without extensive modifications.
  *
  * Error handling is taken care of using the tardigrade_error_tools library. All
  * functions should return tardigrade_error_tools::node* objects in the case of an error
  * and NULL otherwise. Wrapping functions in try, except can be used to do
  * this efficiently.
  *****************************************************************************
  */

#ifndef TARDIGRADE_CONSTITUTIVE_TOOLS_H
#define TARDIGRADE_CONSTITUTIVE_TOOLS_H

#define USE_EIGEN
#include<tardigrade_vector_tools.h>
#include<tardigrade_error_tools.h>

namespace tardigradeConstitutiveTools{

    typedef tardigradeErrorTools::Node errorNode; //!< Redefinition for the error node
    typedef errorNode* errorOut; //!< Redefinition for a pointer to the error node
    typedef double floatType; //!< Define the float values type.
    typedef std::vector< floatType > floatVector; //!< Define a vector of floats
    typedef std::vector< std::vector< floatType > > floatMatrix; //!< Define a matrix of floats

    floatType deltaDirac(const unsigned int i, const unsigned int j);

    errorOut rotateMatrix(const floatVector &A, const floatVector &Q, floatVector &rotatedA);

    errorOut computeRightCauchyGreen( const floatVector &deformationGradient, floatVector &C );

    errorOut computeRightCauchyGreen( const floatVector &deformationGradient, floatVector &C, floatMatrix &dCdF );

    errorOut computeGreenLagrangeStrain(const floatVector &deformationGradient, floatVector &E);

    errorOut computeGreenLagrangeStrain(const floatVector &deformationGradient, floatVector &E, floatMatrix &dEdF);

    errorOut computeDGreenLagrangeStrainDF(const floatVector &deformationGradient, floatMatrix &dEdF);

    errorOut decomposeGreenLagrangeStrain(const floatVector &E, floatVector &Ebar, floatType &J);

    errorOut decomposeGreenLagrangeStrain(const floatVector &E, floatVector &Ebar, floatType &J,
                                          floatMatrix &dEbardE, floatVector &dJdE);

    errorOut mapPK2toCauchy(const floatVector &PK2Stress, const floatVector &deformationGradient, floatVector &cauchyStress);

    errorOut WLF(const floatType &temperature, const floatVector &WLFParameters, floatType &factor);

    errorOut WLF(const floatType &temperature, const floatVector &WLFParameters, floatType &factor, floatType &dfactordT);

    errorOut computeDFDt(const floatVector &velocityGradient, const floatVector &deformationGradient, floatVector &DFDt);

    errorOut computeDFDt(const floatVector &velocityGradient, const floatVector &deformationGradient, floatVector &DFDt,
                         floatMatrix &dDFDtdL, floatMatrix &dDFDtdF);

    errorOut midpointEvolution(const floatType &Dt, const floatVector &Ap, const floatVector &DApDt, const floatVector &DADt,
                               floatVector &dA, floatVector &A, const floatVector &alpha);

    errorOut midpointEvolution(const floatType &Dt, const floatVector &Ap, const floatVector &DApDt, const floatVector &DADt,
                               floatVector &dA, floatVector &A, floatMatrix &DADADt, const floatVector &alpha);

    errorOut midpointEvolution(const floatType &Dt, const floatVector &Ap, const floatVector &DApDt, const floatVector &DADt,
                               floatVector &dA, floatVector &A, floatMatrix &DADADt, floatMatrix &DADADtp,
                               const floatVector &alpha);

    errorOut midpointEvolution(const floatType &Dt, const floatVector &Ap, const floatVector &DApDt, const floatVector &DADt,
                               floatVector &dA, floatVector &A, const floatType alpha=0.5);

    errorOut midpointEvolution(const floatType &Dt, const floatVector &Ap, const floatVector &DApDt, const floatVector &DADt,
                               floatVector &dA, floatVector &A, floatMatrix &DADADt, const floatType alpha=0.5);

    errorOut midpointEvolution(const floatType &Dt, const floatVector &Ap, const floatVector &DApDt, const floatVector &DADt,
                               floatVector &dA, floatVector &A, floatMatrix &DADADt, floatMatrix &DADADtp,
                               const floatType alpha=0.5);

    errorOut evolveF(const floatType &Dt, const floatVector &previousDeformationGradient, const floatVector &Lp, const floatVector &L,
                     floatVector &dF, floatVector &deformationGradient, const floatType alpha=0.5, const unsigned int mode = 1);

    errorOut evolveF(const floatType &Dt, const floatVector &previousDeformationGradient, const floatVector &Lp, const floatVector &L,
                     floatVector &dF, floatVector &deformationGradient, floatMatrix &dFdL, const floatType alpha=0.5, const unsigned int mode = 1);

    errorOut evolveF(const floatType &Dt, const floatVector &previousDeformationGradient, const floatVector &Lp, const floatVector &L,
                     floatVector &deformationGradient, floatType alpha=0.5, const unsigned int mode = 1);

    errorOut evolveF(const floatType &Dt, const floatVector &previousDeformationGradient, const floatVector &Lp, const floatVector &L,
                     floatVector &deformationGradient, floatMatrix &dFdL, const floatType alpha=0.5, const unsigned int mode = 1);

    errorOut evolveF(const floatType &Dt, const floatVector &previousDeformationGradient, const floatVector &Lp, const floatVector &L,
                     floatVector &deformationGradient, floatMatrix &dFdL, floatMatrix &dFdFp, floatMatrix &dFdLp, const floatType alpha=0.5, const unsigned int mode = 1);

    errorOut evolveF(const floatType &Dt, const floatVector &previousDeformationGradient, const floatVector &Lp, const floatVector &L,
                     floatVector &dF, floatVector &deformationGradient, floatMatrix &dFdL, floatMatrix &ddFdFp, floatMatrix &dFdFp, floatMatrix &dFdLp, const floatType alpha=0.5, const unsigned int mode = 1);

    floatType mac(const floatType &x);

    floatType mac(const floatType &x, floatType &dmacdx);

    errorOut computeUnitNormal(const floatVector &A, floatVector &Anorm);

    errorOut computeUnitNormal(const floatVector &A, floatVector &Anorm, floatMatrix &dAnormdA);

    errorOut pullBackVelocityGradient(const floatVector &velocityGradient, const floatVector &deformationGradient,
                                      floatVector &pullBackVelocityGradient);

    errorOut pullBackVelocityGradient(const floatVector &velocityGradient, const floatVector &deformationGradient,
                                      floatVector &pullBackVelocityGradient, floatMatrix &dPullBackLdL,
                                      floatMatrix &dPullBackLdF);

    errorOut quadraticThermalExpansion(const floatType &temperature, const floatType &referenceTemperature,
                                       const floatVector &linearParameters, const floatVector &quadraticParameters,
                                       floatVector &thermalExpansion);

    errorOut quadraticThermalExpansion(const floatType &temperature, const floatType &referenceTemperature,
                                       const floatVector &linearParameters, const floatVector &quadraticParameters,
                                       floatVector &thermalExpansion, floatVector &thermalExpansionJacobian);

    errorOut pushForwardGreenLagrangeStrain(const floatVector &greenLagrangeStrain, const floatVector &deformationGradient,
                                            floatVector &almansiStrain);

    errorOut pushForwardGreenLagrangeStrain(const floatVector &greenLagrangeStrain, const floatVector &deformationGradient,
                                            floatVector &almansiStrain, floatMatrix &dalmansiStraindE, floatMatrix &dalmansiStraindF);

    errorOut pullBackAlmansiStrain( const floatVector &almansiStrain, const floatVector &deformationGradient,
                                    floatVector &greenLagrangeStrain );

    errorOut pullBackAlmansiStrain( const floatVector &almansiStrain, const floatVector &deformationGradient,
                                    floatVector &greenLagrangeStrain, floatMatrix &dEde, floatMatrix &dEdF );

    errorOut computeSymmetricPart( const floatVector &A, floatVector &symmA, unsigned int &dim );

    errorOut computeSymmetricPart( const floatVector &A, floatVector &symmA );

    errorOut computeSymmetricPart( const floatVector &A, floatVector &symmA, floatMatrix &dSymmAdA );

    errorOut pushForwardPK2Stress( const floatVector &PK2, const floatVector &F, floatVector &cauchyStress );

    errorOut pushForwardPK2Stress( const floatVector &PK2, const floatVector &F, floatVector &cauchyStress,
                                   floatMatrix &dCauchyStressdPK2, floatMatrix &dCauchyStressdF );

    errorOut pullBackCauchyStress( const floatVector &cauchyStress, const floatVector &F, floatVector &PK2 );

    errorOut pullBackCauchyStress( const floatVector &cauchyStress, const floatVector &F, floatVector &PK2,
                                   floatMatrix &dPK2dCauchyStress, floatMatrix &dPK2dF );

}

#endif
