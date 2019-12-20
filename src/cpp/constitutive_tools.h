/*
===============================================================================
|                             constitutive_tools.h                            |
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
|                                                                             |
| Error handling is taken care of using the error_tools library. All          |
| functions should return error_tools::node* objects in the case of an error  |
| and NULL otherwise. Wrapping functions in try, except can be used to do     |
| this efficiently.                                                           |
===============================================================================
*/

#ifndef CONSTITUTIVE_TOOLS_H
#define CONSTITUTIVE_TOOLS_H

#define USE_EIGEN
#include<vector_tools.h>
#include<error_tools.h>

namespace constitutiveTools{

    typedef errorTools::Node errorNode; //!Redefinition for the error node
    typedef errorNode* errorOut; //!Redefinition for a pointer to the error node
    typedef double floatType; //!Define the float values type.
    typedef std::vector< floatType > floatVector; //! Define a vector of floats
    typedef std::vector< std::vector< floatType > > floatMatrix; //!Define a matrix of floats

    floatType deltaDirac(const unsigned int i, const unsigned int j);

    errorOut rotateMatrix(const floatVector &A, const floatVector &Q, floatVector &rotatedA);

    errorOut computeGreenLagrangeStrain(const floatVector &deformationGradient, floatVector &E);

    errorOut computeDGreenLagrangeStrainDF(const floatVector &F, floatMatrix &dEdF);

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
                               floatVector &A, const floatType alpha=0.5);

    errorOut midpointEvolution(const floatType &Dt, const floatVector &Ap, const floatVector &DApDt, const floatVector &DADt,
                               floatVector &A, floatMatrix &DADADt, const floatType alpha);

    errorOut evolveF(const floatType &Dt, const floatVector &Fp, const floatVector &Lp, const floatVector &L, 
                     floatVector &F, const floatType alpha=0.5);

    errorOut evolveF(const floatType &Dt, const floatVector &Fp, const floatVector &Lp, const floatVector &L,
                     floatVector &F, floatMatrix &dFdL, const floatType alpha=0.5);

    floatType mac(const floatType &x);

    floatType mac(const floatType &x, floatType &dmacdx);

    errorOut computeUnitNormal(const floatVector &A, floatVector &Anorm);

    errorOut computeUnitNormal(const floatVector &A, floatVector &Anorm, floatMatrix &dAnormdA);

}

#endif
