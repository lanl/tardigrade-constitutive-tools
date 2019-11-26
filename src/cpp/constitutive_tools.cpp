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
}
