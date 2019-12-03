//Tests for constitutive_tools

#include<constitutive_tools.h>
#include<sstream>
#include<fstream>
#include<iostream>

struct cout_redirect{
    cout_redirect( std::streambuf * new_buffer)
        : old( std::cout.rdbuf( new_buffer ) )
    { }

    ~cout_redirect( ) {
        std::cout.rdbuf( old );
    }

    private:
        std::streambuf * old;
};

struct cerr_redirect{
    cerr_redirect( std::streambuf * new_buffer)
        : old( std::cerr.rdbuf( new_buffer ) )
    { }

    ~cerr_redirect( ) {
        std::cerr.rdbuf( old );
    }

    private:
        std::streambuf * old;
};

int testDeltaDirac(std::ofstream &results){
    /*!
     * Test the deltaDirac function in constitutive tools
     * 
     * :param std::ofstream &results: The output file
     */

    if (constitutiveTools::deltaDirac(1, 2) != 0){
        results << "deltaDirac & False\n";
        return 1;
    }

    if (constitutiveTools::deltaDirac(1, 1) != 1){
        results << "deltaDirac & False\n";
        return 1;
    }

    results << "deltaDirac & True\n";
    return 0;
    
}

int testComputeGreenLagrangeStrain(std::ofstream &results){
    /*!
     * Test the computation of the Green-Lagrange strain
     * 
     * :param std::ofstream &results: The output file
     */

    constitutiveTools::floatVector F = {1, 0, 0, 0, 1, 0, 0, 0, 1};
    constitutiveTools::floatVector E;

    constitutiveTools::errorOut ret;
    ret = constitutiveTools::computeGreenLagrangeStrain(F, E);

    if (ret){
        results << "testComputeGreenLagrangeStrain (test 1) & False\n";
        return 1;
    }
    
    if (! vectorTools::fuzzyEquals(E, {0, 0, 0, 0, 0, 0, 0, 0, 0})){
        results << "testComputeGreenLagrangeStrain (test 1) & False\n";
        return 1;
    }

    F = {0.69646919, 0.28613933, 0.22685145,
         0.55131477, 0.71946897, 0.42310646,
         0.98076420, 0.68482974, 0.4809319};

    ret = constitutiveTools::computeGreenLagrangeStrain(F, E);

    if (ret){
        results << "testComputeGreenLagrangeStrain (test 2) & False\n";
        return 1;
    }

    if (! vectorTools::fuzzyEquals(E, {0.37545786,  0.63379879,  0.43147034,
                                       0.63379879,  0.03425154,  0.34933978,
                                       0.43147034,  0.34933978, -0.26911192})){
        results << "testComputeGreenLagrangeStrain (test 2) & False\n";
        return 1;
    }

    results << "testComputeGreenLagrangeSTrain & True\n";
    return 0;
}

int main(){
    /*!
    The main loop which runs the tests defined in the 
    accompanying functions. Each function should output
    the function name followed by & followed by True or False 
    if the test passes or fails respectively.
    */

    //Open the results file
    std::ofstream results;
    results.open("results.tex");

    //Run the tests
    testDeltaDirac(results);
    testComputeGreenLagrangeStrain(results);

    //Close the results file
    results.close();

    return 0;
}
