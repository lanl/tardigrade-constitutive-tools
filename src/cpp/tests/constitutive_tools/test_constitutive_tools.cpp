//Tests for constitutive_tools

#include<constitutive_tools.h>
#include<sstream>
#include<fstream>
#include<iostream>

typedef constitutiveTools::errorOut errorOut;
typedef constitutiveTools::floatType floatType;
typedef constitutiveTools::floatVector floatVector;
typedef constitutiveTools::floatMatrix floatMatrix;

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

int testRotateMatrix(std::ofstream &results){
    /*!
     * Test the rotation of a matrix by an orthogonal rotation matrix..
     * 
     * :param std::ofstream &results: The output file
     */


    floatVector Q = {-0.44956296, -0.88488713, -0.12193405,
                     -0.37866166,  0.31242661, -0.87120891,
                      0.80901699, -0.3454915 , -0.47552826};

    floatVector A = {1, 2, 3, 4, 5, 6, 7, 8, 9};
    
    floatVector rotatedA;
    errorOut ret = constitutiveTools::rotateMatrix(A, Q, rotatedA);

    if (ret){
        results << "testRotatedMatrix (test 1) & False\n";
        return 1;
    }

    if (! vectorTools::fuzzyEquals( rotatedA, {-0.09485264, -3.38815017, -5.39748037,
                                               -1.09823916,  2.23262233,  4.68884658,
                                               -1.68701666,  6.92240128, 12.8622303})){
        results << "testRotatedMatrix (test 1) & False\n";
        return 1;
    }

    //Test rotation back to original frame
    
    floatVector QT(Q.size(), 0);
    for (unsigned int i=0; i<3; i++){
        for (unsigned int j=0; j<3; j++){
            QT[3*j + i] = Q[3*i + j];
        }
    }

    floatVector App;
    ret = constitutiveTools::rotateMatrix(rotatedA, QT, App);

    if (ret){
        results << "testRotateMatrix (test 2) & False\n";
        return 1;
    }

    if (! vectorTools::fuzzyEquals(A, App)){
        results << "testRotateMatrix (test 2) & False\n";
        return 1;
    }

    results << "testRotatedMatrix & True\n";
    return 0;
}

int testComputeGreenLagrangeStrain(std::ofstream &results){
    /*!
     * Test the computation of the Green-Lagrange strain
     * 
     * :param std::ofstream &results: The output file
     */

    floatVector F = {1, 0, 0, 0, 1, 0, 0, 0, 1};
    floatVector E;

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

int testDecomposeGreenLagrangeStrain(std::ofstream &results){
    /*!
     * Test the decomposition of the Green-Lagrange strain into isochoric and 
     * volumetric parts.
     */

    floatVector F = {0.69646919, 0.28613933, 0.22685145,
                     0.55131477, 0.71946897, 0.42310646,
                     0.98076420, 0.68482974, 0.4809319};

    floatType J = vectorTools::determinant(F, 3, 3);
    floatVector Fbar = F/pow(J, 1./3);

    floatVector E, Ebar;
    errorOut ret = constitutiveTools::computeGreenLagrangeStrain(Fbar, Ebar);

    if (ret){
        ret->print();
        results << "testDecomposeGreenLagrangeStrain & False\n";
        return 1;
    }

    ret = constitutiveTools::computeGreenLagrangeStrain(F, E);

    if (ret){
        ret->print();
        results << "testDecomposeGreenLagrangeStrain & False\n";
        return 1;
    }

    floatType JOut;
    floatVector EbarOut;
    ret = constitutiveTools::decomposeGreenLagrangeStrain(E, EbarOut, JOut);

    if (ret){
        ret->print();
        results << "testDecomposeGreenLagrangeStrain & False\n";
        return 1;
    }

    if (!vectorTools::fuzzyEquals(J, JOut)){
        results << "testDecomposeGreenLagrangeStrain (test 1) & False\n";
        return 1;
    }

    if (!vectorTools::fuzzyEquals(EbarOut, Ebar)){
        results << "testDecomposeGreenLagrangeStrain (test 2) & False\n";
        return 1;
    }

    results << "testDecomposeGreenLagrangeStrain & True\n";
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
    testRotateMatrix(results);
    testComputeGreenLagrangeStrain(results);
    testDecomposeGreenLagrangeStrain(results);

    //Close the results file
    results.close();

    return 0;
}
