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

    floatVector EbarOut2;
    floatType JOut2;
    floatMatrix dEbardE;
    floatVector dJdE;
    ret = constitutiveTools::decomposeGreenLagrangeStrain(E, EbarOut2, JOut2, dEbardE, dJdE);

    if (ret){
        ret->print();
        results << "testDecomposeGreenLagrangeStrain & False\n";
        return 1;
    }

    if (!vectorTools::fuzzyEquals(EbarOut, EbarOut2)){
        results << "testDecomposeGreenLagrangeStrain (test 3) & False\n";
        return 1;
    }

    if (!vectorTools::fuzzyEquals(JOut, JOut2)){
        results << "testDecomposeGreenLagrangeStrain (test 4) & False\n";
        return 1;
    }

    floatType eps = 1e-8;
    for (unsigned int i=0; i<E.size(); i++){
        floatVector delta(E.size(), 0);
        delta[i] =  fabs(eps*E[i]);
     
        ret = constitutiveTools::decomposeGreenLagrangeStrain(E + delta, EbarOut2, JOut2);

        if (ret){
            ret->print();
            results << "testDecomposeGreenLagrangeStrain (test 5) & False\n";
            return 1;
        }

        if (!vectorTools::fuzzyEquals((JOut2 - JOut)/delta[i], dJdE[i], 1e-4, 1e-4)){
            results << "testDecomposeGreenLagrangeStrain (test 6) & False\n";
            return 1;
        }
    }

    for (unsigned int i=0; i<E.size(); i++){
        floatVector delta(E.size(), 0);
        delta[i] = fabs(eps*E[i]);
     
        ret = constitutiveTools::decomposeGreenLagrangeStrain(E + delta, EbarOut2, JOut2);

        if (ret){
            ret->print();
            results << "testDecomposeGreenLagrangeStrain (test 5) & False\n";
            return 1;
        }

        floatVector gradCol = (EbarOut2 - EbarOut)/delta[i];

        for (unsigned int j=0; j<gradCol.size(); j++){
            if (!vectorTools::fuzzyEquals(gradCol[j], dEbardE[j][i], 1e-4, 1e-4)){
                results << "testDecomposeGreenLagrangeStrain (test 7) & False\n";
                return 1;
            }
        }
    }

    floatVector badE = {-1, 0, 0, 0, 1, 0, 0, 0, 1};

    ret = constitutiveTools::decomposeGreenLagrangeStrain(badE, EbarOut, JOut);

    if (!ret){
        results << "testDecomposeGreenLagrangeStrain (test 8) & False\n";
        return 1;
    }

    results << "testDecomposeGreenLagrangeStrain & True\n";
    return 0;
}

int testMapPK2toCauchy(std::ofstream &results){
    /*!
     * Test the mapping of the PK2 stress from the reference 
     * configuration to the current configuration.
     * 
     * :param std::ofstream &results: The output file
     */
    
    floatVector F = {1.96469186, -2.13860665, -2.73148546,
                     0.51314769,  2.1946897,  -0.7689354,
                     4.80764198,  1.84829739, -0.19068099};

    floatVector PK2 = {-1.07882482, -1.56821984,  2.29049707,
                       -0.61427755, -4.40322103, -1.01955745,
                        2.37995406, -3.1750827,  -3.24548244};

    floatVector cauchy;

    errorOut error = constitutiveTools::mapPK2toCauchy(PK2, F, cauchy);

    if (error){
        error->print();
        results << "testMapPK2toCauchy & False\n";
        return 1;
    }

    if (!vectorTools::fuzzyEquals(cauchy, {-2.47696057,  0.48015011, -0.28838671,
                                            0.16490963, -0.57481137, -0.92071407,
                                           -0.21450698, -1.22714923, -1.73532173})){
        results << "testMapPK2toCauchy (test 1) & False\n";
        return 1;
    }

    results << "testMapPK2toCauchy & True\n";
    return 0;
}

int testWLF(std::ofstream &results){
    /*!
     * Test the computation of the WLF function.
     * 
     * :param std::ofstream &results: The output file
     */

    floatType T = 145.;
    floatType Tr = 27.5;
    floatType C1 = 18.2;
    floatType C2 = 282.7;
    
    floatType factor, dfactordT;

    floatVector WLFParameters {Tr, C1, C2};

    constitutiveTools::WLF(T, WLFParameters, factor);

    if (!vectorTools::fuzzyEquals(factor, pow(10, -C1*(T - Tr)/(C2 + (T - Tr))))){
        results << "testWLF (test 1) & False\n";
        return 1;
    }

    floatType factor2;
    constitutiveTools::WLF(T, WLFParameters, factor2, dfactordT);
    
    if (!vectorTools::fuzzyEquals(factor, factor2)){
        results << "testWLF (test 2) & False\n";
        return 1;
    }

    floatType delta = fabs(1e-6*T);
    constitutiveTools::WLF(T + delta, WLFParameters, factor2);

    if (!vectorTools::fuzzyEquals(dfactordT, (factor2 - factor)/delta)){
        results << "testWLF (test 3) & False\n";
        return 1;
    }

    results << "testWLF & True\n";
    return 0;
}

int testComputeDGreenLagrangeStrainDF(std::ofstream &results){
    /*!
     * Test the computation of the gradient of the Green-Lagrange 
     * strain w.r.t. the deformation gradient.
     * 
     * :param std::ofstream &results: The output file.
     */

    floatVector F = {0.69646919, 0.28613933, 0.22685145,
                     0.55131477, 0.71946897, 0.42310646,
                     0.98076420, 0.68482974, 0.4809319};

    floatMatrix dEdF;

    errorOut error = constitutiveTools::computeDGreenLagrangeStrainDF(F, dEdF);

    if (error){
        error->print();
        results << "testComputeDGreenLagrangeStrainDF & False\n";
        return 1;
    }

    floatVector E, E2;
    error = constitutiveTools::computeGreenLagrangeStrain(F, E);

    if (error){
        error->print();
        results << "testComputeDGreenLagrangeStrainDF & False\n";
        return 1;
    }

    floatType eps = 1e-6;
    for (unsigned int i=0; i<F.size(); i++){
        floatVector delta(F.size(), 0);

        delta[i] = fabs(eps*F[i]);

        error = constitutiveTools::computeGreenLagrangeStrain(F + delta, E2);

        floatVector gradCol = (E2 - E)/delta[i];

        for (unsigned int j=0; j<gradCol.size(); j++){
            if (!vectorTools::fuzzyEquals(gradCol[j], dEdF[j][i])){
                results << "testComputeDGreenLagrangeStrainDF (test 1) & False\n";
                return 1;
            }
        }
    }
    results << "testComputeDGreenLagrangeStrainDF & True\n";
    return 0;
}

int testMidpointEvolution(std::ofstream &results){
    /*!
     * Test the midpoint evolution algorithm.
     * 
     * :param std::ofstream &results: The output file
     */

    floatType Dt = 2.5;
    floatVector Ap    = {9, 10, 11, 12};
    floatVector DApDt = {1, 2, 3, 4};
    floatVector DADt  = {5, 6, 7, 8};
    floatVector A;

    //Test implicit integration
    errorOut error = constitutiveTools::midpointEvolution(Dt, Ap, DApDt, DADt, A, 0);

    if (error){
        error->print();
        results << "testMidpointEvolution & False\n";
        return 1;
    }

    if (!vectorTools::fuzzyEquals(A, Ap + Dt*DADt)){
        results << "testMidpointEvolution (test 1) & False\n";
        return 1;
    }

    //Test explicit integration
    error = constitutiveTools::midpointEvolution(Dt, Ap, DApDt, DADt, A, 1);

    if (error){
        error->print();
        results << "testMidpointEvolution & False\n";
        return 1;
    }

    if (!vectorTools::fuzzyEquals(A, Ap + Dt*DApDt)){
        results << "testMidpointEvolution (test 2) & False\n";
        return 1;
    }

    //Test midpoint integration
    error = constitutiveTools::midpointEvolution(Dt, Ap, DApDt, DADt, A);

    if (error){
        error->print();
        results << "testMidpointEvolution & False\n";
        return 1;
    }

    if (!vectorTools::fuzzyEquals(A, Ap + Dt*0.5*(DApDt + DADt))){
        results << "testMidpointEvolution (test 3) & False\n";
        return 1;
    }

    //Add test for the jacobian
    floatType alpha = .37;
    floatType eps = 1e-6;
    floatVector A0, Ai;
    floatMatrix DADADt;

    error = constitutiveTools::midpointEvolution(Dt, Ap, DApDt, DADt, A, alpha);
    error = constitutiveTools::midpointEvolution(Dt, Ap, DApDt, DADt, A0, DADADt, alpha);

    if (error){
        error->print();
        results << "testMidpointEvolution & False\n";
        return 1;
    }

    if (!vectorTools::fuzzyEquals(A0, A)){
        results << "testMidpointEvolution (test 4) & False\n";
        return 1;
    }

    for (unsigned int i=0; i<DADt.size(); i++){
        floatVector delta = floatVector(DADt.size(), 0);
        delta[i] = eps*(DADt[i]) + eps;

        error = constitutiveTools::midpointEvolution(Dt, Ap, DApDt, DADt+delta, Ai, alpha);

        floatVector gradCol = (Ai - A0)/delta[i];

        for (unsigned int j=0; j<gradCol.size(); j++){
            if (!vectorTools::fuzzyEquals(DADADt[j][i], gradCol[j])){
                results << "testMidpointEvolution (test 5) & False\n";
                return 1;
            }
        }
        
    }
    
    results << "testMidpointEvolution & True\n";
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
    testComputeDGreenLagrangeStrainDF(results);
    testDecomposeGreenLagrangeStrain(results);
    testMapPK2toCauchy(results);
    testWLF(results);
    testMidpointEvolution(results);

    //Close the results file
    results.close();

    return 0;
}
