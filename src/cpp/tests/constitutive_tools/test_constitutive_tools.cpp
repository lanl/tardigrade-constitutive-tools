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
    floatVector alphaVec = {0.1, 0.2, 0.3, 0.4};
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

    error = constitutiveTools::midpointEvolution(Dt, Ap, DApDt, DADt, A, alphaVec);

    if (error){
        error->print();
        results << "testMidpointEvolution & False\n";
        return 1;
    }

    if (!vectorTools::fuzzyEquals(A, {20.5, 23. , 25.5, 28.})){
        results << "testMidpointEvolution (test 4) & False\n";
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
        results << "testMidpointEvolution (test 5) & False\n";
        return 1;
    }

    for (unsigned int i=0; i<DADt.size(); i++){
        floatVector delta = floatVector(DADt.size(), 0);
        delta[i] = eps*(DADt[i]) + eps;

        error = constitutiveTools::midpointEvolution(Dt, Ap, DApDt, DADt+delta, Ai, alpha);

        floatVector gradCol = (Ai - A0)/delta[i];

        for (unsigned int j=0; j<gradCol.size(); j++){
            if (!vectorTools::fuzzyEquals(DADADt[j][i], gradCol[j])){
                results << "testMidpointEvolution (test 6) & False\n";
                return 1;
            }
        }
        
    }
    
    results << "testMidpointEvolution & True\n";
    return 0;
}

int testComputeDFDt(std::ofstream &results){
    /*!
     * Test the computation of the total time derivative of the 
     * deformation gradient.
     * 
     * :param std::ofstream &results: The output file
     */

    floatVector F = {0.69646919, 0.28613933, 0.22685145,
                     0.55131477, 0.71946897, 0.42310646,
                     0.98076420, 0.68482974, 0.4809319};

    floatVector L = {0.57821272, 0.27720263, 0.45555826,
                     0.82144027, 0.83961342, 0.95322334,
                     0.4768852 , 0.93771539, 0.1056616};

    floatVector answer = {1.00232848, 0.67686793, 0.46754712,
                          1.96988645, 1.49191786, 1.00002629,
                          0.95274131, 0.88347295, 0.55575157};

    floatVector DFDt;

    errorOut error = constitutiveTools::computeDFDt(L, F, DFDt);

    if (error){
        error->print();
        results << "testComputeDFDt & False\n";
        return 1;
    }

    if (!vectorTools::fuzzyEquals(DFDt, answer)){
        results << "testComputeDFDt (test 1) & False\n";
        return 1;
    }

    //Test the jacobians
    floatVector DFDtJ;
    floatMatrix dDFDtdL, dDFDtdF;
    error = constitutiveTools::computeDFDt(L, F, DFDtJ, dDFDtdL, dDFDtdF);

    if (error){
        error->print();
        results << "testComputeDFDt & False\n";
        return 1;
    }

    if (!vectorTools::fuzzyEquals(DFDt, DFDtJ)){
        results << "testComputeDFDt (test 2) & False\n";
        return 1;
    }

    //Use finite differences to estimate the jacobian
    floatType eps = 1e-6;
    for (unsigned int i=0; i<F.size(); i++){

        //Compute finite difference gradient w.r.t. L
        floatVector delta(L.size(), 0);
        delta[i] = eps*fabs(L[i]) + eps;

        error = constitutiveTools::computeDFDt(L + delta, F, DFDtJ);

        if (error){
            error->print();
            results << "testComputeDFDt & False\n";
            return 1;
        }

        floatVector gradCol = (DFDtJ - DFDt)/delta[i];

        for (unsigned int j=0; j<gradCol.size(); j++){
            if (!vectorTools::fuzzyEquals(dDFDtdL[j][i], gradCol[j])){
                results << "testComputeDFDt (test 3) & False\n";
                return 1;
            }
        }

        //Compute finite difference gradient w.r.t. F
        delta = floatVector(F.size(), 0);
        delta[i] = eps*fabs(F[i]) + eps;

        error = constitutiveTools::computeDFDt(L, F + delta, DFDtJ);

        if (error){
            error->print();
            results << "testComputeDFDt & False\n";
            return 1;
        }

        gradCol = (DFDtJ - DFDt)/delta[i];

        for (unsigned int j=0; j<gradCol.size(); j++){
            if (!vectorTools::fuzzyEquals(dDFDtdF[j][i], gradCol[j])){
                results << "testComputeDFDt (test 4) & False\n";
                return 1;
            }
        }


    }

    results << "testComputeDFDt & True\n";
    return 0;
}

int testEvolveF(std::ofstream &results){
    /*!
     * Test the evolution of the deformation gradient.
     * 
     * :param std::ofstream &results: The output file
     */
    
    floatType Dt = 2.7;

    floatVector Fp = {0.69646919, 0.28613933, 0.22685145,
                      0.55131477, 0.71946897, 0.42310646,
                      0.98076420, 0.68482974, 0.4809319};

    floatVector Lp = {0.69006282, 0.0462321 , 0.88086378,
                      0.8153887 , 0.54987134, 0.72085876, 
                      0.66559485, 0.63708462, 0.54378588};

    floatVector L = {0.57821272, 0.27720263, 0.45555826,
                     0.82144027, 0.83961342, 0.95322334,
                     0.4768852 , 0.93771539, 0.1056616};

    //Test 1 (fully explicit)
    floatVector F;
    errorOut error = constitutiveTools::evolveF(Dt, Fp, Lp, L, F, 1);

    if (error){
        error->print();
        results << "testEvolveF & False\n";
        return 1;
    }

    floatVector answer = {4.39551129, 2.53782698, 1.84614498,
                          4.81201673, 3.75047725, 2.48674399,
                          4.62070491, 3.44211354, 2.32252023};

    if (!vectorTools::fuzzyEquals(answer, F)){
        results << "testEvolveF (test 1) & False\n";
        return 1;
    }

    //Test 2 (fully implicit)
    error = constitutiveTools::evolveF(Dt, Fp, Lp, L, F, 0);

    if (error){
        error->print();
        results << "testEvolveF & False\n";
        return 1;
    }

    answer = {0.63522182, -0.1712192 , -0.00846781,
             -0.81250979, -0.19375022, -0.20193394,
             -0.36163914, -0.03662069, -0.05769288};

    if (!vectorTools::fuzzyEquals(answer, F)){
        results << "testEvolveF (test 2) & False\n";
        return 1;
    }

    //Test 3 (midpoint rule)
    error = constitutiveTools::evolveF(Dt, Fp, Lp, L, F, 0.5);

    if (error){
        error->print();
        results << "testEvolveF & False\n";
        return 1;
    }

    answer = {0.20004929, -0.4409338 , -0.18955924,
             -3.59005736, -2.17210401, -1.55661536,
             -1.88391214, -1.13150095, -0.80579654};

    if (!vectorTools::fuzzyEquals(answer, F)){
        results << "testEvolveF (test 3) & False\n";
        return 1;
    }

    //Tests 4 and 5 (jacobian)
    floatVector FJ;
    floatMatrix dFdL;
    error = constitutiveTools::evolveF(Dt, Fp, Lp, L, FJ, dFdL, 0.5);

    if (error){
        error->print(); 
        results << "testEvolveF & False\n";
        return 1;
    }

    if (!vectorTools::fuzzyEquals(F, FJ)){
        results << "testEvolveF (test 4) False\n";
        return 1;
    }

    floatType eps = 1e-6;
    for (unsigned int i=0; i<L.size(); i++){

        floatVector delta(L.size(), 0);
        delta[i] = eps*fabs(L[i]) + eps;

        error = constitutiveTools::evolveF(Dt, Fp, Lp, L + delta, FJ, dFdL, 0.5);

        if (error){
            error->print();
            results << "testEvolveF & False\n";
            return 1;
        }

        floatVector gradCol = (FJ - F)/delta[i];

        for (unsigned int j=0; j<gradCol.size(); j++){
            if (!vectorTools::fuzzyEquals(gradCol[j], dFdL[j][i], 1e-5, 1e-5)){
                results << "testEvolveF (test 5) & False\n";
                return 1;
            }
        }

    }
    
    results << "testEvolve & True\n";
    return 0;
}

int testMac(std::ofstream &results){
    /*!
     * Test the computation of the Macullay brackets.
     * 
     * :param std::ofstream &results: The output file
     */

    floatType x = 1;
    if (!vectorTools::fuzzyEquals(constitutiveTools::mac(x), x)){
        results << "testMac (test 1) & False\n";
        return 1;
    }

    x = -1;
    if (!vectorTools::fuzzyEquals(constitutiveTools::mac(x), 0.)){
        results << "testMac (test 2) & False\n";
        return 1;
    }

    floatType xJ = 2;
    floatType dmacdx;
    if (!vectorTools::fuzzyEquals(constitutiveTools::mac(xJ), constitutiveTools::mac(xJ, dmacdx))){
        results << "testMac (test 3) & False\n";
        return 1;
    }

    if (!vectorTools::fuzzyEquals(dmacdx, 1.)){
        results << "testMac (test 4) & False\n";
        return 1;
    }

    xJ = -2;
    if (!vectorTools::fuzzyEquals(constitutiveTools::mac(xJ), constitutiveTools::mac(xJ, dmacdx))){
        results << "testMac (test 5) & False\n";
        return 1;
    }

    if (!vectorTools::fuzzyEquals(dmacdx, 0.)){
        results << "testMac (test 6) & False\n";
        return 1;
    }

    results << "testMac & True\n";
    return 0;
}

int testComputeUnitNormal(std::ofstream &results){
    /*!
     * Test the computation of the unit normal.
     * 
     * :param std::ofstream &results: The output file
     */

    floatVector A = {1, 2, 3, 4, 5, 6, 7, 8, 9};
    floatVector Anorm;

    errorOut error = constitutiveTools::computeUnitNormal(A, Anorm);

    if (error){
        error->print();
        results << "testComputeUnitNormal & False\n";
        return 1;
    }

    if (!vectorTools::fuzzyEquals(vectorTools::inner(Anorm, Anorm), 1.)){
        results << "testComputeUnitNormal (test 1) & False\n";
        return 1;
    }

    //Check the jacobian
    floatVector AnormJ;
    floatMatrix dAnormdA;

    error = constitutiveTools::computeUnitNormal(A, AnormJ, dAnormdA);

    if (error){
        error->print();
        results << "testComputeUnitNormal & False\n";
        return 1;
    }

    //Check the normalized value
    if (!vectorTools::fuzzyEquals(AnormJ, Anorm)){
        results << "testComputeUnitNormal (test 2) & False\n";
        return 1;
    }

    //Check the gradient
    floatType eps = 1e-6;
    for (unsigned int i=0; i<A.size(); i++){
        floatVector delta(A.size(), 0);
        delta[i] = eps*fabs(A[i]) + eps;

        error = constitutiveTools::computeUnitNormal(A + delta, AnormJ, dAnormdA);

        if (error){
            error->print();
            results << "testComputeUnitNormal & False\n";
            return 1;
        }

        floatVector gradCol = (AnormJ - Anorm)/delta[i];

        for (unsigned int j=0; j<gradCol.size(); j++){
            if (!vectorTools::fuzzyEquals(dAnormdA[j][i], gradCol[j])){
                results << "testComputeUnitNormal (test 3) & False\n";
                return 1;
            }
        }
    }

    results << "testComputeUnitNormal & True\n";
    return 0;
}

int testPullBackVelocityGradient(std::ofstream &results){
    /*!
     * Test the pull back operation on the velocity gradient.
     * 
     * :param std::ofstream &results: The output file.
     */

    floatVector velocityGradient = {0.69006282, 0.0462321 , 0.88086378,
                                    0.8153887 , 0.54987134, 0.72085876,
                                    0.66559485, 0.63708462, 0.54378588};

    floatVector deformationGradient = {0.69646919, 0.28613933, 0.22685145,
                                       0.55131477, 0.71946897, 0.42310646,
                                       0.98076420, 0.68482974, 0.4809319};

    floatVector pullBackL;    
    floatVector expectedPullBackL = {6.32482111,   3.11877752,   2.43195977,
                                    20.19439192,  10.22175689,   7.88052809,
                                   -38.85113898, -18.79212468, -14.76285795};

    errorOut error = constitutiveTools::pullBackVelocityGradient(velocityGradient, deformationGradient, pullBackL);

    if (error){
        error->print();
        results << "testPullBackVelocityGradient & False\n";
        return 1;
    }

    if (!vectorTools::fuzzyEquals(pullBackL, expectedPullBackL)){
        results << "testPullBackVelocityGradient (test 1) & False\n";
        return 1;
    }

    floatVector pullBackLJ;
    floatMatrix dpbLdL, dpbLdF;
    
    //Test of the jacobian
    error = constitutiveTools::pullBackVelocityGradient(velocityGradient, deformationGradient, pullBackLJ, 
                                                        dpbLdL, dpbLdF);

    if (error){
        error->print();
        results << "testPullBackVelocityGradient & False\n";
        return 1;
    }

    if (!vectorTools::fuzzyEquals(pullBackL, pullBackLJ)){
        results << "testPullBackVelocityGradient (test 2) & False\n";
        return 1;
    }

    //Check dpbLdL
    floatType eps = 1e-6;
    for (unsigned int i=0; i<velocityGradient.size(); i++){
        floatVector delta(velocityGradient.size(), 0);
        delta[i] = eps*fabs(velocityGradient[i]) + eps;

        error = constitutiveTools::pullBackVelocityGradient(velocityGradient + delta, deformationGradient, pullBackLJ);

        if (error){
            error->print();
            results << "testPullBackVelocityGradient & False\n";
            return 1;
        }

        floatVector gradCol = (pullBackLJ - pullBackL)/delta[i];

        for (unsigned int j=0; j<gradCol.size(); j++){
            if (!vectorTools::fuzzyEquals(gradCol[j], dpbLdL[j][i])){
                results << "testPullBackVelocityGradient (test 3) & False\n";
                return 1;
            }
        }
    }

    //Check dpbLdF
    for (unsigned int i=0; i<deformationGradient.size(); i++){
        floatVector delta(deformationGradient.size(), 0);
        delta[i] = eps*fabs(deformationGradient[i]) + eps;

        error = constitutiveTools::pullBackVelocityGradient(velocityGradient, deformationGradient + delta, pullBackLJ);

        if (error){
            error->print();
            results << "testPullBackVelocityGradient & False\n";
            return 1;
        }

        floatVector gradCol = (pullBackLJ - pullBackL)/delta[i];

        for (unsigned int j=0; j<gradCol.size(); j++){
            if (!vectorTools::fuzzyEquals(gradCol[j], dpbLdF[j][i], 1e-4)){
                results << "testPullBackVelocityGradient (test 4) & False\n";
                return 1;
            }
        }
    }

    results << "testPullBackVelocityGradient & True\n";
    return 0;
}

int testQuadraticThermalExpansion(std::ofstream &results){
    /*!
     * Test the computation of the thermal expansion using a 
     * quadratic form.
     * 
     * :param std::ofstream &results: The output file.
     */

    floatType temperature = 283.15;
    floatType referenceTemperature = 273.15;

    floatVector linearParameters = {1, 2, 3, 4};
    floatVector quadraticParameters = {5, 6, 7, 8};

    floatVector thermalExpansion;
    errorOut error = constitutiveTools::quadraticThermalExpansion(     temperature, referenceTemperature, 
                                                                  linearParameters,  quadraticParameters, 
                                                                  thermalExpansion);

    if (error){
        error->print();
        results << "testQuadraticThermalExpansion & False\n";
        return 1;
    }

    if (!vectorTools::fuzzyEquals(thermalExpansion, {510., 620., 730., 840.})){
        results << "testQuadraticThermalExpansion (test 1) & False\n";
        return 1;
    }

    floatVector thermalExpansionJ, thermalExpansionJacobian;
    floatType eps = 1e-6;
    floatType delta = eps*temperature + eps;

    error = constitutiveTools::quadraticThermalExpansion(      temperature,     referenceTemperature,
                                                          linearParameters,      quadraticParameters,
                                                         thermalExpansionJ, thermalExpansionJacobian);

    if (error){
        error->print();
        results << "testQuadraticThermalExpansion & False\n";
        return 1;
    }

    if (!vectorTools::fuzzyEquals(thermalExpansion, thermalExpansionJ)){
        results << "testQuadraticThermalExpansion (test 2) & False\n";
        return 1;
    }

    error = constitutiveTools::quadraticThermalExpansion(      temperature + delta,     referenceTemperature,
                                                                  linearParameters,      quadraticParameters,
                                                                 thermalExpansionJ);

    if (error){
        error->print();
        results << "testQuadraticThermalExpansion & False\n";
        return 1;
    }
    
    if (!vectorTools::fuzzyEquals(thermalExpansionJacobian, (thermalExpansionJ - thermalExpansion)/delta, 1e-4)){
        results << "testQuadraticThermalExpansion (test 3) & False\n";
        return 1;
    }

    results << "testQuadraticThermalExpansion & True\n";
    return 0;
}

int testPushForwardGreenLagrangeStrain(std::ofstream &results){
    /*!
     * Test the push-forward operation on the Green-Lagrange strain.
     * 
     * :param std::ofstream &results: The output file.
     */

    floatVector deformationGradient = {0.30027935, -0.72811411,  0.26475099,
                                       1.2285819 ,  0.57663593,  1.43113814,
                                      -0.45871432,  0.2175795 ,  0.54013937};

    floatVector greenLagrangeStrain;
    errorOut error = constitutiveTools::computeGreenLagrangeStrain(deformationGradient, greenLagrangeStrain);

    if (error){
        error->print();
        results << "testPushForwardGreenLagrangeStrain & False\n";
        return 1;
    }

    floatVector almansiStrain = {-0.33393717,  0.0953188 , -0.29053383,
                                  0.0953188 ,  0.35345526,  0.11588247,
                                 -0.29053383,  0.11588247, -0.56150741};

    floatVector result;
    error = constitutiveTools::pushForwardGreenLagrangeStrain(greenLagrangeStrain, deformationGradient, 
                                                              result);

    if (error){
        error->print();
        results << "testPushForwardGreenLagrangeStrain & False\n";
        return 1;
    }

    if (!vectorTools::fuzzyEquals(result, almansiStrain)){
        results << "testPushForwardGreenLagrangeStrain (test 1) & False\n";
        return 1;
    }

    //Test the jacobian
    floatVector resultJ;
    floatMatrix dedE, dedF;
    error = constitutiveTools::pushForwardGreenLagrangeStrain(greenLagrangeStrain, deformationGradient, 
                                                              resultJ, dedE, dedF);

    if (error){
        error->print();
        results << "testPushForwardGreenLagrangeStrain & False\n";
        return 1;
    }

    if (!vectorTools::fuzzyEquals(result, resultJ)){
        results << "testPushForwardGreenLagrangeStrain (test 2) & False\n";
        return 1;
    }

    //Check dedE
    floatType eps = 1e-6;
    for (unsigned int i=0; i<greenLagrangeStrain.size(); i++){
        floatVector delta(greenLagrangeStrain.size(), 0);
        delta[i] = eps*fabs(greenLagrangeStrain[i]) + eps;

        error = constitutiveTools::pushForwardGreenLagrangeStrain(greenLagrangeStrain + delta, deformationGradient, 
                                                                  resultJ);
    
        if (error){
            error->print();
            results << "testPushForwardGreenLagrangeStrain & False\n";
            return 1;
        }

        floatVector grad = (resultJ - result)/delta[i];

        for (unsigned int j=0; j<grad.size(); j++){
            if (!vectorTools::fuzzyEquals(grad[j], dedE[j][i])){
                results << "testPushForwardGreenLagrangeStrain (test 3) & False\n";
                return 1;
            }
        }
    }

    //Check dedF
    for (unsigned int i=0; i<deformationGradient.size(); i++){
        floatVector delta(deformationGradient.size(), 0);
        delta[i] = eps*fabs(deformationGradient[i]) + eps;

        error = constitutiveTools::pushForwardGreenLagrangeStrain(greenLagrangeStrain, deformationGradient + delta, 
                                                                  resultJ);
    
        if (error){
            error->print();
            results << "testPushForwardGreenLagrangeStrain & False\n";
            return 1;
        }

        floatVector grad = (resultJ - result)/delta[i];

        for (unsigned int j=0; j<grad.size(); j++){
            if (!vectorTools::fuzzyEquals(grad[j], dedF[j][i], 1e-5)){
                results << "testPushForwardGreenLagrangeStrain (test 4) & False\n";
                return 1;
            }
        }
    }

    results << "testPushForwardGreenLagrangeStrain & True\n";
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
    testComputeDFDt(results);
    testEvolveF(results);
    testMac(results);
    testComputeUnitNormal(results);
    testPullBackVelocityGradient(results);
    testQuadraticThermalExpansion(results);
    testPushForwardGreenLagrangeStrain(results);

    //Close the results file
    results.close();

    return 0;
}
