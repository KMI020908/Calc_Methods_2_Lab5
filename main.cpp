#include<vector>
#include"methods.cpp"
#include"ioData.cpp"
#include<iomanip>
#include"testFuncs.h"
#include<algorithm>
#include<cstdio>

// Процедура для заполнения названий файла
void getFileNames(std::size_t numOfEq, std::string &SOLUTION_FILE, std::string &DATA_FILE, const std::string &folderName){
    SOLUTION_FILE = "D:\\Calc_Methods_2\\Lab5\\" + folderName + "\\solution" + std::to_string(numOfEq) + ".txt";
    DATA_FILE = "D:\\Calc_Methods_2\\Lab5\\" + folderName + "\\data" + std::to_string(numOfEq) + ".txt";
}

// Проверка решений методом квадратур
template<typename Type>
void checkQuadMethod(std::size_t numOfEq, std::size_t numOfXIntervals, Type a, Type b, Type lambda, Type (*K)(Type, Type), Type (*f)(Type), 
void(*fillSysMatrix)(std::vector<std::vector<Type>>&, Type, Type, std::size_t, Type, Type(*)(Type, Type)), SYSTEM_FLAG sysMethod){
    std::vector<Type> solution;
    std::string SOLUTION_FILE;
    std::string DATA_FILE;
    getFileNames(numOfEq, SOLUTION_FILE, DATA_FILE, "QuadMethod");
    std::vector<Type> dataVec = {a, b, lambda};
    writeVectorFile(dataVec, DATA_FILE);
    writeScalarFile(numOfXIntervals, DATA_FILE, true);
    getSecondFredholmIntegral_QM(solution, numOfXIntervals, a, b, lambda, K, f, fillSysMatrix, sysMethod);
    writeVectorFile(solution, SOLUTION_FILE);
}

// Оценка порядка сходимости метода квадпратур
template<typename Type>
void getErrEstimateQuad(const std::string &fileName, Type (*realSol)(Type x), std::size_t numOfIt, std::size_t numOfXIntervals, Type a, Type b, Type lambda, Type (*K)(Type, Type), Type (*f)(Type), 
void(*fillSysMatrix)(std::vector<std::vector<Type>>&, Type, Type, std::size_t, Type, Type(*)(Type, Type)), SYSTEM_FLAG sysMethod){ 
    std::string SOLUTION_FILE = "D:\\Calc_Methods_2\\Lab5\\errEstQuad\\err" + fileName + ".txt";
    errEstimateQuadMethod(SOLUTION_FILE, realSol, numOfIt, numOfXIntervals, a, b, lambda, K, f, fillSysMatrix, sysMethod);
}


// Проверка решений методом простой итерации
template<typename Type>
void checkSimpleIt(std::size_t numOfEq, std::size_t numOfXIntervals, Type a, Type b, Type lambda, Type (*K)(Type, Type), Type (*f)(Type), Type (*U0)(Type),
Type quadMethod(const std::vector<Type>&, Type, Type, std::size_t, Type, Type (*)(Type, Type)), Type eps = 1e-6, std::size_t stopIt = 50000){
    std::vector<Type> solution;
    std::string SOLUTION_FILE;
    std::string DATA_FILE;
    getFileNames(numOfEq, SOLUTION_FILE, DATA_FILE, "SimpleItMethod");
    std::vector<Type> dataVec = {a, b, lambda};
    writeVectorFile(dataVec, DATA_FILE);
    writeScalarFile(numOfXIntervals, DATA_FILE, true);
    std::cout << "Количество итераций уравнения " << numOfEq << ":\t" << getSecondFredholmIntegral_SIt(solution, U0, numOfXIntervals, a, b, lambda, K, f, quadMethod, eps, stopIt) << '\n';
    writeVectorFile(solution, SOLUTION_FILE);
}

// Зависимость точности от количесвто итераций в методе простой итерации
template<typename Type>
void simpleIterationsAnalys(const std::string &fileName, Type (*realSol)(Type x), std::size_t numOfXIntervals, Type a, Type b, Type lambda, Type (*K)(Type, Type), Type (*f)(Type), Type (*U0)(Type),
Type quadMethod(const std::vector<Type>&, Type, Type, std::size_t, Type, Type (*)(Type, Type)), std::size_t it0, std::size_t itStep, std::size_t numOfErrors, Type eps = 1e-6){
    std::string SOLUTION_FILE = "D:\\Calc_Methods_2\\Lab5\\simpleIterationsAnalys\\it" + fileName + ".txt";

    std::vector<Type> numSol;
    std::vector<Type> realSolVec(numOfXIntervals + 1);
    Type h = (b - a) / numOfXIntervals;
    for (std::size_t i = 0; i < numOfXIntervals + 1; i++){
        realSolVec[i] = realSol(a + i * h);
    }
    
    std::vector<Type> errVec;
    std::size_t tempIt = it0;
    for (std::size_t i = 0; i < numOfErrors; i++){
        getSecondFredholmIntegral_SIt(numSol, U0, numOfXIntervals, a, b, lambda, K, f, quadMethod, eps, tempIt);
        errVec.push_back(normOfVector(numSol - realSolVec, INFINITY));
        tempIt += itStep;
    }
    writeVectorFile(errVec, SOLUTION_FILE);
}


// Проверка решений с вырожденным ядром
template<typename Type>
void checkDegenerateKernel(std::size_t numOfEq,  std::size_t numOfXIntervals, Type a, Type b, Type lambda, 
const std::vector<Type(*)(Type)> &phiVec, const std::vector<Type(*)(Type)> &psiVec, Type (*f)(Type), 
Type (*quadMethod)(Type, std::size_t, Type, Type(*)(Type), Type(*)(Type)), SYSTEM_FLAG sysMethod){
    std::vector<Type> solution;
    std::string SOLUTION_FILE;
    std::string DATA_FILE;
    getFileNames(numOfEq, SOLUTION_FILE, DATA_FILE, "DegKernel");
    std::vector<Type> dataVec = {a, b, lambda};
    writeVectorFile(dataVec, DATA_FILE);
    writeScalarFile(numOfXIntervals, DATA_FILE, true);
    getSecondFredholmIntegral_DegKernel(solution, numOfXIntervals, a, b, lambda, phiVec, psiVec, f, quadMethod, sysMethod);
    writeVectorFile(solution, SOLUTION_FILE);
}

// Зависимость точности от числа слагаемых в разложении ряда тейлора
template<typename Type>
std::size_t degKernelAnalys(const std::string &fileName, Type (*realSol)(Type x),  std::size_t numOfXIntervals, Type a, Type b, Type lambda, 
const std::vector<Type(*)(Type)> &allPhiVec, const std::vector<Type(*)(Type)> &allPsiVec, Type (*f)(Type), 
Type (*quadMethod)(Type, std::size_t, Type, Type(*)(Type), Type(*)(Type)), SYSTEM_FLAG sysMethod){
    std::size_t m = allPhiVec.size();
    if (m != allPsiVec.size()){
        return 0;
    }
    std::string SOLUTION_FILE = "D:\\Calc_Methods_2\\Lab5\\KernelAnalys\\ker" + fileName + ".txt";
    std::vector<Type(*)(Type)> phiVec = {allPhiVec[0]};
    std::vector<Type(*)(Type)> psiVec = {allPsiVec[0]};

    std::vector<Type> numSol;
    std::vector<Type> realSolVec(numOfXIntervals + 1);
    Type h = (b - a) / numOfXIntervals;
    for (std::size_t i = 0; i < numOfXIntervals + 1; i++){
        realSolVec[i] = realSol(a + i * h);
    }
    
    std::vector<Type> errVec;
    for (std::size_t i = 1; i < m; i++){
        phiVec.push_back(allPhiVec[i]);
        psiVec.push_back(allPsiVec[i]);
        getSecondFredholmIntegral_DegKernel(numSol, numOfXIntervals, a, b, lambda, phiVec, psiVec, f, quadMethod, sysMethod);
        errVec.push_back(normOfVector(numSol - realSolVec, INFINITY));
    }
    
    writeVectorFile(errVec, SOLUTION_FILE);
    return m;
}


// Проверка решений с сингулярным ядром 
template<typename Type>
void checkSingularIntegral(std::size_t numOfEq, Type(*f)(Type, Type), std::size_t numOfFinElems){
    std::vector<Type> solution;
    std::string SOLUTION_FILE;
    std::string DATA_FILE;
    getFileNames(numOfEq, SOLUTION_FILE, DATA_FILE, "singularIntegr");
    writeScalarFile(numOfFinElems, DATA_FILE);
    solveSingularIntegralEq(f, numOfFinElems, solution);
    writeVectorFile(solution, SOLUTION_FILE);
}

// Зависимость R от количества разбиений N
template<typename Type>
void checkRCoverage(const std::string &fileName, Type(*f)(Type, Type), std::size_t numOfFinElems0, std::size_t step, std::size_t lastNumOfFinElems){
    std::string SOLUTION_FILE = "D:\\Calc_Methods_2\\Lab5\\RAnalys\\R" + fileName + ".txt";
    std::string DATA_FILE = "D:\\Calc_Methods_2\\Lab5\\RAnalys\\dataR" + fileName + ".txt";
    std::vector<std::size_t> RDataList = {numOfFinElems0, step};
    writeVectorFile(RDataList, DATA_FILE);
    std::vector<Type> RList;
    std::vector<Type> solution;
    for (std::size_t m = numOfFinElems0; m < lastNumOfFinElems; m += step){
        RList.push_back(solveSingularIntegralEq(f, m, solution));
    }
    writeVectorFile(RList, SOLUTION_FILE);
}


// Сделать, чтобы системы решались одним методом
template<typename Type>
void makeSameQ(bool isSame, Type &var, const Type &sameVar){
    if (isSame){
        var = sameVar;
    }
}

// Сделать, чтобы системы были построены одной квадратурой
template<typename Type>
void makeSameSysQuadQ(bool isSameSysQuad, void(**fillSysMatrix)(std::vector<std::vector<Type>>&, Type, Type, std::size_t, Type, Type(*)(Type, Type)), void(**fillSameSysMatrix)(std::vector<std::vector<Type>>&, Type, Type, std::size_t, Type, Type(*)(Type, Type))){
    if (isSameSysQuad){
        *fillSysMatrix = *fillSameSysMatrix;
    }
}

// Сделать, чтобы квадратуры в методе простой итерации были одинаковые
template<typename Type>
void makeSameQuadQ(bool isSameQuad, Type (**quadMethod)(const std::vector<Type>&, Type, Type, std::size_t, Type, Type (*)(Type, Type)), Type (**sameQuadMethod)(const std::vector<Type>&, Type, Type, std::size_t, Type, Type (*)(Type, Type))){
    if (isSameQuad){
        *quadMethod = *sameQuadMethod;
    }
}

template<typename Type>
void temp_main(){
    std::size_t numOfEq;
    Type a, b, lambda, eps;
    std::size_t numOfXIntervals, stopIt, numOfIt, numOfTeylorSums, it0, itStep, numOfErrs;
    Type (*K)(Type x, Type y) = nullptr;
    std::vector<Type(*)(Type)> phiVec = {};
    std::vector<Type(*)(Type)> psiVec = {};
    Type (*f)(Type x) = nullptr;
    Type (*F)(Type x, Type y) = nullptr; 
    Type (*realSol)(Type x) = nullptr;
    void(*fillSysMatrix)(std::vector<std::vector<Type>>&, Type, Type, std::size_t, Type, Type(*)(Type, Type)) = nullptr;
    Type (*quadMethod)(const std::vector<Type>&, Type, Type, std::size_t, Type, Type (*)(Type, Type)) = nullptr;
    Type (*quadMethodMulty)(Type, std::size_t, Type, Type(*)(Type), Type(*)(Type)) = nullptr;
    SYSTEM_FLAG sysMethod;

    // Сделать все методы решения СЛАУ одинаковыми для метода квадратур
    bool isSameSysMethod = true;
    SYSTEM_FLAG sameSysMethod = GM;

    // Сделать все квадратуры в методе квадратур одинаковыми 
    bool isSameSysQuad = true;
    void(*fillSameSysMatrix)(std::vector<std::vector<Type>>&, Type, Type, std::size_t, Type, Type(*)(Type, Type)) = fillSysMatrixSimpson;

    // Сделать все квадратуры в методе простой итерации одинаковыми 
    bool isSameQuadMethodSimpleIt = true;
    Type (*sameQuadMethod)(const std::vector<Type>&, Type, Type, std::size_t, Type, Type (*)(Type, Type)) = trapezoidQuad;

    // Сделать все точности одинаковыми
    bool isSameEps = true;
    Type sameEps = 1e-6;

    // Векторы и функции для вырожденного ядра из методички
    std::vector<Type(*)(Type)> allPhiVecK1 = {
        [](Type x){return (1.0 - x) / 2.0;},
        [](Type x){return std::pow(x, 3.0) / 4.0;},
        [](Type x){return -std::pow(x, 5.0) / 48.0;}, 
        [](Type x){return std::pow(x, 7.0) / 1440.0;}, 
        [](Type x){return -std::pow(x, 9.0) / 80640.0;}, 
        [](Type x){return std::pow(x, 11.0) / 7257600.0;}
    };
    std::vector<Type(*)(Type)> allPsiVecK1 = {
        [](Type s){return 1.0;},
        [](Type s){return std::pow(s, 2.0);}, 
        [](Type s){return std::pow(s, 4.0);}, 
        [](Type s){return std::pow(s, 6.0);}, 
        [](Type s){return std::pow(s, 8.0);}, 
        [](Type s){return std::pow(s, 10.0);}
    };

    // Первый тест из методички при [a, b] = [0, 1]
    numOfEq = 1;
    a = 0.0;
    b = 1.0;
    K = K1;
    f = f1;
    lambda = 1.0;
    numOfXIntervals = 100;
    realSol = [](Type x){return 1.0;};

    sysMethod = GM;
    makeSameQ(isSameSysMethod, sysMethod, sameSysMethod);
    fillSysMatrix = fillSysMatrixTrapezoid;
    makeSameSysQuadQ(isSameSysQuad, &fillSysMatrix, &fillSameSysMatrix);
    //checkQuadMethod(numOfEq, numOfXIntervals, a, b, lambda, K, f, fillSysMatrix, sysMethod);

    eps = 1e-6;
    quadMethod = trapezoidQuad;
    makeSameQuadQ(isSameQuadMethodSimpleIt, &quadMethod, &sameQuadMethod);
    makeSameQ(isSameEps, eps, sameEps);
    //checkSimpleIt(numOfEq, numOfXIntervals, a, b, lambda, K, f, f, quadMethod, eps);
    it0 = 1;
    itStep = 1;
    numOfErrs = 15;
    //simpleIterationsAnalys(std::to_string(numOfEq), realSol, numOfXIntervals, a, b, lambda, K, f, f, quadMethod, it0, itStep, numOfErrs, eps);

    numOfTeylorSums = 5;
    phiVec.resize(numOfTeylorSums);
    psiVec.resize(numOfTeylorSums);
    for (std::size_t i = 0; i < numOfTeylorSums; i++){
        phiVec[i] = allPhiVecK1[i];
        psiVec[i] = allPsiVecK1[i];
    }
    quadMethodMulty = trapezoidQuadMulty;
    sysMethod = GM;
    //checkDegenerateKernel(numOfEq, numOfXIntervals, a, b, lambda, phiVec, psiVec, f, quadMethodMulty, sysMethod);
    //degKernelAnalys(std::to_string(1), realSol, numOfXIntervals, a, b, lambda, allPhiVecK1, allPsiVecK1, f, quadMethodMulty, sysMethod);

    // Первый тест из методички при [a, b] = [0.1, 1]
    numOfEq = 2;
    a = 0.1;
    b = 1.0;
    K = K1;
    f = f1;
    lambda = 1.0;
    numOfXIntervals = 100;

    sysMethod = GM;
    makeSameQ(isSameSysMethod, sysMethod, sameSysMethod);
    fillSysMatrix = fillSysMatrixTrapezoid;
    makeSameSysQuadQ(isSameSysQuad, &fillSysMatrix, &fillSameSysMatrix);
    //checkQuadMethod(numOfEq, numOfXIntervals, a, b, lambda, K, f, fillSysMatrix, sysMethod);

    eps = 1e-6;
    quadMethod = trapezoidQuad;
    makeSameQuadQ(isSameQuadMethodSimpleIt, &quadMethod, &sameQuadMethod);
    makeSameQ(isSameEps, eps, sameEps);
    //checkSimpleIt(numOfEq, numOfXIntervals, a, b, lambda, K, f, f, quadMethod, eps);

    numOfTeylorSums = 5;
    phiVec.resize(numOfTeylorSums);
    psiVec.resize(numOfTeylorSums);
    for (std::size_t i = 0; i < numOfTeylorSums; i++){
        phiVec[i] = allPhiVecK1[i];
        psiVec[i] = allPsiVecK1[i];
    }
    quadMethodMulty = trapezoidQuadMulty;
    sysMethod = GM;
    //checkDegenerateKernel(numOfEq, numOfXIntervals, a, b, lambda, phiVec, psiVec, f, quadMethodMulty, sysMethod);


    // Второй тест из методички при [a, b] = [0, 1]
    numOfEq = 3;
    a = 0.0;
    b = 1.0;
    K = K2;
    f = f2;
    lambda = 1.0;
    numOfXIntervals = 100;

    sysMethod = GM;
    makeSameQ(isSameSysMethod, sysMethod, sameSysMethod);
    fillSysMatrix = fillSysMatrixTrapezoid;
    makeSameSysQuadQ(isSameSysQuad, &fillSysMatrix, &fillSameSysMatrix);
    //checkQuadMethod(numOfEq, numOfXIntervals, a, b, lambda, K, f, fillSysMatrix, sysMethod);

    eps = 1e-6;
    quadMethod = trapezoidQuad;
    makeSameQuadQ(isSameQuadMethodSimpleIt, &quadMethod, &sameQuadMethod);
    makeSameQ(isSameEps, eps, sameEps);
    //checkSimpleIt(numOfEq, numOfXIntervals, a, b, lambda, K, f, f, quadMethod, eps);

    numOfTeylorSums = 5;
    phiVec.resize(numOfTeylorSums);
    psiVec.resize(numOfTeylorSums);
    for (std::size_t i = 0; i < numOfTeylorSums; i++){
        phiVec[i] = allPhiVecK1[i];
        psiVec[i] = allPsiVecK1[i];
    }
    quadMethodMulty = trapezoidQuadMulty;
    sysMethod = GM;
    //checkDegenerateKernel(numOfEq, numOfXIntervals, a, b, lambda, phiVec, psiVec, f, quadMethodMulty, sysMethod);


    // Второй тест из методички при [a, b] = [0.1, 1]
    numOfEq = 4;
    a = 0.1;
    b = 1.0;
    K = K2;
    f = f2;
    lambda = 1.0;
    numOfXIntervals = 100;

    sysMethod = GM;
    makeSameQ(isSameSysMethod, sysMethod, sameSysMethod);
    fillSysMatrix = fillSysMatrixTrapezoid;
    makeSameSysQuadQ(isSameSysQuad, &fillSysMatrix, &fillSameSysMatrix);
    //checkQuadMethod(numOfEq, numOfXIntervals, a, b, lambda, K, f, fillSysMatrix, sysMethod);

    eps = 1e-6;
    quadMethod = trapezoidQuad;
    makeSameQuadQ(isSameQuadMethodSimpleIt, &quadMethod, &sameQuadMethod);
    makeSameQ(isSameEps, eps, sameEps);
    //checkSimpleIt(numOfEq, numOfXIntervals, a, b, lambda, K, f, f, quadMethod, eps);

    numOfTeylorSums = 5;
    phiVec.resize(numOfTeylorSums);
    psiVec.resize(numOfTeylorSums);
    for (std::size_t i = 0; i < numOfTeylorSums; i++){
        phiVec[i] = allPhiVecK1[i];
        psiVec[i] = allPsiVecK1[i];
    }
    quadMethodMulty = trapezoidQuadMulty;
    sysMethod = GM;
    //checkDegenerateKernel(numOfEq, numOfXIntervals, a, b, lambda, phiVec, psiVec, f, quadMethodMulty, sysMethod);


    // Тест для проверки аналитического решения
    numOfEq = 5;
    a = 0.0;
    b = 1.0;
    K = K3;
    f = f3;
    lambda = 1.0;
    numOfXIntervals = 100;

    sysMethod = GM;
    makeSameQ(isSameSysMethod, sysMethod, sameSysMethod);
    fillSysMatrix = fillSysMatrixTrapezoid;
    makeSameSysQuadQ(isSameSysQuad, &fillSysMatrix, &fillSameSysMatrix);
    //checkQuadMethod(numOfEq, numOfXIntervals, a, b, lambda, K, f, fillSysMatrix, sysMethod);
    
    realSol = [](Type x){return std::sin(x);};
    numOfIt = 7;
    fillSysMatrix = fillSysMatrixTrapezoid;
    //getErrEstimateQuad("Trapezoid", realSol, numOfIt, numOfXIntervals, a, b, lambda, K, f, fillSysMatrix, sysMethod);
    fillSysMatrix = fillSysMatrixSimpson;
    //getErrEstimateQuad("Simpson", realSol, numOfIt, numOfXIntervals, a, b, lambda, K, f, fillSysMatrix, sysMethod);

    eps = 1e-6;
    quadMethod = trapezoidQuad;
    makeSameQuadQ(isSameQuadMethodSimpleIt, &quadMethod, &sameQuadMethod);
    makeSameQ(isSameEps, eps, sameEps);
    checkSimpleIt(numOfEq, numOfXIntervals, a, b, lambda, K, f, f, quadMethod, eps);
    it0 = 1;
    itStep = 1;
    numOfErrs = 15;
    //simpleIterationsAnalys(std::to_string(numOfEq), realSol, numOfXIntervals, a, b, lambda, K, f, f, quadMethod, it0, itStep, numOfErrs, eps);

    numOfTeylorSums = 1;
    phiVec.resize(numOfTeylorSums);
    psiVec.resize(numOfTeylorSums);
    phiVec = {[](Type x){return x;}};
    psiVec = {[](Type s){return s;}};
    quadMethodMulty = trapezoidQuadMulty;
    sysMethod = GM;
    //checkDegenerateKernel(numOfEq, numOfXIntervals, a, b, lambda, phiVec, psiVec, f, quadMethodMulty, sysMethod);


    // Второй тест для проверки аналитического решения
    numOfEq = 6;
    a = 0.0;
    b = M_PI_2;
    K = K4;
    f = f4;
    lambda = 1.0;
    numOfXIntervals = 100;

    sysMethod = GM;
    makeSameQ(isSameSysMethod, sysMethod, sameSysMethod);
    fillSysMatrix = fillSysMatrixTrapezoid;
    makeSameSysQuadQ(isSameSysQuad, &fillSysMatrix, &fillSameSysMatrix);
    //checkQuadMethod(numOfEq, numOfXIntervals, a, b, lambda, K, f, fillSysMatrix, sysMethod);
    
    realSol = [](Type x){return std::sin(x);};
    numOfIt = 7;
    fillSysMatrix = fillSysMatrixTrapezoid;
    //getErrEstimateQuad("Trapezoid2", realSol, numOfIt, numOfXIntervals, a, b, lambda, K, f, fillSysMatrix, sysMethod);
    fillSysMatrix = fillSysMatrixSimpson;
    //getErrEstimateQuad("Simpson2", realSol, numOfIt, numOfXIntervals, a, b, lambda, K, f, fillSysMatrix, sysMethod);

    eps = 1e-6;
    quadMethod = trapezoidQuad;
    makeSameQuadQ(isSameQuadMethodSimpleIt, &quadMethod, &sameQuadMethod);
    makeSameQ(isSameEps, eps, sameEps);
    //checkSimpleIt(numOfEq, numOfXIntervals, a, b, lambda, K, f, f, quadMethod, eps);

    std::vector<Type(*)(Type)> allPhiVecK4 = {
        [](Type x){return 1.0;},
        [](Type x){return x;},
        [](Type x){return std::pow(x, 2.0) / 2.0;},
        [](Type x){return std::pow(x, 3.0) / 6.0;}, 
        [](Type x){return std::pow(x, 4.0) / 24.0;}, 
        [](Type x){return std::pow(x, 5.0) / 120.0;}, 
        [](Type x){return std::pow(x, 6.0) / 720;},
        [](Type x){return std::pow(x, 7.0) / 5040;}, 
        [](Type x){return std::pow(x, 8.0) / 40320;}, 
        [](Type x){return std::pow(x, 9.0) / 362880.0;},
        [](Type x){return std::pow(x, 10.0) / 3628800.0;},
        [](Type x){return std::pow(x, 11.0) / 39916800.0;},
        [](Type x){return std::pow(x, 12.0) / 479001600.0;},
        [](Type x){return std::pow(x, 13.0) / 6227020800.0;}
    };
    std::vector<Type(*)(Type)> allPsiVecK4 = {
        [](Type s){return 1.0;},
        [](Type s){return s;},
        [](Type s){return std::pow(s, 2.0);},
        [](Type s){return std::pow(s, 3.0);}, 
        [](Type s){return std::pow(s, 4.0);}, 
        [](Type s){return std::pow(s, 5.0);}, 
        [](Type s){return std::pow(s, 6.0);},
        [](Type s){return std::pow(s, 7.0);}, 
        [](Type s){return std::pow(s, 8.0);}, 
        [](Type s){return std::pow(s, 9.0);},
        [](Type s){return std::pow(s, 10.0);},
        [](Type s){return std::pow(s, 11.0);},
        [](Type s){return std::pow(s, 12.0);},
        [](Type s){return std::pow(s, 13.0);}
    };
    numOfTeylorSums = 10;
    phiVec.resize(numOfTeylorSums);
    psiVec.resize(numOfTeylorSums);
    for (std::size_t i = 0; i < numOfTeylorSums; i++){
        phiVec[i] = allPhiVecK4[i];
        psiVec[i] = allPsiVecK4[i];
    }
    quadMethodMulty = trapezoidQuadMulty;
    sysMethod = GM;
    //checkDegenerateKernel(numOfEq, numOfXIntervals, a, b, lambda, phiVec, psiVec, f, quadMethodMulty, sysMethod);
    //degKernelAnalys(std::to_string(6), realSol, numOfXIntervals, a, b, lambda, allPhiVecK4, allPsiVecK4, f, quadMethodMulty, sysMethod);


    // Сингулярные уравнения

    // Сингулярное уравнение первый вариант
    numOfEq = 1;
    F = fSing1;
    numOfXIntervals = 200;
    //checkSingularIntegral(numOfEq, F, numOfXIntervals);
    //checkRCoverage("1", F, 15, 5, 210);

    // Сингулярное уравнение девятый вариант
    numOfEq = 9;
    F = fSing9;
    numOfXIntervals = 200;
    //checkSingularIntegral(numOfEq, F, numOfXIntervals);
    //checkRCoverage("9", F, 15, 5, 210);

    // Сингулярное уравнение девятый вариант
    numOfEq = 8;
    F = fSing8;
    numOfXIntervals = 200;
    //checkSingularIntegral(numOfEq, F, numOfXIntervals);
    //checkRCoverage("8", F, 15, 5, 210);

    // Сингулярное уравнение девятый вариант
    numOfEq = 20;
    F = fSing20;
    numOfXIntervals = 200;
    //checkSingularIntegral(numOfEq, F, numOfXIntervals);
    //checkRCoverage("20", F, 15, 5, 210);
}

int main(){
    temp_main<double>();
    return 0;
}