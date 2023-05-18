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
    std::size_t numOfXIntervals, stopIt, numOfIt, numOfTeylorSums;
    Type (*K)(Type x, Type y) = nullptr;
    Type (*f)(Type x) = nullptr;
    Type (*U0)(Type x) = nullptr;
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
    std::vector<Type(*)(Type)> phiVec = {};
    std::vector<Type(*)(Type)> psiVec = {};


    // Первый тест из методички при [a, b] = [0, 1]
    numOfEq = 1;
    a = 0.0;
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
    U0 = [](Type x){return 0.0;};
    quadMethod = trapezoidQuad;
    makeSameQuadQ(isSameQuadMethodSimpleIt, &quadMethod, &sameQuadMethod);
    makeSameQ(isSameEps, eps, sameEps);
    //checkSimpleIt(numOfEq, numOfXIntervals, a, b, lambda, K, f, U0, quadMethod, eps);

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
    U0 = [](Type x){return 0.0;};
    quadMethod = trapezoidQuad;
    makeSameQuadQ(isSameQuadMethodSimpleIt, &quadMethod, &sameQuadMethod);
    makeSameQ(isSameEps, eps, sameEps);
    //checkSimpleIt(numOfEq, numOfXIntervals, a, b, lambda, K, f, U0, quadMethod, eps);

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
    U0 = [](Type x){return 0.0;};
    quadMethod = trapezoidQuad;
    makeSameQuadQ(isSameQuadMethodSimpleIt, &quadMethod, &sameQuadMethod);
    makeSameQ(isSameEps, eps, sameEps);
    //checkSimpleIt(numOfEq, numOfXIntervals, a, b, lambda, K, f, U0, quadMethod, eps);

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
    U0 = [](Type x){return 0.0;};
    quadMethod = trapezoidQuad;
    makeSameQuadQ(isSameQuadMethodSimpleIt, &quadMethod, &sameQuadMethod);
    makeSameQ(isSameEps, eps, sameEps);
    //checkSimpleIt(numOfEq, numOfXIntervals, a, b, lambda, K, f, U0, quadMethod, eps);

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
    //U0 = [](Type x){return 0.0;};
    quadMethod = trapezoidQuad;
    makeSameQuadQ(isSameQuadMethodSimpleIt, &quadMethod, &sameQuadMethod);
    makeSameQ(isSameEps, eps, sameEps);
    //checkSimpleIt(numOfEq, numOfXIntervals, a, b, lambda, K, f, f, quadMethod, eps);

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
    //U0 = [](Type x){return 0.0;};
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
        [](Type x){return std::pow(x, 10.0) / 3628800.0;}
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
        [](Type s){return std::pow(s, 10.0);}
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
    checkDegenerateKernel(numOfEq, numOfXIntervals, a, b, lambda, phiVec, psiVec, f, quadMethodMulty, sysMethod);
}

int main(){
    temp_main<double>();
    return 0;
}