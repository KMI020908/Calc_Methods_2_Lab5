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

// Оценка порядка сходимости метода квадпратур
template<typename Type>
void getErrEstimateQuad(const std::string &fileName, Type (*realSol)(Type x), std::size_t numOfIt, std::size_t numOfXIntervals, Type a, Type b, Type lambda, Type (*K)(Type, Type), Type (*f)(Type), 
void(*fillSysMatrix)(std::vector<std::vector<Type>>&, Type, Type, std::size_t, Type, Type(*)(Type, Type)), SYSTEM_FLAG sysMethod){ 
    std::string SOLUTION_FILE = "D:\\Calc_Methods_2\\Lab5\\errEstQuad\\err" + fileName + ".txt";
    errEstimateQuadMethod(SOLUTION_FILE, realSol, numOfIt, numOfXIntervals, a, b, lambda, K, f, fillSysMatrix, sysMethod);
}
/*
// Решенние уравнения Пуассона
template<typename Type>
void getPoisson2DEquationSolution(std::size_t numOfEq, Type L1, Type L2, Type tau,
std::size_t numOfXIntervals, std::size_t numOfYIntervals, CONDS_FLAG condsX, CONDS_FLAG condsY, 
Type(*U0)(Type x, Type y), Type(*xi)(Type x, Type y), Type(*psi)(Type x, Type y), Type(*f)(Type x, Type y), Type eps){ 
    std::string SOLUTION_FILE;
    std::string DATA_FILE;
    std::string INTERVAL_FILE;
    std::vector<Type> dataVec = {L1, L2, tau, eps};
    std::vector<std::size_t> numOfIntervalsVec = {numOfXIntervals, numOfYIntervals};
    
    getFileNames(numOfEq, SOLUTION_FILE, DATA_FILE, INTERVAL_FILE);
    std::cout << "Количество итераций уравнения " << numOfEq << ":\t" << solve2DStationaryPoissonEquation(SOLUTION_FILE, L1, L2, tau, numOfXIntervals, numOfYIntervals, condsX, condsY, U0, xi, psi, f, eps) << '\n' << '\n';
    writeVectorFile(dataVec, DATA_FILE);
    writeVectorFile(numOfIntervalsVec, INTERVAL_FILE);
}

template<typename Type>
void getCoverageEstimateTau(std::size_t numOfEq, Type L1, Type L2, Type tau0, Type tauEnd, std::size_t tauIntervals,
std::size_t numOfXIntervals, std::size_t numOfYIntervals, CONDS_FLAG condsX, CONDS_FLAG condsY, 
Type(*U0)(Type x, Type y), Type(*T)(Type x, Type y), Type(*Q)(Type x, Type y), Type(*f)(Type x, Type y), Type eps){
    if (tauEnd > tau0){
        std::string SOLUTION_FILE = "D:\\Calc_Methods_2\\Lab4\\coverageTimeEst\\coverage" + std::to_string(numOfEq) + ".txt";
        std::string TAU_FILE = "D:\\Calc_Methods_2\\Lab4\\coverageTimeEst\\tau" + std::to_string(numOfEq) + ".txt";
        std::vector<Type> tauVec;
        std::vector<Type> errVec;
        Type tauStep = (tauEnd - tau0) / tauIntervals;
        for (std::size_t i = 0; i < tauIntervals; i++){
            tauVec.push_back(tau0 + i * tauStep);
            errVec.push_back(solve2DStationaryPoissonEquation(SOLUTION_FILE, L1, L2, tauVec[i], numOfXIntervals, numOfYIntervals, condsX, condsY, U0, T, Q, f, eps));
        }
        writeVectorFile(tauVec, TAU_FILE);
        writeVectorFile(errVec, SOLUTION_FILE);
    }
}
*/

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
    std::size_t numOfXIntervals, stopIt, numOfIt;
    Type (*K)(Type x, Type y) = nullptr;
    Type (*f)(Type x) = nullptr;
    Type (*U0)(Type x) = nullptr;
    Type (*realSol)(Type x) = nullptr;
    void(*fillSysMatrix)(std::vector<std::vector<Type>>&, Type, Type, std::size_t, Type, Type(*)(Type, Type)) = nullptr;
    Type (*quadMethod)(const std::vector<Type>&, Type, Type, std::size_t, Type, Type (*)(Type, Type)) = nullptr;
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


    // Первый тест из методички при [a, b] = [0, 1]
    numOfEq = 1;
    a = 0.0;
    b = 1.0;
    K = K1;
    f = f1;
    lambda = 1.0;
    numOfXIntervals = 10;

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


    // Специальный тест для проверки порядка сходимости метода квадратур
    numOfEq = 5;
    a = 0.0;
    b = 1.0;
    K = K3;
    f = f3;
    lambda = 1.0;
    numOfXIntervals = 10;

    sysMethod = GM;
    makeSameQ(isSameSysMethod, sysMethod, sameSysMethod);
    fillSysMatrix = fillSysMatrixTrapezoid;
    makeSameSysQuadQ(isSameSysQuad, &fillSysMatrix, &fillSameSysMatrix);
    //checkQuadMethod(numOfEq, numOfXIntervals, a, b, lambda, K, f, fillSysMatrix, sysMethod);
    
    realSol = [](Type x){return std::sin(x);};
    numOfIt = 7;
    fillSysMatrix = fillSysMatrixTrapezoid;
    getErrEstimateQuad("Trapezoid", realSol, numOfIt, numOfXIntervals, a, b, lambda, K, f, fillSysMatrix, sysMethod);
    fillSysMatrix = fillSysMatrixSimpson;
    getErrEstimateQuad("Simpson", realSol, numOfIt, numOfXIntervals, a, b, lambda, K, f, fillSysMatrix, sysMethod);

    eps = 1e-6;
    //U0 = [](Type x){return 0.0;};
    quadMethod = trapezoidQuad;
    makeSameQuadQ(isSameQuadMethodSimpleIt, &quadMethod, &sameQuadMethod);
    makeSameQ(isSameEps, eps, sameEps);
    //checkSimpleIt(numOfEq, numOfXIntervals, a, b, lambda, K, f, f, quadMethod, eps);
}

int main(){
    temp_main<double>();
    return 0;
}