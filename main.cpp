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

// Оценка порядка сходимости при разных sigma при известном аналитическом решении
template<typename Type>
void getCoverageEstimate(std::size_t numOfEq, Type(*realSol)(Type t, Type x, Type y), std::size_t numOfIt, Type L1, Type L2, Type timeEnd,
std::size_t numOfXIntervals, std::size_t numOfYIntervals, std::size_t numOfTIntervals, CONDS_FLAG condsX, CONDS_FLAG condsY, 
Type(*U0)(Type x, Type y), Type(*T)(Type x, Type y), Type(*Q)(Type x, Type y), Type(*f)(Type t, Type x, Type y)){ 
    std::string SOLUTION_FILE = "D:\\Calc_Methods_2\\Lab4\\coverageEst\\coverage" + std::to_string(numOfEq) + ".txt";
    getRealSolEstimatePoisson2DEq(SOLUTION_FILE, realSol, numOfIt, L1, L2, timeEnd, numOfXIntervals, numOfYIntervals, numOfTIntervals, condsX, condsY, U0, T, Q, f);
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
void makeSameSysQ(bool isSameSysMethod, SYSTEM_FLAG &sysMethod, const SYSTEM_FLAG &sameSysMethod){
    if (isSameSysMethod){
        sysMethod = sameSysMethod;
    }
}

// Сделать, чтобы системы были построены одной квадратурой
template<typename Type>
void makeSameQuadQ(bool isSameSysQuad, void(**fillSysMatrix)(std::vector<std::vector<Type>>&, Type, Type, std::size_t, Type, Type(*)(Type, Type)), void(**fillSameSysMatrix)(std::vector<std::vector<Type>>&, Type, Type, std::size_t, Type, Type(*)(Type, Type))){
    if (isSameSysQuad){
        *fillSysMatrix = *fillSameSysMatrix;
    }
}

template<typename Type>
void temp_main(){
    std::size_t numOfEq;
    Type a, b, lambda;
    std::size_t numOfXIntervals;
    Type (*K)(Type x, Type y) = nullptr;
    Type (*f)(Type x) = nullptr;
    void(*fillSysMatrix)(std::vector<std::vector<Type>>&, Type, Type, std::size_t, Type, Type(*)(Type, Type)) = nullptr;
    SYSTEM_FLAG sysMethod;

    bool isSameSysMethod = true;
    SYSTEM_FLAG sameSysMethod = GM;

    bool isSameSysQuad = true;
    void(*fillSameSysMatrix)(std::vector<std::vector<Type>>&, Type, Type, std::size_t, Type, Type(*)(Type, Type)) = fillSysMatrixSimpson;

    // Первый тест из методички при [a, b] = [0, 1]
    numOfEq = 1;
    a = 0.0;
    b = 1.0;
    K = K1;
    f = f1;
    lambda = 1.0;
    numOfXIntervals = 100;
    sysMethod = GM;
    makeSameSysQ(isSameSysMethod, sysMethod, sameSysMethod);
    fillSysMatrix = fillSysMatrixTrapezoid;
    makeSameQuadQ(isSameSysQuad, &fillSysMatrix, &fillSameSysMatrix);
    checkQuadMethod(numOfEq, numOfXIntervals, a, b, lambda, K, f, fillSysMatrix, sysMethod);

    // Первый тест из методички при [a, b] = [0.1, 1]
    numOfEq = 2;
    a = 0.1;
    b = 1.0;
    K = K1;
    f = f1;
    lambda = 1.0;
    numOfXIntervals = 100;
    sysMethod = GM;
    makeSameSysQ(isSameSysMethod, sysMethod, sameSysMethod);
    fillSysMatrix = fillSysMatrixTrapezoid;
    makeSameQuadQ(isSameSysQuad, &fillSysMatrix, &fillSameSysMatrix);
    checkQuadMethod(numOfEq, numOfXIntervals, a, b, lambda, K, f, fillSysMatrix, sysMethod);

    // Второй тест из методички при [a, b] = [0, 1]
    numOfEq = 3;
    a = 0.0;
    b = 1.0;
    K = K2;
    f = f2;
    lambda = 1.0;
    numOfXIntervals = 100;
    sysMethod = GM;
    makeSameSysQ(isSameSysMethod, sysMethod, sameSysMethod);
    fillSysMatrix = fillSysMatrixTrapezoid;
    makeSameQuadQ(isSameSysQuad, &fillSysMatrix, &fillSameSysMatrix);
    checkQuadMethod(numOfEq, numOfXIntervals, a, b, lambda, K, f, fillSysMatrix, sysMethod);

    // Второй тест из методички при [a, b] = [0.1, 1]
    numOfEq = 4;
    a = 0.1;
    b = 1.0;
    K = K2;
    f = f2;
    lambda = 1.0;
    numOfXIntervals = 100;
    sysMethod = GM;
    makeSameSysQ(isSameSysMethod, sysMethod, sameSysMethod);
    fillSysMatrix = fillSysMatrixTrapezoid;
    makeSameQuadQ(isSameSysQuad, &fillSysMatrix, &fillSameSysMatrix);
    checkQuadMethod(numOfEq, numOfXIntervals, a, b, lambda, K, f, fillSysMatrix, sysMethod);
}

int main(){
    temp_main<double>();

    /*
    double(*K)(double x, double y) = [](double x, double y){return 1.0;};
    double h = 1.0;
    double lambda = 1.0;
    std::size_t nx = 6;
    double a = 0.0;
    std::vector<std::vector<double>> matrix;

    fillSysMatrixSimpson(matrix, a, h, nx, lambda, K);
    std::cout << matrix;
    */
    return 0;
}