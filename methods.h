#ifndef METHODS_H
#define METHODS_H

#include<vector>
#include<cmath>
#include"ioData.h"
#include<iostream>
#include <functional>

// Методы 5-й семестр
template<typename Type>
SOLUTION_FLAG gaussMethod(std::vector<std::vector<Type>> &lCoefs, std::vector<Type> &rCoefs, std::vector<Type> &solution, 
Type accuracy = 1e-7);

template<typename Type>
SOLUTION_FLAG gaussMethodFull(std::vector<std::vector<Type>> &lCoefs, std::vector<Type> &rCoefs, std::vector<Type> &solution, 
Type accuracy = 1e-7); // Метод Гаусса с полным выбором главного элемента

template<typename Type>
SOLUTION_FLAG qrMethod(std::vector<std::vector<Type>> &lCoefs, std::vector<Type> &rCoefs, std::vector<Type> &solution, 
Type accuracy = 1e-7);

template<typename Type>
SOLUTION_FLAG tridiagonalAlgoritm(const std::vector<std::vector<Type>> &lCoefs, const std::vector<Type> &rCoefs, std::vector<Type> &solution);

template<typename Type>
SOLUTION_FLAG tridiagonalAlgoritm(const std::vector<Type> &diag, std::vector<Type> &lDiag, std::vector<Type> &uDiag, const std::vector<Type> &rVec, std::vector<Type> &solution);

template<typename Type>
SOLUTION_FLAG subTridiagonalAlgoritm(const std::vector<Type> &diag, std::vector<Type> &lDiag, std::vector<Type> &uDiag, const std::vector<Type> &rVec, 
std::vector<Type> &solution, std::size_t startIndex = 0);

template<typename Type>
SOLUTION_FLAG uniDimTridiagonalAlgoritm(const std::vector<Type> &diag, std::vector<Type> &lDiag, std::vector<Type> &uDiag, const std::vector<Type> &rVec, 
std::vector<Type> &solution, std::size_t dim);

template<typename Type>
std::size_t simpleItMethod(const std::vector<std::vector<Type>> &lCoefs, const std::vector<Type> &rCoefs, const std::vector<Type> 
&firstVec, std::vector<Type> &solution, Type tao, Type accuracy = 1e-7, double p = 2.0, Type epsilon_0 = 1e-4, std::size_t stopIt = 100000);

template<typename Type>
std::size_t JacobiMethod(const std::vector<std::vector<Type>> &lCoefs, const std::vector<Type> &rCoefs, 
const std::vector<Type> &firstVec, std::vector<Type> &solution, Type accuracy = 1e-7, double p = 2.0, Type epsilon_0 = 1e-7, std::size_t stopIt = 100000);

template<typename Type>
std::size_t relaxationMethod(const std::vector<std::vector<Type>> &lCoefs, const std::vector<Type> &rCoefs, 
const std::vector<Type> &firstVec, std::vector<Type> &solution, Type accuracy = 1e-7, Type omega = 1.0, double p = 2.0, Type epsilon_0 = 1e-7, std::size_t stopIt = 100000);

template<typename Type>
std::size_t relaxationMethodFor3Diag(const std::vector<Type> &a, const std::vector<Type> &b, const std::vector<Type> & c, const std::vector<Type> &d, 
const std::vector<Type> &firstVec, std::vector<Type> &solution, Type accuracy  = 1e-7, Type omega = 1.0, double p = 2.0, Type epsilon_0 = 1e-7, std::size_t stopIt = 100000);

// Вспомогательные функции
template<typename Type>
QUADRATIC_FLAG findQMatrix(std::vector<std::vector<Type>> &lCoefs, std::vector<std::vector<Type>> &Q, Type accuracy = 1e-6);

template<typename Type>
QUADRATIC_FLAG findQBlock(std::vector<std::vector<Type>> &matrix, std::vector<std::vector<Type>> &Q, std::size_t rows, Type accuracy = 1e-6);

template<typename Type>
QUADRATIC_FLAG findQMatrix3Diag(std::vector<std::vector<Type>> &lCoefs, std::vector<std::vector<Type>> &Q, Type accuracy = 1e-6);

template<typename Type>
QUADRATIC_FLAG findQBlock3Diag(std::vector<std::vector<Type>> &matrix, std::vector<std::vector<Type>> &Q, std::size_t rows, Type accuracy = 1e-6);

template<typename Type>
QUADRATIC_FLAG findQMatrixHess(std::vector<std::vector<Type>> &lCoefs, std::vector<std::vector<Type>> &Q, Type accuracy = 1e-6);

template<typename Type>
QUADRATIC_FLAG findQBlockHess(std::vector<std::vector<Type>> &lCoefs, std::vector<std::vector<Type>> &Q, std::size_t rows, Type accuracy = 1e-6);

template<typename Type>
Type findResidual(const std::vector<std::vector<Type>> &lCoefs, const std::vector<Type> &rCoefs, const std::vector<Type> &solution); // Найти невязку

template<typename Type>
Type findMatrixNorm1(const std::vector<std::vector<Type>> &matrix); // Октаэдрическая норма матрицы

template<typename Type>
Type findCond_1(const std::vector<std::vector<Type>> &matrix); // Число обусловленности с октаэдоической метрикой

template<typename Type>
Type findMatrixNormInf(const std::vector<std::vector<Type>> &matrix); // Кубическая норма матрицы

template<typename Type>
Type normOfMatrix(const std::vector<std::vector<Type>> &matrix, double p = 1.0); 

template<typename Type>
Type findCond_inf(const std::vector<std::vector<Type>> &matrix); // Число обусловленности с кубической метрикой

template<typename Type>
Type norm1OfVector(const std::vector<Type> &vector); // Октаэдрическая норма вектора

template<typename Type>
Type norm2OfVector(const std::vector<Type> &vector); // Квадратичная норма вектора

template<typename Type>
Type normOfVector(const std::vector<Type> &vector, double p = 2.0);

template<typename Type>
Type normInfOfVector(const std::vector<Type> &vector); // Кубическая норма вектора

template<typename Type>
INVERTIBLE_FLAG invertMatrix(const std::vector<std::vector<Type>> &inputMatrix, std::vector<std::vector<Type>> &resMatrix,
    SOLUTION_FLAG (*method)(std::vector<std::vector<Type>>& matrix, std::vector<Type>& rVec, std::vector<Type>& sol, Type accuracy) = qrMethod); // Обратная матрица

template<typename Type>
std::size_t transposeMatrix(std::vector<std::vector<Type>> &matrix);

template<typename Type>
Type findLowerBoundOfcond1(std::vector<std::vector<Type>> &lCoefs, std::vector<Type> &rCoefs, Type delta1 = 1e-1, Type delta2 = 1e-3, Type delta3 = 1e-6);

template<typename Type>
Type findLowerBoundOfcondInf(std::vector<std::vector<Type>> &lCoefs, std::vector<Type> &rCoefs, Type delta1 = 1e-1, Type delta2 = 1e-3, Type delta3 = 1e-6);

template<typename Type>
MULTIPLIED_FLAG multiplyMatrix(std::vector<std::vector<Type>> &matrix1, const std::vector<std::vector<Type>> &matrix2, std::vector<std::vector<Type>> &result);

template<typename Type>
MULTIPLIED_FLAG multiplyMatrix(const std::vector<std::vector<Type>> &matrix, const std::vector<Type> &vec, std::vector<Type> &result);

// Перегрузки
template<typename Type>
std::ostream& operator<<(std::ostream &os, const std::vector<std::vector<Type>> &matrix);

template<typename Type>
std::ostream& operator<<(std::ostream &os, const std::vector<Type> &vector);

template<typename Type>
std::vector<Type> operator+(const std::vector<Type>& vec1, const std::vector<Type>& vec2);

template<typename Type>
std::vector<Type> operator-(const std::vector<Type>& vec1, const std::vector<Type>& vec2);

template<typename Type>
std::vector<Type> operator*(const std::vector<std::vector<Type>> &matrix, const std::vector<Type> &vec);

template<typename Type>
std::vector<Type> operator*(Type num, const std::vector<Type> &vec);

template<typename Type>
Type dot(const std::vector<Type> &v1, const std::vector<Type> &v2);

template<typename Type>
QUADRATIC_FLAG findCanonicalFormSimpleIt(const std::vector<std::vector<Type>> &lCoefs, 
const std::vector<Type> &rCoefs, std::vector<std::vector<Type>> &C, std::vector<Type> &y, Type tao);

template<typename Type>
QUADRATIC_FLAG findCanonicalFormJacobi(const std::vector<std::vector<Type>> &lCoefs, 
const std::vector<Type> &rCoefs, std::vector<std::vector<Type>> &C, std::vector<Type> &y);

template<typename Type>
QUADRATIC_FLAG findCanonicalFormRelaxation(const std::vector<std::vector<Type>> &lCoefs, 
const std::vector<Type> &rCoefs, std::vector<std::vector<Type>> &C, std::vector<Type> &y, Type omega = 1.0);

template<typename Type>
QUADRATIC_FLAG findCanonicalFormRelaxation2(const std::vector<std::vector<Type>> &lCoefs, 
const std::vector<Type> &rCoefs, std::vector<std::vector<Type>> &C, std::vector<Type> &y, Type omega = 1.0);

template<typename Type>
Type findLowerBoundOfIterations(const std::vector<std::vector<Type>> &lCoefs, 
const std::vector<Type> &rCoefs, const std::vector<Type> &firstVec, Type accuracy, ITERATION_METHOD_FLAG method, 
Type tao, Type omega = 1.0, double p = 2.0);

template<typename Type>
std::size_t findExactItersSimpleItMethod(const std::vector<std::vector<Type>> &lCoefs, const std::vector<Type> &rCoefs, 
const std::vector<Type> &firstVec, std::vector<Type> &solution, std::vector<Type> &rightSolution, Type tao, Type accuracy, double p, std::size_t stopIt = SIZE_MAX);

template<typename Type>
std::size_t findExactItersJacobiMethod(const std::vector<std::vector<Type>> &lCoefs, const std::vector<Type> &rCoefs, 
const std::vector<Type> &firstVec, std::vector<Type> &solution, std::vector<Type> &rightSolution, Type accuracy, double p, std::size_t stopIt = SIZE_MAX);

template<typename Type>
std::size_t findExactRelaxationMethod(const std::vector<std::vector<Type>> &lCoefs, const std::vector<Type> &rCoefs, 
const std::vector<Type> &firstVec, std::vector<Type> &solution, std::vector<Type> &rightSolution, Type accuracy, Type omega, double p, std::size_t stopIt = SIZE_MAX);

template<typename Type>
Type findNormOfErrAfterEstIt_SIT(const std::vector<std::vector<Type>> &lCoefs, const std::vector<Type> &rCoefs, 
const std::vector<Type> &firstVec, std::vector<Type> &solution, std::vector<Type> &rightSolution, Type bound, Type tao, Type accuracy, double p, std::size_t stopIt = SIZE_MAX);

template<typename Type>
Type findNormOfErrAfterEstIt_JAC(const std::vector<std::vector<Type>> &lCoefs, const std::vector<Type> &rCoefs, 
const std::vector<Type> &firstVec, std::vector<Type> &solution, std::vector<Type> &rightSolution, Type bound, Type accuracy, double p, std::size_t stopIt = SIZE_MAX);

template<typename Type>
Type findNormOfErrAfterEstIt_REL(const std::vector<std::vector<Type>> &lCoefs, const std::vector<Type> &rCoefs, 
const std::vector<Type> &firstVec, std::vector<Type> &solution, std::vector<Type> &rightSolution, Type bound, Type tao, Type accuracy, double p, std::size_t stopIt = SIZE_MAX);

template <typename Type>
std::size_t findEigenNumsQRMethod(std::vector<std::vector<Type>> &matrix, std::vector<Type> &eigenList, Type accuracy = 1e-6, bool hasShift = false, bool is3Diag = false);

template <typename Type>
QUADRATIC_FLAG getHessenbergMatrix(std::vector<std::vector<Type>> &matrix, Type accuracy = 1e-6, bool isSymmetric = false);

template<typename Type>
std::size_t findEigenNumsQRMethodHessenberg(std::vector<std::vector<Type>> &matrix, std::vector<Type> &eigenList, Type accuracy = 1e-6, bool hasShift = true, bool isSymmetric = false);

template<typename Type>
std::size_t invertItersMethod(const std::vector<std::vector<Type>> &matrix, std::vector<std::vector<Type>> &eigenMatrix, const std::vector<Type> &startEigenList,
Type accuracy = 1e-6, bool is3Diag = false);

template<typename Type>
Type invertItersMethodRayleigh(const std::vector<std::vector<Type>> &matrix, std::vector<Type> &startVec, 
std::vector<Type> &eigenVec, Type accuracy = 1e-6, bool is3Diag = false);

template<typename Type>
FILE_FLAG writeRayleighSwPool(const std::vector<std::vector<Type>> &matrix, Type step, const std::string& OUT_FILE_PATH, 
Type accuracy = 1e-6, bool is3Diag = false);

// Лаба 3
template<typename Type>
std::size_t getUniformGrid(std::vector<Type> &xGrid, std::vector<Type> &fGrid, 
Type (*f)(Type x), Type firstX, Type lastX, std::size_t numOfFinitElems);

template<typename Type>
std::size_t getChebyshevGrid(std::vector<Type> &xGrid, std::vector<Type> &fGrid, 
Type (*f)(Type x), Type firstX, Type lastX, std::size_t numOfFinitElems);

template<typename Type>
Type c(Type x, const std::vector<Type> &xGrid, std::size_t k);

template<typename Type>
Type LagrangePolynom(Type x, const std::vector<Type> &xGrid, const std::vector<Type> &fGrid);

template<typename Type>
SOLUTION_FLAG getLagrangeInterpolation(const std::vector<Type> &xVec, std::vector<Type> &fVec, const std::vector<Type> &xGrid, const std::vector<Type> &fGrid);

template<typename Type>
SOLUTION_FLAG findSplineCoefs(const std::vector<Type> &xGrid, const std::vector<Type> &fGrid, 
std::vector<Type> &a, std::vector<Type> &b, std::vector<Type> &c, std::vector<Type> &d);

template<typename Type>
SOLUTION_FLAG getCubeSplineInterpolation(const std::vector<Type> &xVec, std::vector<Type> &fVec, const std::vector<Type> &xGrid, const std::vector<Type> &fGrid, Type accuracy = 1e-10);

template<typename Type>
std::size_t getInterpolationErrorsLagrangeUniform(Type (*f)(Type x), Type firstX, Type lastX, std::size_t numOfErr,
std::vector<Type> &uniError, Type accuracy = 1e-10);

template<typename Type>
std::size_t getInterpolationErrorsLagrangeChebyshev(Type (*f)(Type x), Type firstX, Type lastX, std::size_t numOfErr,
std::vector<Type> &chebError, Type accuracy = 1e-10);

template<typename Type>
std::size_t getInterpolationErrorsSplineUniform(Type (*f)(Type x), Type firstX, Type lastX, std::size_t numOfErr,
std::vector<Type> &uniError, Type accuracy = 1e-10);

template<typename Type>
std::size_t getInterpolationErrorsSplineChebyshev(Type (*f)(Type x), Type firstX, Type lastX, std::size_t numOfErr,
std::vector<Type> &chebError, Type accuracy = 1e-10);

template<typename Type>
std::size_t getSpeedEstimateInPoint(Type (*f)(Type x), Type firstX, Type lastX, Type xi, std::size_t numOfFinEl0, 
std::vector<Type> &stepVec, std::vector<Type> &errResult, std::vector<Type> &speedResult, std::size_t stopIt = 2, Type accuracy = 1e-10);

template<typename Type>
std::size_t getSpeedEstimate(Type (*f)(Type x), Type firstX, Type lastX, std::size_t numOfFinEl0, std::vector<Type> &xiVec,
std::vector<Type> &err1, std::vector<Type> &err2, std::vector<Type> &speedResult, Type accuracy = 1e-10);

// Лаба 5
template<typename Type>
Type getEquationSolutionBissection(Type (*f)(Type x), Type firstX, Type lastX, Type accuracy = 1e-6, std::size_t stopIteration = 10000);

template<typename Type>
std::size_t getIterationsBissection(Type (*f)(Type x), Type firstX, Type lastX, Type accuracy = 1e-6, std::size_t stopIteration =10000);

template<typename Type>
Type diff(Type (*f)(Type x), Type x, Type h = 1e-4);

template<typename Type>
Type partialDiff2D(Type (*f)(Type x, Type y), Type x, Type y, VARIABLE_FLAG flag, Type h = 1e-4);

template<typename Type>
void getJacobiMatrix2D(std::vector<std::vector<Type>> &matrix, Type (*f1)(Type x, Type y), Type (*f2)(Type x, Type y), Type x, Type y, Type h = 1e-4);

template<typename Type>
Type getEquationSolutionNewthon(Type (*f)(Type x), Type firstX, Type h = 1e-4, Type accuracy = 1e-6, std::size_t stopIteration = 10000);

template<typename Type>
Type getEquationSolutionNewthonModified(Type (*f)(Type x), Type firstX, Type h = 1e-4, Type accuracy = 1e-6, std::size_t stopIteration = 10000);

template<typename Type>
std::size_t getIterationsNewthon(Type (*f)(Type x), Type firstX, Type h = 1e-4, Type accuracy = 1e-6, std::size_t stopIteration = 10000);

template<typename Type>
std::size_t getSystemSolutionNewthon(std::vector<Type> &solution, Type (*f1)(Type x, Type y), Type (*f2)(Type x, Type y),
Type firstX, Type firstY, Type accuracy = 1e-6, Type h = 1e-4, Type p = 2.0, std::size_t stopIteration = 10000);

template<typename Type>
std::size_t getSystemSolutionNewthonAnalytic(std::vector<Type> &solution, std::vector<Type> (*getJacobiMatrixElems)(Type x, Type y), 
Type (*f1)(Type x, Type y), Type (*f2)(Type x, Type y), Type firstX, Type firstY, Type accuracy, Type p, std::size_t stopIteration = 10000);

template<typename Type>
std::size_t locoliseRoots(Type (*f)(Type x), Type firstX, Type lastX, std::size_t numOfElems, 
std::vector<std::vector<Type>> &segmentMatrix);

template<typename Type>
FILE_FLAG writeNewthonSwPool(Type (*reF)(Type x, Type y), Type (*imF)(Type x, Type y), Type R, std::size_t n, 
Type h, Type accuracy, const std::string &OUT_FILE_PATH, std::size_t stopIteration = 10000);

template<typename Type>
Type getConvergEstimateNewthon(Type (*f)(Type x), Type x_0, Type realX, Type h = 1e-4, Type accuracy = 1e-6, std::size_t stopIteration = 10000);

// Методы 6-й семестр

// Лаб 1
template<typename Type>
Type getUniformGrid(Type a, Type b, std::size_t numOfFinElems, std::vector<Type> &xGrid);

template<typename Type>
Type partialDiff(Type (*f)(Type, std::vector<Type>&), std::size_t varPosition, Type t, const std::vector<Type> &x, Type h = 1e-4);

template<typename Type>
Type partialDiff(std::vector<Type> (*fSys)(Type t, const std::vector<Type>& x), std::size_t eqPosition, std::size_t varPosition, Type t, const std::vector<Type> &x, Type h = 1e-4);

template<typename Type>
std::size_t forwardEulerMethod(std::vector<Type>(*f)(Type t, const std::vector<Type> &U), Type t0, Type T, const std::vector<Type> &U0, std::size_t numOfTimeInterv,
std::vector<std::vector<Type>> &solution);

template<typename Type>
std::size_t backwardEulerMethod(std::vector<Type>(*f)(Type t, const std::vector<Type> &U), Type t0, Type T, const std::vector<Type> &U0, std::size_t numOfTimeInterv,
std::vector<std::vector<Type>> &solution, Type h = 1e-4, Type eps = 1e-6, std::size_t iterParam = 1);

template<typename Type>
std::size_t symmetricScheme(std::vector<Type>(*f)(Type t, const std::vector<Type> &U), Type t0, Type T, const std::vector<Type> &U0, std::size_t numOfTimeInterv,
std::vector<std::vector<Type>> &solution, Type h = 1e-4, Type eps = 1e-6, std::size_t iterParam = 1);

template<typename Type>
std::size_t RungeKuttaMethod2(std::vector<Type>(*f)(Type t, const std::vector<Type> &U), Type t0, Type T, const std::vector<Type> &U0, std::size_t numOfTimeInterv,
std::vector<std::vector<Type>> &solution, bool autoStep = true, Type eps = 1e-4, Type lowEps = 1e-6);

template<typename Type>
std::size_t RungeKuttaMethod4(std::vector<Type>(*f)(Type t, const std::vector<Type> &U), Type t0, Type T, const std::vector<Type> &U0, std::size_t numOfTimeInterv,
std::vector<std::vector<Type>> &solution, bool autoStep = true, Type eps = 1e-6, Type lowEps = 1e-8);

template<typename Type>
std::size_t iterationOfRungeKutta4(std::vector<Type>(*f)(Type t, const std::vector<Type> &U), Type t, Type tau,
const std::vector<Type> &y0, std::vector<Type> &y);

template<typename Type>
std::size_t AdamsMethod(std::vector<Type>(*f)(Type t, const std::vector<Type> &U), Type t0, Type T, const std::vector<Type> &U0, std::size_t numOfTimeInterv,
std::vector<std::vector<Type>> &solution);

template<typename Type>
std::size_t predicCorrect(std::vector<Type>(*f)(Type t, const std::vector<Type> &U), Type t0, Type T, const std::vector<Type> &U0, std::size_t numOfTimeInterv,
std::vector<std::vector<Type>> &solution);

template<typename Type>
std::size_t getPhaseTraces(std::vector<Type>(*f)(Type t, const std::vector<Type> &U), Type t0, Type T, std::size_t numOfTimeInterv, 
DIFF_METHOD_FLAG flag, Type L, std::size_t N, std::vector<std::vector<Type>> &dataMatrix, Type h = 1e-4, Type eps = 1e-6, std::size_t iterParam = 1);

template<typename Type>
Type getSpeedEstimateDiffSystem(std::vector<Type>(*f)(Type t, const std::vector<Type> &U), Type t0, Type T, const std::vector<Type> &U0, std::size_t numOfTimeInterv, 
DIFF_METHOD_FLAG flag, std::vector<Type> &speedResult, Type h = 1e-4, Type eps = 1e-6, std::size_t iterParam = 1);

template<typename Type>
Type getSpeedEstimateDiffSystem(std::vector<Type>(*f)(Type t, const std::vector<Type> &U), Type(*realSolution)(Type t), Type t0, Type T, const std::vector<Type> &U0, std::size_t numOfTimeInterv, 
DIFF_METHOD_FLAG flag, std::vector<Type> &speedResult, Type h = 1e-4, Type eps = 1e-6, std::size_t iterParam = 1);

template<typename Type>
std::size_t RungeKuttaMethodStepAnalys2(std::vector<Type>(*f)(Type t, const std::vector<Type> &U), Type t0, Type T, const std::vector<Type> &U0, std::size_t numOfTimeInterv,
std::vector<std::vector<Type>> &dataMatrix, Type eps = 1e-6, Type lowEps = 1e-8);

template<typename Type>
std::size_t RungeKuttaMethodStepAnalys4(std::vector<Type>(*f)(Type t, const std::vector<Type> &U), Type t0, Type T, const std::vector<Type> &U0, std::size_t numOfTimeInterv,
std::vector<std::vector<Type>> &dataMatrix, Type eps = 1e-6, Type lowEps = 1e-8);

// Лаб 2
template<typename Type>
Type aCoef(Type(*K)(Type x), Type h, std::size_t i);

template<typename Type>
FILE_FLAG solveHeatEquation(const std::string &solutionFile, Type rho, Type c, Type(*K)(Type x), Type L, Type timeEnd,
std::size_t numOfXIntervals, std::size_t numOfTimeIntervals, Type sigma, CONDS_FLAG flag, Type(*T0)(Type x), Type(*bound1)(Type t), Type(*bound2)(Type t));

template<typename Type>
FILE_FLAG solveHeatQuasilinearEquation(const std::string &solutionFile, Type rho, Type c, Type alpha, Type beta, Type gamma, Type L, Type timeEnd,
std::size_t numOfXIntervals, std::size_t numOfTimeIntervals, CONDS_FLAG flag, Type(*T0)(Type x), Type(*bound1)(Type t), Type(*bound2)(Type t), std::size_t numOfIters = 4);

template<typename Type>
FILE_FLAG getSpeedEstimateHeatEq(const std::string &speedFile, Type rho, Type c, Type(*K)(Type x), Type L, Type timeEnd,
std::size_t numOfXIntervals, std::size_t numOfTimeIntervals, Type sigma, CONDS_FLAG flag, Type(*T0)(Type x), Type(*bound1)(Type t), Type(*bound2)(Type t));

template<typename Type>
FILE_FLAG getSpeedEstimateQuasilinearHeatEq(const std::string &solutionFile, Type rho, Type c, Type alpha, Type beta, Type gamma, Type L, Type timeEnd,
std::size_t numOfXIntervals, std::size_t numOfTimeIntervals, CONDS_FLAG flag, Type(*T0)(Type x), Type(*bound1)(Type t), Type(*bound2)(Type t), std::size_t numOfIters = 4);

template<typename Type>
FILE_FLAG getRealSpeedEstimateHeatEq(const std::string &speedFile, Type (*realSol)(Type t, Type x), Type rho, Type c, Type(*K)(Type x), Type L, Type timeEnd,
std::size_t numOfXIntervals, std::size_t numOfTimeIntervals, Type sigma, CONDS_FLAG flag, Type(*T0)(Type x), Type(*bound1)(Type t), Type(*bound2)(Type t));

// Лаб 3
template<typename Type>
FILE_FLAG solveWaveEquation(const std::string &solutionFile, Type a, Type L, Type timeEnd,
std::size_t numOfXIntervals, std::size_t numOfTimeIntervals, Type(*U0)(Type x), Type(*Ut0)(Type x), Type(*bound1)(Type t), Type(*bound2)(Type t), Type x0 = 0.0);

template<typename Type>
FILE_FLAG getRealSpeedEstimateWaveEq(const std::string &speedFile, Type (*realSol)(Type t, Type x), Type a, Type L, Type timeEnd,
std::size_t numOfXIntervals, std::size_t numOfTimeIntervals, Type(*U0)(Type x), Type(*Ut0)(Type x), Type(*bound1)(Type t), Type(*bound2)(Type t));

// Лаб 4
template<typename Type>
Type maxVecElem(std::vector<Type> vec);

template<typename Type>
Type normC2Ddiff(const std::vector<std::vector<Type>> &m1, const std::vector<std::vector<Type>> &m2);

// Сеточная вторая частная производная по X для РЕГУЛЯРНЫХ точек
template<typename Type>
Type secondPartialDiffXRegular(const std::vector<std::vector<Type>> &fMatrix, std::size_t i, std::size_t j, Type hX, Type hY, Type(*boundFlux)(Type, Type));

// Сеточная вторая частная производная по X ВПЕРЕД для НЕРЕГУЛЯРНЫХ точек
template<typename Type>
Type secondPartialDiffXForward(const std::vector<std::vector<Type>> &fMatrix, std::size_t i, std::size_t j, Type hX, Type hY, Type(*boundFlux)(Type, Type));

// Сеточная вторая частная производная по X НАЗАД для НЕРЕГУЛЯРНЫХ точек
template<typename Type>
Type secondPartialDiffXBackward(const std::vector<std::vector<Type>> &fMatrix, std::size_t i, std::size_t j, Type hX, Type hY, Type(*boundFlux)(Type, Type));

// Сеточная вторая частная производная по Y для РЕГУЛЯРНЫХ точек
template<typename Type>
Type secondPartialDiffYRegular(const std::vector<std::vector<Type>> &fMatrix, std::size_t i, std::size_t j, Type hX, Type hY, Type(*boundFlux)(Type, Type));

// Сеточная вторая частная производная по Y ВПЕРЕД для НЕРЕГУЛЯРНЫХ точек
template<typename Type>
Type secondPartialDiffYForward(const std::vector<std::vector<Type>> &fMatrix, std::size_t i, std::size_t j, Type hX, Type hY, Type(*boundFlux)(Type, Type));

// Сеточная вторая частная производная по Y НАЗАД для НЕРЕГУЛЯРНЫХ точек
template<typename Type>
Type secondPartialDiffYBackward(const std::vector<std::vector<Type>> &fMatrix, std::size_t i, std::size_t j, Type hX, Type hY, Type(*boundFlux)(Type, Type));

// Заполняем элементы матрицы на границе, если задана температура
template<typename Type>
void fillBoundMatrixElems(std::vector<std::vector<Type>> &matrix, Type (*T)(Type, Type), const std::vector<BOUND_FLAG> &condsX, Type h1, const std::vector<BOUND_FLAG> &condsY, Type h2);

// Нахождение температуры вдоль Ox при f = f(x, y)
template<typename Type>
void getXTemperature(std::vector<Type> &solution, std::vector<Type> &lowDiag, std::vector<Type> &mainDiag, std::vector<Type> &upDiag, std::vector<Type> &fVec,
const std::vector<std::vector<Type>> &matrix, std::size_t j, Type hX, Type hY, Type tau, const std::vector<BOUND_FLAG> &condsY, 
Type(*T)(Type, Type), Type(*Q)(Type, Type), Type(*f)(Type, Type), Type(*secondPartialDiffY)(const std::vector<std::vector<Type>>&, std::size_t, std::size_t, Type, Type, Type(*)(Type, Type)));

// Нахождение температуры вдоль Oy при f = f(x, y)
template<typename Type>
void getYTemperature(std::vector<Type> &solution, std::vector<Type> &lowDiag, std::vector<Type> &mainDiag, std::vector<Type> &upDiag, std::vector<Type> &fVec,
const std::vector<std::vector<Type>> &matrix, std::size_t i, Type hX, Type hY, Type tau, const std::vector<BOUND_FLAG> &condsX, 
Type(*T)(Type, Type), Type(*Q)(Type, Type), Type(*f)(Type, Type), Type(*secondPartialDiffX)(const std::vector<std::vector<Type>>&, std::size_t, std::size_t, Type, Type, Type(*)(Type, Type)));

void fillCondsVec(std::vector<BOUND_FLAG> &condsVec, CONDS_FLAG conds);

template<typename Type>
std::size_t solve2DStationaryPoissonEquation(const std::string &solutionFile, Type L1, Type L2, Type tau, std::size_t numOfXIntervals, std::size_t numOfYIntervals, 
CONDS_FLAG condsX, CONDS_FLAG condsY, Type(*U0)(Type, Type), Type (*T)(Type, Type), Type (*Q)(Type, Type), Type(*f)(Type, Type), Type eps);

// Нахождение температуры вдоль Ox при f = f(t, x, y)
template<typename Type>
void getXTemperature(std::vector<Type> &solution, std::vector<Type> &lowDiag, std::vector<Type> &mainDiag, std::vector<Type> &upDiag, std::vector<Type> &fVec,
const std::vector<std::vector<Type>> &matrix, std::size_t j, std::size_t k, Type hX, Type hY, Type tau, const std::vector<BOUND_FLAG> &condsY, 
Type(*T)(Type, Type), Type(*Q)(Type, Type), Type(*f)(Type, Type, Type), Type(*secondPartialDiffY)(const std::vector<std::vector<Type>>&, std::size_t, std::size_t, Type, Type, Type(*)(Type, Type)));

// Нахождение температуры вдоль Oy при f = f(t, x, y)
template<typename Type>
void getYTemperature(std::vector<Type> &solution, std::vector<Type> &lowDiag, std::vector<Type> &mainDiag, std::vector<Type> &upDiag, std::vector<Type> &fVec,
const std::vector<std::vector<Type>> &matrix, std::size_t i, std::size_t k, Type hX, Type hY, Type tau, const std::vector<BOUND_FLAG> &condsX, 
Type(*T)(Type, Type), Type(*Q)(Type, Type), Type(*f)(Type, Type, Type), Type(*secondPartialDiffX)(const std::vector<std::vector<Type>>&, std::size_t, std::size_t, Type, Type, Type(*)(Type, Type)));

// Нахождение нормы невязки двумерного уравнения теплопроводности
template<typename Type>
Type get2DHeatEqNormOfResidual(Type (*realSol)(Type t, Type x, Type y), Type L1, Type L2, Type timeEnd, std::size_t numOfXIntervals, std::size_t numOfYIntervals, std::size_t numOfTIntervals,
CONDS_FLAG condsX, CONDS_FLAG condsY, Type(*U0)(Type, Type), Type (*T)(Type, Type), Type (*Q)(Type, Type), Type(*f)(Type, Type, Type));

// Нахождение норм при уменьшении шагов
template<typename Type>
FILE_FLAG getRealSolEstimatePoisson2DEq(const std::string &speedFile, Type (*realSol)(Type t, Type x, Type y), std::size_t numOfIt, Type L1, Type L2, Type timeEnd, std::size_t numOfXIntervals, std::size_t numOfYIntervals, std::size_t numOfTIntervals,
CONDS_FLAG condsX, CONDS_FLAG condsY, Type(*U0)(Type, Type), Type (*T)(Type, Type), Type (*Q)(Type, Type), Type(*f)(Type, Type, Type));

// Лаб 5
// Матрица соотвествующая методу трапеций, O(h^2)
template<typename Type>
void fillSysMatrixTrapezoid(std::vector<std::vector<Type>> &sysMatrix, Type a, Type h, std::size_t numOfXIntervals, Type lambda, Type (*K)(Type, Type));

// Матрица, соответсвующая методу Симпсона, O(h^4) 
template<typename Type>
void fillSysMatrixSimpson(std::vector<std::vector<Type>> &sysMatrix, Type a, Type h, std::size_t numOfXIntervals, Type lambda, Type (*K)(Type, Type));

// Решение интегрального уравнения методом квадратур
template<typename Type>
Type getSecondFredholmIntegral_QM(std::vector<Type> &solution, std::size_t numOfXIntervals, Type a, Type b, Type lambda, Type (*K)(Type, Type), Type (*f)(Type), 
void(*fillSysMatrix)(std::vector<std::vector<Type>>&, Type, Type, std::size_t, Type, Type(*)(Type, Type)), SYSTEM_FLAG sysMethod);

// Оценка порядка точности
template<typename Type>
FILE_FLAG errEstimateQuadMethod(const std::string &speedFile, Type (*realSol)(Type x), std::size_t numOfIt, std::size_t numOfXIntervals, Type a, Type b, Type lambda, Type (*K)(Type, Type), Type (*f)(Type), 
void(*fillSysMatrix)(std::vector<std::vector<Type>>&, Type, Type, std::size_t, Type, Type(*)(Type, Type)), SYSTEM_FLAG sysMethod);


// Квадратурные формулы для вычислений интеграла в уравнении Фредгольма в точке x
template<typename Type>
Type trapezoidQuad(const std::vector<Type> &UVec, Type x, Type a, std::size_t numOfXIntervals, Type h, Type (*K)(Type, Type));

// Решние интегрального уравнения методом простой итерации
template<typename Type>
std::size_t getSecondFredholmIntegral_SIt(std::vector<Type> &solution, Type (*U0)(Type), std::size_t numOfXIntervals, Type a, Type b, Type lambda, Type (*K)(Type, Type), Type (*f)(Type), 
Type (*quadMethod)(const std::vector<Type>&, Type, Type, std::size_t, Type, Type (*)(Type, Type)), Type eps = 1e-6, std::size_t stopIt = 50000);


// Квадратурные формулы для вычислений интеграла произведения функций в уравнении Фредгольма с вырожденным ядром
template<typename Type>
Type trapezoidQuadMulty(Type a, std::size_t numOfXIntervals, Type h, Type(*f1)(Type), Type(*f2)(Type));

// Решение интегрального уравнения Фредгольма с вырожденным ядром, т.е. K(x, s) = Sum[phi_k(x) * psi_k(s), {k, 1, m}]
template<typename Type>
Type getSecondFredholmIntegral_DegKernel(std::vector<Type> &solution, std::size_t numOfXIntervals, Type a, Type b, Type lambda, 
const std::vector<Type(*)(Type)> &phiVec, const std::vector<Type(*)(Type)> &psiVec, Type (*f)(Type), 
Type (*quadMethod)(Type, std::size_t, Type, Type(*)(Type), Type(*)(Type)), SYSTEM_FLAG sysMethod);


// Вычисление специального интеграла с особенностью
template<typename Type>
Type solveSingularIntegralEq(Type (*f)(Type, Type), std::size_t numOfFinElems, std::vector<Type> &solution);

#endif
