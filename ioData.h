#ifndef IO_DATA_H
#define IO_DATA_H

#include<string>
#include<fstream> 
#include<iomanip> 
#include<vector>
#include"PRNG.h"
#include"Flags.h"

template<typename Type>
FILE_FLAG generateRandomTest(size_t equationDim, Type minValue, Type maxValue, const std::string& FILE_PATH);

template<typename Type>
FILE_FLAG readData(std::vector<std::vector<Type>> &lCoefs, std::vector<Type> &rCoefs,  const std::string& IN_FILE_PATH);

template<typename Type>
FILE_FLAG writeData(const std::vector<Type> &solution, const std::vector<Type> &startPoint, Type accuracy, const std::string& OUT_FILE_PATH, std::size_t numOfIt, 
Type tao = 0.0, Type omega = 0.0);

template<typename Type>
FILE_FLAG readMatrix(std::vector<std::vector<Type>> &matrix, const std::string& IN_FILE_PATH);

template<typename Type>
FILE_FLAG addData(const std::vector<Type> &solution, const std::vector<Type> &startPoint, Type accuracy, const std::string& OUT_FILE_PATH, std::size_t numOfIt, 
Type tao = 0.0, Type omega = 0.0);

template<typename Type>
FILE_FLAG writeBoundOfIterations(Type bound, const std::string& OUT_FILE_PATH);

template<typename Type>
FILE_FLAG writeCanonicalForm(const std::vector<std::vector<Type>> &C, const std::vector<Type> &y, const std::string& OUT_FILE_PATH);

template<typename Type>
FILE_FLAG writeQRMatrix(const std::vector<std::vector<Type>> &Q, const std::vector<std::vector<Type>> &R, const std::string& OUT_FILE_PATH);

template<typename Type>
FILE_FLAG writeMatrixFile(const std::vector<std::vector<Type>> &matrix, const std::string& OUT_FILE_PATH, bool add = false);

template<typename Type>
FILE_FLAG writeVectorFile(const std::vector<Type> &vector, const std::string& OUT_FILE_PATH, bool add = false);

template<typename Type>
FILE_FLAG writeScalarFile(Type num, const std::string& OUT_FILE_PATH, bool add = false);

template<typename Type>
FILE_FLAG writeResidual(Type residual, const std::string& OUT_FILE_PATH);

template<typename Type>
FILE_FLAG addPerturbation(const std::vector<Type> &solution, const std::string& OUT_FILE_PATH, Type perturbation, SOLUTION_FLAG FLAG = HAS_SOLUTION);

template<typename Type>
FILE_FLAG writeConds(Type cond_1, Type cond_inf, const std::string& OUT_FILE_PATH);

template<typename Type>
FILE_FLAG writeLowerBounds(Type lowerBound1, Type lowerBoundInf, const std::string& OUT_FILE_PATH);

template<typename Type>
FILE_FLAG writePointsOfSimpleItMethod(const std::vector<std::vector<Type>> &lCoefs, const std::vector<Type> &rCoefs, const std::vector<Type> 
&firstVec, std::vector<Type> &solution, Type tao, const std::string& OUT_FILE_PATH, Type accuracy = 1e-7, double p = 2.0, Type epsilon_0 = 1e-4);

template<typename Type>
FILE_FLAG writePointsOfJacobiMethod(const std::vector<std::vector<Type>> &lCoefs, const std::vector<Type> &rCoefs, 
const std::vector<Type> &firstVec, std::vector<Type> &solution, const std::string& OUT_FILE_PATH, Type accuracy = 1e-7, double p = 2.0, Type epsilon_0 = 1e-4);

template<typename Type>
FILE_FLAG writePointsOfRelaxationMethod(const std::vector<std::vector<Type>> &lCoefs, const std::vector<Type> &rCoefs, 
const std::vector<Type> &firstVec, std::vector<Type> &solution,  const std::string& OUT_FILE_PATH, Type accuracy = 1e-7, Type omega = 1.0, double p = 2.0, Type epsilon_0 = 1e-4);

template<typename Type>
FILE_FLAG writeNormOfError(std::vector<Type> &solution, std::vector<Type> &realSolution, const std::string& OUT_FILE_PATH, double p = 2.0);

template<typename Type>
FILE_FLAG writeNormErrAfterEstIt(Type normErr, const std::string& OUT_FILE_PATH);

FILE_FLAG writeExactIters(std::size_t exactIterations, const std::string& OUT_FILE_PATH);

template<typename Type>
FILE_FLAG writeEigenData(std::size_t numOfIterations, std::vector<Type> &eigList, const std::string& OUT_FILE_PATH, 
bool hasShift = false, bool makeHessenberg = false, bool add = false);

template<typename Type>
FILE_FLAG writeEigenVec(std::size_t numOfIters, std::vector<std::vector<Type>> &eigMatrix, std::vector<Type> &eigList, const std::string& OUT_FILE_PATH, 
bool add = false);

// Лаба 3
template<typename Type>
FILE_FLAG readValueTable(const std::vector<Type> &xVec, const std::vector<Type> &fVec, const std::string &FILE_PATH);

template<typename Type>
FILE_FLAG writeEqResBissection(Type sol, Type firstX, Type lastX, Type accuracy, const std::string& OUT_FILE_PATH, bool add = false);

template<typename Type>
FILE_FLAG writeEqResNewthon(Type sol, Type x0, Type accuracy, const std::string& OUT_FILE_PATH, bool add = false);

template<typename Type>
FILE_FLAG writeSysResNewthon(const std::vector<Type> &sol, Type x0, Type y0, Type accuracy, const std::string& OUT_FILE_PATH, bool add = false);

FILE_FLAG writeIters(std::size_t exactIterations, const std::string& OUT_FILE_PATH);

template<typename Type>
std::string getFileNameHeatEq(std::size_t numOfEq, CONDS_FLAG flag, Type sigma, const std::string &name, const std::string &folder = "LinearHeatEq");

std::string getFileNameQuasiHeatEq(std::size_t numOfEq, CONDS_FLAG flag, const std::string &name, const std::string &folder = "QuasilinearHeatEq");

template<typename Type>
FILE_FLAG readMatrix(std::vector<std::vector<Type>> &matrix, std::size_t rows, std::size_t cols, const std::string& IN_FILE_PATH);

#endif