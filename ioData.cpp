#include "ioData.h"

template<typename Type>
FILE_FLAG generateRandomTest(size_t equationDim, Type minValue, Type maxValue, const std::string& FILE_PATH){
    std::ofstream file;
    file.open(FILE_PATH);
    if (!file.is_open())
        exit(NOT_OPEN);
    file << equationDim << '\n';
    PRNG generator;
    for (std::size_t i = 0; i < equationDim; i++){
        for (std::size_t j = 0; j < equationDim + 1; j++)
            file << getRandomNumber(generator, minValue, maxValue) << '\t';
        file << '\n';
    }
    file.close();
    return IS_CLOSED;
}

template<typename Type>
FILE_FLAG readData(std::vector<std::vector<Type>> &lCoefs, std::vector<Type> &rCoefs,  const std::string& IN_FILE_PATH) {
	std::ifstream file;
	file.open(IN_FILE_PATH);
	if (!file.is_open())
		exit(NOT_OPEN);
	std::size_t n;
	file >> n;
    std::vector<Type> hVec;
    hVec.reserve(n);
    Type hValue = 0;
    lCoefs.clear();
    rCoefs.clear();
    for (std::size_t i = 0; i < n; i++){
        for (std::size_t j = 0; j < n; j++){ 
            file >> hValue;
            hVec.push_back(hValue);    
        }
        file >> hValue;
        rCoefs.push_back(hValue);
        lCoefs.push_back(hVec);
        hVec.clear();
    }
	file.close();
	return IS_CLOSED;
}

template<typename Type>
FILE_FLAG readMatrix(std::vector<std::vector<Type>> &matrix, const std::string& IN_FILE_PATH){
	std::ifstream file;
	file.open(IN_FILE_PATH);
	if (!file.is_open())
		exit(NOT_OPEN);
	std::size_t n;
	file >> n;
    std::vector<Type> hVec;
    hVec.reserve(n);
    Type hValue = 0;
    matrix.clear();
    for (std::size_t i = 0; i < n; i++){
        for (std::size_t j = 0; j < n; j++){ 
            file >> hValue;
            hVec.push_back(hValue);    
        }
        matrix.push_back(hVec);
        hVec.clear();
    }
	file.close();
	return IS_CLOSED;
}

template<typename Type>
FILE_FLAG readMatrix(std::vector<std::vector<Type>> &matrix, std::size_t rows, std::size_t cols, const std::string& IN_FILE_PATH){
	std::ifstream file;
	file.open(IN_FILE_PATH);
	if (!file.is_open())
		exit(NOT_OPEN);
    std::vector<Type> hVec;
    hVec.reserve(cols);
    Type hValue = 0;
    matrix.clear();
    for (std::size_t i = 0; i < rows; i++){
        for (std::size_t j = 0; j < cols; j++){ 
            file >> hValue;
            hVec.push_back(hValue);    
        }
        matrix.push_back(hVec);
        hVec.clear();
    }
	file.close();
	return IS_CLOSED;
}

template<typename Type>
FILE_FLAG writeData(const std::vector<Type> &solution, const std::vector<Type> &startPoint, Type accuracy, const std::string& OUT_FILE_PATH, std::size_t numOfIt, 
Type tao, Type omega){
	std::ofstream file;
	file.open(OUT_FILE_PATH);
	if (!file.is_open())
		exit(NOT_OPEN);
    if (!numOfIt){
        file << "The solution:" << '\n';
	    file << "No solution.";
        file.close();
        return IS_CLOSED;
    }
    file << "Start point: " << startPoint << '\n';
    if (omega != 0.0){
        file << "omega = " << omega << '\n';
    }
    if (tao != 0.0){
        file << "tao = " << tao << '\n';
    }
    file << "The solution:" << '\n';
	file << "X = " << std::setprecision(std::abs(std::log10(accuracy))) << solution << '\n';
    file << "Number of iterations: " << numOfIt;
	file.close();
	return IS_CLOSED;
}

template<typename Type>
FILE_FLAG addData(const std::vector<Type> &solution, const std::vector<Type> &startPoint, Type accuracy, const std::string& OUT_FILE_PATH, std::size_t numOfIt, 
Type tao, Type omega){
	std::ofstream file;
	file.open(OUT_FILE_PATH, std::ios::app);
	if (!file.is_open())
		exit(NOT_OPEN);
    if (!numOfIt){
        file << "The solution:" << '\n';
	    file << "No solution.";
        file.close();
        return IS_CLOSED;
    }
    file << '\n' << '\n';
    file << "Start point: " << startPoint << '\n';
    if (omega != 0.0){
        file << "omega = " << omega << '\n';
    }
    if (tao != 0.0){
        file << "tao = " << tao << '\n';
    }
    file << "The solution:" << '\n';
	file << "X = " << std::setprecision(std::abs(std::log10(accuracy))) << solution << '\n';
    file << "Number of iterations: " << numOfIt;
	file.close();
	return IS_CLOSED;
}

template<typename Type>
FILE_FLAG writeBoundOfIterations(Type bound, const std::string& OUT_FILE_PATH){
    std::ofstream file;
	file.open(OUT_FILE_PATH, std::ios::app);
	if (!file.is_open())
		exit(NOT_OPEN);
    file << '\n';
    file << "Number of iterations > " << bound;
    file.close();
    return IS_CLOSED;
}

template<typename Type>
FILE_FLAG writeCanonicalForm(const std::vector<std::vector<Type>> &C, const std::vector<Type> &y, const std::string& OUT_FILE_PATH){
    std::ofstream file;
	file.open(OUT_FILE_PATH, std::ios::app);
	if (!file.is_open())
		exit(NOT_OPEN);
    std::size_t dimMatrix = C.size(); 
    file << '\n' << '\n';
    file << "C matrix:" << '\n';
    for (std::size_t i = 0; i < dimMatrix; i++){
        for (std::size_t j = 0; j < dimMatrix; j++){
            file << C[i][j] << '\t';    
        }
        file << '\n';
    }
    file << '\n';
    file << "Norm 1 of C: " << normOfMatrix(C, 1.0) << '\n';
    file << "Norm Inf of C: " << normOfMatrix(C, INFINITY) << '\n';
    file << '\n';
    file << "y vector:" << '\n';
    file << y;
    file.close();
	return IS_CLOSED;
}

template<typename Type>
FILE_FLAG writeQRMatrix(const std::vector<std::vector<Type>> &Q, const std::vector<std::vector<Type>> &R, const std::string& OUT_FILE_PATH){
	std::ofstream file;
	file.open(OUT_FILE_PATH, std::ios::app);
	if (!file.is_open())
		exit(NOT_OPEN);
    std::size_t dimMatrix = Q.size(); 
    file << '\n' << '\n';
    file << "Q matrix:" << '\n';
    for (std::size_t i = 0; i < dimMatrix; i++){
        for (std::size_t j = 0; j < dimMatrix; j++){
            file << Q[i][j] << '\t';    
        }
        file << '\n';
    }
    file << '\n';
    file << "R matrix:" << '\n';
    for (std::size_t i = 0; i < dimMatrix; i++){
        for (std::size_t j = 0; j < dimMatrix; j++){
            file << R[i][j] << '\t';    
        }
        file << '\n';
    }
	file.close();
	return IS_CLOSED;
}

template<typename Type>
FILE_FLAG writeScalarFile(Type num, const std::string& OUT_FILE_PATH, bool add){
	std::ofstream file;
	if (add){
        file.open(OUT_FILE_PATH, std::ios::app);
        file << '\n' << '\n';
    }
    else{
        file.open(OUT_FILE_PATH);
    }
	if (!file.is_open())
		exit(NOT_OPEN);
    file << num;
	file.close();
	return IS_CLOSED;
}

template<typename Type>
FILE_FLAG writeMatrixFile(const std::vector<std::vector<Type>> &matrix, const std::string& OUT_FILE_PATH, bool add){
    std::ofstream file;
	if (add){
        file.open(OUT_FILE_PATH, std::ios::app);
        file << '\n' << '\n';
    }
    else{
        file.open(OUT_FILE_PATH);
    }
	if (!file.is_open())
		exit(NOT_OPEN);
    std::size_t rows = matrix.size();
    for (std::size_t i = 0; i < rows; i++){
        for (std::size_t j = 0; j < matrix[i].size(); j++){
            file << matrix[i][j] << '\t';
        }
        file << '\n';
    }
    file.close();
    return IS_CLOSED;
} 

template<typename Type>
FILE_FLAG writeVectorFile(const std::vector<Type> &vector, const std::string& OUT_FILE_PATH, bool add){
    std::ofstream file;
	if (add){
        file.open(OUT_FILE_PATH, std::ios::app);
        file << '\n' << '\n';
    }
    else{
        file.open(OUT_FILE_PATH);
    }
	if (!file.is_open())
		exit(NOT_OPEN);
    std::size_t rows = vector.size();
    if (rows == 0){
        file.close();
        return IS_CLOSED;
    }
    for (std::size_t i = 0; i < rows - 1; i++){
        file << vector[i] << '\n';
    }
     file << vector[rows - 1];
    file.close();
    return IS_CLOSED;
} 

template<typename Type>
FILE_FLAG writeResidual(Type residual, const std::string& OUT_FILE_PATH){
	std::ofstream file;
	file.open(OUT_FILE_PATH, std::ios::app);
	if (!file.is_open())
		exit(NOT_OPEN); 
    file << '\n';
    file << "Residual = " << residual;
	file.close();
	return IS_CLOSED;
}

template<typename Type>
FILE_FLAG addPerturbation(const std::vector<Type> &solution, const std::string& OUT_FILE_PATH, Type perturbation, SOLUTION_FLAG FLAG){
	std::ofstream file;
	file.open(OUT_FILE_PATH, std::ios::app);
	if (!file.is_open())
		exit(NOT_OPEN);
    if (FLAG == NO_SOLUTION){
        file << '\n' << '\n';
        file << "The solution:" << '\n';
	    file << "No solution.";
        file.close();
        return IS_CLOSED;
    }
    file << '\n' << '\n';
    file << "The solution after perturbation = " << perturbation << '\n';
	file << "X1 = " << "{ ";
    for (std::size_t i = 0; i < solution.size() - 1; i++)
        file << solution[i] << ", ";
    file << solution[solution.size() - 1] << ' ';
    file << '}';
	file.close();
	return IS_CLOSED;
}

template<typename Type>
FILE_FLAG writeConds(Type cond_1, Type cond_inf, const std::string& OUT_FILE_PATH){
    std::ofstream file;
	file.open(OUT_FILE_PATH, std::ios::app);
	if (!file.is_open())
		exit(NOT_OPEN);
    file << '\n' << '\n';
    file << "cond_1 A = " << cond_1;
    file << '\n' << '\n';
    file << "cond_inf A = " << cond_inf;
    file.close();
    return IS_CLOSED;
} 

template<typename Type>
FILE_FLAG writeLowerBounds(Type lowerBound1, Type lowerBoundInf, const std::string& OUT_FILE_PATH){
    std::ofstream file;
	file.open(OUT_FILE_PATH, std::ios::app);
	if (!file.is_open())
		exit(NOT_OPEN);
    file << '\n' << '\n';
    file << lowerBound1 << " ≤ cond_1 A";
    file << '\n' << '\n';
    file << lowerBoundInf << " ≤ cond_inf A";
    file.close();
    return IS_CLOSED;
}

template<typename Type>
FILE_FLAG writePointsOfSimpleItMethod(const std::vector<std::vector<Type>> &lCoefs, const std::vector<Type> &rCoefs, const std::vector<Type> 
&firstVec, std::vector<Type> &solution, Type tao, const std::string& OUT_FILE_PATH, Type accuracy, double p, Type epsilon_0){
    std::ofstream file;
	file.open(OUT_FILE_PATH);
	if (!file.is_open())
		exit(NOT_OPEN);
    std::size_t rows = lCoefs.size(); // Количество строк в СЛАУ
    solution.resize(rows); // Искомое решение
    std::size_t cols = 0;
    if (rows != 0)
        cols = lCoefs[0].size();
    else{
        file.close();
        return IS_CLOSED;
    }
    std::vector<Type> prev_solution = firstVec;
    solution = prev_solution;
    file << prev_solution << '\n';
    for (std::size_t i = 0; i < rows; i++){
        Type sum = 0.0;
        for (std::size_t j = 0; j < cols; j++){
            sum += lCoefs[i][j] * prev_solution[j];
        }
        solution[i] = prev_solution[i] + tao * (rCoefs[i] - sum);
        file << solution << '\n';
    }
    std::vector<Type> diffVec = solution - prev_solution;
    Type diffNorm = normOfVector(diffVec, p);
    while (diffNorm / (normOfVector(prev_solution, p) + epsilon_0) > accuracy || diffNorm > accuracy){
        prev_solution = solution;
        for (std::size_t i = 0; i < rows; i++){
            Type sum = 0.0;
            for (std::size_t j = 0; j < cols; j++){
                sum += lCoefs[i][j] * prev_solution[j];
            }
            solution[i] = prev_solution[i] + tao * (rCoefs[i] - sum);
            file << solution << '\n';
        }
        diffVec = solution - prev_solution;
        diffNorm = normOfVector(diffVec, p);
    }
    file.close();
    return IS_CLOSED;
}

template<typename Type>
FILE_FLAG writePointsOfJacobiMethod(const std::vector<std::vector<Type>> &lCoefs, const std::vector<Type> &rCoefs, 
const std::vector<Type> &firstVec, std::vector<Type> &solution, const std::string& OUT_FILE_PATH, Type accuracy, double p, Type epsilon_0){
    std::ofstream file;
	file.open(OUT_FILE_PATH);
	if (!file.is_open())
		exit(NOT_OPEN);
    std::size_t rows = lCoefs.size(); // Количество строк в СЛАУ
    solution.resize(rows); // Искомое решение
    std::size_t cols = 0;
    if (rows != 0)
        cols = lCoefs[0].size();
    else{
        file.close();
        return IS_CLOSED;
    }
    std::vector<Type> prev_solution = firstVec;
    solution = prev_solution;
    file << prev_solution << '\n';
    for (std::size_t i = 0; i < rows; i++){
        Type sum = 0.0;
        for (std::size_t j = 0; j < cols; j++){
            if (i != j){
                sum += lCoefs[i][j] * prev_solution[j];
            }
        }
        solution[i] = (1/lCoefs[i][i]) * (rCoefs[i] - sum);
        file << solution << '\n';
    }
    std::vector<Type> diffVec = solution - prev_solution;
    Type diffNorm = normOfVector(diffVec, p);
    while (diffNorm / (normOfVector(prev_solution, p) + epsilon_0) > accuracy){
        prev_solution = solution;
        for (std::size_t i = 0; i < rows; i++){
            Type sum = 0.0;
            for (std::size_t j = 0; j < cols; j++){
                if (i != j){
                    sum += lCoefs[i][j] * prev_solution[j];
                }
            }
            solution[i] = (1/lCoefs[i][i]) * (rCoefs[i] - sum);
            file << solution << '\n';
        }   
        diffVec = solution - prev_solution;
        diffNorm = normOfVector(diffVec, p);
    }
    file.close();
    return IS_CLOSED;
}

template<typename Type>
FILE_FLAG writePointsOfRelaxationMethod(const std::vector<std::vector<Type>> &lCoefs, const std::vector<Type> &rCoefs, 
const std::vector<Type> &firstVec, std::vector<Type> &solution,  const std::string& OUT_FILE_PATH, Type accuracy, Type omega, double p, Type epsilon_0){
    std::ofstream file;
	file.open(OUT_FILE_PATH);
	if (!file.is_open())
		exit(NOT_OPEN);
    std::size_t rows = lCoefs.size(); // Количество строк в СЛАУ
    solution.resize(rows); // Искомое решение
    std::size_t cols = 0;
    if (rows != 0)
        cols = lCoefs[0].size();
    else{
        file.close();
        return IS_CLOSED;
    }
    std::vector<Type> prev_solution = firstVec;
    solution = prev_solution;
    file << prev_solution << '\n';
    for (std::size_t i = 0; i < rows; i++){
        Type sum1 = 0.0;
        for (std::size_t j = 0; j < i; j++){
            sum1 += lCoefs[i][j] * solution[j];
        }
        Type sum2 = 0.0;
        for (std::size_t j = i + 1; j < cols; j++){
            sum2 += lCoefs[i][j] * solution[j];
        }
        solution[i] = (1 - omega) * prev_solution[i] - (omega / lCoefs[i][i]) * (sum1 + sum2 - rCoefs[i]);
        file << solution << '\n';
    }
    std::vector<Type> diffVec = solution - prev_solution;
    Type diffNorm = normOfVector(diffVec, p);
    while (diffNorm / (normOfVector(prev_solution, p) + epsilon_0) > accuracy || diffNorm > accuracy){
        prev_solution = solution;
        for (std::size_t i = 0; i < rows; i++){
            Type sum1 = 0.0;
            for (std::size_t j = 0; j < i; j++){
                sum1 += lCoefs[i][j] * solution[j];
            }
            Type sum2 = 0.0;
            for (std::size_t j = i + 1; j < cols; j++){
                sum2 += lCoefs[i][j] * solution[j];
            }
            solution[i] = (1 - omega) * prev_solution[i] - (omega / lCoefs[i][i]) * (sum1 + sum2 - rCoefs[i]);
            file << solution << '\n';
        }
        diffVec = solution - prev_solution;
        diffNorm = normOfVector(diffVec, p);
    }
    file.close();
    return IS_CLOSED;
}

template<typename Type>
FILE_FLAG writeNormOfError(std::vector<Type> &solution, std::vector<Type> &realSolution, const std::string& OUT_FILE_PATH, double p){
    std::ofstream file;
	file.open(OUT_FILE_PATH, std::ios::app);
	if (!file.is_open())
		exit(NOT_OPEN);
    file << '\n';
    file << "Norm of error: " << normOfVector(solution - realSolution, p);
    file.close();
    return IS_CLOSED;
}

FILE_FLAG writeExactIters(std::size_t exactIterations, const std::string& OUT_FILE_PATH){
	std::ofstream file;
	file.open(OUT_FILE_PATH, std::ios::app);
	if (!file.is_open())
		exit(NOT_OPEN); 
    file << '\n';
    file << "Exact iterations = " << exactIterations;
	file.close();
	return IS_CLOSED;
}

template<typename Type>
FILE_FLAG writeNormErrAfterEstIt(Type normErr, const std::string& OUT_FILE_PATH){
	std::ofstream file;
	file.open(OUT_FILE_PATH, std::ios::app);
	if (!file.is_open())
		exit(NOT_OPEN); 
    file << '\n';
    file << "Norm of error after ets. iterations = " << normErr;
	file.close();
	return IS_CLOSED;
}

// Лаба 3
template<typename Type>
FILE_FLAG writeEigenData(std::size_t numOfIterations, std::vector<Type> &eigList, const std::string& OUT_FILE_PATH, 
bool hasShift, bool makeHessenberg, bool add){
	std::ofstream file;
    if (add){
        file.open(OUT_FILE_PATH, std::ios::app);
        file << '\n' << '\n';
    }
    else{
        file.open(OUT_FILE_PATH);
    }
	if (!file.is_open())
		exit(NOT_OPEN); 
    if (!hasShift && !makeHessenberg){
        file << "QR method without shifts and without Hessenberg transform. :" << '\n';
        file << "Eigen values of matrix: " << eigList << '\n';
        file << "Number of iterations: " << numOfIterations;
    }
    if (hasShift && !makeHessenberg){
        file << "QR method with shifts and without Hessenberg transform. :" << '\n';
        file << "Eigen values of matrix: " << eigList << '\n';
        file << "Number of iterations: " << numOfIterations;
    }
    if (!hasShift && makeHessenberg){
        file << "QR method without shifts and with Hessenberg transform. :" << '\n';
        file << "Eigen values of matrix: " << eigList << '\n';
        file << "Number of iterations: " << numOfIterations;
    }
    if (hasShift && makeHessenberg){
        file << "QR method with shifts and with Hessenberg transform. :" << '\n';
        file << "Eigen values of matrix: " << eigList << '\n';
        file << "Number of iterations: " << numOfIterations;
    }
	file.close();
	return IS_CLOSED;
}

template<typename Type>
FILE_FLAG writeEigenVec(std::size_t numOfIters, std::vector<std::vector<Type>> &eigMatrix, std::vector<Type> &eigList, const std::string& OUT_FILE_PATH, 
bool add){
	std::ofstream file;
    if (add){
        file.open(OUT_FILE_PATH, std::ios::app);
        file << '\n' << '\n';
    }
    else{
        file.open(OUT_FILE_PATH);
    }
	if (!file.is_open())
		exit(NOT_OPEN);
    file << "Number of iterations: " << numOfIters << '\n' << '\n'; 
    std::size_t rows = eigMatrix.size();
    for (std::size_t i = 0; i < rows; i++){
        file << "Eigen number: " << eigList[i] << " ----> " << "Eigen vector: " << eigMatrix[i] << '\n' << '\n';
    }
	file.close();
	return IS_CLOSED;
}

// Лаба 3

template<typename Type>
FILE_FLAG writeValueTable(const std::vector<Type> &xVec, const std::vector<Type> &fVec, const std::string &FILE_PATH){
    std::ofstream file;
    file.open(FILE_PATH);
    if (!file.is_open()){
        return(NOT_OPEN);
    }
    
    std::size_t numOfUnits = xVec.size();
    if (fVec.size() != numOfUnits){
        file.close();
        return IS_CLOSED;
    }
    for (std::size_t i = 0; i < numOfUnits; i++){
        file << xVec[i] << '\t' << '\t' << fVec[i] << '\n';
    }
    file.close();
    return IS_CLOSED;
}

template<typename Type>
FILE_FLAG writeEqResBissection(Type sol, Type firstX, Type lastX, Type accuracy, const std::string& OUT_FILE_PATH, bool add){
	std::ofstream file;
	if (add){
        file.open(OUT_FILE_PATH, std::ios::app);
        file << '\n' << '\n';
    }
    else{
        file.open(OUT_FILE_PATH);
    }
	if (!file.is_open())
		exit(NOT_OPEN); 
    file << "accuracy = " << accuracy;
    file << '\n';
    file << "[a, b] = " << '[' << firstX << ", " << lastX << ']';
    file << '\n';
    file << "x = " << sol;
	file.close();
	return IS_CLOSED;
}

template<typename Type>
FILE_FLAG writeEqResNewthon(Type sol, Type x0, Type accuracy, const std::string& OUT_FILE_PATH, bool add){
	std::ofstream file;
	if (add){
        file.open(OUT_FILE_PATH, std::ios::app);
        file << '\n' << '\n';
    }
    else{
        file.open(OUT_FILE_PATH);
    }
	if (!file.is_open())
		exit(NOT_OPEN); 
    file << "accuracy = " << accuracy;
    file << '\n';
    file << "x0 = " << x0;
    file << '\n';
    file << "x = " << sol;
	return IS_CLOSED;
}

template<typename Type>
FILE_FLAG writeSysResNewthon(const std::vector<Type> &sol, Type x0, Type y0, Type accuracy, const std::string& OUT_FILE_PATH, bool add){
	std::ofstream file;
	if (add){
        file.open(OUT_FILE_PATH, std::ios::app);
        file << '\n' << '\n';
    }
    else{
        file.open(OUT_FILE_PATH);
    }
	if (!file.is_open())
		exit(NOT_OPEN); 
    file << "accuracy = " << accuracy;
    file << '\n';
    file << "{x0, y0} = " << '{' << x0 << ", " << y0 << '}';
    file << '\n';
    file << "x = " << sol;
	return IS_CLOSED;
}

FILE_FLAG writeIters(std::size_t exactIterations, const std::string& OUT_FILE_PATH){
	std::ofstream file;
	file.open(OUT_FILE_PATH, std::ios::app);
	if (!file.is_open())
		exit(NOT_OPEN); 
    file << '\n';
    file << "Number of iterations = " << exactIterations;
	file.close();
	return IS_CLOSED;
}

template<typename Type>
std::string getFileNameHeatEq(std::size_t numOfEq, CONDS_FLAG flag, Type sigma, const std::string &name, const std::string &folder){
    std::string condsFolder;
    switch (flag){
        case LT_RT:
            condsFolder = "FirstConds\\";
            break;
        case LT_RQ:
            condsFolder = "MixedCondsLeftT\\";
            break;
        case LQ_RT:
            condsFolder = "MixedCondsRightT\\";
            break;
        case LQ_RQ:
            condsFolder = "SecondConds\\";
            break;
        
        default:
            condsFolder = "FirstConds\\";
            break;
    }
    std::string sigmaFolder;
    if (sigma == 0.0){
        sigmaFolder = "Forward\\";
    }
    else{
        if (sigma == 1.0){
           sigmaFolder = "Backward\\"; 
        }
        else{
            if (sigma == 0.5){
                sigmaFolder = "Symmetric\\";  
            }
            else{
                sigmaFolder = "Another\\";
            }
        }
    }
    return "D:\\Calc_Methods_2\\Lab2\\" + folder + "\\" + condsFolder + sigmaFolder + name + std::to_string(numOfEq) + ".txt";
}

std::string getFileNameQuasiHeatEq(std::size_t numOfEq, CONDS_FLAG flag, const std::string &name, const std::string &folder){
    std::string condsFolder;
    switch (flag){
        case LT_RT:
            condsFolder = "FirstConds\\";
            break;
        case LT_RQ:
            condsFolder = "MixedCondsLeftT\\";
            break;
        case LQ_RT:
            condsFolder = "MixedCondsRightT\\";
            break;
        case LQ_RQ:
            condsFolder = "SecondConds\\";
            break;
        
        default:
            condsFolder = "FirstConds\\";
            break;
    }
    return "D:\\Calc_Methods_2\\Lab2\\" + folder + "\\" + condsFolder + name + std::to_string(numOfEq) + ".txt";
}