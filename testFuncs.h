// Файл, содержащий функции для тестов
#ifndef TEST_FUNC_H
#define TEST_FUNC_H
#include<cmath>

// Первый тест из методички
template<typename Type>
Type f1(Type x){
    return 0.5 * (1.0 + std::sin(x));
}

template<typename Type>
Type K1(Type x, Type y){
    return 0.5 * (1.0 - x * std::cos(x * y));
}

// Втрой тест из методички
template<typename Type>
Type f2(Type x){
    return std::pow(x, 2.0) + std::sqrt(x);
}

template<typename Type>
Type K2(Type x, Type y){
    return 0.5 * (1.0 - x * std::cos(x * y));
}

// Тест с известным аналитическим решением
template<typename Type>
Type f3(Type x){
    return x * (std::cos(1.0) - std::sin(1.0)) + std::sin(x);
}

template<typename Type>
Type K3(Type x, Type y){
    return x * y;
}

// Второй тест с известным аналитическим решением
template<typename Type>
Type f4(Type x){
    return std::sin(x) - (1.0 + x * std::exp(M_PI_2 * x)) / (1.0 + std::pow(x, 2.0));
}

template<typename Type>
Type K4(Type x, Type y){
    return std::exp(x * y);
}


#endif