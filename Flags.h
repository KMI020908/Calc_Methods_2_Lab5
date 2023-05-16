#ifndef FLAGS_H
#define FLAGS_H

enum SOLUTION_FLAG{
    NO_SOLUTION, // 0
    HAS_SOLUTION // 1 
};

enum FILE_FLAG{
    NOT_OPEN, // 0
    IS_CLOSED // 1
};

enum INVERTIBLE_FLAG{
    NOT_INVERTIBLE, // 0
    IS_INVERTIBLE // 1
};

enum QUADRATIC_FLAG{
    NOT_QUADRATIC, // 0
    IS_QUADRATIC // 1
};

enum MULTIPLIED_FLAG{
    NOT_MULTIPLIED, // 0
    IS_MULTIPLIED // 1
};

enum CORRECT_INPUT_FLAG{
    ERROR, // 0
    CORRECT // 1
};

enum ITERATION_METHOD_FLAG{
    SIMPLE_IT, // 0
    JACOBI, // 1
    RELAXATION // 2
};

enum GRID_FLAG{
    UNIFORM, // 0
    CHEBYSHEV // 1
};

enum VARIABLE_FLAG{
    X, // 0
    Y  // 1  
};

enum DIFF_METHOD_FLAG{
    FW_EULER, // 0
    BW_EULER, // 1
    SYM_SCHEME, // 2
    TWICE_RG, // 3
    FOURTH_RG, // 4
    FOURTH_AD, // 5
    PREDICT_CORRECT // 6
};

enum CONDS_FLAG{
    LT_RT,  // 0 - Температура на обоих концах
    LT_RQ,  // 1 - Температура на левом конце, а поток на правом
    LQ_RT,  // 2 - Поток на левом конце, а температура на правом
    LQ_RQ   // 3 - Поток на обоих концах
};

enum AXIS_FLAG{
    OX, // 0 - проход вдоль Ox
    OY // 1 - проход вдоль Oy
};

enum BOUND_FLAG{
    Temp, // 0 - на границе задано распределение температуры
    Flux // 1 - на границе задан тепловой поток 
};

enum SYSTEM_FLAG{
    GM, // 0 - метод Гаусса
    QR, // 1 - QR метод
    SM, // 2 - Простая итерация
    JM, // 3 - Метод Якоби
    RM  // 4 - Метод Релаксации
};

#endif 