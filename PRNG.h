#ifndef PRNG_H
#define PRNG_H

#include <random>
#include <cassert>

struct PRNG
{
    std::mt19937 engine;
};

void initGenerator(PRNG& generator)
{
    std::random_device device;
    generator.engine.seed(device());
}

unsigned random(PRNG& generator, unsigned minValue, unsigned maxValue)
{
    assert(minValue < maxValue);
    std::uniform_int_distribution<unsigned> distribution(minValue, maxValue);
    return distribution(generator.engine);
}

template <typename Type>
Type getRandomNumber(PRNG& generator, Type minValue, Type maxValue)
{
    assert(minValue < maxValue);
    std::uniform_real_distribution<Type> distribution(minValue, maxValue);
    return distribution(generator.engine);
}
#endif