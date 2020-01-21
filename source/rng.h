//
// Created by ilaria on 2020-01-09.
//

#ifndef RNG_H
#define RNG_H
#include <random>
#include <iostream>
#include "pcg/pcg_random.hpp"


namespace rn{
    extern pcg32 rng;
    extern void seed(long n);
    extern int uniform_integer_box(const int min, const int max);
    extern double uniform_real_box(const double min, const double max);
}


#endif
