#ifndef SIMPLE_PIC_PARAM_HPP
#define SIMPLE_PIC_PARAM_HPP

#include <cstdio>
#include <string>

// define

// static
namespace
{
    const std::string PATH("data/");

    constexpr int MAX_TIME_STEP = 2000;
    constexpr int OUTPUT_STEP   = 1;
    constexpr int PARTICLE_SORT = 10;

    // system
    constexpr int LX0 = 100;
    constexpr int LY0 = 100;
    constexpr int LZ0 = 25;
    constexpr int LX = LX0 + 5;
    constexpr int LY = LY0 + 5;
    constexpr int LZ = LZ0 + 5;

    // x = [2 (LX-3))
    constexpr int X0 = 2;
    constexpr int Y0 = 2;
    constexpr int Z0 = 2;
    constexpr int X1 = LX - 3;
    constexpr int Y1 = LY - 3;
    constexpr int Z1 = LZ - 3;

    // sort
    constexpr int PARTICLE_SORT_PREC = 10;

    constexpr double C = 0.5;
    constexpr double C2 = C * C;
};

#endif

