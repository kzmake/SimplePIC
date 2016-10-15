#ifndef SIMPLE_PIC_PARAM_HPP
#define SIMPLE_PIC_PARAM_HPP

#include <cstdio>
#include <string>

// define

// static
namespace
{
    const std::string PATH("data/");

    const int MAX_TIME_STEP = 100;
    const int OUTPUT_STEP   = 1;
    const int PARTICLE_SORT = 100;

    // system
    const int LX0 = 50;
    const int LY0 = 50;
    const int LZ0 = 50;
    const int LX = LX0 + 5;
    const int LY = LY0 + 5;
    const int LZ = LZ0 + 5;

    // x = [2 (LX-3))
    const int X0 = 2;
    const int Y0 = 2;
    const int Z0 = 2;
    const int X1 = LX - 3;
    const int Y1 = LY - 3;
    const int Z1 = LZ - 3;

    // sort
    const int PARTICLE_SORT_PREC = 10;

    const double C = 0.5;
    const double C2 = C * C;
};

#endif

