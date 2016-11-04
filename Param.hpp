#ifndef SIMPLE_PIC_PARAM_HPP
#define SIMPLE_PIC_PARAM_HPP

#include <cstdio>
#include <string>

namespace
{
    constexpr int NUM_DENS      = 100;

    constexpr double MASS_RATIO = 100.0;
    constexpr double ELE_MASS   = 0.0625;
    constexpr double ION_MASS   = ELE_MASS * MASS_RATIO;

    constexpr double ELE_WPE    = 0.50;
    constexpr double ION_WPE    = sqrt(ELE_MASS / ION_MASS) * ELE_WPE;

    constexpr double ELE_Q      = - ELE_WPE * sqrt(ELE_MASS / (double)NUM_DENS);
    constexpr double ION_Q      = - ELE_Q;

    constexpr int RANDOM_SEED   = 1000;

    constexpr double C          = 0.5;
    constexpr double C2         = C * C;

    constexpr double ELE_VTH    = 0.1 * C;
    constexpr double ION_VTH    = sqrt(ELE_MASS / ION_MASS) * ELE_VTH;

    const std::string PATH("data/");

    constexpr int MAX_TIME_STEP = 1000;
    constexpr int OUTPUT_STEP   = 50;
    constexpr int SORT_STEP = 1;

    // system
    constexpr int LX0 = 5;
    constexpr int LY0 = 50;
    constexpr int LZ0 = 50;
    constexpr int LX  = LX0 + 5;
    constexpr int LY  = LY0 + 5;
    constexpr int LZ  = LZ0 + 5;

    // x = [2 (LX-3))
    constexpr int X0  = 2;
    constexpr int Y0  = 2;
    constexpr int Z0  = 2;
    constexpr int X1  = LX - 3;
    constexpr int Y1  = LY - 3;
    constexpr int Z1  = LZ - 3;

    // sort
    //constexpr int SORT_PREC = 10;
};

#endif

