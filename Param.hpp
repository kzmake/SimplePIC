#ifndef SIMPLE_PIC_PARAM_HPP
#define SIMPLE_PIC_PARAM_HPP

#include <cstdio>
#include <string>
#include <cmath>

//#define PIC_PML

namespace
{
    constexpr int NUM_DENS      = 1;

    constexpr double MASS_RATIO = 100.0;
    constexpr double ELE_MASS   = 0.0625;
    constexpr double ION_MASS   = ELE_MASS * MASS_RATIO;

    constexpr double ELE_WPE    = 0.05; // wpe * delt <= 0.1
    constexpr double ION_WPE    = std::sqrt(ELE_MASS / ION_MASS) * ELE_WPE;

    constexpr double ELE_Q      = - ELE_WPE * std::sqrt(ELE_MASS / (double)NUM_DENS);
    constexpr double ION_Q      = - ELE_Q;

    constexpr int RANDOM_SEED   = 1000;

    constexpr double C          = 0.5;
    constexpr double C2         = C * C;

    constexpr double ELE_VTH    = 0.1 * C;
    constexpr double ION_VTH    = std::sqrt(ELE_MASS / ION_MASS) * ELE_VTH;

    const std::string PATH("data/");

    constexpr int MAX_TIME_STEP = 10;
    constexpr int OUTPUT_STEP   = 10;
    constexpr int SORT_STEP = 50;

#ifdef PIC_PML
    constexpr int L  = 32;
    constexpr int M  = 4;
    constexpr double R0 = 1.0e-10;

    constexpr double SIGMA_MAX = - std::log(R0) * (M + 1.0) * C / (2.0 * L);
#endif

    // system
    constexpr int LX0 = 5;
    constexpr int LY0 = 32;
    constexpr int LZ0 = 32;

#ifdef PIC_PML
    constexpr int LX  = LX0 + 5;
    constexpr int LY  = LY0 + 2 * (L + 1);
    constexpr int LZ  = LZ0 + 2 * (L + 1);

    // x = [2 (LX-3))
    constexpr int X0  = 2;
    constexpr int Y0  = (L + 1);
    constexpr int Z0  = (L + 1);
    constexpr int X1  = LX - 3;
    constexpr int Y1  = LY - (L + 1);
    constexpr int Z1  = LZ - (L + 1);
#else
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
#endif
};

#endif
