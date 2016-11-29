#ifndef SIMPLE_PIC_PARAM_HPP
#define SIMPLE_PIC_PARAM_HPP

#include <cstdio>
#include <string>
#include <cmath>
#include <type_traits>

#include "Boundary.hpp"
#include "ShapeFactor.hpp"
#include "Vector.hpp"
#include "Particle.hpp"

template<Boundary> class FieldFrame;
template<Boundary> class PlasmaFrame;
template<Boundary, Shape> class SolverFrame;

//using Field = FieldFrame<PML>;
using Field = FieldFrame<Boundary::GradZero>;
using Plasma = PlasmaFrame<Boundary::Periodic>;
using Solver = SolverFrame<Boundary::Periodic, Shape::TSC>;

namespace
{
    constexpr int NUM_DENS = 20;

    constexpr double MASS_RATIO = 100.0;
    constexpr double ELE_MASS = 0.0625;
    constexpr double ION_MASS = ELE_MASS * MASS_RATIO;

    constexpr double ELE_WPE = 0.05; // wpe * delt <= 0.1
    constexpr double ION_WPE = std::sqrt(1.0 / MASS_RATIO) * ELE_WPE;

    constexpr double ELE_Q = - ELE_WPE * std::sqrt(ELE_MASS / (double)NUM_DENS);
    constexpr double ION_Q = - ELE_Q;

    constexpr int RANDOM_SEED = 1000;

    constexpr double C = 0.5;
    constexpr double C2 = C * C;

    constexpr double ELE_VTH = 0.1 * C;
    constexpr double ION_VTH = std::sqrt(1.0 / MASS_RATIO) * ELE_VTH;

    const std::string PATH("data/");

    constexpr int MAX_TIME_STEP = 4000;
    constexpr int OUTPUT_STEP = 10;
    constexpr int SORT_STEP = 50;

    // Boundaryconditions
    constexpr bool enabledPML = std::is_same<Field, FieldFrame<Boundary::PML>>::value; // PML
    constexpr bool enabledABS = std::is_same<Field, FieldFrame<Boundary::Absorbing>>::value; // Absorbing
    constexpr bool enabledMSK = std::is_same<Field, FieldFrame<Boundary::Masking>>::value; // Masking
    constexpr bool enabledPRI = std::is_same<Field, FieldFrame<Boundary::Periodic>>::value; // Periodic

    constexpr int ERROR_VALUE = -1.0;

    // PML and Absorbing
    constexpr int L = (enabledPML || enabledABS) ? 32 : ERROR_VALUE;
    constexpr int M = (enabledPML || enabledABS) ? 4 : ERROR_VALUE;
    constexpr double R0 = (enabledPML || enabledABS) ? 1.0e-10 : ERROR_VALUE;
    constexpr double SIGMA_MAX = (enabledPML || enabledABS) ? (- std::log(R0) * (M + 1.0) * C / (2.0 * L)) : ERROR_VALUE;

    //Masking
    constexpr int MASKING_L = (enabledMSK) ? 100 : ERROR_VALUE;

    // system
    constexpr int LX0 = 5;
    constexpr int LY0 = 440;
    constexpr int LZ0 = 440;

    constexpr int LX = LX0 + 5;
    constexpr int LY = (enabledPML) ? LY0 + 2 * (L + 1)   // PML
                                    : LY0 + 5;            // Periodic and Masking and Absorbing
    constexpr int LZ = (enabledPML) ? LZ0 + 2 * (L + 1)   // PML
                     : (enabledABS) ? LZ0 + 2 * (L + 1)   // Absorbing
                     : (enabledMSK) ? LZ0 + 2 * MASKING_L // Masking
                                    : LZ0 + 5;            // Periodic

    // [X0 X1), [Y0 Y1), [Z0 Z1)
    constexpr int X0 = 2;
    constexpr int Y0 = (enabledPML) ? (L + 1): // PML
                                      2;       // Periodic and Masking and Absorbing

    constexpr int Z0 = (enabledPML) ? (L + 1)   // PML
                     : (enabledABS) ? (L + 1)   // Absorbing
                     : (enabledMSK) ? MASKING_L // Masking
                                    : 2;        // Periodic

    constexpr int X1 = LX - 3;
    constexpr int Y1 = (enabledPML) ? LY - (L + 1) // PML
                                    : LY - 3;      // Periodic and Masking and Absorbing
    constexpr int Z1 = (enabledPML) ? LZ - (L + 1)   // PML
                     : (enabledABS) ? LZ - (L + 1)   // Absorbing
                     : (enabledMSK) ? LZ - MASKING_L // Masking
                                    : LZ - 3;        // Periodic
};

namespace Filament
{
    constexpr int R = 20;
    constexpr int R2 = R * R;
};

#endif

