#ifndef SIMPLE_PIC_BOUNDARY_HPP
#define SIMPLE_PIC_BOUNDARY_HPP

enum class Boundary : int
{
    // Boundary
    Periodic = 1,
    PML = 2, // x: periodic, z: PML
    GradZero = 3, // x: periodic, z: grad(p) == 0

    Absorbing = -1, // x: periodic, y: PML
    Masking = -2, // x: periodic, y: masking
};

#endif

