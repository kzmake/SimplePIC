#ifndef SIMPLE_PIC_PARTICLE_HPP
#define SIMPLE_PIC_PARTICLE_HPP

#include "Vector.hpp"

class Particle
{
  public:
    long int id;
    Vector r, v;
  public:
    Particle() : id(0) {}
    Particle(const long int id) : id(id) {}
    Particle(const long int id, const Vector& r, const Vector& v) : id(id), r(r), v(v) {}
    Particle(const Vector& r, const Vector& v) : r(r), v(v) {}
};

#endif

