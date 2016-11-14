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

static_assert(std::is_trivially_copyable<Particle>::value == true, "trivially_copyableじゃないよ");
static_assert(std::is_standard_layout<Particle>::value == true, "standard_layoutじゃないよ");


static_assert(sizeof(double) == 8, "doubleが８bytesじゃーないよ");
static_assert(sizeof(Particle) == sizeof(double) * 7, "paddingされてるよ");

static_assert(offsetof(Particle, id)  == sizeof(double) * 0, "メンバー変数がずれてるよ");
static_assert(offsetof(Particle, r.x) == sizeof(double) * 1, "メンバー変数がずれてるよ");
static_assert(offsetof(Particle, r.y) == sizeof(double) * 2, "メンバー変数がずれてるよ");
static_assert(offsetof(Particle, r.z) == sizeof(double) * 3, "メンバー変数がずれてるよ");
static_assert(offsetof(Particle, v.x) == sizeof(double) * 4, "メンバー変数がずれてるよ");
static_assert(offsetof(Particle, v.y) == sizeof(double) * 5, "メンバー変数がずれてるよ");
static_assert(offsetof(Particle, v.z) == sizeof(double) * 6, "メンバー変数がずれてるよ");
#endif

