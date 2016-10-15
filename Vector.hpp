#ifndef SIMPLE_PIC_VECTOR_HPP
#define SIMPLE_PIC_VECTOR_HPP

#include <cstdio>
#include <cmath>

class Vector
{
  public:
    union
    {
        double v[3];
        struct
        {
            double x, y, z;
        };
    };
  public:
    Vector() : v() {}
    Vector(const double x, const double y, const double z) : v{x, y, z} {}

    void Set(const double x, const double y, const double z) { v[0] = x; v[1] = y; v[2] = z; }
    void Zero() { v[0] = 0.0; v[1] = 0.0; v[2] = 0.0; }

    const Vector operator+(const Vector& b) const { return Vector(v[0] + b.v[0], v[1] + b.v[1], v[2] + b.v[2]); }
    const Vector operator-(const Vector& b) const { return Vector(v[0] - b.v[0], v[1] - b.v[1], v[2] - b.v[2]); }
    Vector& operator+=(const Vector& b) { v[0] += b.v[0]; v[1] += b.v[1]; v[2] += b.v[2]; return *this; }
    Vector& operator-=(const Vector& b) { v[0] -= b.v[0]; v[1] -= b.v[1]; v[2] -= b.v[2]; return *this; }

    const Vector operator*(const double b) const { return Vector(v[0] * b, v[1] * b, v[2] * b); }
    const Vector operator/(const double b) const { return Vector(v[0] / b, v[1] / b, v[2] / b); }
    Vector& operator*=(const double b) { v[0] *= b; v[1] *= b; v[2] *= b; return *this; }
    Vector& operator/=(const double b) { v[0] /= b; v[1] /= b; v[2] /= b; return *this; }

    Vector CrossProduct(const Vector& b) const { return Vector(v[1] * b.v[2] - v[2] * b.v[1], v[2] * b.v[0] - v[0] * b.v[2], v[0] * b.v[1] - v[1] * b.v[0]); }
    double DotProduct(const Vector& b) const { return (v[0] * b.v[0] + v[1] * b.v[1] + v[2] * b.v[2]); }
    double Mag() const { return sqrt(Mag2()); }
    double Mag2() const { return (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]); }
};

static const Vector operator*(const double a, const Vector& b) { return b * a; }

#if 0
static_assert(std::is_pod<Vector>::value == true, "PODじゃないよ");
static_assert(std::is_trivial<Vector>::value == true, "trivialじゃないよ");
#endif
static_assert(std::is_trivially_copyable<Vector>::value == true, "trivially_copyableじゃないよ");
static_assert(std::is_standard_layout<Vector>::value == true, "standard_layoutじゃないよ");


static_assert(sizeof(double) == 8, "doubleが８bytesじゃーないよ");
static_assert(sizeof(Vector) == sizeof(double) * 3, "paddingされてるよ");

static_assert(offsetof(Vector, v[0]) == sizeof(double) * 0, "メンバー変数がずれてるよ");
static_assert(offsetof(Vector, v[1]) == sizeof(double) * 1, "メンバー変数がずれてるよ");
static_assert(offsetof(Vector, v[2]) == sizeof(double) * 2, "メンバー変数がずれてるよ");
static_assert(offsetof(Vector, x) == offsetof(Vector, v[0]), "メンバー変数がずれてるよ");
static_assert(offsetof(Vector, y) == offsetof(Vector, v[1]), "メンバー変数がずれてるよ");
static_assert(offsetof(Vector, z) == offsetof(Vector, v[2]), "メンバー変数がずれてるよ");

#endif
