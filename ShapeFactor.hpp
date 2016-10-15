#ifndef SIMPLE_PIC_SHAPEFACTOR_HPP
#define SIMPLE_PIC_SHAPEFACTOR_HPP

#include <cmath>

enum Shape
{
    // Shape
    NGP = 0,  // 0th
    CIC = 1,  // 1st
    TSC = 2,  // 2nd

    Spline = TSC,
};

template<Shape, int N>
class ShapeFactor
{
  public:
    void operator()(double (&sf)[N], const double rx);
    void operator()(double (&sf)[N], const double rx, const int shift);

};

template<>inline void ShapeFactor<NGP, 3>::operator()(double (&sf)[3], const double rx)
{
    // fix me
}

template<> inline void ShapeFactor<NGP, 5>::operator()(double (&sf)[5], const double rx)
{
    // fix me
}

template<> inline void ShapeFactor<CIC, 3>::operator()(double (&sf)[3], const double rx)
{
    int i = int(rx);
    double dx = rx - (i + .5);
    double ax = fabs(dx);

    sf[0] = 0.0;
    sf[1] = 1.0 - ax;
    sf[2] = 0.0;

    sf[1+(int)copysign(1, dx)] = ax;
}

template<> inline void ShapeFactor<CIC, 5>::operator()(double (&sf)[5], const double rx)
{
    int i = int(rx);
    double dx = rx - (i + .5);
    double ax = fabs(dx);

    sf[0] = 0.0;
    sf[1] = 0.0;
    sf[2] = 1.0 - ax;
    sf[3] = 0.0;
    sf[4] = 0.0;

    sf[2+(int)copysign(1, dx)] = ax;
}

template<> inline void ShapeFactor<CIC, 5>::operator()(double (&sf)[5], const double rx, const int shift)
{
    int i = int(rx);
    double dx = rx - (i + .5);
    double ax = fabs(dx);

    sf[0] = 0.0;
    sf[1] = 0.0;
    sf[2] = 0.0;
    sf[3] = 0.0;
    sf[4] = 0.0;

    sf[2 + shift] = 1.0 - ax;
    sf[2 + shift + (int)copysign(1, dx)] = ax;
}

template<> inline void ShapeFactor<TSC, 3>::operator()(double (&sf)[3], const double rx)
{
    int i = int(rx);
    double dx = rx - (i + .5);

    sf[0] = 0.5 * (0.5-dx)*(0.5-dx);
    sf[1] = 0.75 - dx*dx;
    sf[2] = 0.5 * (0.5+dx)*(0.5+dx);
}

template<> inline void ShapeFactor<TSC, 5>::operator()(double (&sf)[5], const double rx)
{
    int i = int(rx);
    double dx = rx - (i + .5);

    sf[0] = 0.0;
    sf[1] = 0.5 * (0.5-dx)*(0.5-dx);
    sf[2] = 0.75 - dx*dx;
    sf[3] = 0.5 * (0.5+dx)*(0.5+dx);
    sf[4] = 0.0;
}


template<> inline void ShapeFactor<TSC, 5>::operator()(double (&sf)[5], const double rx, const int shift)
{
    int i = int(rx);
    double dx = rx - (i + .5);

    sf[0] = 0.0;
    sf[4] = 0.0;

    sf[2 - shift] = 0.0;

    sf[1 + shift] = 0.5 * (0.5-dx)*(0.5-dx);
    sf[2 + shift] = 0.75 - dx*dx;
    sf[3 + shift] = 0.5 * (0.5+dx)*(0.5+dx);
}


#endif
