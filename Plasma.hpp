#ifndef SIMPLE_PIC_PLASMA_HPP
#define SIMPLE_PIC_PLASMA_HPP

#include <vector>

#include "Param.hpp"

#include "Particle.hpp"
#include "Vector.hpp"

class Plasma
{
  public:
    double m;
    double q;
    std::vector<Particle> p;

  public:
    Plasma() {}
    Plasma(double m, double q) : m(m), q(q) {}
    Plasma(double m, double q, std::vector<Particle> p) : m(m), q(q) { this->p = std::move(p); }

    void SetParticle(std::vector<Particle> p) { this->p = std::move(p); }

    void UpdateR()
    {
        for (long int n = 0; n < p.size(); ++n) { p[n].r += p[n].v; }
    }
    
    void BoundaryR()
    {
        for (long int n = 0; n < p.size(); ++n)
        {
            if (p[n].r.x >= X1) { p[n].r.x -= LX0; }
	        if (p[n].r.x <  X0) { p[n].r.x += LX0; }
	        if (p[n].r.y >= Y1) { p[n].r.y -= LY0; }
	        if (p[n].r.y <  Y0) { p[n].r.y += LY0; }
            if (p[n].r.z >= Z1) { p[n].r.z -= LZ0; }
	        if (p[n].r.z <  Z0) { p[n].r.z += LZ0; }
        }
    }
};

#endif

