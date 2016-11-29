#ifndef SIMPLE_PIC_PLASMA_BASE_HPP
#define SIMPLE_PIC_PLASMA_BASE_HPP

#include <vector>
#include <algorithm>

#include "Param.hpp"

class PlasmaBase
{
  public:
    double m;
    double q;
    std::vector<Particle> p;

  protected:
    PlasmaBase() {}
    
  public:
    ~PlasmaBase() {}
    
    void Sort()
    {
        std::sort(p.begin(), p.end(),
                [](const Particle& a, const Particle& b)
                {
                    return int(a.r.z) + int(a.r.y)*LZ + int(a.r.x)*LY*LZ
                         < int(b.r.z) + int(b.r.y)*LZ + int(b.r.x)*LY*LZ;
                });
    }

    void UpdateR()
    {
        for (unsigned long int n = 0; n < p.size(); ++n)
        {
            p[n].r += p[n].v;
        }
    }
};

#endif

