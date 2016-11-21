#ifndef SIMPLE_PIC_PLASMA_HPP
#define SIMPLE_PIC_PLASMA_HPP

#include <vector>
#include <algorithm>

#include "Param.hpp"

#include "Particle.hpp"
#include "Vector.hpp"

class Plasma
{
  public:
    double m;
    double q;
    std::vector<Particle> p;

    std::vector<Particle> mpiParticlesX0;
    std::vector<Particle> mpiParticlesX1;

    std::vector<double> mpiSendBuf;
    std::vector<double> mpiRecvBuf;

  public:
    Plasma() {}

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
        for (unsigned long int n = 0; n < p.size(); ++n) { p[n].r += p[n].v; }
    }
    
    void BoundaryR()
    {
        for (unsigned long int n = 0; n < p.size(); ++n)
        {
	        if (p[n].r.y >= Y1) { p[n].r.y -= LY0; }
	        if (p[n].r.y <  Y0) { p[n].r.y += LY0; }
            if (p[n].r.z >= Z1) { p[n].r.z -= LZ0; }
	        if (p[n].r.z <  Z0) { p[n].r.z += LZ0; }
        }

        auto MPIMoveParticle = [this](std::vector<Particle>& p, const bool reverse = false)
        {
            long sendSize = p.size();
            long recvSize = 0;

            int destRank, srcRank;

            if (reverse != true)
                MPI_Cart_shift(comm, 0,  1, &srcRank, &destRank);
            else
                MPI_Cart_shift(comm, 0, -1, &srcRank, &destRank);
            
            MPI_Status status;
            MPI_Sendrecv(&sendSize, 1, MPI_LONG, destRank, 201,
                         &recvSize, 1, MPI_LONG,  srcRank, 201,
                         comm, &status);

            mpiRecvBuf.reserve(recvSize * 7);
            mpiRecvBuf.resize(recvSize * 7);

            MPI_Sendrecv(&p[0], sendSize * 7, MPI_DOUBLE, destRank, 202,
                &mpiRecvBuf[0], recvSize * 7, MPI_DOUBLE,  srcRank, 202,
                comm, &status);

            p.reserve(recvSize);
            p.resize(recvSize);

            memcpy(&p[0], &mpiRecvBuf[0], sizeof(Particle) * recvSize);
        };

        unsigned long int pSize = p.size();
        for (unsigned long int n = 0; n < pSize; ++n)
        {
            if (p[n].r.x >= X1)
            {
                p[n].r.x -= LX0;
                mpiParticlesX1.push_back(p[n]);
                p[n].id = -1;
            }
	        if (p[n].r.x <  X0)
            {
                p[n].r.x += LX0;
                mpiParticlesX0.push_back(p[n]);
                p[n].id = -1;
            }
        }

        MPIMoveParticle(mpiParticlesX1, false);
        MPIMoveParticle(mpiParticlesX0,  true);

        p.erase(std::remove_if(p.begin(), p.end(), [](const Particle& p){ return p.id == -1; }), p.end());

        pSize = mpiParticlesX0.size();
        for (unsigned long int n = 0; n < pSize; ++n)
        {
            p.push_back(mpiParticlesX0[n]);
        }

        pSize = mpiParticlesX1.size();
        for (unsigned long int n = 0; n < pSize; ++n)
        {
            p.push_back(mpiParticlesX1[n]);
        }


        mpiParticlesX1.clear();
        mpiParticlesX0.clear();
    }
};

#endif

