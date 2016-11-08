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
            unsigned long int sendSize = p.size();
            unsigned long int recvSize = 0;

            mpiSendBuf.resize(7 * sendSize);

            unsigned long int idx;
            
            idx = 0;
            for (unsigned long int n = 0; n < sendSize; ++n)
            {
                memcpy(&mpiSendBuf[idx++], &p[n].id, sizeof(double));
                mpiSendBuf[idx++] = p[n].r.x;
                mpiSendBuf[idx++] = p[n].r.y;
                mpiSendBuf[idx++] = p[n].r.z;
                mpiSendBuf[idx++] = p[n].v.x;
                mpiSendBuf[idx++] = p[n].v.y;
                mpiSendBuf[idx++] = p[n].v.z;
            }

            int forward = (MPI::COMM_WORLD.Get_rank() + 1) % MPI::COMM_WORLD.Get_size();
            int backward = MPI::COMM_WORLD.Get_rank() - 1;
            if (backward < 0) backward = MPI::COMM_WORLD.Get_size() - 1;

            MPI::Status status;
            
            if (reverse != true)
                MPI::COMM_WORLD.Sendrecv(
                        &sendSize, 1, MPI::INT,  forward, 1,
                        &recvSize, 1, MPI::INT, backward, 1, status);
            else
                MPI::COMM_WORLD.Sendrecv(
                        &sendSize, 1, MPI::INT, backward, 1,
                        &recvSize, 1, MPI::INT,  forward, 1, status);


            mpiRecvBuf.resize(7 * recvSize);

            if (reverse != true)
                MPI::COMM_WORLD.Sendrecv(
                        &mpiSendBuf[0], sendSize * 7, MPI::DOUBLE,  forward, 1,
                        &mpiRecvBuf[0], recvSize * 7, MPI::DOUBLE, backward, 1, status);
            else
                MPI::COMM_WORLD.Sendrecv(
                        &mpiSendBuf[0], sendSize * 7, MPI::DOUBLE, backward, 1,
                        &mpiRecvBuf[0], recvSize * 7, MPI::DOUBLE,  forward, 1, status);

            p.resize(recvSize);

            idx = 0;
            for (unsigned long int n = 0; n < recvSize; ++n)
            {
                memcpy(&p[n].id, &mpiRecvBuf[idx++], sizeof(double));
                p[n].r.x = mpiRecvBuf[idx++];
                p[n].r.y = mpiRecvBuf[idx++];
                p[n].r.z = mpiRecvBuf[idx++];
                p[n].v.x = mpiRecvBuf[idx++];
                p[n].v.y = mpiRecvBuf[idx++];
                p[n].v.z = mpiRecvBuf[idx++];
            }
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

