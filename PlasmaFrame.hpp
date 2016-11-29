#ifndef SIMPLE_PIC_PLASMA_FRAME_HPP
#define SIMPLE_PIC_PLASMA_FRAME_HPP

#include <vector>
#include <algorithm>

#include "Param.hpp"
#include "PlasmaBase.hpp"

template<Boundary>
class PlasmaFrame : public PlasmaBase
{
  private:
    std::vector<Particle> pBuf;
    std::vector<Particle> pX0Buf;
    std::vector<Particle> pX1Buf;

  public:
    PlasmaFrame() {}
    ~PlasmaFrame() {}

    void BoundaryR();

  private:
    void MPIMoveParticle(std::vector<Particle>& p, const bool reverse = false)
    {        
        long sendSize = p.size();
        long recvSize = -1;

        int destRank, srcRank;

        if (reverse != true)
            MPI_Cart_shift(comm, 0,  1, &srcRank, &destRank);
        else
            MPI_Cart_shift(comm, 0, -1, &srcRank, &destRank);

        MPI_Status status;
        MPI_Sendrecv(&sendSize, 1, MPI_LONG, destRank, 301,
                     &recvSize, 1, MPI_LONG,  srcRank, 301,
                     comm, &status);

        pBuf.resize(recvSize);

        MPI_Sendrecv(&p[0], sendSize * 7, MPI_DOUBLE, destRank, 302,
                  &pBuf[0], recvSize * 7, MPI_DOUBLE,  srcRank, 302,
                  comm, &status);

        //p.swap(pBuf);
        
        p.resize(recvSize);
        memcpy(&p[0], &pBuf[0], sizeof(Particle) * recvSize);
    }
};

template<> void PlasmaFrame<Boundary::Periodic>::BoundaryR()
{
    unsigned long pX0, pX1;
    auto pSize = p.size();

    pX0 = 0;
    pX1 = 0;

    for (unsigned long n = 0; n < pSize; ++n)
    {
        if (p[n].r.x >= X1) ++pX1;
        if (p[n].r.x <  X0) ++pX0;
    }

    pX1Buf.resize(pX1);
    pX0Buf.resize(pX0);

    pX0 = 0;
    pX1 = 0;

    for (unsigned long n = 0; n < pSize; ++n)
    {
        if (p[n].r.x >= X1)
        {
            p[n].r.x -= LX0;
            pX1Buf[pX1++] = p[n];
            p[n].id = -1;
        }
        if (p[n].r.x <  X0)
        {
            p[n].r.x += LX0;
            pX0Buf[pX0++] = p[n];
            p[n].id = -1;
        }
    }

    p.erase(
        std::remove_if(p.begin(), p.end(),
            [](const Particle& p)
            {
                return p.id == -1;
            })
        , p.end());

    MPIMoveParticle(pX1Buf, false);
    MPIMoveParticle(pX0Buf,  true);
    
    if (pX1Buf.size() > 0)
    {
        pSize = p.size();
        p.resize(pSize + pX1Buf.size());
        memcpy(&p[pSize], &pX1Buf[0], sizeof(Particle) * pX1Buf.size());
    }

    if (pX0Buf.size() > 0)
    {
        pSize = p.size();
        p.resize(pSize + pX0Buf.size());
        memcpy(&p[pSize], &pX0Buf[0], sizeof(Particle) * pX0Buf.size());
    }
    
    pX1Buf.clear();
    pX0Buf.clear();

    pSize = p.size();
    for (unsigned long n = 0; n < pSize; ++n)
    {
        if (p[n].r.y >= Y1) { p[n].r.y -= LY0; }
        if (p[n].r.y <  Y0) { p[n].r.y += LY0; }
        if (p[n].r.z >= Z1) { p[n].r.z -= LZ0; }
        if (p[n].r.z <  Z0) { p[n].r.z += LZ0; }
    }
}

#endif

