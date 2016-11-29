#ifndef SIMPLE_PIC_SOLVER_FRAME_HPP
#define SIMPLE_PIC_SOLVER_FRAME_HPP

#include <vector>
#include <cstdlib>

#include "Param.hpp"
#include "SolverBase.hpp"

template<Boundary, Shape S>
class SolverFrame : public SolverBase<S>
{
  private:
    std::vector<Vector> fBuf;

  public:
    SolverFrame()
    {
        try
        {
            fBuf.resize(LY * LZ);
        }
        catch(std::bad_alloc)
        {
            fprintf(stderr, "Error: bad_alloc : %s\n", __FILE__); 
            abort();
        }
    }
    
    ~SolverFrame()
    {
    }

    void BoundaryJ();

  private:
    void MPIAddField(Vector***& v, const int dstX, const int srcX, const bool reverse = false)
    {
        int destRank, srcRank;

        if (reverse != true)
            MPI_Cart_shift(comm, 0,  1, &srcRank, &destRank);
        else
            MPI_Cart_shift(comm, 0, -1, &srcRank, &destRank);

        MPI_Status status;
        MPI_Sendrecv(&v[srcX][0][0], LY * LZ * 3, MPI_DOUBLE, destRank, 201,
                           &fBuf[0], LY * LZ * 3, MPI_DOUBLE,  srcRank, 201,
                           comm, &status);

        for (int j = 0; j < LY; ++j)
        for (int k = 0; k < LZ; ++k)
        {
            v[dstX][j][k] += fBuf[j * LZ + k];
        }
    }
};

template<> void SolverFrame<Boundary::Periodic, Shape::TSC>::BoundaryJ()
{
    // X
    {
        MPIAddField(J, X0+0, X1+0, false);
        MPIAddField(J, X0+1, X1+1, false);
        MPIAddField(J, X0+2, X1+2, false);
        MPIAddField(J, X1-1, X0-1, true);
        MPIAddField(J, X1-2, X0-2, true);
    }
    // Y
    for (int i = 0; i < LX; ++i)
    for (int k = 0; k < LZ; ++k)
    {
        J[i][Y0+0][k] += J[i][Y1+0][k];
        J[i][Y0+1][k] += J[i][Y1+1][k];
        J[i][Y0+2][k] += J[i][Y1+2][k];
        J[i][Y1-1][k] += J[i][Y0-1][k];
        J[i][Y1-2][k] += J[i][Y0-2][k];
    }
    // Z
    for (int i = 0; i < LX; ++i)
    for (int j = 0; j < LY; ++j)
    {
        J[i][j][Z0+0] += J[i][j][Z1+0];
        J[i][j][Z0+1] += J[i][j][Z1+1];
        J[i][j][Z0+2] += J[i][j][Z1+2];
        J[i][j][Z1-1] += J[i][j][Z0-1];
        J[i][j][Z1-2] += J[i][j][Z0-2];
    }
}

template<> void SolverFrame<Boundary::Periodic, Shape::CIC>::BoundaryJ()
{
    // X
    {
        MPIAddField(J, X0+0, X1+0, false);
        MPIAddField(J, X0+1, X1+1, false);
        MPIAddField(J, X0+2, X1+2, false);
        MPIAddField(J, X1-1, X0-1, true);
        MPIAddField(J, X1-2, X0-2, true); // ?
    }
    // Y
    for (int i = 0; i < LX; ++i)
    for (int k = 0; k < LZ; ++k)
    {
        J[i][Y0+0][k] += J[i][Y1+0][k];
        J[i][Y0+1][k] += J[i][Y1+1][k];
        J[i][Y0+2][k] += J[i][Y1+2][k];
        J[i][Y1-1][k] += J[i][Y0-1][k];
        J[i][Y1-2][k] += J[i][Y0-2][k]; // ?
    }
    // Z
    for (int i = 0; i < LX; ++i)
    for (int j = 0; j < LY; ++j)
    {
        J[i][j][Z0+0] += J[i][j][Z1+0];
        J[i][j][Z0+1] += J[i][j][Z1+1];
        J[i][j][Z0+2] += J[i][j][Z1+2];
        J[i][j][Z1-1] += J[i][j][Z0-1];
        J[i][j][Z1-2] += J[i][j][Z0-2]; // ?
    }
}

template<> void SolverFrame<Boundary::Periodic, Shape::NGP>::BoundaryJ()
{
    // X
    {
        MPIAddField(J, X0+0, X1+0, false);
        MPIAddField(J, X0+1, X1+1, false);
        MPIAddField(J, X0+2, X1+2, false);
        MPIAddField(J, X1-1, X0-1, true);
        MPIAddField(J, X1-2, X0-2, true); // ?
    }
    // Y
    for (int i = 0; i < LX; ++i)
    for (int k = 0; k < LZ; ++k)
    {
        J[i][Y0+0][k] += J[i][Y1+0][k];
        J[i][Y0+1][k] += J[i][Y1+1][k];
        J[i][Y0+2][k] += J[i][Y1+2][k];
        J[i][Y1-1][k] += J[i][Y0-1][k];
        J[i][Y1-2][k] += J[i][Y0-2][k]; // ?
    }
    // Z
    for (int i = 0; i < LX; ++i)
    for (int j = 0; j < LY; ++j)
    {
        J[i][j][Z0+0] += J[i][j][Z1+0];
        J[i][j][Z0+1] += J[i][j][Z1+1];
        J[i][j][Z0+2] += J[i][j][Z1+2];
        J[i][j][Z1-1] += J[i][j][Z0-1];
        J[i][j][Z1-2] += J[i][j][Z0-2]; // ?
    }
}

#endif

