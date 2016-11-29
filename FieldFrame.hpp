#ifndef SIMPLE_PIC_FIELD_FRAME_HPP
#define SIMPLE_PIC_FIELD_FRAME_HPP

#include <cstdlib>

#include "Param.hpp"
#include "FieldBase.hpp"

template<Boundary>
class FieldFrame : public FieldBase
{
  public:

  public:
    FieldFrame() {}
    ~FieldFrame() {}

    void BoundaryB();
    void BoundaryE();

  private:    
    void MPICopyField(Vector***& v, const int dstX, const int srcX, const bool reverse = false)
    {
        int destRank, srcRank;

        if (reverse != true)
            MPI_Cart_shift(comm, 0,  1, &srcRank, &destRank);
        else
            MPI_Cart_shift(comm, 0, -1, &srcRank, &destRank);

        MPI_Status status;
        MPI_Sendrecv(&v[srcX][0][0], LY * LZ * 3, MPI_DOUBLE, destRank, 101,
                     &v[dstX][0][0], LY * LZ * 3, MPI_DOUBLE,  srcRank, 101,
                           comm, &status);
    }
};

template<> void FieldFrame<Boundary::Periodic>::BoundaryB()
{
    // X
    {
        MPICopyField(B, X1+0, X0+0, true);
        MPICopyField(B, X1+1, X0+1, true);
        MPICopyField(B, X1+2, X0+2, true);
        MPICopyField(B, X0-1, X1-1, false);
        MPICopyField(B, X0-2, X1-2, false);
    }
    // Y
    for (int i = 0; i < LX; ++i)
    for (int k = 0; k < LZ; ++k)
    {
        B[i][Y1+0][k] = B[i][Y0+0][k];
        B[i][Y1+1][k] = B[i][Y0+1][k];
        B[i][Y1+2][k] = B[i][Y0+2][k];
        B[i][Y0-1][k] = B[i][Y1-1][k];
        B[i][Y0-2][k] = B[i][Y1-2][k];
    }
    // Z
    for (int i = 0; i < LX; ++i)
    for (int j = 0; j < LY; ++j)
    {
        B[i][j][Z1+0] = B[i][j][Z0+0];
        B[i][j][Z1+1] = B[i][j][Z0+1];
        B[i][j][Z1+2] = B[i][j][Z0+2];
        B[i][j][Z0-1] = B[i][j][Z1-1];
        B[i][j][Z0-2] = B[i][j][Z1-2];
    }
}

template<> void FieldFrame<Boundary::Periodic>::BoundaryE()
{
    // X
    {
        MPICopyField(E, X1+0, X0+0, true);
        MPICopyField(E, X1+1, X0+1, true);
        MPICopyField(E, X1+2, X0+2, true);
        MPICopyField(E, X0-1, X1-1, false);
        MPICopyField(E, X0-2, X1-2, false);
    }
    // Y
    for (int i = 0; i < LX; ++i)
    for (int k = 0; k < LZ; ++k)
    {
        E[i][Y1+0][k] = E[i][Y0+0][k];
        E[i][Y1+1][k] = E[i][Y0+1][k];
        E[i][Y1+2][k] = E[i][Y0+2][k];
        E[i][Y0-1][k] = E[i][Y1-1][k];
        E[i][Y0-2][k] = E[i][Y1-2][k];
    }
    // Z
    for (int i = 0; i < LX; ++i)
    for (int j = 0; j < LY; ++j)
    {
        E[i][j][Z1+0] = E[i][j][Z0+0];
        E[i][j][Z1+1] = E[i][j][Z0+1];
        E[i][j][Z1+2] = E[i][j][Z0+2];
        E[i][j][Z0-1] = E[i][j][Z1-1];
        E[i][j][Z0-2] = E[i][j][Z1-2];
    }
}

template<>
class FieldFrame<Boundary::PML> : public FieldBase
{
  public:
    double* sigmaL0E;
    double* sigmaL1E;
    double* sigmaL0E_;
    double* sigmaL1E_;

    Vector*** pmlY0Ex;
    Vector*** pmlY0Ey;
    Vector*** pmlY0Ez;
    Vector*** pmlY1Ex;
    Vector*** pmlY1Ey;
    Vector*** pmlY1Ez;
    
    Vector*** pmlZ0Ex;
    Vector*** pmlZ0Ey;
    Vector*** pmlZ0Ez;
    Vector*** pmlZ1Ex;
    Vector*** pmlZ1Ey;
    Vector*** pmlZ1Ez;

    double* sigmaL0B;
    double* sigmaL1B;
    double* sigmaL0B_;
    double* sigmaL1B_;

    Vector*** pmlY0Bx;
    Vector*** pmlY0By;
    Vector*** pmlY0Bz;
    Vector*** pmlY1Bx;
    Vector*** pmlY1By;
    Vector*** pmlY1Bz;
    
    Vector*** pmlZ0Bx;
    Vector*** pmlZ0By;
    Vector*** pmlZ0Bz;
    Vector*** pmlZ1Bx;
    Vector*** pmlZ1By;
    Vector*** pmlZ1Bz;

  public:
    FieldFrame()
    {
        try
        {            
            sigmaL0E  = new double [L];
            sigmaL0B  = new double [L];
            sigmaL0E_ = new double [L];
            sigmaL0B_ = new double [L];

            sigmaL1E  = new double [L];
            sigmaL1B  = new double [L];
            sigmaL1E_ = new double [L];
            sigmaL1B_ = new double [L];

            pmlY0Ex = CreateArray(LX, L, LZ);
            pmlY0Ey = CreateArray(LX, L, LZ);
            pmlY0Ez = CreateArray(LX, L, LZ);
            pmlY1Ex = CreateArray(LX, L, LZ);
            pmlY1Ey = CreateArray(LX, L, LZ);
            pmlY1Ez = CreateArray(LX, L, LZ);

            pmlZ0Ex = CreateArray(LX, LY, L);
            pmlZ0Ey = CreateArray(LX, LY, L);
            pmlZ0Ez = CreateArray(LX, LY, L);
            pmlZ1Ex = CreateArray(LX, LY, L);
            pmlZ1Ey = CreateArray(LX, LY, L);
            pmlZ1Ez = CreateArray(LX, LY, L);

            pmlY0Bx = CreateArray(LX, L, LZ);
            pmlY0By = CreateArray(LX, L, LZ);
            pmlY0Bz = CreateArray(LX, L, LZ);
            pmlY1Bx = CreateArray(LX, L, LZ);
            pmlY1By = CreateArray(LX, L, LZ);
            pmlY1Bz = CreateArray(LX, L, LZ);
            
            pmlZ0Bx = CreateArray(LX, LY, L);
            pmlZ0By = CreateArray(LX, LY, L);
            pmlZ0Bz = CreateArray(LX, LY, L);
            pmlZ1Bx = CreateArray(LX, LY, L);
            pmlZ1By = CreateArray(LX, LY, L);
            pmlZ1Bz = CreateArray(LX, LY, L);
        }
        catch(std::bad_alloc)
        {
            fprintf(stderr, "Errer: bad_alloc : %s\n", __FILE__);            

            delete[] sigmaL0E;
            delete[] sigmaL0B;
            delete[] sigmaL0E_;
            delete[] sigmaL0B_;

            delete[] sigmaL1E;
            delete[] sigmaL1B;
            delete[] sigmaL1E_;
            delete[] sigmaL1B_;

            DeleteArray(pmlY0Ex);
            DeleteArray(pmlY0Ey);
            DeleteArray(pmlY0Ez);
            DeleteArray(pmlY1Ex);
            DeleteArray(pmlY1Ey);
            DeleteArray(pmlY1Ez);

            DeleteArray(pmlZ0Ex);
            DeleteArray(pmlZ0Ey);
            DeleteArray(pmlZ0Ez);
            DeleteArray(pmlZ1Ex);
            DeleteArray(pmlZ1Ey);
            DeleteArray(pmlZ1Ez);

            DeleteArray(pmlY0Bx);
            DeleteArray(pmlY0By);
            DeleteArray(pmlY0Bz);
            DeleteArray(pmlY1Bx);
            DeleteArray(pmlY1By);
            DeleteArray(pmlY1Bz);

            DeleteArray(pmlZ0Bx);
            DeleteArray(pmlZ0By);
            DeleteArray(pmlZ0Bz);
            DeleteArray(pmlZ1Bx);
            DeleteArray(pmlZ1By);
            DeleteArray(pmlZ1Bz);

            abort();
        }

        for (int i = 0; i < LX; ++i)
	    for (int j = 0; j < L ; ++j)
	    for (int k = 0; k < LZ; ++k)
        {
            pmlY0Ex[i][j][k].Zero();
            pmlY0Ey[i][j][k].Zero();
            pmlY0Ez[i][j][k].Zero();
            pmlY1Ex[i][j][k].Zero();
            pmlY1Ey[i][j][k].Zero();
            pmlY1Ez[i][j][k].Zero();

            pmlY0Bx[i][j][k].Zero();
            pmlY0By[i][j][k].Zero();
            pmlY0Bz[i][j][k].Zero();
            pmlY1Bx[i][j][k].Zero();
            pmlY1By[i][j][k].Zero();
            pmlY1Bz[i][j][k].Zero();
        }

        for (int i = 0; i < LX; ++i)
	    for (int j = 0; j < LY; ++j)
	    for (int k = 0; k < L ; ++k)
        {
            pmlZ0Ex[i][j][k].Zero();
            pmlZ0Ey[i][j][k].Zero();
            pmlZ0Ez[i][j][k].Zero();
            pmlZ1Ex[i][j][k].Zero();
            pmlZ1Ey[i][j][k].Zero();
            pmlZ1Ez[i][j][k].Zero();

            pmlZ0Bx[i][j][k].Zero();
            pmlZ0By[i][j][k].Zero();
            pmlZ0Bz[i][j][k].Zero();
            pmlZ1Bx[i][j][k].Zero();
            pmlZ1By[i][j][k].Zero();
            pmlZ1Bz[i][j][k].Zero();
        }

        for (int n = 0; n < L; ++n)
        {
            double sigL0E = SIGMA_MAX * std::pow((n + 0.5) / L, M);
            double sigL0B = SIGMA_MAX * std::pow((n + 1.0) / L, M);

            double sigL1E = SIGMA_MAX * std::pow((n + 0.5) / L, M);
            double sigL1B = SIGMA_MAX * std::pow((n + 0.0) / L, M);

                        
            sigmaL0E[n]  = 1.0 / (1.0 + sigL0E * 0.5);
            sigmaL0E_[n] = (1.0 - sigL0E * 0.5);
            sigmaL1E[n]  = 1.0 / (1.0 + sigL1E * 0.5);
            sigmaL1E_[n] = (1.0 - sigL1E * 0.5);

            sigmaL0B[n]  = 1.0 / (1.0 + sigL0B * 0.25);
            sigmaL0B_[n] = (1.0 - sigL0B * 0.25);
            sigmaL1B[n]  = 1.0 / (1.0 + sigL1B * 0.25);
            sigmaL1B_[n] = (1.0 - sigL1B * 0.25);
        }
    }
    
    ~FieldFrame()
    {
        delete[] sigmaL0E;
        delete[] sigmaL0B;
        delete[] sigmaL0E_;
        delete[] sigmaL0B_;

        delete[] sigmaL1E;
        delete[] sigmaL1B;
        delete[] sigmaL1E_;
        delete[] sigmaL1B_;

        DeleteArray(pmlY0Ex);
        DeleteArray(pmlY0Ey);
        DeleteArray(pmlY0Ez);
        DeleteArray(pmlY1Ex);
        DeleteArray(pmlY1Ey);
        DeleteArray(pmlY1Ez);

        DeleteArray(pmlZ0Ex);
        DeleteArray(pmlZ0Ey);
        DeleteArray(pmlZ0Ez);
        DeleteArray(pmlZ1Ex);
        DeleteArray(pmlZ1Ey);
        DeleteArray(pmlZ1Ez);

        DeleteArray(pmlY0Bx);
        DeleteArray(pmlY0By);
        DeleteArray(pmlY0Bz);
        DeleteArray(pmlY1Bx);
        DeleteArray(pmlY1By);
        DeleteArray(pmlY1Bz);

        DeleteArray(pmlZ0Bx);
        DeleteArray(pmlZ0By);
        DeleteArray(pmlZ0Bz);
        DeleteArray(pmlZ1Bx);
        DeleteArray(pmlZ1By);
        DeleteArray(pmlZ1Bz);
    }

    void BoundaryB()
    {
        // PML - Y0 [0 Y0)
        for (int j = 1; j < Y0; ++j)
        for (int i = X0; i < X1; ++i)
        for (int k = Z0; k < Z1; ++k)
        {
            int l = L - j;
            
            pmlY0By[i][l][k].x = sigmaL0B[l] * (sigmaL0B_[l] * pmlY0By[i][l][k].x - 0.5 * (E[i][j][k].z - E[i][j-1][k].z));
            pmlY0Bz[i][l][k].x =               (               pmlY0Bz[i][l][k].x + 0.5 * (E[i][j][k].y - E[i][j][k-1].y));

            pmlY0Bz[i][l][k].y =               (               pmlY0Bz[i][l][k].y - 0.5 * (E[i][j][k].x - E[i][j][k-1].x));
            pmlY0Bx[i][l][k].y =               (               pmlY0Bx[i][l][k].y + 0.5 * (E[i][j][k].z - E[i-1][j][k].z));
        
            pmlY0Bx[i][l][k].z =               (               pmlY0Bx[i][l][k].z - 0.5 * (E[i][j][k].y - E[i-1][j][k].y));
            pmlY0By[i][l][k].z = sigmaL0B[l] * (sigmaL0B_[l] * pmlY0By[i][l][k].z + 0.5 * (E[i][j][k].x - E[i][j-1][k].x));

            B[i][j][k].x = pmlY0By[i][l][k].x + pmlY0Bz[i][l][k].x;
            B[i][j][k].y = pmlY0Bz[i][l][k].y + pmlY0Bx[i][l][k].y;
            B[i][j][k].z = pmlY0Bx[i][l][k].z + pmlY0By[i][l][k].z;
        }

        // PML - Y1 [Y1 LY)
        for (int j = Y1; j < LY-1; ++j)
        for (int i = X0; i < X1; ++i)
        for (int k = Z0; k < Z1; ++k)
        {
            int l = j - Y1;
            
            pmlY1By[i][l][k].x = sigmaL1B[l] * (sigmaL1B_[l] * pmlY1By[i][l][k].x - 0.5 * (E[i][j][k].z - E[i][j-1][k].z));
            pmlY1Bz[i][l][k].x =               (               pmlY1Bz[i][l][k].x + 0.5 * (E[i][j][k].y - E[i][j][k-1].y));

            pmlY1Bz[i][l][k].y =               (               pmlY1Bz[i][l][k].y - 0.5 * (E[i][j][k].x - E[i][j][k-1].x));
            pmlY1Bx[i][l][k].y =               (               pmlY1Bx[i][l][k].y + 0.5 * (E[i][j][k].z - E[i-1][j][k].z));
        
            pmlY1Bx[i][l][k].z =               (               pmlY1Bx[i][l][k].z - 0.5 * (E[i][j][k].y - E[i-1][j][k].y));
            pmlY1By[i][l][k].z = sigmaL1B[l] * (sigmaL1B_[l] * pmlY1By[i][l][k].z + 0.5 * (E[i][j][k].x - E[i][j-1][k].x));

            B[i][j][k].x = pmlY1By[i][l][k].x + pmlY1Bz[i][l][k].x;
            B[i][j][k].y = pmlY1Bz[i][l][k].y + pmlY1Bx[i][l][k].y;
            B[i][j][k].z = pmlY1Bx[i][l][k].z + pmlY1By[i][l][k].z;
        }

        // PML - Z0 [1 Z0) ... 
        // 
        // x: |0|1|2|3|4||X0~
        //    |+|*|*|*|*||S ~  (+: perfect, *: PML, S: simulate box)
        // l: |4|3|2|1|0||  ~
        for (int k = 1; k < Z0; ++k)
        for (int i = X0; i < X1; ++i)
        for (int j = Y0; j < Y1; ++j)
        {
            int l = L - k;
            
            pmlZ0By[i][j][l].x =               (               pmlZ0By[i][j][l].x - 0.5 * (E[i][j][k].z - E[i][j-1][k].z));
            pmlZ0Bz[i][j][l].x = sigmaL0B[l] * (sigmaL0B_[l] * pmlZ0Bz[i][j][l].x + 0.5 * (E[i][j][k].y - E[i][j][k-1].y));

            pmlZ0Bz[i][j][l].y = sigmaL0B[l] * (sigmaL0B_[l] * pmlZ0Bz[i][j][l].y - 0.5 * (E[i][j][k].x - E[i][j][k-1].x));
            pmlZ0Bx[i][j][l].y =               (               pmlZ0Bx[i][j][l].y + 0.5 * (E[i][j][k].z - E[i-1][j][k].z));
        
            pmlZ0Bx[i][j][l].z =               (               pmlZ0Bx[i][j][l].z - 0.5 * (E[i][j][k].y - E[i-1][j][k].y));
            pmlZ0By[i][j][l].z =               (               pmlZ0By[i][j][l].z + 0.5 * (E[i][j][k].x - E[i][j-1][k].x));

            B[i][j][k].x = pmlZ0By[i][j][l].x + pmlZ0Bz[i][j][l].x;
            B[i][j][k].y = pmlZ0Bz[i][j][l].y + pmlZ0Bx[i][j][l].y;
            B[i][j][k].z = pmlZ0Bx[i][j][l].z + pmlZ0By[i][j][l].z;
        }

        // PML - Z1 [Z1 LZ)
        //             
        // X1 + x: ~ X1-1||0|1|2|3|4|
        //         ~    S||+|*|*|*|*|  (+: perfect, *: PML, S: simulate box)
        //      l: ~     ||0|1|2|3|4|

        for (int k = Z1; k < LZ-1; ++k)
        for (int i = X0; i < X1; ++i)
        for (int j = Y0; j < Y1; ++j)
        {
            int l = k - Z1;
            
            pmlZ1By[i][j][l].x =               (               pmlZ1By[i][j][l].x - 0.5 * (E[i][j][k].z - E[i][j-1][k].z));
            pmlZ1Bz[i][j][l].x = sigmaL1B[l] * (sigmaL1B_[l] * pmlZ1Bz[i][j][l].x + 0.5 * (E[i][j][k].y - E[i][j][k-1].y));

            pmlZ1Bz[i][j][l].y = sigmaL1B[l] * (sigmaL1B_[l] * pmlZ1Bz[i][j][l].y - 0.5 * (E[i][j][k].x - E[i][j][k-1].x));
            pmlZ1Bx[i][j][l].y =               (               pmlZ1Bx[i][j][l].y + 0.5 * (E[i][j][k].z - E[i-1][j][k].z));
        
            pmlZ1Bx[i][j][l].z =               (               pmlZ1Bx[i][j][l].z - 0.5 * (E[i][j][k].y - E[i-1][j][k].y));
            pmlZ1By[i][j][l].z =               (               pmlZ1By[i][j][l].z + 0.5 * (E[i][j][k].x - E[i][j-1][k].x));

            B[i][j][k].x = pmlZ1By[i][j][l].x + pmlZ1Bz[i][j][l].x;
            B[i][j][k].y = pmlZ1Bz[i][j][l].y + pmlZ1Bx[i][j][l].y;
            B[i][j][k].z = pmlZ1Bx[i][j][l].z + pmlZ1By[i][j][l].z;
        }

        // PML - R1 [0 Y0) [0 Z0)
        for (int i = X0; i < X1; ++i)
        for (int j = 1; j < Y0; ++j)
        for (int k = 1; k < Z0; ++k)
        {
            int lj = L - j;
            int lk = L - k;
            
            pmlY0By[i][lj][k].x = sigmaL0B[lj] * (sigmaL0B_[lj] * pmlY0By[i][lj][k].x - 0.5 * (E[i][j][k].z - E[i][j-1][k].z));
            pmlY0Bz[i][lj][k].x = sigmaL0B[lk] * (sigmaL0B_[lk] * pmlY0Bz[i][lj][k].x + 0.5 * (E[i][j][k].y - E[i][j][k-1].y));

            pmlY0Bz[i][lj][k].y = sigmaL0B[lk] * (sigmaL0B_[lk] * pmlY0Bz[i][lj][k].y - 0.5 * (E[i][j][k].x - E[i][j][k-1].x));
            pmlY0Bx[i][lj][k].y =                (                pmlY0Bx[i][lj][k].y + 0.5 * (E[i][j][k].z - E[i-1][j][k].z));
        
            pmlY0Bx[i][lj][k].z =                (                pmlY0Bx[i][lj][k].z - 0.5 * (E[i][j][k].y - E[i-1][j][k].y));
            pmlY0By[i][lj][k].z = sigmaL0B[lj] * (sigmaL0B_[lj] * pmlY0By[i][lj][k].z + 0.5 * (E[i][j][k].x - E[i][j-1][k].x));

            B[i][j][k].x = pmlY0By[i][lj][k].x + pmlY0Bz[i][lj][k].x;
            B[i][j][k].y = pmlY0Bz[i][lj][k].y + pmlY0Bx[i][lj][k].y;
            B[i][j][k].z = pmlY0Bx[i][lj][k].z + pmlY0By[i][lj][k].z;
        }

        // PML - R2 [0 Y0) [Z1 LZ)
        for (int i = X0; i < X1; ++i)
        for (int j = 1; j < Y0; ++j)
        for (int k = Z1; k < LZ-1; ++k)
        {
            int lj = L - j;
            int lk = k - Z1;
            
            pmlY0By[i][lj][k].x = sigmaL0B[lj] * (sigmaL0B_[lj] * pmlY0By[i][lj][k].x - 0.5 * (E[i][j][k].z - E[i][j-1][k].z));
            pmlY0Bz[i][lj][k].x = sigmaL1B[lk] * (sigmaL1B_[lk] * pmlY0Bz[i][lj][k].x + 0.5 * (E[i][j][k].y - E[i][j][k-1].y));

            pmlY0Bz[i][lj][k].y = sigmaL1B[lk] * (sigmaL1B_[lk] * pmlY0Bz[i][lj][k].y - 0.5 * (E[i][j][k].x - E[i][j][k-1].x));
            pmlY0Bx[i][lj][k].y =                (                pmlY0Bx[i][lj][k].y + 0.5 * (E[i][j][k].z - E[i-1][j][k].z));
        
            pmlY0Bx[i][lj][k].z =                (                pmlY0Bx[i][lj][k].z - 0.5 * (E[i][j][k].y - E[i-1][j][k].y));
            pmlY0By[i][lj][k].z = sigmaL0B[lj] * (sigmaL0B_[lj] * pmlY0By[i][lj][k].z + 0.5 * (E[i][j][k].x - E[i][j-1][k].x));

            B[i][j][k].x = pmlY0By[i][lj][k].x + pmlY0Bz[i][lj][k].x;
            B[i][j][k].y = pmlY0Bz[i][lj][k].y + pmlY0Bx[i][lj][k].y;
            B[i][j][k].z = pmlY0Bx[i][lj][k].z + pmlY0By[i][lj][k].z;
        }

        // PML - R3 [Y1 LY) [0 Z0)
        for (int i = X0; i < X1; ++i)
        for (int j = Y1; j < LY-1; ++j)
        for (int k = 1; k < Z0; ++k)
        {
            int lj = j - Y1;
            int lk = L - k;

            pmlY1By[i][lj][k].x = sigmaL1B[lj] * (sigmaL1B_[lj] * pmlY1By[i][lj][k].x - 0.5 * (E[i][j][k].z - E[i][j-1][k].z));
            pmlY1Bz[i][lj][k].x = sigmaL0B[lk] * (sigmaL0B_[lk] * pmlY1Bz[i][lj][k].x + 0.5 * (E[i][j][k].y - E[i][j][k-1].y));

            pmlY1Bz[i][lj][k].y = sigmaL0B[lk] * (sigmaL0B_[lk] * pmlY1Bz[i][lj][k].y - 0.5 * (E[i][j][k].x - E[i][j][k-1].x));
            pmlY1Bx[i][lj][k].y =                (                pmlY1Bx[i][lj][k].y + 0.5 * (E[i][j][k].z - E[i-1][j][k].z));
        
            pmlY1Bx[i][lj][k].z =                (                pmlY1Bx[i][lj][k].z - 0.5 * (E[i][j][k].y - E[i-1][j][k].y));
            pmlY1By[i][lj][k].z = sigmaL1B[lj] * (sigmaL1B_[lj] * pmlY1By[i][lj][k].z + 0.5 * (E[i][j][k].x - E[i][j-1][k].x));

            B[i][j][k].x = pmlY1By[i][lj][k].x + pmlY1Bz[i][lj][k].x;
            B[i][j][k].y = pmlY1Bz[i][lj][k].y + pmlY1Bx[i][lj][k].y;
            B[i][j][k].z = pmlY1Bx[i][lj][k].z + pmlY1By[i][lj][k].z;
        }

        // PML - R4 [Y1 LY) [0 Z0)
        for (int i = X0; i < X1; ++i)
        for (int j = Y1; j < LY-1; ++j)
        for (int k = Z1; k < LZ-1; ++k)
        {
            int lj = j - Y1;
            int lk = k - Z1;

            pmlY1By[i][lj][k].x = sigmaL1B[lj] * (sigmaL1B_[lj] * pmlY1By[i][lj][k].x - 0.5 * (E[i][j][k].z - E[i][j-1][k].z));
            pmlY1Bz[i][lj][k].x = sigmaL1B[lk] * (sigmaL1B_[lk] * pmlY1Bz[i][lj][k].x + 0.5 * (E[i][j][k].y - E[i][j][k-1].y));

            pmlY1Bz[i][lj][k].y = sigmaL1B[lk] * (sigmaL1B_[lk] * pmlY1Bz[i][lj][k].y - 0.5 * (E[i][j][k].x - E[i][j][k-1].x));
            pmlY1Bx[i][lj][k].y =                (                pmlY1Bx[i][lj][k].y + 0.5 * (E[i][j][k].z - E[i-1][j][k].z));
        
            pmlY1Bx[i][lj][k].z =                (                pmlY1Bx[i][lj][k].z - 0.5 * (E[i][j][k].y - E[i-1][j][k].y));
            pmlY1By[i][lj][k].z = sigmaL1B[lj] * (sigmaL1B_[lj] * pmlY1By[i][lj][k].z + 0.5 * (E[i][j][k].x - E[i][j-1][k].x));

            B[i][j][k].x = pmlY1By[i][lj][k].x + pmlY1Bz[i][lj][k].x;
            B[i][j][k].y = pmlY1Bz[i][lj][k].y + pmlY1Bx[i][lj][k].y;
            B[i][j][k].z = pmlY1Bx[i][lj][k].z + pmlY1By[i][lj][k].z;
        }

        // X
        {
            MPICopyField(B, X1+0, X0+0, true);
            MPICopyField(B, X1+1, X0+1, true);
            MPICopyField(B, X1+2, X0+2, true);
            MPICopyField(B, X0-1, X1-1, false);
            MPICopyField(B, X0-2, X1-2, false);
        }
    }

    void BoundaryE()
    {
        // PML - Y0 [0 Y0)
        for (int i = X0; i < X1; ++i)
        for (int j =  1; j < Y0; ++j)
        for (int k = Z0; k < Z1; ++k)
        {
            int l = L - j;
            
            pmlY0Ey[i][l][k].x = sigmaL0E[l] * (sigmaL0E_[l] * pmlY0Ey[i][l][k].x + C2 * (B[i][j+1][k].z - B[i][j][k].z));
            pmlY0Ez[i][l][k].x =               (               pmlY0Ez[i][l][k].x - C2 * (B[i][j][k+1].y - B[i][j][k].y));

            pmlY0Ez[i][l][k].y =               (               pmlY0Ez[i][l][k].y + C2 * (B[i][j][k+1].x - B[i][j][k].x));
            pmlY0Ex[i][l][k].y =               (               pmlY0Ex[i][l][k].y - C2 * (B[i+1][j][k].z - B[i][j][k].z));

            pmlY0Ex[i][l][k].z =               (               pmlY0Ex[i][l][k].z + C2 * (B[i+1][j][k].y - B[i][j][k].y));
            pmlY0Ey[i][l][k].z = sigmaL0E[l] * (sigmaL0E_[l] * pmlY0Ey[i][l][k].z - C2 * (B[i][j+1][k].x - B[i][j][k].x));

            E[i][j][k].x = pmlY0Ey[i][l][k].x + pmlY0Ez[i][l][k].x;
            E[i][j][k].y = pmlY0Ez[i][l][k].y + pmlY0Ex[i][l][k].y;
            E[i][j][k].z = pmlY0Ex[i][l][k].z + pmlY0Ey[i][l][k].z;
        }

        // PML - Y1 [Y1 LY)
        for (int i = X0; i < X1; ++i)
        for (int j = Y1; j < LY-1; ++j)
        for (int k = Z0; k < Z1; ++k)
        {
            int l = j - Y1;
            
            pmlY1Ey[i][l][k].x = sigmaL1E[l] * (sigmaL1E_[l] * pmlY1Ey[i][l][k].x + C2 * (B[i][j+1][k].z - B[i][j][k].z));
            pmlY1Ez[i][l][k].x =               (               pmlY1Ez[i][l][k].x - C2 * (B[i][j][k+1].y - B[i][j][k].y));

            pmlY1Ez[i][l][k].y =               (               pmlY1Ez[i][l][k].y + C2 * (B[i][j][k+1].x - B[i][j][k].x));
            pmlY1Ex[i][l][k].y =               (               pmlY1Ex[i][l][k].y - C2 * (B[i+1][j][k].z - B[i][j][k].z));
        
            pmlY1Ex[i][l][k].z =               (               pmlY1Ex[i][l][k].z + C2 * (B[i+1][j][k].y - B[i][j][k].y));
            pmlY1Ey[i][l][k].z = sigmaL1E[l] * (sigmaL1E_[l] * pmlY1Ey[i][l][k].z - C2 * (B[i][j+1][k].x - B[i][j][k].x));

            E[i][j][k].x = pmlY1Ey[i][l][k].x + pmlY1Ez[i][l][k].x;
            E[i][j][k].y = pmlY1Ez[i][l][k].y + pmlY1Ex[i][l][k].y;
            E[i][j][k].z = pmlY1Ex[i][l][k].z + pmlY1Ey[i][l][k].z;
        }

        // PML - Z0 [0 Z0)
        for (int i = X0; i < X1; ++i)
        for (int j = Y0; j < Y1; ++j)
        for (int k = 1; k < Z0; ++k)
        {
            int l = L - k;
            
            pmlZ0Ey[i][j][l].x =               (               pmlZ0Ey[i][j][l].x + C2 * (B[i][j+1][k].z - B[i][j][k].z));
            pmlZ0Ez[i][j][l].x = sigmaL0E[l] * (sigmaL0E_[l] * pmlZ0Ez[i][j][l].x - C2 * (B[i][j][k+1].y - B[i][j][k].y));

            pmlZ0Ez[i][j][l].y = sigmaL0E[l] * (sigmaL0E_[l] * pmlZ0Ez[i][j][l].y + C2 * (B[i][j][k+1].x - B[i][j][k].x));
            pmlZ0Ex[i][j][l].y =               (               pmlZ0Ex[i][j][l].y - C2 * (B[i+1][j][k].z - B[i][j][k].z));
        
            pmlZ0Ex[i][j][l].z =               (               pmlZ0Ex[i][j][l].z + C2 * (B[i+1][j][k].y - B[i][j][k].y));
            pmlZ0Ey[i][j][l].z =               (               pmlZ0Ey[i][j][l].z - C2 * (B[i][j+1][k].x - B[i][j][k].x));

            E[i][j][k].x = pmlZ0Ey[i][j][l].x + pmlZ0Ez[i][j][l].x;
            E[i][j][k].y = pmlZ0Ez[i][j][l].y + pmlZ0Ex[i][j][l].y;
            E[i][j][k].z = pmlZ0Ex[i][j][l].z + pmlZ0Ey[i][j][l].z;
        }

        // PML - Z1 [Z1 LZ)
        for (int i = X0; i < X1; ++i)
        for (int j = Y0; j < Y1; ++j)
        for (int k = Z1; k < LZ-1; ++k)
        {
            int l = k - Z1;
            
            pmlZ1Ey[i][j][l].x =               (               pmlZ1Ey[i][j][l].x + C2 * (B[i][j+1][k].z - B[i][j][k].z));
            pmlZ1Ez[i][j][l].x = sigmaL1E[l] * (sigmaL1E_[l] * pmlZ1Ez[i][j][l].x - C2 * (B[i][j][k+1].y - B[i][j][k].y));

            pmlZ1Ez[i][j][l].y = sigmaL1E[l] * (sigmaL1E_[l] * pmlZ1Ez[i][j][l].y + C2 * (B[i][j][k+1].x - B[i][j][k].x));
            pmlZ1Ex[i][j][l].y =               (               pmlZ1Ex[i][j][l].y - C2 * (B[i+1][j][k].z - B[i][j][k].z));
        
            pmlZ1Ex[i][j][l].z =               (               pmlZ1Ex[i][j][l].z + C2 * (B[i+1][j][k].y - B[i][j][k].y));
            pmlZ1Ey[i][j][l].z =               (               pmlZ1Ey[i][j][l].z - C2 * (B[i][j+1][k].x - B[i][j][k].x));

            E[i][j][k].x = pmlZ1Ey[i][j][l].x + pmlZ1Ez[i][j][l].x;
            E[i][j][k].y = pmlZ1Ez[i][j][l].y + pmlZ1Ex[i][j][l].y;
            E[i][j][k].z = pmlZ1Ex[i][j][l].z + pmlZ1Ey[i][j][l].z;
        }

        // PML - R1 [0 Y0) [0 Z0)
        for (int i = X0; i < X1; ++i)
        for (int j = 1; j < Y0; ++j)
        for (int k = 1; k < Z0; ++k)
        {
            int lj = L - j;
            int lk = L - k;
            
            pmlY0Ey[i][lj][k].x = sigmaL0E[lj] * (sigmaL0E_[lj] * pmlY0Ey[i][lj][k].x + C2 * (B[i][j+1][k].z - B[i][j][k].z));
            pmlY0Ez[i][lj][k].x = sigmaL0E[lk] * (sigmaL0E_[lk] * pmlY0Ez[i][lj][k].x - C2 * (B[i][j][k+1].y - B[i][j][k].y));

            pmlY0Ez[i][lj][k].y = sigmaL0E[lk] * (sigmaL0E_[lk] * pmlY0Ez[i][lj][k].y + C2 * (B[i][j][k+1].x - B[i][j][k].x));
            pmlY0Ex[i][lj][k].y =                (                pmlY0Ex[i][lj][k].y - C2 * (B[i+1][j][k].z - B[i][j][k].z));

            pmlY0Ex[i][lj][k].z =                (                pmlY0Ex[i][lj][k].z + C2 * (B[i+1][j][k].y - B[i][j][k].y));
            pmlY0Ey[i][lj][k].z = sigmaL0E[lj] * (sigmaL0E_[lj] * pmlY0Ey[i][lj][k].z - C2 * (B[i][j+1][k].x - B[i][j][k].x));

            E[i][j][k].x = pmlY0Ey[i][lj][k].x + pmlY0Ez[i][lj][k].x;
            E[i][j][k].y = pmlY0Ez[i][lj][k].y + pmlY0Ex[i][lj][k].y;
            E[i][j][k].z = pmlY0Ex[i][lj][k].z + pmlY0Ey[i][lj][k].z;
        }

        // PML - R2 [0 Y0) [Z1 LZ)
        for (int i = X0; i < X1; ++i)
        for (int j = 1; j < Y0; ++j)
        for (int k = Z1; k < LZ-1; ++k)
        {
            int lj = L - j;
            int lk = k - Z1;
            
            pmlY0Ey[i][lj][k].x = sigmaL0E[lj] * (sigmaL0E_[lj] * pmlY0Ey[i][lj][k].x + C2 * (B[i][j+1][k].z - B[i][j][k].z));
            pmlY0Ez[i][lj][k].x = sigmaL1E[lk] * (sigmaL1E_[lk] * pmlY0Ez[i][lj][k].x - C2 * (B[i][j][k+1].y - B[i][j][k].y));

            pmlY0Ez[i][lj][k].y = sigmaL1E[lk] * (sigmaL1E_[lk] * pmlY0Ez[i][lj][k].y + C2 * (B[i][j][k+1].x - B[i][j][k].x));
            pmlY0Ex[i][lj][k].y =                (                pmlY0Ex[i][lj][k].y - C2 * (B[i+1][j][k].z - B[i][j][k].z));

            pmlY0Ex[i][lj][k].z =                (                pmlY0Ex[i][lj][k].z + C2 * (B[i+1][j][k].y - B[i][j][k].y));
            pmlY0Ey[i][lj][k].z = sigmaL0E[lj] * (sigmaL0E_[lj] * pmlY0Ey[i][lj][k].z - C2 * (B[i][j+1][k].x - B[i][j][k].x));

            E[i][j][k].x = pmlY0Ey[i][lj][k].x + pmlY0Ez[i][lj][k].x;
            E[i][j][k].y = pmlY0Ez[i][lj][k].y + pmlY0Ex[i][lj][k].y;
            E[i][j][k].z = pmlY0Ex[i][lj][k].z + pmlY0Ey[i][lj][k].z;
        }

        // PML - R3 [Y1 LY) [0 Z0)
        for (int i = X0; i < X1; ++i)
        for (int j = Y1; j < LY-1; ++j)
        for (int k = 1; k < Z0; ++k)
        {
            int lj = j - Y1;
            int lk = L - k;

            pmlY1Ey[i][lj][k].x = sigmaL1E[lj] * (sigmaL1E_[lj] * pmlY1Ey[i][lj][k].x + C2 * (B[i][j+1][k].z - B[i][j][k].z));
            pmlY1Ez[i][lj][k].x = sigmaL0E[lk] * (sigmaL0E_[lk] * pmlY1Ez[i][lj][k].x - C2 * (B[i][j][k+1].y - B[i][j][k].y));

            pmlY1Ez[i][lj][k].y = sigmaL0E[lk] * (sigmaL0E_[lk] * pmlY1Ez[i][lj][k].y + C2 * (B[i][j][k+1].x - B[i][j][k].x));
            pmlY1Ex[i][lj][k].y =                (                pmlY1Ex[i][lj][k].y - C2 * (B[i+1][j][k].z - B[i][j][k].z));

            pmlY1Ex[i][lj][k].z =                (                pmlY1Ex[i][lj][k].z + C2 * (B[i+1][j][k].y - B[i][j][k].y));
            pmlY1Ey[i][lj][k].z = sigmaL1E[lj] * (sigmaL1E_[lj] * pmlY1Ey[i][lj][k].z - C2 * (B[i][j+1][k].x - B[i][j][k].x));

            E[i][j][k].x = pmlY1Ey[i][lj][k].x + pmlY1Ez[i][lj][k].x;
            E[i][j][k].y = pmlY1Ez[i][lj][k].y + pmlY1Ex[i][lj][k].y;
            E[i][j][k].z = pmlY1Ex[i][lj][k].z + pmlY1Ey[i][lj][k].z;
        }

        // PML - R4 [Y1 LY) [0 Z0)
        for (int i = X0; i < X1; ++i)
        for (int j = Y1; j < LY-1; ++j)
        for (int k = Z1; k < LZ-1; ++k)
        {
            int lj = j - Y1;
            int lk = k - Z1;

            pmlY1Ey[i][lj][k].x = sigmaL1E[lj] * (sigmaL1E_[lj] * pmlY1Ey[i][lj][k].x + C2 * (B[i][j+1][k].z - B[i][j][k].z));
            pmlY1Ez[i][lj][k].x = sigmaL1E[lk] * (sigmaL1E_[lk] * pmlY1Ez[i][lj][k].x - C2 * (B[i][j][k+1].y - B[i][j][k].y));

            pmlY1Ez[i][lj][k].y = sigmaL1E[lk] * (sigmaL1E_[lk] * pmlY1Ez[i][lj][k].y + C2 * (B[i][j][k+1].x - B[i][j][k].x));
            pmlY1Ex[i][lj][k].y =                (                pmlY1Ex[i][lj][k].y - C2 * (B[i+1][j][k].z - B[i][j][k].z));

            pmlY1Ex[i][lj][k].z =                (                pmlY1Ex[i][lj][k].z + C2 * (B[i+1][j][k].y - B[i][j][k].y));
            pmlY1Ey[i][lj][k].z = sigmaL1E[lj] * (sigmaL1E_[lj] * pmlY1Ey[i][lj][k].z - C2 * (B[i][j+1][k].x - B[i][j][k].x));

            E[i][j][k].x = pmlY1Ey[i][lj][k].x + pmlY1Ez[i][lj][k].x;
            E[i][j][k].y = pmlY1Ez[i][lj][k].y + pmlY1Ex[i][lj][k].y;
            E[i][j][k].z = pmlY1Ex[i][lj][k].z + pmlY1Ey[i][lj][k].z;
        }

        // X
        {
            MPICopyField(E, X1+0, X0+0, true);
            MPICopyField(E, X1+1, X0+1, true);
            MPICopyField(E, X1+2, X0+2, true);
            MPICopyField(E, X0-1, X1-1, false);
            MPICopyField(E, X0-2, X1-2, false);
        }
    }

  private:    
    void MPICopyField(Vector***& v, const int dstX, const int srcX, const bool reverse = false)
    {
        int destRank, srcRank;

        if (reverse != true)
            MPI_Cart_shift(comm, 0,  1, &srcRank, &destRank);
        else
            MPI_Cart_shift(comm, 0, -1, &srcRank, &destRank);

        MPI_Status status;
        MPI_Sendrecv(&v[srcX][0][0], LY * LZ * 3, MPI_DOUBLE, destRank, 101,
                     &v[dstX][0][0], LY * LZ * 3, MPI_DOUBLE,  srcRank, 101,
                     comm, &status);
    }
};

template<> void FieldFrame<Boundary::GradZero>::BoundaryB()
{
    // Y
    for (int i = 0; i < LX; ++i)
    for (int k = 0; k < LZ; ++k)
    {
        B[i][Y0+0][k].x = B[i][Y0+1][k].x;
        B[i][Y0+0][k].y = B[i][Y0+1][k].y + 0.5 * (B[i][Y0+1][k].y - B[i][Y0+2][k].y);
        B[i][Y0+0][k].z = B[i][Y0+1][k].z;

        B[i][Y0-1][k].x = B[i][Y0+0][k].x;
        B[i][Y0-1][k].y = B[i][Y0+0][k].y;
        B[i][Y0-1][k].z = B[i][Y0+0][k].z;

        
        B[i][Y1+0][k].x = B[i][Y1-1][k].x;
        B[i][Y1-1][k].y = B[i][Y1-2][k].y + 0.5 * (B[i][Y1-2][k].y - B[i][Y1-3][k].y);
        B[i][Y1+0][k].z = B[i][Y1-1][k].z;

        B[i][Y1+1][k].x = B[i][Y1+0][k].x;
        B[i][Y1+0][k].y = B[i][Y1-1][k].y;
        B[i][Y1+1][k].z = B[i][Y1+0][k].z;
    }
    // Z
    for (int i = 0; i < LX; ++i)
    for (int j = 0; j < LY; ++j)
    {
        B[i][j][Z0+0].x = B[i][j][Z0+1].x;
        B[i][j][Z0+0].y = B[i][j][Z0+1].y;
        B[i][j][Z0+0].z = B[i][j][Z0+1].z + 0.5 * (B[i][j][Z0+1].z - B[i][j][Z0+2].z);

        B[i][j][Z0-1].x = B[i][j][Z0+0].x;
        B[i][j][Z0-1].y = B[i][j][Z0+0].y;
        B[i][j][Z0-1].z = B[i][j][Z0+0].z;

        
        B[i][j][Z1+0].x = B[i][j][Z1-1].x;
        B[i][j][Z1+0].y = B[i][j][Z1-1].y;
        B[i][j][Z1-1].z = B[i][j][Z1-2].z + 0.5 * (B[i][j][Z1-2].z - B[i][j][Z1-3].z);

        B[i][j][Z1+1].x = B[i][j][Z1+0].x;
        B[i][j][Z1+1].y = B[i][j][Z1+0].y;
        B[i][j][Z1+0].z = B[i][j][Z1-1].z;
    }
    // X
    {
        MPICopyField(B, X1+0, X0+0, true);
        MPICopyField(B, X1+1, X0+1, true);
        MPICopyField(B, X1+2, X0+2, true);
        MPICopyField(B, X0-1, X1-1, false);
        MPICopyField(B, X0-2, X1-2, false);
    }
}

template<> void FieldFrame<Boundary::GradZero>::BoundaryE()
{
    // Y
    for (int i = 0; i < LX; ++i)
    for (int k = 0; k < LZ; ++k)
    {
        E[i][Y0+0][k].x = E[i][Y0+1][k].x + 0.5 * (E[i][Y0+1][k].x - E[i][Y0+2][k].x);
        E[i][Y0+0][k].y = E[i][Y0+1][k].y;
        E[i][Y0+0][k].z = E[i][Y0+1][k].z + 0.5 * (E[i][Y0+1][k].z - E[i][Y0+2][k].z);

        E[i][Y0-1][k].x = E[i][Y0+0][k].x;
        E[i][Y0-1][k].y = E[i][Y0+0][k].y;
        E[i][Y0-1][k].z = E[i][Y0+0][k].z;
        

        E[i][Y1-1][k].x = E[i][Y1-2][k].x + 0.5 * (E[i][Y1-2][k].x - E[i][Y1-3][k].x);
        E[i][Y1+0][k].y = E[i][Y1-1][k].y;
        E[i][Y1-1][k].z = E[i][Y1-2][k].z + 0.5 * (E[i][Y1-2][k].z - E[i][Y1-3][k].z);

        E[i][Y1+0][k].x = E[i][Y1-1][k].x;
        E[i][Y1+1][k].y = E[i][Y1+0][k].y;
        E[i][Y1+0][k].z = E[i][Y1-1][k].z;
    }
    // Z
    for (int i = 0; i < LX; ++i)
    for (int j = 0; j < LY; ++j)
    {
        E[i][j][Z0+0].x = E[i][j][Z0+1].x + 0.5 * (E[i][j][Z0+1].x - E[i][j][Z0+2].x);
        E[i][j][Z0+0].y = E[i][j][Z0+1].y + 0.5 * (E[i][j][Z0+1].y - E[i][j][Z0+2].y);
        E[i][j][Z0+0].z = E[i][j][Z0+1].z;

        E[i][j][Z0-1].x = E[i][j][Z0+0].x;
        E[i][j][Z0-1].y = E[i][j][Z0+0].y;
        E[i][j][Z0-1].z = E[i][j][Z0+0].z;
        

        E[i][j][Z1-1].x = E[i][j][Z1-2].x + 0.5 * (E[i][j][Z1-2].x - E[i][j][Z1-3].x);
        E[i][j][Z1-1].y = E[i][j][Z1-2].y + 0.5 * (E[i][j][Z1-2].y - E[i][j][Z1-3].y);
        E[i][j][Z1+0].z = E[i][j][Z1-1].z;

        E[i][j][Z1+0].x = E[i][j][Z1-1].x;
        E[i][j][Z1+0].y = E[i][j][Z1-1].y;
        E[i][j][Z1+1].z = E[i][j][Z1+0].z;
    }
    // X
    {
        MPICopyField(E, X1+0, X0+0, true);
        MPICopyField(E, X1+1, X0+1, true);
        MPICopyField(E, X1+2, X0+2, true);
        MPICopyField(E, X0-1, X1-1, false);
        MPICopyField(E, X0-2, X1-2, false);
    }
}

template<>
class FieldFrame<Boundary::Absorbing> : public FieldBase
{
  public:
    double* sigmaL0E;
    double* sigmaL1E;
    double* sigmaL0E_;
    double* sigmaL1E_;
    
    Vector*** pmlZ0Ex;
    Vector*** pmlZ0Ey;
    Vector*** pmlZ0Ez;
    Vector*** pmlZ1Ex;
    Vector*** pmlZ1Ey;
    Vector*** pmlZ1Ez;

    double* sigmaL0B;
    double* sigmaL1B;
    double* sigmaL0B_;
    double* sigmaL1B_;

    Vector*** pmlZ0Bx;
    Vector*** pmlZ0By;
    Vector*** pmlZ0Bz;
    Vector*** pmlZ1Bx;
    Vector*** pmlZ1By;
    Vector*** pmlZ1Bz;

  public:
    FieldFrame()
    {
        try
        {            
            sigmaL0E  = new double [L];
            sigmaL0B  = new double [L];
            sigmaL0E_ = new double [L];
            sigmaL0B_ = new double [L];

            sigmaL1E  = new double [L];
            sigmaL1B  = new double [L];
            sigmaL1E_ = new double [L];
            sigmaL1B_ = new double [L];

            pmlZ0Ex = CreateArray(LX, LY, L);
            pmlZ0Ey = CreateArray(LX, LY, L);
            pmlZ0Ez = CreateArray(LX, LY, L);
            pmlZ1Ex = CreateArray(LX, LY, L);
            pmlZ1Ey = CreateArray(LX, LY, L);
            pmlZ1Ez = CreateArray(LX, LY, L);
            
            pmlZ0Bx = CreateArray(LX, LY, L);
            pmlZ0By = CreateArray(LX, LY, L);
            pmlZ0Bz = CreateArray(LX, LY, L);
            pmlZ1Bx = CreateArray(LX, LY, L);
            pmlZ1By = CreateArray(LX, LY, L);
            pmlZ1Bz = CreateArray(LX, LY, L);
        }
        catch(std::bad_alloc)
        {
            fprintf(stderr, "Errer: bad_alloc : %s\n", __FILE__);            

            delete[] sigmaL0E;
            delete[] sigmaL0B;
            delete[] sigmaL0E_;
            delete[] sigmaL0B_;

            delete[] sigmaL1E;
            delete[] sigmaL1B;
            delete[] sigmaL1E_;
            delete[] sigmaL1B_;

            DeleteArray(pmlZ0Ex);
            DeleteArray(pmlZ0Ey);
            DeleteArray(pmlZ0Ez);
            DeleteArray(pmlZ1Ex);
            DeleteArray(pmlZ1Ey);
            DeleteArray(pmlZ1Ez);

            DeleteArray(pmlZ0Bx);
            DeleteArray(pmlZ0By);
            DeleteArray(pmlZ0Bz);
            DeleteArray(pmlZ1Bx);
            DeleteArray(pmlZ1By);
            DeleteArray(pmlZ1Bz);

            abort();
        }

        for (int i = 0; i < LX; ++i)
	    for (int j = 0; j < LY; ++j)
	    for (int k = 0; k < L ; ++k)
        {
            pmlZ0Ex[i][j][k].Zero();
            pmlZ0Ey[i][j][k].Zero();
            pmlZ0Ez[i][j][k].Zero();
            pmlZ1Ex[i][j][k].Zero();
            pmlZ1Ey[i][j][k].Zero();
            pmlZ1Ez[i][j][k].Zero();

            pmlZ0Bx[i][j][k].Zero();
            pmlZ0By[i][j][k].Zero();
            pmlZ0Bz[i][j][k].Zero();
            pmlZ1Bx[i][j][k].Zero();
            pmlZ1By[i][j][k].Zero();
            pmlZ1Bz[i][j][k].Zero();
        }

        for (int n = 0; n < L; ++n)
        {
            double sigL0E = SIGMA_MAX * std::pow((n + 0.5) / L, M);
            double sigL0B = SIGMA_MAX * std::pow((n + 1.0) / L, M);
            double sigL1E = SIGMA_MAX * std::pow((n + 0.5) / L, M);
            double sigL1B = SIGMA_MAX * std::pow((n + 0.0) / L, M);

            sigmaL0E[n]  = 1.0 / (1.0 + sigL0E * 0.5);
            sigmaL0E_[n] = (1.0 - sigL0E * 0.5);
            sigmaL1E[n]  = 1.0 / (1.0 + sigL1E * 0.5);
            sigmaL1E_[n] = (1.0 - sigL1E * 0.5);

            sigmaL0B[n]  = 1.0 / (1.0 + sigL0B * 0.25);
            sigmaL0B_[n] = (1.0 - sigL0B * 0.25);
            sigmaL1B[n]  = 1.0 / (1.0 + sigL1B * 0.25);
            sigmaL1B_[n] = (1.0 - sigL1B * 0.25);
        }
    }
    
    ~FieldFrame()
    {
        delete[] sigmaL0E;
        delete[] sigmaL0B;
        delete[] sigmaL0E_;
        delete[] sigmaL0B_;

        delete[] sigmaL1E;
        delete[] sigmaL1B;
        delete[] sigmaL1E_;
        delete[] sigmaL1B_;

        DeleteArray(pmlZ0Ex);
        DeleteArray(pmlZ0Ey);
        DeleteArray(pmlZ0Ez);
        DeleteArray(pmlZ1Ex);
        DeleteArray(pmlZ1Ey);
        DeleteArray(pmlZ1Ez);

        DeleteArray(pmlZ0Bx);
        DeleteArray(pmlZ0By);
        DeleteArray(pmlZ0Bz);
        DeleteArray(pmlZ1Bx);
        DeleteArray(pmlZ1By);
        DeleteArray(pmlZ1Bz);
    }

    void BoundaryB()
    {
        // PML - Z0 [1 Z0) ... 
        for (int i = X0; i < X1; ++i)
        for (int j = 1; j < LY-1; ++j)
        for (int k = 1; k < Z0; ++k)
        {
            int l = L - k;
            
            pmlZ0By[i][j][l].x =               (               pmlZ0By[i][j][l].x - 0.5 * (E[i][j][k].z - E[i][j-1][k].z));
            pmlZ0Bz[i][j][l].x = sigmaL0B[l] * (sigmaL0B_[l] * pmlZ0Bz[i][j][l].x + 0.5 * (E[i][j][k].y - E[i][j][k-1].y));

            pmlZ0Bz[i][j][l].y = sigmaL0B[l] * (sigmaL0B_[l] * pmlZ0Bz[i][j][l].y - 0.5 * (E[i][j][k].x - E[i][j][k-1].x));
            pmlZ0Bx[i][j][l].y =               (               pmlZ0Bx[i][j][l].y + 0.5 * (E[i][j][k].z - E[i-1][j][k].z));
        
            pmlZ0Bx[i][j][l].z =               (               pmlZ0Bx[i][j][l].z - 0.5 * (E[i][j][k].y - E[i-1][j][k].y));
            pmlZ0By[i][j][l].z =               (               pmlZ0By[i][j][l].z + 0.5 * (E[i][j][k].x - E[i][j-1][k].x));

            B[i][j][k].x = pmlZ0By[i][j][l].x + pmlZ0Bz[i][j][l].x;
            B[i][j][k].y = pmlZ0Bz[i][j][l].y + pmlZ0Bx[i][j][l].y;
            B[i][j][k].z = pmlZ0Bx[i][j][l].z + pmlZ0By[i][j][l].z;
        }

        // PML - Z1 [Z1 LZ)
        for (int i = X0; i < X1; ++i)
        for (int j = 1; j < LY-1; ++j)
        for (int k = Z1; k < LZ-1; ++k)
        {
            int l = k - Z1;
            
            pmlZ1By[i][j][l].x =               (               pmlZ1By[i][j][l].x - 0.5 * (E[i][j][k].z - E[i][j-1][k].z));
            pmlZ1Bz[i][j][l].x = sigmaL1B[l] * (sigmaL1B_[l] * pmlZ1Bz[i][j][l].x + 0.5 * (E[i][j][k].y - E[i][j][k-1].y));

            pmlZ1Bz[i][j][l].y = sigmaL1B[l] * (sigmaL1B_[l] * pmlZ1Bz[i][j][l].y - 0.5 * (E[i][j][k].x - E[i][j][k-1].x));
            pmlZ1Bx[i][j][l].y =               (               pmlZ1Bx[i][j][l].y + 0.5 * (E[i][j][k].z - E[i-1][j][k].z));
        
            pmlZ1Bx[i][j][l].z =               (               pmlZ1Bx[i][j][l].z - 0.5 * (E[i][j][k].y - E[i-1][j][k].y));
            pmlZ1By[i][j][l].z =               (               pmlZ1By[i][j][l].z + 0.5 * (E[i][j][k].x - E[i][j-1][k].x));

            B[i][j][k].x = pmlZ1By[i][j][l].x + pmlZ1Bz[i][j][l].x;
            B[i][j][k].y = pmlZ1Bz[i][j][l].y + pmlZ1Bx[i][j][l].y;
            B[i][j][k].z = pmlZ1Bx[i][j][l].z + pmlZ1By[i][j][l].z;
        }
        // X
        {
            MPICopyField(B, X1+0, X0+0, true);
            MPICopyField(B, X1+1, X0+1, true);
            MPICopyField(B, X1+2, X0+2, true);
            MPICopyField(B, X0-1, X1-1, false);
            MPICopyField(B, X0-2, X1-2, false);
        }
        // Y
        for (int i = 0; i < LX; ++i)
        for (int k = 0; k < LZ; ++k)
        {
            B[i][Y1+0][k] = B[i][Y0+0][k];
            B[i][Y1+1][k] = B[i][Y0+1][k];
            B[i][Y1+2][k] = B[i][Y0+2][k];
            B[i][Y0-1][k] = B[i][Y1-1][k];
            B[i][Y0-2][k] = B[i][Y1-2][k];
        }
    }

    void BoundaryE()
    {
        // PML - Z0 [0 Z0)
        for (int i = X0; i < X1; ++i)
        for (int j = 1; j < LY-1; ++j)
        for (int k = 1; k < Z0; ++k)
        {
            int l = L - k;
            
            pmlZ0Ey[i][j][l].x =               (               pmlZ0Ey[i][j][l].x + C2 * (B[i][j+1][k].z - B[i][j][k].z));
            pmlZ0Ez[i][j][l].x = sigmaL0E[l] * (sigmaL0E_[l] * pmlZ0Ez[i][j][l].x - C2 * (B[i][j][k+1].y - B[i][j][k].y));

            pmlZ0Ez[i][j][l].y = sigmaL0E[l] * (sigmaL0E_[l] * pmlZ0Ez[i][j][l].y + C2 * (B[i][j][k+1].x - B[i][j][k].x));
            pmlZ0Ex[i][j][l].y =               (               pmlZ0Ex[i][j][l].y - C2 * (B[i+1][j][k].z - B[i][j][k].z));
        
            pmlZ0Ex[i][j][l].z =               (               pmlZ0Ex[i][j][l].z + C2 * (B[i+1][j][k].y - B[i][j][k].y));
            pmlZ0Ey[i][j][l].z =               (               pmlZ0Ey[i][j][l].z - C2 * (B[i][j+1][k].x - B[i][j][k].x));

            E[i][j][k].x = pmlZ0Ey[i][j][l].x + pmlZ0Ez[i][j][l].x;
            E[i][j][k].y = pmlZ0Ez[i][j][l].y + pmlZ0Ex[i][j][l].y;
            E[i][j][k].z = pmlZ0Ex[i][j][l].z + pmlZ0Ey[i][j][l].z;
        }

        // PML - Z1 [Z1 LZ)
        for (int i = X0; i < X1; ++i)
        for (int j = 1; j < LY-1; ++j)
        for (int k = Z1; k < LZ-1; ++k)
        {
            int l = k - Z1;
            
            pmlZ1Ey[i][j][l].x =               (               pmlZ1Ey[i][j][l].x + C2 * (B[i][j+1][k].z - B[i][j][k].z));
            pmlZ1Ez[i][j][l].x = sigmaL1E[l] * (sigmaL1E_[l] * pmlZ1Ez[i][j][l].x - C2 * (B[i][j][k+1].y - B[i][j][k].y));

            pmlZ1Ez[i][j][l].y = sigmaL1E[l] * (sigmaL1E_[l] * pmlZ1Ez[i][j][l].y + C2 * (B[i][j][k+1].x - B[i][j][k].x));
            pmlZ1Ex[i][j][l].y =               (               pmlZ1Ex[i][j][l].y - C2 * (B[i+1][j][k].z - B[i][j][k].z));
        
            pmlZ1Ex[i][j][l].z =               (               pmlZ1Ex[i][j][l].z + C2 * (B[i+1][j][k].y - B[i][j][k].y));
            pmlZ1Ey[i][j][l].z =               (               pmlZ1Ey[i][j][l].z - C2 * (B[i][j+1][k].x - B[i][j][k].x));

            E[i][j][k].x = pmlZ1Ey[i][j][l].x + pmlZ1Ez[i][j][l].x;
            E[i][j][k].y = pmlZ1Ez[i][j][l].y + pmlZ1Ex[i][j][l].y;
            E[i][j][k].z = pmlZ1Ex[i][j][l].z + pmlZ1Ey[i][j][l].z;
        }

        // X
        {
            MPICopyField(E, X1+0, X0+0, true);
            MPICopyField(E, X1+1, X0+1, true);
            MPICopyField(E, X1+2, X0+2, true);
            MPICopyField(E, X0-1, X1-1, false);
            MPICopyField(E, X0-2, X1-2, false);
        }
        // Y
        for (int i = 0; i < LX; ++i)
        for (int k = 0; k < LZ; ++k)
        {
            E[i][Y1+0][k] = E[i][Y0+0][k];
            E[i][Y1+1][k] = E[i][Y0+1][k];
            E[i][Y1+2][k] = E[i][Y0+2][k];
            E[i][Y0-1][k] = E[i][Y1-1][k];
            E[i][Y0-2][k] = E[i][Y1-2][k];
        }
    }

  private:    
    void MPICopyField(Vector***& v, const int dstX, const int srcX, const bool reverse = false)
    {
        int destRank, srcRank;

        if (reverse != true)
            MPI_Cart_shift(comm, 0,  1, &srcRank, &destRank);
        else
            MPI_Cart_shift(comm, 0, -1, &srcRank, &destRank);

        MPI_Status status;
        MPI_Sendrecv(&v[srcX][0][0], LY * LZ * 3, MPI_DOUBLE, destRank, 101,
                     &v[dstX][0][0], LY * LZ * 3, MPI_DOUBLE,  srcRank, 101,
                     comm, &status);
    }
};

#endif

