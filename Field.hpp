#ifndef SIMPLE_PIC_FIELD_HPP
#define SIMPLE_PIC_FIELD_HPP

#include <cstdlib>

#include "Param.hpp"
#include "Vector.hpp"

class Field
{
  public:
    Vector*** E;
    Vector*** B;
#ifdef MPI_PIC
    double* mpiBuf;
#endif
  public:
    Field()
    {
        try
        {
            E = CreateArray(LX, LY, LZ);
            B = CreateArray(LX, LY, LZ);
#ifdef MPI_PIC        
            mpiBuf = new double [3 * LY * LZ];
#endif
        }
        catch(std::bad_alloc)
        {
            fprintf(stderr, "エラー: bad_alloc 動的メモリ確保エラー: %s\n", __FILE__);
            DeleteArray(E);
            DeleteArray(B);
#ifdef MPI_PIC
            delete[] mpiBuf;
#endif
            abort();
        }

        for (int i = 0; i < LX; ++i)
	    for (int j = 0; j < LY; ++j)
	    for (int k = 0; k < LZ; ++k)
        {
            B[i][j][k].Zero();
        }

        for (int i = 0; i < LX; ++i)
	    for (int j = 0; j < LY; ++j)
	    for (int k = 0; k < LZ; ++k)
        {
            E[i][j][k].Zero();
        }

    }
    ~Field()
    {
        DeleteArray(E);
        DeleteArray(B);
#ifdef MPI_PIC
        delete[] mpiBuf;
#endif
    }

    Vector*** CreateArray(const int x, const int y, const int z)
    {
        Vector*** v  = new Vector** [x];
        v[0]         = new Vector*  [x * y];
        v[0][0]      = new Vector   [x * y * z]();

        for (int i = 0; i < x; ++i)
        {
            v[i] = v[0] + i * y;

            for (int j = 0; j < y; ++j)
            {
                v[i][j] = v[0][0] + i * y * z + j * z;
            }
        }

        return v;
    }

    void DeleteArray(Vector***& v)
    {
        delete[] v[0][0];
        delete[] v[0];
        delete[] v;
        v = nullptr;
    }

#if 1
    void UpdateB()
    {
        for (int i = X0; i < X1; ++i)
        for (int j = Y0; j < Y1; ++j)
        for (int k = Z0; k < Z1; ++k)
        {
            B[i][j][k].x -= 0.5 * (E[i][j][k].z - E[i][j-1][k].z - E[i][j][k].y + E[i][j][k-1].y);
            B[i][j][k].y -= 0.5 * (E[i][j][k].x - E[i][j][k-1].x - E[i][j][k].z + E[i-1][j][k].z);
            B[i][j][k].z -= 0.5 * (E[i][j][k].y - E[i-1][j][k].y - E[i][j][k].x + E[i][j-1][k].x);
        }
    }

    void UpdateE()
    {
        for (int i = X0; i < X1; ++i)
        for (int j = Y0; j < Y1; ++j)
        for (int k = Z0; k < Z1; ++k)
        {
            E[i][j][k].x += C2 * (B[i][j+1][k].z - B[i][j][k].z - B[i][j][k+1].y + B[i][j][k].y);
            E[i][j][k].y += C2 * (B[i][j][k+1].x - B[i][j][k].x - B[i+1][j][k].z + B[i][j][k].z);
            E[i][j][k].z += C2 * (B[i+1][j][k].y - B[i][j][k].y - B[i][j+1][k].x + B[i][j][k].x);
        }
    }
#else
    void UpdateB()
    {
        for (int i = X0; i < X1; ++i)
        for (int j = Y0; j < Y1; ++j)
        for (int k = Z0; k < Z1; ++k)
        {
            B[i][j][k].x += 0.5 * C * (E[i][j][k].y - E[i][j][k-1].y - E[i][j][k].z + E[i][j-1][k].z);
            B[i][j][k].y += 0.5 * C * (E[i][j][k].z - E[i-1][j][k].z - E[i][j][k].x + E[i][j][k-1].x);
            B[i][j][k].z += 0.5 * C * (E[i][j][k].x - E[i][j-1][k].x - E[i][j][k].y + E[i-1][j][k].y);
        }
    }

    void UpdateE()
    {
        for (int i = X0; i < X1; ++i)
        for (int j = Y0; j < Y1; ++j)
        for (int k = Z0; k < Z1; ++k)
        {
            E[i][j][k].x -= C * (B[i][j][k+1].y - B[i][j][k].y - B[i][j+1][k].z + B[i][j][k].z);
            E[i][j][k].y -= C * (B[i+1][j][k].z - B[i][j][k].z - B[i][j][k+1].x + B[i][j][k].x);
            E[i][j][k].z -= C * (B[i][j+1][k].x - B[i][j][k].x - B[i+1][j][k].y + B[i][j][k].y);
        }
    }
#endif
    void BoundaryB()
    {
#ifdef MPI_PIC
        auto MPICopyField = [this](Vector***& v, const int dstX, const int srcX, const bool reverse = false)
        {
            for (int j = 0; j < LY; ++j)
            for (int k = 0; k < LZ; ++k)
            {
                mpiBuf[j * LZ * 3 + k * 3 + 0] = v[srcX][j][k].x;
                mpiBuf[j * LZ * 3 + k * 3 + 1] = v[srcX][j][k].y;
                mpiBuf[j * LZ * 3 + k * 3 + 2] = v[srcX][j][k].z;
            }

            int forward = (MPI::COMM_WORLD.Get_rank() + 1) % MPI::COMM_WORLD.Get_size();
            int backward = MPI::COMM_WORLD.Get_rank() - 1;
            if (backward < 0) backward = MPI::COMM_WORLD.Get_size() - 1;

            MPI::Status status;
            
            if (reverse != true)
                MPI::COMM_WORLD.Sendrecv_replace(mpiBuf, 3 * LY * LZ, MPI::DOUBLE,
                    forward, 0, backward, 0, status);
            else
                MPI::COMM_WORLD.Sendrecv_replace(mpiBuf, 3 * LY * LZ, MPI::DOUBLE,
                    backward, 0, forward, 0, status);

            for (int j = 0; j < LY; ++j)
            for (int k = 0; k < LZ; ++k)
            {
                v[dstX][j][k].x = mpiBuf[j * LZ * 3 + k * 3 + 0];
                v[dstX][j][k].y = mpiBuf[j * LZ * 3 + k * 3 + 1];
                v[dstX][j][k].z = mpiBuf[j * LZ * 3 + k * 3 + 2];
            }

        };
#endif
        // periodic - Copy
        // X
#ifdef MPI_PIC
        {
            MPICopyField(B, X1+0, X0+0, true);
            MPICopyField(B, X1+1, X0+1, true);
            MPICopyField(B, X1+2, X0+2, true);
            MPICopyField(B, X0-1, X1-1, false);
            MPICopyField(B, X0-2, X1-2, false);
        }
#else
        for (int j = 0; j < LY; ++j)
        for (int k = 0; k < LZ; ++k)
        {
            B[X1+0][j][k] = B[X0+0][j][k];
            B[X1+1][j][k] = B[X0+1][j][k];
            B[X1+2][j][k] = B[X0+2][j][k];
            B[X0-1][j][k] = B[X1-1][j][k];
            B[X0-2][j][k] = B[X1-2][j][k];
        }
#endif
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

    void BoundaryE()
    {
#ifdef MPI_PIC
        auto MPICopyField = [this](Vector***& v, const int dstX, const int srcX, const bool reverse = false)
        {
            for (int j = 0; j < LY; ++j)
            for (int k = 0; k < LZ; ++k)
            {
                mpiBuf[j * LZ * 3 + k * 3 + 0] = v[srcX][j][k].x;
                mpiBuf[j * LZ * 3 + k * 3 + 1] = v[srcX][j][k].y;
                mpiBuf[j * LZ * 3 + k * 3 + 2] = v[srcX][j][k].z;
            }

            int forward = (MPI::COMM_WORLD.Get_rank() + 1) % MPI::COMM_WORLD.Get_size();
            int backward = MPI::COMM_WORLD.Get_rank() - 1;
            if (backward < 0) backward = MPI::COMM_WORLD.Get_size() - 1;
            
            MPI::Status status; 
            
            if (reverse != true)
                MPI::COMM_WORLD.Sendrecv_replace(mpiBuf, 3 * LY * LZ, MPI::DOUBLE,
                     forward, 0, backward, 0, status);
            else
                MPI::COMM_WORLD.Sendrecv_replace(mpiBuf, 3 * LY * LZ, MPI::DOUBLE,
                    backward, 0, forward, 0, status);

            for (int j = 0; j < LY; ++j)
            for (int k = 0; k < LZ; ++k)
            {
                v[dstX][j][k].x = mpiBuf[j * LZ * 3 + k * 3 + 0];
                v[dstX][j][k].y = mpiBuf[j * LZ * 3 + k * 3 + 1];
                v[dstX][j][k].z = mpiBuf[j * LZ * 3 + k * 3 + 2];
            }

        };
#endif

        // periodic - Copy
        // X
#ifdef MPI_PIC
        {
            MPICopyField(E, X1+0, X0+0, true);
            MPICopyField(E, X1+1, X0+1, true);
            MPICopyField(E, X1+2, X0+2, true);
            MPICopyField(E, X0-1, X1-1, false);
            MPICopyField(E, X0-2, X1-2, false);
        }
#else
        for (int j = 0; j < LY; ++j)
        for (int k = 0; k < LZ; ++k)
        {
            E[X1+0][j][k] = E[X0+0][j][k];
            E[X1+1][j][k] = E[X0+1][j][k];
            E[X1+2][j][k] = E[X0+2][j][k];
            E[X0-1][j][k] = E[X1-1][j][k];
            E[X0-2][j][k] = E[X1-2][j][k];
        }
#endif
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
};

#endif
