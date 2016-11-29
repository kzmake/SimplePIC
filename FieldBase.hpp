#ifndef SIMPLE_PIC_FIELD_BASE_HPP
#define SIMPLE_PIC_FIELD_BASE_HPP

#include <cstdlib>

#include "Param.hpp"

class FieldBase
{
  public:
    Vector*** E;
    Vector*** B;

  protected:
    FieldBase()
    {
        try
        {
            E = CreateArray(LX, LY, LZ);
            B = CreateArray(LX, LY, LZ);
        }
        catch(std::bad_alloc)
        {
            fprintf(stderr, "Errer: bad_alloc : %s\n", __FILE__);
            DeleteArray(E);
            DeleteArray(B);

            abort();
        }

        for (int i = 0; i < LX; ++i)
	    for (int j = 0; j < LY; ++j)
	    for (int k = 0; k < LZ; ++k)
        {
            E[i][j][k].Zero();
            B[i][j][k].Zero();
        }
    }

  public:
    ~FieldBase()
    {
        DeleteArray(E);
        DeleteArray(B);
    }

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

  protected:
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
};

#endif
