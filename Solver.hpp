#ifndef SIMPLE_PIC_SOLVER_HPP
#define SIMPLE_PIC_SOLVER_HPP

#include <vector>
#include <cstdlib>

#include "Param.hpp"

#include "Vector.hpp"
#include "ShapeFactor.hpp"
#include "Plasma.hpp"
#include "Field.hpp"


template<Shape S>
class Solver
{
  public:
    ShapeFactor<S, 3> sf3;
    ShapeFactor<S, 5> sf5;

    Vector*** J;
    Vector*** Bc;
    Vector*** Ec;

    double* mpiBuf;

  public:
    Solver()
    {
        try
        {
            J  = CreateArray(LX, LY, LZ);
            Ec = CreateArray(LX, LY, LZ);
            Bc = CreateArray(LX, LY, LZ);

            mpiBuf = new double [3 * LY * LZ];
        }
        catch(std::bad_alloc)
        {
            fprintf(stderr, "エラー: bad_alloc 動的メモリ確保エラー: %s\n", __FILE__); 
            DeleteArray(J);
            DeleteArray(Ec);

            delete[] mpiBuf;
            
            abort();
        }

        for (int i = 0; i < LX; ++i)
	    for (int j = 0; j < LY; ++j)
	    for (int k = 0; k < LZ; ++k)
	    {
		    J[i][j][k].Zero();
        }

        //muda 
        for (int i = 0; i < LX; ++i)
	    for (int j = 0; j < LY; ++j)
	    for (int k = 0; k < LZ; ++k)
	    {
		    Bc[i][j][k].Zero();
            Ec[i][j][k].Zero();
        }
    }
    ~Solver()
    {
        DeleteArray(J);
        DeleteArray(Ec);
        DeleteArray(Bc);

        delete[] mpiBuf;
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

    void DensityDecomposition(Plasma& plasma, Field& field)
    {
        std::vector<Particle>& p = plasma.p;

        double s0[3][5];
        double ds[3][5];
        double w[5][5][5][3];

        const double w1_3 = (1.0/3.0);

        Vector J0[5][5][5];


        for (unsigned long int n = 0; n < p.size(); ++n)
        {
            Vector r0 = p[n].r;
            Vector r1 = p[n].r + p[n].v;

            sf5(s0[0], r0.x);
            sf5(s0[1], r0.y);
            sf5(s0[2], r0.z);

            sf5(ds[0], r1.x, int(r1.x) - int(r0.x));
            sf5(ds[1], r1.y, int(r1.y) - int(r0.y));
            sf5(ds[2], r1.z, int(r1.z) - int(r0.z));

            for(int i = 0; i < 5; ++i)
            {
                ds[0][i] -= s0[0][i];
                ds[1][i] -= s0[1][i];
                ds[2][i] -= s0[2][i];
            }

            for(int i = 0; i < 5; ++i)
            for(int j = 0; j < 5; ++j)
            for(int k = 0; k < 5; ++k)
            {
                w[i][j][k][0] = ds[0][i] * (s0[1][j]*s0[2][k] + 0.5*ds[1][j]*s0[2][k] + 0.5*s0[1][j]*ds[2][k] + w1_3*ds[1][j]*ds[2][k]);
                w[i][j][k][1] = ds[1][j] * (s0[0][i]*s0[2][k] + 0.5*ds[0][i]*s0[2][k] + 0.5*s0[0][i]*ds[2][k] + w1_3*ds[0][i]*ds[2][k]);
                w[i][j][k][2] = ds[2][k] * (s0[0][i]*s0[1][j] + 0.5*ds[0][i]*s0[1][j] + 0.5*s0[0][i]*ds[1][j] + w1_3*ds[0][i]*ds[1][j]);

            }
            
            for(int i = 0; i < 5; ++i)
            for(int j = 0; j < 5; ++j)
            for(int k = 0; k < 5; ++k)
            {
                J0[i][j][k].Zero();
            }

            for(int i = 0; i < 5-1; ++i)
            for(int j = 0; j < 5; ++j)
            for(int k = 0; k < 5; ++k)
            {
                J0[i+1][j][k].x = J0[i][j][k].x - plasma.q * w[i][j][k][0];
            }

            for(int i = 0; i < 5; ++i)
            for(int j = 0; j < 5-1; ++j)
            for(int k = 0; k < 5; ++k)
            {
                J0[i][j+1][k].y = J0[i][j][k].y - plasma.q * w[i][j][k][1];
            }

            for(int i = 0; i < 5; ++i)
            for(int j = 0; j < 5; ++j)
            for(int k = 0; k < 5-1; ++k)
            {
                J0[i][j][k+1].z = J0[i][j][k].z - plasma.q * w[i][j][k][2];
            }

            int ii = int(r0.x) - 2;
            int jj = int(r0.y) - 2;
            int kk = int(r0.z) - 2;

            for(int i = 0; i < 5; ++i)
            for(int j = 0; j < 5; ++j)
            for(int k = 0; k < 5; ++k)
            {
                J[ii+i][jj+j][kk+k].x += J0[i][j][k].x;
                J[ii+i][jj+j][kk+k].y += J0[i][j][k].y;
                J[ii+i][jj+j][kk+k].z += J0[i][j][k].z;
            } 
        }
    }
    
    void BoundaryJ()
    {
#ifdef MPI_PIC
        auto MPIAddField = [this](Vector***& v, const int dstX, const int srcX, const bool reverse = false)
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
                v[dstX][j][k].x += mpiBuf[j * LZ * 3 + k * 3 + 0];
                v[dstX][j][k].y += mpiBuf[j * LZ * 3 + k * 3 + 1];
                v[dstX][j][k].z += mpiBuf[j * LZ * 3 + k * 3 + 2];
            }

        };
#endif
        // periodic - Add
        // X
#ifdef MPI_PIC
        {
            MPIAddField(J, X0+0, X1+0, false);
            MPIAddField(J, X0+1, X1+1, false);
            MPIAddField(J, X1-1, X0-1, true);
            MPIAddField(J, X1-2, X0-2, true);
        }
#else
        for (int j = 0; j < LY; ++j)
        for (int k = 0; k < LZ; ++k)
        {
            J[X0+0][j][k] += J[X1+0][j][k];
            J[X0+1][j][k] += J[X1+1][j][k];
            J[X1-1][j][k] += J[X0-1][j][k];
            J[X1-2][j][k] += J[X0-2][j][k];
        }
#endif
        // Y
        for (int i = 0; i < LX; ++i)
        for (int k = 0; k < LZ; ++k)
        {
            J[i][Y0+0][k] += J[i][Y1+0][k];
            J[i][Y0+1][k] += J[i][Y1+1][k];
            J[i][Y1-1][k] += J[i][Y0-1][k];
            J[i][Y1-2][k] += J[i][Y0-2][k];
        }
        // Z
        for (int i = 0; i < LX; ++i)
        for (int j = 0; j < LY; ++j)
        {
            J[i][j][Z0+0] += J[i][j][Z1+0];
            J[i][j][Z0+1] += J[i][j][Z1+1];
            J[i][j][Z1-1] += J[i][j][Z0-1];
            J[i][j][Z1-2] += J[i][j][Z0-2];
        }
    }  

    void UpdateEbyJ(Field& field)
    {
        for (int i = X0; i < X1; ++i)
	    for (int j = Y0; j < Y1; ++j)
	    for (int k = Z0; k < Z1; ++k)
	    {
		    field.E[i][j][k] -= J[i][j][k];
        }
    }

    void ClearJ()
    {
        for (int i = 0; i < LX; ++i)
	    for (int j = 0; j < LY; ++j)
	    for (int k = 0; k < LZ; ++k)
	    {
		    J[i][j][k].Zero();
        }
    }

    
    void CalcOnCenter(Field& field)
    {
        for (int i = 0; i < LX-1; ++i)
	    for (int j = 0; j < LY-1; ++j)
	    for (int k = 0; k < LZ-1; ++k)
	    {
		    Ec[i][j][k].x = 0.5 * (field.E[i][j][k].x + field.E[i+1][j][k].x);
            Ec[i][j][k].y = 0.5 * (field.E[i][j][k].y + field.E[i][j+1][k].y);
            Ec[i][j][k].z = 0.5 * (field.E[i][j][k].z + field.E[i][j][k+1].z);

		    Bc[i][j][k].x = 0.25 * (field.B[i][j][k].x + field.B[i][j+1][k].x + field.B[i][j][k+1].x + field.B[i][j+1][k+1].x);
	    	Bc[i][j][k].y = 0.25 * (field.B[i][j][k].y + field.B[i+1][j][k].y + field.B[i][j][k+1].y + field.B[i+1][j][k+1].y);
	    	Bc[i][j][k].z = 0.25 * (field.B[i][j][k].z + field.B[i+1][j][k].z + field.B[i][j+1][k].z + field.B[i+1][j+1][k].z);
	    }
    }
        
    void BunemanBoris(Plasma& plasma)
    {
        Vector u_minus, u_, u_plus, u, t, s;
        Vector e, b;
        double gamma;
        
        std::vector<Particle> &p = plasma.p;
        for (unsigned long int n = 0; n < p.size(); ++n)
        {
            double sx[3] = {};
            double sy[3] = {};
            double sz[3] = {};

            sf3(sx, p[n].r.x);
            sf3(sy, p[n].r.y);
            sf3(sz, p[n].r.z);

            int ii = int(p[n].r.x) - 1;
            int jj = int(p[n].r.y) - 1;
            int kk = int(p[n].r.z) - 1;

            e.Zero();
            b.Zero();

            for(int i = 0; i < 3; ++i)
            for(int j = 0; j < 3; ++j)
            for(int k = 0; k < 3; ++k)
            {
                e += Ec[ii+i][jj+j][kk+k] * sx[i]*sy[j]*sz[k];
                b += Bc[ii+i][jj+j][kk+k] * sx[i]*sy[j]*sz[k];
            }
            
            b *= 0.5 * (plasma.q / plasma.m);     // q/m * B * delt_t/2
            e *= 0.5 * (plasma.q / plasma.m);     // q/m * E * delt_t/2
#if 1            
            gamma = C / sqrt(C2 - p[n].v.Mag2());

            u_minus = gamma * p[n].v + e;

            gamma = C / sqrt(C2 + u_minus.Mag2());

            b *= gamma;

            gamma = 2.0 / (1.0 + b.Mag2());

            u_plus = (u_minus + u_minus.CrossProduct(b)) * gamma;

            u_minus += u_plus.CrossProduct(b) + e;

            gamma = C / sqrt(C2 + u_minus.Mag2());

            p[n].v = u_minus * gamma;
#else
            gamma = C / sqrt(C2 - p[n].v.Mag2());

            u_minus = gamma * p[n].v + e;

            t = b * C / sqrt(C2 + u_minus.Mag2());

            u_ = u_minus + u_minus.CrossProduct(t);

            s = 2.0 / (1.0 + t.Mag2()) * t;

            u_plus = u_minus + u_.CrossProduct(s);

            u = u_plus + e;

            gamma = C / sqrt(C2 + u.Mag2());

            p[n].v = u * gamma;
#endif
        } 
    }
};
#endif
