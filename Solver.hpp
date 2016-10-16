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

  public:
    Solver()
    {
        try
        {
            J  = CreateArray(LX, LY, LZ);
            Ec = CreateArray(LX, LY, LZ);
            Bc = CreateArray(LX, LY, LZ);
        }
        catch(std::bad_alloc)
        {
            fprintf(stderr, "エラー: bad_alloc 動的メモリ確保エラー: %s\n", __FILE__); 
            DeleteArray(J);
            DeleteArray(Ec);
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

        for (unsigned long int n = 0; n < p.size(); ++n)
        {
            Vector r0 = p[n].r;
            Vector r1 = p[n].r + p[n].v;
            double s0[3][5] = {};
#if 1
            sf5(s0[0], r0.x);
            sf5(s0[1], r0.y);
            sf5(s0[2], r0.z);
#else
            sf5(s0[0], r0.x, int(r0.x) - int(r1.x));
            sf5(s0[1], r0.y, int(r0.y) - int(r1.y));
            sf5(s0[2], r0.z, int(r0.z) - int(r1.z));
#endif

            
            double ds[3][5] = {};
#if 1
            sf5(ds[0], r1.x, int(r1.x) - int(r0.x));
            sf5(ds[1], r1.y, int(r1.y) - int(r0.y));
            sf5(ds[2], r1.z, int(r1.z) - int(r0.z));
#else
            sf5(ds[0], r1.x);
            sf5(ds[1], r1.y);
            sf5(ds[2], r1.z);
#endif

            for(int i = 0; i < 5; ++i)
            {
                ds[0][i] -= s0[0][i];
                ds[1][i] -= s0[1][i];
                ds[2][i] -= s0[2][i];
            }

            double w[5][5][5][3] = {};

            const double w1_3 = (1.0/3.0);

            for(int i = 0; i < 5; ++i)
            for(int j = 0; j < 5; ++j)
            for(int k = 0; k < 5; ++k)
            {
                w[i][j][k][0] = ds[0][i] * (s0[1][j]*s0[2][k] + 0.5*ds[1][j]*s0[2][k] + 0.5*s0[1][j]*ds[2][k] + w1_3*ds[1][j]*ds[2][k]);
                w[i][j][k][1] = ds[1][j] * (s0[0][i]*s0[2][k] + 0.5*ds[0][i]*s0[2][k] + 0.5*s0[0][i]*ds[2][k] + w1_3*ds[0][i]*ds[2][k]);
                w[i][j][k][2] = ds[2][k] * (s0[0][i]*s0[1][j] + 0.5*ds[0][i]*s0[1][j] + 0.5*s0[0][i]*ds[1][j] + w1_3*ds[0][i]*ds[1][j]);

            }

            Vector J0[5][5][5] = {};

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
#if 0
        Vector I0;
        I0.Zero();
        for (int i = X0; i < X1; ++i)
	    for (int j = Y0; j < Y1; ++j)
	    for (int k = Z0; k < Z1; ++k)
	    {
		    I0 += J[i][j][k];
        }

        printf("  I( %f, %f, %f)\n", I0.x, I0.y, I0.z);
#endif
        // periodic - Add
        // X
        for (int j = 0; j < LY; ++j)
        for (int k = 0; k < LZ; ++k)
        {
            J[X0+0][j][k] += J[X1+0][j][k];
            J[X0+1][j][k] += J[X1+1][j][k];
            J[X1-1][j][k] += J[X0-1][j][k];
            J[X1-2][j][k] += J[X0-2][j][k];
        }
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
#if 0
        Vector I0;
        I0.Zero();
        for (int i = X0; i < X1; ++i)
	    for (int j = Y0; j < Y1; ++j)
	    for (int k = Z0; k < Z1; ++k)
	    {
		    I0 += J[i][j][k];
        }

        printf("  I( %f, %f, %f)\n", I0.x, I0.y, I0.z);
#endif

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
            if (n == 100)
            {
                //printf("  -> (%ld) V = %f, %f, %f  (%f)\n", n, p[n].v.x, p[n].v.y, p[n].v.z, p[n].v.Mag2());
            }
            
            // shape factor で　rからe, b を取得
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
            
            b *= 0.5 * (plasma.q / plasma.m) / C; // q/m * Bc * delt_t/2
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

    void VillasenorBuneman(Plasma& plasma)
    {
        for (int i = 0; i < LX; ++i)
	    for (int j = 0; j < LY; ++j)
	    for (int k = 0; k < LZ; ++k)
	    {
		    J[i][j][k].Zero();
        }

        Vector r;
        std::vector<Particle> &p = plasma.p;
        //printf("  >plasma.q = %f, plasma.p.size(%ld)\n", plasma.q, p.size());

        for (unsigned long int n = 0; n < p.size(); ++n)
        {
            r = p[n].r + p[n].v;

            SplitX(r, p[n].r, plasma.q);
        }
    }

    void SplitX(const Vector& r, const Vector& r0, const double q)
    {
        if (int(r.x + 0.5) != int(r0.x + 0.5))
        {
            Vector r1;
            r1.x = 0.5 * (int(r.x + 0.5) + int(r0.x + 0.5));
            r1.y = r0.y + (r.y - r0.y) * ((r1.x - r0.x) / (r.x - r0.x));
            r1.z = r0.z + (r.z - r0.z) * ((r1.x - r0.x) / (r.x - r0.x));

            SplitY(r,  r1, q);
            SplitY(r1, r0, q);
        }
        else
        {
            SplitY(r,  r0, q);
        }
    }
    
    void SplitY(const Vector& r, const Vector& r0, const double q)
    {
        if (int(r.y + 0.5) != int(r0.y + 0.5))
        {
            Vector r1;
            r1.y = 0.5 * (int(r.y + 0.5) + int(r0.y + 0.5));
            r1.z = r0.z + (r.z - r0.z) * ((r1.y - r0.y) / (r.y - r0.y));
            r1.x = r0.x + (r.x - r0.x) * ((r1.y - r0.y) / (r.y - r0.y));

            SplitZ(r,  r1, q);
            SplitZ(r1, r0, q);
        }
        else
        {
            SplitZ(r,  r0, q);
        }
    }

    void SplitZ(const Vector& r, const Vector& r0, const double q)
    {
        if (int(r.z + 0.5) != int(r0.z + 0.5))
        {
            Vector r1;
            r1.z = 0.5 * (int(r.z + 0.5) + int(r0.z + 0.5));
            r1.x = r0.x + (r.x - r0.x) * ((r1.z - r0.z) / (r.z - r0.z));
            r1.y = r0.y + (r.y - r0.y) * ((r1.z - r0.z) / (r.z - r0.z));

            Depsit(r,  r1, q);
            Depsit(r1, r0, q);
        }
        else
        {
            Depsit(r,  r0, q);
        }
    }

    void Depsit(const Vector& r, const Vector& r0, const double q)
    {
        Vector mid = 0.5 * (r + r0);

        int gi = int(mid.x + 0.5);
        int gj = int(mid.y + 0.5);
        int gk = int(mid.z + 0.5);

        Vector d, c;
        d.x = mid.x - gi + 0.5;
        d.y = mid.y - gj + 0.5;
        d.z = mid.z - gk + 0.5;

        c.x = 1.0 - d.x;
        c.y = 1.0 - d.y;
        c.z = 1.0 - d.z;

        Vector diff = (r - r0);
        Vector flux = q * diff;
        double delt = q * diff.x * diff.y * diff.z / 12;

        int i = 0;
        int j = 0;
        int k = 0;

        i = gi;
        j = gj;
        k = gk;

        J[i][j  ][k+2].x += flux.x * d.y * d.z + delt;
        J[i][j-1][k+2].x += flux.x * c.y * d.z - delt;
        J[i][j  ][k+1].x += flux.x * d.y * c.z - delt;
        J[i][j-1][k+1].x += flux.x * c.y * c.z + delt;

        J[i  ][j][k+2].y += flux.y * d.z * d.x + delt;
        J[i  ][j][k+1].y += flux.y * c.z * d.x - delt;
        J[i-1][j][k+2].y += flux.y * d.z * c.x - delt;
        J[i-1][j][k+1].y += flux.y * c.z * c.x + delt;

        J[i  ][j  ][k+1].z += flux.z * d.x * d.y + delt;
        J[i-1][j  ][k+1].z += flux.z * c.x * d.y - delt;
        J[i  ][j-1][k+1].z += flux.z * d.x * c.y - delt;
        J[i-1][j-1][k+1].z += flux.z * c.x * c.y + delt;
        
    }
};

#endif
