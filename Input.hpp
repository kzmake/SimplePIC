#ifndef SIMPLE_PIC_INPUT_HPP
#define SIMPLE_PIC_INPUT_HPP

#include <cstdio>
#include <cmath>
#include <vector>
#include <random>

#if 0
void Input(std::vector<Plasma>& p, Field& f)
{
    int seed = RANDOM_SEED + MPI::COMM_WORLD.Get_rank();

    const double wpe = ELE_WPE;

    Plasma ions;
	Plasma eles;

    eles.m = ELE_MASS;
    ions.m = ION_MASS;

    eles.q = - wpe * sqrt(eles.m / (double)NUM_DENS);
    ions.q = - eles.q;

    double vte = ELE_VTH;
    double vti = ION_VTH;
    double vme = sqrt(1.5) * vte;
    double vmi = sqrt(1.5) * vti;
    double std_devi_e = vme / sqrt(3.0);
    double std_devi_i = vmi / sqrt(3.0);

    Vector r, vi, ve;
    std::mt19937 mt;
    mt.seed(seed);
    std::uniform_real_distribution<double> dist(0.0, 1.0);

    // Box-Muller method
    auto randomBoxMuller = [](double a, double b)
    {
        double alpha = 1.0 - a;
        double beta  = b;

        return sqrt(-2.0*log(alpha)) * cos(2.0*M_PI*beta);
    };

    long int id = MPI::COMM_WORLD.Get_rank() * 1000000;
    for (int i = 0; i < LX0; ++i)
	for (int j = 0; j < LY0; ++j)
    for (int k = 0; k < LZ0; ++k)
    for (int n = 0; n < NUM_DENS; ++n)
	{
		r.x = X0 + i + dist(mt);
		r.y = Y0 + j + dist(mt);
        r.z = Z0 + k + dist(mt);

ION_VELO:
        vi.x = randomBoxMuller(dist(mt), dist(mt));
		vi.y = randomBoxMuller(dist(mt), dist(mt));
		vi.z = randomBoxMuller(dist(mt), dist(mt));
        vi *= std_devi_i;
        
        if (vi.Mag2() > C2) goto ION_VELO;

ELE_VELO:
        ve.x = randomBoxMuller(dist(mt), dist(mt));
		ve.y = randomBoxMuller(dist(mt), dist(mt));
		ve.z = randomBoxMuller(dist(mt), dist(mt));
        ve *= std_devi_e;

        if (ve.Mag2() > C2) goto ELE_VELO;
        
        Particle ion(id, r, vi);
		Particle ele(id, r, ve);

		ions.p.push_back(ion);
		eles.p.push_back(ele);

		++id;
	}

    for (int i = X0; i < X1; ++i)
	for (int j = Y0; j < Y1; ++j)
	for (int k = Z0; k < Z1; ++k)
	{
        //f.B[i][j][k].z = 0.1;
	}

    p.push_back(ions);
    p.push_back(eles);
}

#else
void Input(std::vector<Plasma>& p, Field& f)
{
    Plasma ions;
	Plasma eles;

    auto gaussian = [](double x, double sigma)
    {
        return 1.0/sqrt(2.0*M_PI*sigma) * exp(- (x * x) / (2*sigma*sigma));
    };

    for (int i = 0; i < LX; ++i)
	for (int j = 0; j < LY; ++j)
	for (int k = 0; k < LZ; ++k)
    {
        f.E[i][j][k].x = gaussian(j - LY/2, 6) * gaussian(k - LZ/2, 6);
        //f.E[i][j][k].x = gaussian(k - LZ/2);
	}

    p.push_back(ions);
    p.push_back(eles);
}
#endif

#endif

