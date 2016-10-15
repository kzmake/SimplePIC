#ifndef SIMPLE_PIC_INPUT_HPP
#define SIMPLE_PIC_INPUT_HPP

#include <cmath>
#include <vector>
#include <random>

#if 0
void Input(std::vector<Plasma>& p, Field& f)
{
    const long int seed = 1000;

    const int NUM_DENS = 1;
    const int MASS_RATIO = 100;
    const double wpe = 0.50;

    Plasma ions;
	Plasma eles;

    eles.m = 0.0625;
    ions.m = eles.m * MASS_RATIO;

    eles.q = - wpe  * sqrt(eles.m / (double)NUM_DENS);
    ions.q = - eles.q;

    printf("ions.q = %f eles.q = %f\n", ions.q, eles.q);

    double vte = 0.1 * C;
    double vti = sqrt(eles.m / ions.m) * vte;
    double vme = sqrt(1.5) * vte;
    double vmi = sqrt(1.5) * vti;
    double std_devi_e = vme / sqrt(3.0);
    double std_devi_i = vmi / sqrt(3.0);    Vector r, vi, ve;
    std::mt19937 mt;
    mt.seed(seed);
    std::uniform_real_distribution<double> dist(0.0, 1.0);
    long int id;

    for (int k = Z0; k < Z1; ++k)
	for (int j = Y0; j < Y1; ++j)
	for (int i = X0; i < X1; ++i)
	{
        double x = (i - 2 - LX0/2); 
        double y = (j - 2 - LY0/2);
        double R = x * x + y * y;

        f.E[i][j][k].z = 1.0 / exp(R);
        //if(k == Z0 + LZ0/2) printf("%e \t", 1 / exp(R));
	}

    p.push_back(ions);
    p.push_back(eles);
}
#else
void Input(std::vector<Plasma>& p, Field& f)
{
    const long int seed = 1000;

    const int NUM_DENS = 5;
    const int MASS_RATIO = 100;
    const double wpe = 0.50;

    Plasma ions;
	Plasma eles;

    eles.m = 0.0625;
    ions.m = eles.m * MASS_RATIO;

    eles.q = - wpe  * sqrt(eles.m / (double)NUM_DENS);
    ions.q = - eles.q;

    printf("ions.q = %f eles.q = %f\n", ions.q, eles.q);

    double vte = 0.1 * C;
    double vti = sqrt(eles.m / ions.m) * vte;
    double vme = sqrt(1.5) * vte;
    double vmi = sqrt(1.5) * vti;
    double std_devi_e = vme / sqrt(3.0);
    double std_devi_i = vmi / sqrt(3.0);

    Vector r, vi, ve;
    std::mt19937 mt;
    mt.seed(seed);
    std::uniform_real_distribution<double> dist(0.0, 1.0);
    long int id = 0;

    // Box-Muller method
    auto randomBoxMuller = [](double a, double b)
    {
        double alpha = 1.0 - a;
        double beta  = b;

        return sqrt(-2.0*log(alpha)) * cos(2.0*M_PI*beta);
    };
#if 1
    for (int k = 0; k < LZ0; ++k)
	for (int j = 0; j < LY0; ++j)
	for (int i = 0; i < LX0; ++i)
    for (int n = 0; n < NUM_DENS; ++n)
	{
		r.x = X0 + i + dist(mt);
		r.y = Y0 + j + dist(mt);
        r.z = Z0 + k + dist(mt);

        vi.x = randomBoxMuller(dist(mt), dist(mt));
		vi.y = randomBoxMuller(dist(mt), dist(mt));
		vi.z = randomBoxMuller(dist(mt), dist(mt));

        ve.x = randomBoxMuller(dist(mt), dist(mt));
		ve.y = randomBoxMuller(dist(mt), dist(mt));
		ve.z = randomBoxMuller(dist(mt), dist(mt));

        vi *= std_devi_i;
        ve *= std_devi_e;

        Particle ion(id, r, vi);
		Particle ele(id, r, ve);

		ions.p.push_back(ion);
		eles.p.push_back(ele);

		++id;
	}
#else
    {
        r.x = LX/2 + dist(mt);
        r.y = LY/2 + dist(mt);
        r.z = LZ/2 + dist(mt);

        vi.x = randomBoxMuller(dist(mt), dist(mt));
        vi.y = randomBoxMuller(dist(mt), dist(mt));
        vi.z = randomBoxMuller(dist(mt), dist(mt));

        ve.x = randomBoxMuller(dist(mt), dist(mt));
        ve.y = randomBoxMuller(dist(mt), dist(mt));
        ve.z = randomBoxMuller(dist(mt), dist(mt));

        vi *= std_devi_i;
        ve *= std_devi_e;
#if 0
        vi.x = 0.0;//randomBoxMuller(dist(mt), dist(mt));
        vi.y = 0.0;//randomBoxMuller(dist(mt), dist(mt));
        vi.z = 0.1;//randomBoxMuller(dist(mt), dist(mt));

        ve.x = 0.0;//randomBoxMuller(dist(mt), dist(mt));
        ve.y = 0.0;//randomBoxMuller(dist(mt), dist(mt));
        ve.z = -0.1;//randomBoxMuller(dist(mt), dist(mt));
#endif
        if (vi.Mag2() > C) printf("vi > C\n");
        if (ve.Mag2() > C) printf("ve > C\n");

        printf(" -- #(%ld) R (%f, %f, %f)\n", id, r.x, r.y, r.z);
        printf(" --        Vi(%f, %f, %f)\n", vi.x, vi.y, vi.z);
        printf(" --        Ve(%f, %f, %f)\n", ve.x, ve.y, ve.z);

        Particle ion(id, r, vi);
        Particle ele(id, r, ve);

        ions.p.push_back(ion);
		eles.p.push_back(ele);
    }

#endif

    for (int k = Z0; k < Z1; ++k)
	for (int j = Y0; j < Y1; ++j)
	for (int i = X0; i < X1; ++i)
	{
        //f.B[i][j][k].z = 0.1;
	}

    p.push_back(ions);
    p.push_back(eles);
}
#endif

#endif
