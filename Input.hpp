#ifndef SIMPLE_PIC_INPUT_HPP
#define SIMPLE_PIC_INPUT_HPP

#include <cstdio>
#include <cmath>
#include <vector>
#include <random>
#include <limits>

#if 0
template<BoundaryCondition BP, BoundaryCondition BF>
void Input(std::vector<Plasma<BP>>& p, Field<BF>& f)
{
    int rank;
    MPI_Comm_rank(comm, &rank);

    int seed = RANDOM_SEED + rank;

    const double wpe = ELE_WPE;

    Plasma<BP> ions;
	Plasma<BP> eles;

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

    long int id = rank;
    for (int i = X0; i < X1; ++i)
	for (int j = Y0; j < Y1; ++j)
    for (int k = Z0; k < Z1; ++k)
    for (int n = 0; n < NUM_DENS; ++n)
	{
		r.x = i + dist(mt);
		r.y = j + dist(mt);
        r.z = k + dist(mt);

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

		id += 100;
	}

    p.push_back(ions);
    p.push_back(eles);
}

#else
#if 0
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
	}

    p.push_back(ions);
    p.push_back(eles);
}
#else
#if 0
// PML
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
        f.E[i][j][k].x = gaussian(k - LZ/2, 7);
        f.B[i][j][k].y = gaussian(k - LZ/2, 7) / C;
    }

    p.push_back(ions);
    p.push_back(eles);
}
#else
// GradZero
void Input(std::vector<Plasma>& p, Field& f)
{
    int rank;
    MPI_Comm_rank(comm, &rank);
    
    int seed = RANDOM_SEED + rank;

    const double wpe = ELE_WPE;

    Plasma ions;
	Plasma eles;

    eles.m = ELE_MASS;
    ions.m = ION_MASS;

    eles.q = - wpe * sqrt(eles.m / (double)NUM_DENS);
    ions.q = - eles.q;

    // Ampere
    constexpr double j0 = 0.01 * C;
    constexpr double vj = j0 / (NUM_DENS * ELE_WPE * sqrt(ELE_MASS / NUM_DENS)); // 0.0145
    constexpr double vji = 0.5 * vj;
    constexpr double vje = 0.5 * vj;

    constexpr double vRasio = vj;// / ELE_VTH;
    constexpr double skin_de = C / ELE_WPE;
    constexpr double pRasio = 0.25 * (vRasio * vRasio) / (skin_de * skin_de) * Filament::R2;

    printf(" vRasio = %f\n", vRasio);
    printf(" vth(max) = %f\n", sqrt(pRasio));

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

    long int id = rank;
    for (int i = X0; i < X1; ++i)
	for (int j = LY/2 - Filament::R; j <= LY/2 + Filament::R; ++j)
    for (int k = LZ/2 - Filament::R; k <= LZ/2 + Filament::R; ++k)
    for (int n = 0; n < NUM_DENS; ++n)
	{        
		r.x = i + dist(mt);
		r.y = j + dist(mt);
        r.z = k + dist(mt);

        double radius2 = (r.z - LZ/2) * (r.z - LZ/2) + (r.y - LY/2) * (r.y - LY/2);

        if (radius2 > Filament::R2) continue;

        double vte = sqrt(pRasio - 0.25 * (vRasio * vRasio) / (skin_de * skin_de) * radius2);
        double vti = vte * std::sqrt(ELE_MASS / ION_MASS);
        double vme = sqrt(1.5) * vte;
        double vmi = sqrt(1.5) * vti;
        double std_devi_e = vme / sqrt(3.0);
        double std_devi_i = vmi / sqrt(3.0);

ION_VELO:
        vi.x = randomBoxMuller(dist(mt), dist(mt));
		vi.y = randomBoxMuller(dist(mt), dist(mt));
		vi.z = randomBoxMuller(dist(mt), dist(mt));
        vi *= std_devi_i;
        //vi *= 0.0;

        vi.x += vji;
        
        if (vi.Mag2() > C2) goto ION_VELO;

ELE_VELO:
        ve.x = randomBoxMuller(dist(mt), dist(mt));
		ve.y = randomBoxMuller(dist(mt), dist(mt));
		ve.z = randomBoxMuller(dist(mt), dist(mt));
        ve *= std_devi_e;
        //ve *= 0.0;
    
        ve.x -= vje;

        if (ve.Mag2() > C2) goto ELE_VELO;
        
        Particle ion(id, r, vi);
		Particle ele(id, r, ve);

		ions.p.push_back(ion);
		eles.p.push_back(ele);


		id += 100;
	}

    // Ampare B
    for (int i = 0; i < LX; ++i)
	for (int j = 0; j < LY; ++j)
	for (int k = 0; k < LZ; ++k)
    {
        double radius_z = sqrt((k + 0.5 - LZ/2)*(k + 0.5 - LZ/2) + (j - LY/2)*(j - LY/2));
        double radius_y = sqrt((k - LZ/2)*(k - LZ/2) + (j + 0.5 - LY/2)*(j + 0.5 - LY/2));

        double b0z, b0y;

        if (radius_z * radius_z < Filament::R2)
            b0z = 0.5 * j0 / C2 * radius_z;
        else
            b0z = 0.5 * j0 / C2 * Filament::R2 / radius_z;

        if (radius_y * radius_y < Filament::R2)
            b0y = 0.5 * j0 / C2 * radius_y;
        else
            b0y = 0.5 * j0 / C2 * Filament::R2 / radius_y;

        double radz = atan2((j - LY/2), (k + 0.5 - LZ/2));
        double rady = atan2((j + 0.5 - LY/2), (k - LZ/2));

        f.B[i][j][k].z = + b0z * sin(radz);
        f.B[i][j][k].y = - b0y * cos(rady);

        //f.B[i][j][k].z *= 0.5 * (tanh((radius_z + 320)/ 10) - tanh((radius_z - 320) / 10));
        //f.B[i][j][k].y *= 0.5 * (tanh((radius_y + 320)/ 10) - tanh((radius_y - 320) / 10));
    }


    p.push_back(ions);
    p.push_back(eles);
}
#endif
#endif
#endif

#endif

