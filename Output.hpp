#ifndef SIMPLE_PIC_OUTPUT_HPP
#define SIMPLE_PIC_OUTPUT_HPP

static const std::string kFilenameEnergyIonsK = PATH + "energy_ions_k" + ".txt";
static const std::string kFilenameEnergyElesK = PATH + "energy_eles_k" + ".txt";
static const std::string kFilenameEnergyE = PATH + "energy_f_e" + ".txt";
static const std::string kFilenameEnergyB = PATH + "energy_f_b" + ".txt";

FILE* energyIonsK_fp = nullptr;
FILE* energyElesK_fp = nullptr;
FILE* energyE_fp     = nullptr;
FILE* energyB_fp     = nullptr;


void Output(std::vector<Plasma>& plasma, Field& field, const int ts)
{
    Vector energyB, energyE;

    energyB.Zero();
    energyE.Zero();
    for (int i = X0; i < X1; ++i)
    for (int j = Y0; j < Y1; ++j)
    for (int k = Z0; k < Z1; ++k)
    {
        energyB.x += field.B[i][j][k].x * field.B[i][j][k].x;
        energyB.y += field.B[i][j][k].y * field.B[i][j][k].y;
        energyB.z += field.B[i][j][k].z * field.B[i][j][k].z;

        energyE.x += field.E[i][j][k].x * field.E[i][j][k].x;
        energyE.y += field.E[i][j][k].y * field.E[i][j][k].y;
        energyE.z += field.E[i][j][k].z * field.E[i][j][k].z;
    }

    energyB *= 0.5;
    energyE *= 0.5;



    Vector v2;

    Vector energyK[2];
    double energyR[2];
    for(int s = 0; s < plasma.size(); ++s)
    {
        energyK[s].Zero();
        std::vector<Particle> &p = plasma[s].p;
        for (long int n = 0; n < p.size(); ++n)
        {
            v2.x = p[n].v.x * p[n].v.x;
            v2.y = p[n].v.y * p[n].v.y;
            v2.z = p[n].v.z * p[n].v.z;

            energyK[s] += v2;
            energyR[s] += (C / sqrt(C2 - (v2.x + v2.y + v2.z)) - 1.0);
        }

        energyK[s] *= (0.5 * plasma[s].m);
        energyR[s] *= (plasma[s].m * C2);
    }


    auto WriteEnegy = [](FILE*& fp, const std::string& filename, const Vector& v, const int ts)
    {
        if(fp == nullptr) 
            fp = fopen(filename.c_str(), "w+");

        fprintf(fp, "%d \t %e \t %e \t %e \t\n", ts, v.x, v.y, v.z);

        if(ts == MAX_TIME_STEP) fclose(fp);
    };

    WriteEnegy(energyB_fp, kFilenameEnergyB, energyB, ts);
    WriteEnegy(energyE_fp, kFilenameEnergyE, energyE, ts);
    WriteEnegy(energyIonsK_fp, kFilenameEnergyIonsK, energyK[0], ts);
    WriteEnegy(energyElesK_fp, kFilenameEnergyElesK, energyK[1], ts);



    auto WriteField = [](const std::string s, Vector*** V, const int ts)
    {
        char cts[8];
        sprintf(cts, "%06d", ts);
        std::string filename("");
        filename += PATH + s + cts + ".txt";
        FILE *fp;
        fp = fopen(filename.c_str(), "w+");

        for (int i = X0; i < X1; ++i)
        for (int j = X0; j < Y1; ++j)
        //for (int k = 0; k < LZ; ++k)
        {
            int k = Z0 + LZ0/2;
            fprintf(fp, "%e \t %e \t %e \t\n", V[i][j][k].x, V[i][j][k].y, V[i][j][k].z);
        }
        
        fclose(fp);
    };

    WriteField("f_e", field.E, ts);
    WriteField("f_b", field.B, ts);


    double ions = energyK[0].x + energyK[0].y + energyK[0].z;
    double eles = energyK[1].x + energyK[1].y + energyK[1].z;
    double b = energyB.x + energyB.y + energyB.z;
    double e = energyE.x + energyE.y + energyE.z;
    //printf("=====\n");
    printf("E:%f  B:%f  ions:%f  eles%f  (%f)\n", e, b, ions, eles, ions + eles + e + b);
    //printf("Total:%f\t\n", ions + eles + e + b);

    
}

#endif
