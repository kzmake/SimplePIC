#ifndef SIMPLE_PIC_OUTPUT_HPP
#define SIMPLE_PIC_OUTPUT_HPP

#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/json_parser.hpp>
#include <boost/foreach.hpp>
#include <boost/optional.hpp>

template<Shape SF>
void OutputProfile(std::vector<Plasma>& plasma, Field& field, Solver<SF>& solver)
{
    int rank, size;
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &size);
    
    if (rank == 0)
    {
        std::string filename("");
        filename += PATH + "profile.json";

        // JSON
        boost::property_tree::ptree pt;
        {
            pt.put("MPI.size", size);

            pt.put("PIC.density",    NUM_DENS);
            pt.put("PIC.C"      ,           C);
            pt.put("PIC.rseed"  , RANDOM_SEED);

            pt.put("PIC.shapefactor", SF);

            pt.put("PIC.LX0", LX0);
            pt.put("PIC.LY0", LY0);
            pt.put("PIC.LZ0", LZ0);
            
            pt.put("PIC.LX", LX);
            pt.put("PIC.LY", LY);
            pt.put("PIC.LZ", LZ);

#ifdef PIC_PML
            pt.put("PIC.PML.L",    L);
            pt.put("PIC.PML.M",    M);
            pt.put("PIC.PML.R0",  R0);

            pt.put("PIC.PML.SIGMA_MAX", SIGMA_MAX);
#endif

            pt.put("PIC.timestep.max" , MAX_TIME_STEP);
            pt.put("PIC.timestep.step",   OUTPUT_STEP);

            boost::property_tree::ptree particles;
            {
                boost::property_tree::ptree ele;
                ele.put("mass", ELE_MASS);
                ele.put("wpe" ,  ELE_WPE);
                ele.put("q"   ,    ELE_Q);
                ele.put("vth" ,  ELE_VTH);

                particles.push_back(std::make_pair("ele", ele));
            }
            {
                boost::property_tree::ptree ion;
                ion.put("mass", ION_MASS);
                ion.put("wpe" ,  ION_WPE);
                ion.put("q"   ,    ION_Q);
                ion.put("vth" ,  ION_VTH);

                particles.push_back(std::make_pair("ion", ion));
            }

            pt.add_child("PIC.particle", particles);
            pt.put("PIC.particle.size", plasma[0].p.size());
        }
        
        boost::property_tree::write_json(filename, pt);
    }
}

template<Shape SF>
void Output(std::vector<Plasma>& plasma, Field& field, Solver<SF>& solver, Timer& t, const int ts)
{
    int rank;
    MPI_Comm_rank(comm, &rank);
    
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

    Vector energyK[2];
    double energyR[2];
    for(unsigned int s = 0; s < plasma.size(); ++s)
    {
        Vector v2;

        energyK[s].Zero();
        energyR[s] = 0.0;
        std::vector<Particle> &p = plasma[s].p;
        for (unsigned long int n = 0; n < p.size(); ++n)
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

    auto MPIReduceSumEnergy = [](Vector& energy)
    {
        Vector totalEnergy;

        MPI_Allreduce(&energy, &totalEnergy, 3, MPI_DOUBLE, MPI_SUM, comm);

        memcpy(&energy, &totalEnergy, sizeof(Vector));
    };

    auto MPIReduceSumEnergyR = [](double& energy)
    {
        double totalEnergy;

        MPI_Allreduce(&energy, &totalEnergy, 1, MPI_DOUBLE, MPI_SUM, comm);

        memcpy(&energy, &totalEnergy, sizeof(double));
    };

    MPIReduceSumEnergy(energyB);
    MPIReduceSumEnergy(energyE);

    for(unsigned int s = 0; s < plasma.size(); ++s)
    {
        MPIReduceSumEnergy(energyK[s]);
        MPIReduceSumEnergyR(energyR[s]);
    }

    auto WriteEnegy = [](FILE*& fp, const std::string& filename, const Vector& v, const int ts)
    {
        if(fp == nullptr) fp = fopen(filename.c_str(), "w+");

        fprintf(fp, "%d\t%e\t%e\t%e\n", ts, v.x, v.y, v.z);

        if(ts == MAX_TIME_STEP) fclose(fp);
    };

    auto WriteEnegyR = [](FILE*& fp, const std::string& filename, const double& v, const int ts)
    {
        if(fp == nullptr) fp = fopen(filename.c_str(), "w+");

        fprintf(fp, "%d\t%e\n", ts, v);

        if(ts == MAX_TIME_STEP) fclose(fp);
    };

    auto WriteField = [](const std::string s, Vector*** V, const int ts, const int rank)
    {
        char cts[8], crank[4];
        snprintf(cts, sizeof(cts), "%06d", ts);
        snprintf(crank, sizeof(crank), "%02d", rank);
        std::string filename("");
        filename += PATH + s + "/" + s + cts + "_r" + crank;
        FILE *fp;
        fp = fopen(filename.c_str(), "w+");

        fwrite(&V[0][0][0], sizeof(double) * 3, LX*LY*LZ, fp);
        
        fclose(fp);
    };

    if (rank == 0)
    {
        static const std::string kFilenameEnergyIonsK = PATH + "energy_ions_k" + ".txt";
        static const std::string kFilenameEnergyElesK = PATH + "energy_eles_k" + ".txt";
        static const std::string kFilenameEnergyIonsR = PATH + "energy_ions_r" + ".txt";
        static const std::string kFilenameEnergyElesR = PATH + "energy_eles_r" + ".txt";
        static const std::string kFilenameEnergyE = PATH + "energy_f_e" + ".txt";
        static const std::string kFilenameEnergyB = PATH + "energy_f_b" + ".txt";

        static FILE* energyIonsK_fp = nullptr;
        static FILE* energyElesK_fp = nullptr;
        static FILE* energyIonsR_fp = nullptr;
        static FILE* energyElesR_fp = nullptr;
        static FILE* energyE_fp     = nullptr;
        static FILE* energyB_fp     = nullptr;


        WriteEnegy(energyB_fp, kFilenameEnergyB, energyB, ts);
        WriteEnegy(energyE_fp, kFilenameEnergyE, energyE, ts);
        WriteEnegy(energyIonsK_fp, kFilenameEnergyIonsK, energyK[0], ts);
        WriteEnegy(energyElesK_fp, kFilenameEnergyElesK, energyK[1], ts);
        WriteEnegyR(energyIonsR_fp, kFilenameEnergyIonsR, energyR[0], ts);
        WriteEnegyR(energyElesR_fp, kFilenameEnergyElesR, energyR[1], ts);
    }

    if (ts % OUTPUT_STEP == 0)
    {
        WriteField("f_e",  field.E, ts, rank);
        WriteField("f_b",  field.B, ts, rank);
        WriteField("f_j", solver.J, ts, rank);
    }

    auto WriteTimer = [](FILE*& fp, const std::string& filename, const std::string t, const int ts)
    {
        if(fp == nullptr) fp = fopen(filename.c_str(), "w+");

        fprintf(fp, "%d\t%s\n", ts, t.c_str());

        if(ts == MAX_TIME_STEP) fclose(fp);
    };

    if (rank == 0)
    {
        static const std::string kTimer = PATH + "timer" + ".txt";
        static FILE* timer_fp = nullptr;

        WriteTimer(timer_fp, kTimer,  t[0].format(3, "%w") + t[1].format(3, "\t%w"), ts);
    }
}

#endif
