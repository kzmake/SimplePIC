#include <cstdio>
#include <cmath>
#include <vector>

#include <boost/timer/timer.hpp>

#include <mpi.h>

using Timer = std::vector<boost::timer::cpu_timer>;

#include "Param.hpp"

#include "Plasma.hpp"
#include "Field.hpp"
#include "Solver.hpp"
#include "ShapeFactor.hpp"

#include "Input.hpp"
#include "Output.hpp"

int main(int argc, char **argv)
{

    MPI::Init(argc, argv);

    Field f;
    std::vector<Plasma> p;
    Solver<TSC> s;
    //Solver<CIC> s;
    
    Timer t(6);
    
    {
        t[0].start();
        Input(p, f);
        t[0].stop();

        OutputProfile(p, f, s);
        Output(p, f, s, t, 0);
    }

    for(int ts = 1; ts <= MAX_TIME_STEP; ++ts)
    {
        MPI::COMM_WORLD.Barrier();

        t[0].start();

        if (MPI::COMM_WORLD.Get_rank() == 0)
        {
            printf("%d\n", ts);
        }
        
        t[3].start();

        f.UpdateB();
        f.BoundaryB();

        t[3].stop();
        t[2].start();

        s.CalcOnCenter(f);
        for(unsigned int n = 0; n < p.size(); ++n)
        {
            s.BunemanBoris(p[n]);
        }

        t[2].stop();
        t[3].resume();

        f.UpdateB();
        
#ifdef PIC_PML
        f.BoundaryB_PML();
#endif
        f.BoundaryB();
        
        f.UpdateE();

        t[3].stop();
        t[1].start();

        s.ClearJ();
        for(unsigned int n = 0; n < p.size(); ++n)
        {
            s.DensityDecomposition(p[n], f);
        }
        
        t[1].stop();
        t[3].resume();

        s.BoundaryJ();
        s.UpdateEbyJ(f);

#ifdef PIC_PML
        f.BoundaryE_PML();
#endif  
        f.BoundaryE();
        t[3].stop();
        t[2].resume();

        for(unsigned int n = 0; n < p.size(); ++n)
        {  
            p[n].UpdateR();
            p[n].BoundaryR();
        }
        
        t[2].stop();

        MPI::COMM_WORLD.Barrier();

        t[0].stop();
        t[4].start();
        
        if (ts % SORT_STEP == 0)
        for(unsigned int n = 0; n < p.size(); ++n)
        {  
            p[n].Sort();
        }
        t[4].stop();

        if (MPI::COMM_WORLD.Get_rank() == 0)
        {
            std::string result = t[0].format(3, "total：%ws")
                + t[1].format(3, " | DensityDecomposition：%ws")
                + t[2].format(3, " | Mover & Pusher：%ws")
                + t[3].format(3, " | FDTD：%ws")
                + t[4].format(3, " | Sort：%ws");
            printf("  %s\n", result.c_str());
        }
        
        t[5].start();
        Output(p, f, s, t, ts);
        t[5].stop();

        if (MPI::COMM_WORLD.Get_rank() == 0)
        {
            std::string result = t[5].format(3, "output：%ws");
            printf("  %s\n", result.c_str());
        }


    }

    MPI::Finalize();

    return 0;
}

