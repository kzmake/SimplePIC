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
    
    Timer t(2);
    
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
        
        f.UpdateB();
        f.BoundaryB();

        s.CalcOnCenter(f);
        for(unsigned int n = 0; n < p.size(); ++n)
        {
            s.BunemanBoris(p[n]);
        }

        f.UpdateB();
        
#ifdef PIC_PML
        f.BoundaryB_PML();
#endif
        f.BoundaryB();
        
        f.UpdateE();

        s.ClearJ();

        t[1].start();
        for(unsigned int n = 0; n < p.size(); ++n)
        {
            printf(" p.size (%ld)\n", p[n].p.size());
            s.DensityDecomposition(p[n], f);
        }
        t[1].stop();

        s.BoundaryJ();
        s.UpdateEbyJ(f);

#ifdef PIC_PML
        f.BoundaryE_PML();
#endif  
        f.BoundaryE();

        for(unsigned int n = 0; n < p.size(); ++n)
        {  
            p[n].UpdateR();
            p[n].BoundaryR();
        }

        MPI::COMM_WORLD.Barrier();

        t[0].stop();

        if (ts % SORT_STEP == 0)
        for(unsigned int n = 0; n < p.size(); ++n)
        {  
            p[n].Sort();
        }

        std::string result = t[0].format(3, "total：%ws") + t[1].format(3, " | DensityDecomposition：%ws");

        if (MPI::COMM_WORLD.Get_rank() == 0)
        {
            printf("  %s\n", result.c_str());
        }
        
        Output(p, f, s, t, ts);
    }

    MPI::Finalize();

    return 0;
}
