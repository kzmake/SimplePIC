#include <cstdio>
#include <cmath>

#include <boost/timer/timer.hpp>

#ifdef MPI_PIC
#include <mpi.h>
#endif

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

    {
        boost::timer::cpu_timer timer;
        
        Input(p, f);

        OutputProfile(p, f, s);
        timer.stop();

        Output(p, f, s, timer, 0);
    }

    for(int ts = 1; ts <= MAX_TIME_STEP; ++ts)
    {
        MPI::COMM_WORLD.Barrier();
        boost::timer::cpu_timer timer;

        if (MPI::COMM_WORLD.Get_rank() == 0) printf("%d\n", ts);
        
        f.UpdateB();
        //f.BoundaryB();

        s.CalcOnCenter(f);
        for(unsigned int n = 0; n < p.size(); ++n)
        {
            s.BunemanBoris(p[n]);
        }

        f.UpdateB();
        f.BoundaryB();

        f.UpdateE();

        s.ClearJ();
        for(unsigned int n = 0; n < p.size(); ++n)
        {
            s.DensityDecomposition(p[n], f);
        }

        s.BoundaryJ();
        s.UpdateEbyJ(f);
        
        f.BoundaryE();

        for(unsigned int n = 0; n < p.size(); ++n)
        {  
            p[n].UpdateR();
            p[n].BoundaryR();
        }

        MPI::COMM_WORLD.Barrier();
        timer.stop();

        if (ts % SORT_STEP == 0)
        for(unsigned int n = 0; n < p.size(); ++n)
        {  
            p[n].Sort();
        }

        std::string result = timer.format(3, "total：%ws | user：%us | system：%ss (CPU: %p)");
        if (MPI::COMM_WORLD.Get_rank() == 0) printf("  %s\n", result.c_str());
        
        Output(p, f, s, timer, ts);
    }

    MPI::Finalize();

    return 0;
}
