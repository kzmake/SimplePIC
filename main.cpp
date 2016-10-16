#include <cstdio>
#include <cmath>

#include "Param.hpp"

#include "Plasma.hpp"
#include "Field.hpp"
#include "Solver.hpp"
#include "ShapeFactor.hpp"

#include "Input.hpp"
#include "Output.hpp"

int main()
{
    Field f;
    std::vector<Plasma> p;
    Solver<TSC> s;
    //Solver<CIC> s;

    Input(p, f);

    Output(p, f, 0);

    for(int ts = 1; ts <= MAX_TIME_STEP; ++ts)
    {
        printf("%d\n", ts);

        f.UpdateB();
        f.BoundaryB();
        
        s.CalcOnCenter(f);
        for(unsigned int n = 0; n < p.size(); ++n)
        {
            s.BunemanBoris(p[n]);
        }

        f.UpdateB();
        f.BoundaryB();

        f.UpdateE();


        for(unsigned int n = 0; n < p.size(); ++n)
        {
            //s.VillasenorBuneman(p[n]);
            s.DensityDecomposition(p[n], f);
        }

#if 0
        Vector I1;
        I1.Zero();
        for(int n = 0; n < p.size(); ++n)
        {
            for(int k = 0; k < p[n].p.size(); ++k)
            {
                I1 += p[n].q * p[n].p[k].v;
            }
        }
        printf(" qv( %f, %f, %f)\n", I1.x, I1.y, I1.z);
#endif
        s.BoundaryJ();
        s.UpdateEbyJ(f);

        f.BoundaryE();

        for(unsigned int n = 0; n < p.size(); ++n)
        {  
            p[n].UpdateR();
            p[n].BoundaryR();
        }

        if (ts % OUTPUT_STEP == 0) Output(p, f, ts);
    } 

    return 0;
}
