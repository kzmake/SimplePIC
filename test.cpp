#include <cstdio>
#include <cmath>
#include "ShapeFactor.hpp"

int main(int argc, char** argv)
{
    ShapeFactor<TSC, 5> sf;

    double s[3][5];

    double rx_ = 0.25;
    double rx = 1.25;

    sf(s[1], rx, int(rx) - int(rx_));

    double sum = .0;
    for(int j = 0; j < 3; ++j)
    {
        for(int i = 0; i < 5; ++i)
        {
            sum += s[j][i];
            printf("%0.5f  ", s[j][i]);
        }
        printf("\n");
    }

    printf("\n");

    printf("sum = %0.5f\n", sum);
    
    return 0;
}
