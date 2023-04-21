/* Some utility functions (not part of the algorithm) */

# include <stdio.h>
# include <stdlib.h>
# include <math.h>

# include "nsga2.h"
# include "rand.h"

/* Function to return the maximum of two variables */
double nsga2_optimization::maximum (double a, double b)
{
    if (a>b)
    {
        return(a);
    }
    return (b);
}

/* Function to return the minimum of two variables */
double nsga2_optimization::minimum (double a, double b)
{
    if (a<b)
    {
        return (a);
    }
    return (b);
}

void nsga2_optimization::swap_int(int* i, int* j)
{
  int helper = *i;
  *i = *j;
  *j = helper;
}
