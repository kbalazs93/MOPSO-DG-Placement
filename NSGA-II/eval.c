/* Routine for evaluating population members  */

# include <stdio.h>
# include <stdlib.h>
# include <math.h>

# include "nsga2.h"
# include "rand.h"

/* Routine to evaluate objective function values and constraints for a population */
void nsga2_optimization::evaluate_pop (population *pop, void *inp, void *out)
{
    int i;
    for (i=0; i<nsga2Params.popsize; i++)
    {
        evaluate_ind (&(pop->ind[i]), inp, out);
    }
    return;
}

/* Routine to evaluate objective function values and constraints for an individual */
void nsga2_optimization::evaluate_ind (individual *ind, void *inp, void *out)
{
    int j;
    test_problem (ind->xreal, ind->xbin, ind->gene, ind->obj, ind->constr);

    if (nsga2Params.ncon==0)
    {
        ind->constr_violation = 0.0;
    }
    else
    {
        ind->constr_violation = 0.0;
        for (j=0; j<nsga2Params.ncon; j++)
        {
            if (ind->constr[j]<0.0)
            {
                ind->constr_violation += ind->constr[j];
            }
        }
    }
    return;
}
