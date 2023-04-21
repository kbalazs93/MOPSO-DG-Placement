/* Routines for storing population data into files */

# include <stdio.h>
# include <stdlib.h>
# include <math.h>
# include <iostream>
# include <iomanip>

# include "../NetworkModel/network.h"
# include "nsga2.h"
# include "rand.h"

using namespace std;
/* Balazs fuggveny: Ez amelyik visszaadja a Hypervolume metrikat*/
/* EZ a fuggveny CSAK AZ ELSO 2 objektumot irassa ki. Ha  osszetettebbet akarsz,
Modositsd */

double nsga2_optimization::calculate2Dhypervolume( population *pop) // MSB - activeLosses, middle - VSI, LSB - vD
{
    double first_ref = obj_empty[0];
    double second_ref = obj_empty[1];

    int feasnum = 0;
    for(int i=0; i<nsga2Params.popsize; i++)
      if (pop->ind[i].constr_violation == 0.0 &&  pop->ind[i].rank==1 &&
        pop->ind[i].obj[0] < first_ref && pop->ind[i].obj[1] < second_ref)
      {
        feasnum++;
      }

    double first_func[feasnum];
    double second_func[feasnum];


    int j = 0;
    for (int i=0; i<nsga2Params.popsize; i++)
    {
      if (pop->ind[i].constr_violation == 0.0 && pop->ind[i].rank==1 &&
      pop->ind[i].obj[0] < first_ref && pop->ind[i].obj[1] < second_ref)
      {
        first_func[j] = pop->ind[i].obj[0];
        second_func[j++] = pop->ind[i].obj[1];
      }
    }

    for(int i=0; i<feasnum; i++)
    {
      for(int j=i; j>0 && first_func[j-1] > first_func[j]; j--)
      {
        double cont;
        cont = first_func[j];
        first_func[j] = first_func[j-1];
        first_func[j-1] = cont;
        cont = second_func[j];
        second_func[j] = second_func[j-1];
        second_func[j-1] = cont;
      }
    }

    double hypervolume;
    hypervolume = ((first_ref - first_func[0])/first_ref) * ((second_ref - second_func[0])/second_ref) ;
    for (int i=1; i<feasnum; ++i)
    {
      hypervolume += ((first_ref - first_func[i])/first_ref) * ((second_func[i-1] - second_func[i])/second_ref) ;
    }

    return hypervolume;
}

/* Balazs fuggveny: Ez amelyik visszaadja a spacing metrikat*/
/* EZ a fuggveny CSAK AZ ELSO 2 objektumot irassa ki. Ha  osszetettebbet akarsz,
Modositsd */
double nsga2_optimization::calculateSpacing( population *pop)
{
  double first_ref = obj_empty[0];
  double second_ref = obj_empty[1];

  int feasnum = 0;
  for(int i=0; i<nsga2Params.popsize; i++)
    if (pop->ind[i].constr_violation == 0.0  &&  pop->ind[i].rank==1 &&
      pop->ind[i].obj[0] < first_ref && pop->ind[i].obj[1] < second_ref)
      {
        feasnum++;
      }

  double first_func[feasnum];
  double second_func[feasnum];


  int j = 0;
  for (int i=0; i<nsga2Params.popsize; i++)
  {
    if (pop->ind[i].constr_violation == 0.0 && pop->ind[i].rank==1 &&
    pop->ind[i].obj[0] < first_ref && pop->ind[i].obj[1] < second_ref)
    {
      first_func[j] = pop->ind[i].obj[0];
      second_func[j++] = pop->ind[i].obj[1];
    }
  }

  for(int i=0; i<feasnum; i++)
  {
    for(int j=i; j>0 && first_func[j-1] > first_func[j]; j--)
    {
      double cont;
      cont = first_func[j];
      first_func[j] = first_func[j-1];
      first_func[j-1] = cont;
      cont = second_func[j];
      second_func[j] = second_func[j-1];
      second_func[j-1] = cont;
    }
  }

  double d[feasnum];
  double d_below[feasnum];

  d[0] = abs((first_func[1] - first_func[0])/first_ref) + abs((second_func[1] - second_func[0])/second_ref);
  d_below[0] = d[0];
  for (int i=1; i<feasnum-1; i++)
  {
    d_below[i] = abs((first_func[i+1] - first_func[i])/first_ref) + abs((second_func[i+1] - second_func[i])/second_ref);
    d[i] = fmin(d_below[i], d_below[i-1]);
  }
  d[feasnum - 1] = d_below[feasnum - 2];

  double d_mean = 0;
  for(int i=0; i<feasnum; i++)
  {
    d_mean += d[i];
  }
  d_mean /= feasnum;

  double devSum = 0;
  for(int i=0; i<feasnum; i++)
  {
    devSum += pow(d_mean - d[i], 2);
  }

  return sqrt((1.0/(feasnum-1)) * devSum);
}

int nsga2_optimization::circle_choice( population *pop)
{
  double first_base = obj_empty[0];
  double second_base = obj_empty[1];

  double base_point[2];
  base_point[0] = first_base;
  base_point[1] = second_base;

  int fitness_vals = 2;

  double ref_point[fitness_vals];
  for (int i=0; i<fitness_vals; i++)
  {
    if (pop->ind[i].constr_violation == 0.0 && pop->ind[i].rank==1)
    {
      ref_point[i] = pop->ind[0].obj[i];
      for(int j=0; j<nsga2Params.popsize; j++)
      {
        if( pop->ind[j].obj[i] < ref_point[i] )
          ref_point[i] = pop->ind[j].obj[i];
      }
    }
  }

  int min_ind = 0;
  double minimum_radius_square=-1;

  for (int i=0; i<nsga2Params.popsize; ++i)
  {
    if (pop->ind[i].constr_violation == 0.0 && pop->ind[i].rank==1 )
    {
      double radius_squere = 0;
      for (int j=0; j<fitness_vals; j++)
        radius_squere += pow((pop->ind[i].obj[j] - ref_point[j]) / base_point[j],2);
      if (minimum_radius_square == -1)
        minimum_radius_square = radius_squere;
      else
        if (radius_squere < minimum_radius_square)
        {
          minimum_radius_square = radius_squere;
          min_ind = i;
        }
      }
  }
  return min_ind;
}


/* BALAZS fuggveny, amelyik kivalassza a fuzzy best reszecsket es kiirja az eredmenyeit*/
/* EZ a fuggveny CSAK AZ ELSO 2 objektumot irassa ki. Ha  osszetettebbet akarsz,
Modositsd */
int nsga2_optimization::fuzzy_choice(  population *pop)
{
  double first_ref = obj_empty[0];
  double second_ref = obj_empty[1];

  int fitness_vals = 2;
  double f_min[2], f_max[2];

  int feasnum = 0;
  for(int i=0; i<nsga2Params.popsize; i++)
    if (pop->ind[i].constr_violation == 0.0 && pop->ind[i].rank==1 &&
      pop->ind[i].obj[0] < first_ref && pop->ind[i].obj[1] < second_ref)
      {
        feasnum++;
      }

  double obj1[feasnum], obj2[feasnum];


  int j = 0;
  for (int i=0; i<nsga2Params.popsize; i++)
  {
    if (pop->ind[i].constr_violation == 0.0 && pop->ind[i].rank==1 &&
      pop->ind[i].obj[0] < first_ref && pop->ind[i].obj[1] < second_ref)
    {
        obj1[j] = pop->ind[i].obj[0];
        obj2[j++] = pop->ind[i].obj[1];
    }
  }

  for(int i=0; i<fitness_vals; i++)
  {
    double* cur_obj;
    if (i==0)
      cur_obj = obj1;
    else if (i==1)
      cur_obj = obj2;

    f_min[i] = cur_obj[0];
    f_max[i] = cur_obj[0];
    for(int j=0; j<feasnum; j++)
    {
      if (cur_obj[j] < f_min[i])
        f_min[i] = cur_obj[j];
      if (cur_obj[j] > f_max[i])
        f_max[i] = cur_obj[j];
    }
  }

  double alphamax=0;
  int alphamax_ind=0;
  double alpha[feasnum];
  for(int i=0; i<feasnum; i++)
  {
    alpha[i] = 0;
    for(int j=0; j<fitness_vals; j++)
    {
      double* cur_obj;
      if (j==0)
        cur_obj = obj1;
      else if (j==1)
        cur_obj = obj2;

      if (cur_obj[i] == f_min[j])
        alpha[i] += 1.0;
      else if (cur_obj[i] == f_max[j])
        continue;
      else
        alpha[i] += ( f_max[j] - cur_obj[i] ) / (f_max[j] - f_min[j]) ;
      if (alpha[i] > alphamax)
      {
        alphamax = alpha[i];
        alphamax_ind = i;
      }
    }
  }
  return alphamax_ind;
}

int nsga2_optimization::fuzzy_choice_within_range( population *pop)
{
  double first_ref = obj_empty[0];
  double second_ref = obj_empty[1];

  int fitness_vals = 2;
  double f_min[2], f_max[2];

  int feasnum = 0;
  for(int i=0; i<nsga2Params.popsize; i++)
    if (pop->ind[i].constr_violation == 0.0 && pop->ind[i].rank==1 &&
      pop->ind[i].obj[0] < first_ref && pop->ind[i].obj[1] < second_ref)
      {
        feasnum++;
      }

  double obj1[feasnum], obj2[feasnum];
  int original_ind[feasnum];


  int j = 0;
  for (int i=0; i<nsga2Params.popsize; i++)
  {
    if (pop->ind[i].constr_violation == 0.0 && pop->ind[i].rank==1
      && pop->ind[i].obj[0] < first_ref && pop->ind[i].obj[1] < second_ref)
    {
        obj1[j] = pop->ind[i].obj[0];
        obj2[j] = pop->ind[i].obj[1];
        original_ind[j++] = i;
    }
  }

  for(int i=0; i<fitness_vals; i++)
  {
    double* cur_obj;
    if (i==0)
      cur_obj = obj1;
    else if (i==1)
      cur_obj = obj2;

    f_min[i] = cur_obj[0];
    f_max[i] = cur_obj[0];
    for(int j=0; j<feasnum; j++)
    {
      if (cur_obj[j] < f_min[i])
        f_min[i] = cur_obj[j];
      if (cur_obj[j] > f_max[i])
        f_max[i] = cur_obj[j];
    }
  }

  double alphamax=0;
  int alphamax_ind=0;
  double alpha[feasnum];
  for(int i=0; i<feasnum; i++)
  {
    alpha[i] = 0;
    for(int j=0; j<fitness_vals; j++)
    {
      double* cur_obj;
      if (j==0)
        cur_obj = obj1;
      else if (j==1)
        cur_obj = obj2;

      if (cur_obj[i] == f_min[j])
        alpha[i] += 1.0;
      else if (cur_obj[i] == f_max[j])
        continue;
      else
        alpha[i] += ( f_max[j] - cur_obj[i] ) / (f_max[j] - f_min[j]) ;
      if (alpha[i] > alphamax)
      {
        alphamax = alpha[i];
        alphamax_ind = i;
      }
    }
  }

  return original_ind[alphamax_ind];
}

int nsga2_optimization::calculateFeasnum(population *pop)
{
  feasnum = 0;
  for(int i=0; i<nsga2Params.popsize; ++i)
  {
    if (pop->ind[i].constr_violation == 0.0 && pop->ind[i].rank==1 )
    {
        feasnum++;
    }
  }
  return feasnum;
}

void nsga2_optimization::calculateMinMax(population *pop)
{
  if (f_min != NULL)
    delete[] f_min;
  if (f_max != NULL)
    delete[] f_max;
  f_min = new double[nsga2Params.nobj];
  f_max = new double[nsga2Params.nobj];

  f_min[0] = pop->ind[0].obj[0];
  f_min[1] = pop->ind[0].obj[1];
  f_max[0] = pop->ind[0].obj[0];
  f_max[1] = pop->ind[0].obj[1];

  for(int i=0; i<nsga2Params.popsize; i++)
  {
    if (pop->ind[i].constr_violation == 0.0 &&  pop->ind[i].rank==1 )
    {
      if (pop->ind[i].obj[0] < f_min[0])
        f_min[0] = pop->ind[i].obj[0];
      if (pop->ind[i].obj[0] > f_max[0])
        f_max[0] = pop->ind[i].obj[0];

      if (pop->ind[i].obj[1] < f_min[1])
        f_min[1] = pop->ind[i].obj[1];
      if (pop->ind[i].obj[1] > f_max[1])
        f_max[1] = pop->ind[i].obj[1];
    }
  }
}

void nsga2_optimization::log_table_form( FILE *fpt)
{
  population *pop = parent_pop;

  network empty_netw; // Needed, becasue of convertIndex_LayerB2Longitudinal()

  int my_best = circle_choice( pop);


  int fitness_vals = 2;
  double f_min[fitness_vals];
  f_min[0] = pop->ind[0].obj[0];
  f_min[1] = pop->ind[0].obj[1];
  double f_max[fitness_vals];
  f_max[0] = pop->ind[0].obj[0];
  f_max[1] = pop->ind[0].obj[1];

  for(int i=0; i<nsga2Params.popsize; i++)
  {
    if (pop->ind[i].constr_violation == 0.0 &&  pop->ind[i].rank==1 )
    {
      if (pop->ind[i].obj[0] < f_min[0])
        f_min[0] = pop->ind[i].obj[0];
      if (pop->ind[i].obj[0] > f_max[0])
        f_max[0] = pop->ind[i].obj[0];

      if (pop->ind[i].obj[1] < f_min[1])
        f_min[1] = pop->ind[i].obj[1];
      if (pop->ind[i].obj[1] > f_max[1])
        f_max[1] = pop->ind[i].obj[1];
    }
  }


  fprintf(fpt, "Circle Best Losses: %lf\n", pop->ind[my_best].obj[0]*1e3 );
  fprintf(fpt, "Circle Best VSI: %lf\n", 1/pop->ind[my_best].obj[1] );

  fprintf(fpt, "\n\n");

  fprintf(fpt, "Particle position:\n");
  fprintf(fpt, "Node notation: study notation (layer based notation)\n");
  for(int i=0; i<nsga2Params.nreal; i++)
  {
    if ((nsga2Params.nreal - 1) % 3 == 0 && i < nsga2Params.nreal - 1)
    {
      if (i % 3 == 0)
        fprintf(fpt, "%d (%lf)  ", empty_netw.convertIndex_LayerB2Longitudinal(round(pop->ind[my_best].xreal[i])), round(pop->ind[my_best].xreal[i]));
      else if (i % 3 == 1)
        fprintf(fpt, "P: %lf   ", pop->ind[my_best].xreal[i]);
      else if (i % 3 == 2)
        fprintf(fpt, "Q: %lf   ", pop->ind[my_best].xreal[i]);
    }
    else if (i < nsga2Params.nreal - 1)
    {
      if (i % 2 == 0)
        fprintf(fpt, "%d (%lf)  ", empty_netw.convertIndex_LayerB2Longitudinal(round(pop->ind[my_best].xreal[i])), round(pop->ind[my_best].xreal[i]));
      else if (i % 2 == 1)
        fprintf(fpt, "%lf,   ", pop->ind[my_best].xreal[i]);
    }
  }
  fprintf(fpt, "%lf\n", round(pop->ind[my_best].xreal[nsga2Params.nreal - 1]));

  fprintf(fpt, "\n\n\n" );

  fprintf(fpt, "Active Losses range: %lf - %lf\n", f_min[0]*1e3, f_max[0]*1e3);
  fprintf(fpt, "VSI range: %lf - %lf\n", 1/f_max[1], 1/f_min[1]);

  fprintf(fpt, "\n" );



  double first_ref = obj_empty[0];
  double second_ref = obj_empty[1];

  int feasnum = 0;
  for(int i=0; i<nsga2Params.popsize; ++i)
  {
    if (pop->ind[i].constr_violation == 0.0 && pop->ind[i].rank==1 )
    {
      if (pop->ind[i].obj[0] < first_ref
        && pop->ind[i].obj[1] < second_ref)
        {
          feasnum++;
        }
    }
  }
  fprintf(fpt, "reference active losses: %lf\n", first_ref*1e3);
  fprintf(fpt, "reference VSI: %lf\n", 1/second_ref);
  fprintf(fpt, "\n");

  fprintf(fpt, "hypervolume: %lf\n", calculate2Dhypervolume( pop));
  fprintf(fpt, "feasable solutions: %d/%d\n", feasnum, nsga2Params.popsize);
  fprintf(fpt, "spacing: %lf\n", calculateSpacing( pop));
}

void nsga2_optimization::log_master_form(ofstream& os)
{
    population* pop = parent_pop ;
    int my_best = circle_choice(pop);


    int dgNum = bin_params[0];
    int dgType = bin_params[1];
    int hasTapChanger = bin_params[3];
    int load_flow_type = bin_params[2];

    int var_num;
    if(dgType == 1 || dgType == 2)
      var_num = 2;
    else if(dgType == 3)
      var_num = 3;

    os<<fixed;
    int sorted_indexes[dgNum];
    for(int i=0; i<dgNum; i++)
      sorted_indexes[i] = i;



    for(int i=0; i<dgNum; i++)
      for(int j=0; j<dgNum-i-1; j++)
      {
        if (int(round(pop->ind[my_best].xreal[sorted_indexes[j]*var_num])) > int(round(pop->ind[my_best].xreal[sorted_indexes[j+1]*var_num])))
          swap_int(sorted_indexes+j, sorted_indexes+j+1);
      }

    if(hasTapChanger)
      os<<-int(round(pop->ind[my_best].xreal[nsga2Params.nreal-1]))<<" ";
    for(int i=0; i<dgNum; i++)
    {
      os<<int(round(pop->ind[my_best].xreal[sorted_indexes[i]*var_num]));
      if(i<dgNum-1)
        os<<",";
    }
    os<<" ";

    if (dgType == 3)
      os<<"P:";

    if(dgType == 1 || dgType == 3)
    {
      for(int i=0; i<dgNum; i++)
      {
        os<<setprecision(3)<<pop->ind[my_best].xreal[sorted_indexes[i]*var_num+1]*1e-6;
        if(i<dgNum-1)
          os<<",";
      }
      os<<" ";
    }

    if (dgType == 3)
      os<<"Q:";

    if(dgType == 2 || dgType == 3)
    {
      int qoffset = 1;
      if(dgType == 3)
        qoffset = 2;
      for(int i=0; i<dgNum; i++)
      {
        os<<setprecision(3)<<pop->ind[my_best].xreal[sorted_indexes[i]*var_num+qoffset]*1e-6;
        if(i<dgNum-1)
          os<<",";
      }
      os<<" ";
    }

    os<<setprecision(2)<<pop->ind[my_best].obj[0]*1e3<<" ";
    os<<setprecision(2)<<1.0/pop->ind[my_best].obj[1]<<" ";

    os<<endl;
}


void nsga2_optimization::createParetoOctaveGraph(FILE *fpt)
{
  population *pop = parent_pop;

  //double Vmin = 0.95;
  //double Vmax = 1.05;
  //int load_flow_type = 2;  // 1 - const power, 2 - const current, 3 - mixed model
  //network netw("./NetworkModel/IEEE33.txt", Vmin, Vmax, load_flow_type);
  //netw.solvePowerFlow("summer", 12);

  double first_ref = obj_empty[0];
  double second_ref = obj_empty[1];

  int feasnum = 0;
  for(int i=0; i<nsga2Params.popsize; i++)
    if (pop->ind[i].constr_violation == 0.0 && pop->ind[i].rank==1)// &&
  //    pop->ind[i].obj[0] < first_ref && pop->ind[i].obj[1] < second_ref)
    {
      feasnum++;
    }

  double first_func[feasnum];
  double second_func[feasnum];


  int j = 0;
  for (int i=0; i<nsga2Params.popsize; i++)
  {
    if (pop->ind[i].constr_violation == 0.0 && pop->ind[i].rank==1 )//&&
  //  pop->ind[i].obj[0] < first_ref && pop->ind[i].obj[1] < second_ref)
    {
      first_func[j] = pop->ind[i].obj[0];
      second_func[j++] = pop->ind[i].obj[1];
    }
  }

  fprintf(fpt, "p_front_nsga = [...\n");
  for (int i=0; i<feasnum-1; i++)
  {
    fprintf(fpt, "\t%lf %lf\n", first_func[i]*1e3, 1/second_func[i]);
  }
  fprintf(fpt, "\t%lf %lf\n", first_func[feasnum-1]*1e3, 1/second_func[feasnum-1]);

  fprintf(fpt, "\n\n\n" );
  fprintf(fpt, "ref_point = [%lf  %lf]", first_ref, second_ref);
}


void nsga2_optimization::report_pop ( FILE *fpt)
{
    fprintf(fpt,"# of objectives = %d, # of constraints = %d, # of real_var = %d, # of bits of bin_var = %d, constr_violation, rank, crowding_distance\n",nsga2Params.nobj,nsga2Params.ncon,nsga2Params.nreal,nsga2Params.bitlength);
    population *pop = parent_pop;

    int i, j, k;
    for (i=0; i<nsga2Params.popsize; i++)
    {
        for (j=0; j<nsga2Params.nobj; j++)
        {
            fprintf(fpt,"%e\t",pop->ind[i].obj[j]);
        }
        if (nsga2Params.ncon!=0)
        {
            for (j=0; j<nsga2Params.ncon; j++)
            {
                fprintf(fpt,"%e\t",pop->ind[i].constr[j]);
            }
        }
        if (nsga2Params.nreal!=0)
        {
            for (j=0; j<nsga2Params.nreal; j++)
            {
                fprintf(fpt,"%e\t",pop->ind[i].xreal[j]);
            }
        }
        if (nsga2Params.nbin!=0)
        {
            for (j=0; j<nsga2Params.nbin; j++)
            {
                for (k=0; k<nsga2Params.nbits[j]; k++)
                {
                    fprintf(fpt,"%d\t",pop->ind[i].gene[j][k]);
                }
            }
        }
        fprintf(fpt,"%e\t",pop->ind[i].constr_violation);
        fprintf(fpt,"%d\t",pop->ind[i].rank);
        fprintf(fpt,"%e\n",pop->ind[i].crowd_dist);
    }
    return;
}

/* Function to print the information of feasible and non-dominated population in a file */
void nsga2_optimization::report_feasible ( FILE *fpt)
{
    fprintf(fpt,"# of objectives = %d, # of constraints = %d, # of real_var = %d, # of bits of bin_var = %d, constr_violation, rank, crowding_distance\n",nsga2Params.nobj,nsga2Params.ncon,nsga2Params.nreal,nsga2Params.bitlength);
    population *pop = parent_pop;
    int i, j, k;
    for (i=0; i<nsga2Params.popsize; i++)
    {
        if (pop->ind[i].constr_violation == 0.0 && pop->ind[i].rank==1)
        {
            for (j=0; j<nsga2Params.nobj; j++)
            {
                fprintf(fpt,"%e\t",pop->ind[i].obj[j]);
            }
            if (nsga2Params.ncon!=0)
            {
                for (j=0; j<nsga2Params.ncon; j++)
                {
                    fprintf(fpt,"%e\t",pop->ind[i].constr[j]);
                }
            }
            if (nsga2Params.nreal!=0)
            {
                for (j=0; j<nsga2Params.nreal; j++)
                {
                    fprintf(fpt,"%e\t",pop->ind[i].xreal[j]);
                }
            }
            if (nsga2Params.nbin!=0)
            {
                for (j=0; j<nsga2Params.nbin; j++)
                {
                    for (k=0; k<nsga2Params.nbits[j]; k++)
                    {
                        fprintf(fpt,"%d\t",pop->ind[i].gene[j][k]);
                    }
                }
            }
            fprintf(fpt,"%e\t",pop->ind[i].constr_violation);
            fprintf(fpt,"%d\t",pop->ind[i].rank);
            fprintf(fpt,"%e\n",pop->ind[i].crowd_dist);
        }
    }
    return;
}

void nsga2_optimization::log_MinMax(ofstream& os, int i)
{
  os<<fixed;
  calculateMinMax(parent_pop);
  double printed_min, printed_max;
  printed_min = f_min[i];
  printed_max = f_max[i];
  if (i==0)
  {
    printed_min *= 1e3;
    printed_max *= 1e3;
  }
  else if (i==1)
  {
    printed_min = 1.0/f_max[i];
    printed_max = 1.0/f_min[i];
  }
  os<<setprecision(2)<<printed_min<<"-"<<setprecision(2)<<printed_max<<endl;
}
