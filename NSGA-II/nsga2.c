/* NSGA-II routine (implementation of the 'main' function) */

# include "nsga2.h"
# include "rand.h"
#include <string>

using namespace std;

nsga2_optimization::nsga2_optimization(int argc, char **argv, double* real_params, int real_params_num, int* bin_params, int bin_params_num)
{
  void *inp = NULL;
  void *out = NULL;

  f_min = NULL;
  f_max = NULL;

  set_real_params(real_params, real_params_num);
  set_bin_params(bin_params, bin_params_num);

  ReadParameters(argc, argv);

  InitNSGA2( inp, out);
}

nsga2_optimization::nsga2_optimization(double seed, int popsize, int ngen,
  int nobj, int ncon, int nreal, double* min_realvar, double* max_realvar,
  double pcross_real, double pmut_real, double distr_ind_cross_real, double distr_ind_mut_real,
   int nbin, int* nbits, double* min_binvar, double* max_binvar, double pcross_bin, double pmut_bin,
   double* real_params, int real_params_num, int* bin_params, int bin_params_num,
   int choice, int obj1, int obj2, int obj3, int angle1, int angle2)
{
  set_real_params(real_params, real_params_num);
  set_bin_params(bin_params, bin_params_num);

  f_min = NULL;
  f_max = NULL;

  nsga2Params.seed = seed;
  if (nsga2Params.seed<=0.0 || nsga2Params.seed>=1.0)
  {
      printf("\n Entered seed value is wrong, seed value must be in (0,1) \n");
      exit(1);
  }

  nsga2Params.popsize = popsize;
  if (nsga2Params.popsize<4 || (nsga2Params.popsize%4)!= 0)
  {
      printf("\n population size read is : %d",nsga2Params.popsize);
      printf("\n Wrong population size entered, hence exiting \n");
      exit (1);
  }

  nsga2Params.ngen = ngen;
  if (nsga2Params.ngen<1)
  {
      printf("\n number of generations read is : %d",nsga2Params.ngen);
      printf("\n Wrong nuber of generations entered, hence exiting \n");
      exit (1);
  }

  nsga2Params.nobj = nobj;
  if (nsga2Params.nobj<1)
  {
      printf("\n number of objectives entered is : %d",nsga2Params.nobj);
      printf("\n Wrong number of objectives entered, hence exiting \n");
      exit (1);
  }

  nsga2Params.ncon = ncon;
  if (nsga2Params.ncon<0)
  {
      printf("\n number of constraints entered is : %d",nsga2Params.ncon);
      printf("\n Wrong number of constraints enetered, hence exiting \n");
      exit (1);
  }

  nsga2Params.nreal = nreal;
  if (nsga2Params.nreal<0)
  {
      printf("\n number of real variables entered is : %d",nsga2Params.nreal);
      printf("\n Wrong number of variables entered, hence exiting \n");
      exit (1);
  }
  if (nsga2Params.nreal != 0)
  {
      nsga2Params.min_realvar = (double *)malloc(nsga2Params.nreal*sizeof(double));
      nsga2Params.max_realvar = (double *)malloc(nsga2Params.nreal*sizeof(double));
      for (int i=0; i<nsga2Params.nreal; i++)
      {
          nsga2Params.min_realvar[i] = min_realvar[i];
          nsga2Params.max_realvar[i] = max_realvar[i];
          if (nsga2Params.max_realvar[i] <= nsga2Params.min_realvar[i])
          {
              printf("\n Wrong limits entered for the min and max bounds of real variable, hence exiting \n");
              exit(1);
          }
      }

      nsga2Params.pcross_real = pcross_real;
      if (nsga2Params.pcross_real<0.0 || nsga2Params.pcross_real>1.0)
      {
          printf("\n Probability of crossover entered is : %e",nsga2Params.pcross_real);
          printf("\n Entered value of probability of crossover of real variables is out of bounds, hence exiting \n");
          exit (1);
      }

      nsga2Params.pmut_real = pmut_real;
      if (nsga2Params.pmut_real<0.0 || nsga2Params.pmut_real>1.0)
      {
          printf("\n Probability of mutation entered is : %e",nsga2Params.pmut_real);
          printf("\n Entered value of probability of mutation of real variables is out of bounds, hence exiting \n");
          exit (1);
      }

      nsga2Params.eta_c = distr_ind_cross_real;
      if (nsga2Params.eta_c<=0)
      {
          printf("\n The value entered is : %e",nsga2Params.eta_c);
          printf("\n Wrong value of distribution index for crossover entered, hence exiting \n");
          exit (1);
      }

      nsga2Params.eta_m = distr_ind_mut_real;
      if (nsga2Params.eta_m<=0)
      {
          printf("\n The value entered is : %e",nsga2Params.eta_m);
          printf("\n Wrong value of distribution index for mutation entered, hence exiting \n");
          exit (1);
      }
  }

  nsga2Params.nbin = nbin;
  if (nsga2Params.nbin<0)
  {
      printf ("\n number of binary variables entered is : %d",nsga2Params.nbin);
      printf ("\n Wrong number of binary variables entered, hence exiting \n");
      exit(1);
  }

  if (nsga2Params.nbin != 0)
  {
      nsga2Params.nbits = (int *)malloc(nsga2Params.nbin*sizeof(int));
      nsga2Params.min_binvar = (double *)malloc(nsga2Params.nbin*sizeof(double));
      nsga2Params.max_binvar = (double *)malloc(nsga2Params.nbin*sizeof(double));
      for (int i=0; i<nsga2Params.nbin; i++)
      {
          nsga2Params.nbits[i] = nbits[i];
          if (nsga2Params.nbits[i] < 1)
          {
              printf("\n Wrong number of bits for binary variable entered, hence exiting");
              exit(1);
          }
          nsga2Params.min_binvar[i] = min_binvar[i];
          nsga2Params.max_binvar[i] = max_binvar[i];
          if (nsga2Params.max_binvar[i] <= nsga2Params.min_binvar[i])
          {
              printf("\n Wrong limits entered for the min and max bounds of binary variable entered, hence exiting \n");
              exit(1);
          }
      }
      nsga2Params.pcross_bin = pcross_bin;
      if (nsga2Params.pcross_bin<0.0 || nsga2Params.pcross_bin>1.0)
      {
          printf("\n Probability of crossover entered is : %e",nsga2Params.pcross_bin);
          printf("\n Entered value of probability of crossover of binary variables is out of bounds, hence exiting \n");
          exit (1);
      }
      nsga2Params.pmut_bin = pmut_bin;
      if (nsga2Params.pmut_bin<0.0 || nsga2Params.pmut_bin>1.0)
      {
          printf("\n Probability of mutation entered is : %e",nsga2Params.pmut_bin);
          printf("\n Entered value of probability  of mutation of binary variables is out of bounds, hence exiting \n");
          exit (1);
      }
  }
  if (nsga2Params.nreal==0 && nsga2Params.nbin==0)
  {
      printf("\n Number of real as well as binary variables, both are zero, hence exiting \n");
      exit(1);
  }

  nsga2Params.choice = choice;
  if (nsga2Params.choice!=0 && nsga2Params.choice!=1)
  {
      printf("\n Entered the wrong nsga2Params.choice, hence exiting, nsga2Params.choice entered was %d\n",nsga2Params.choice);
      exit(1);
  }

  if (nsga2Params.choice==1)
  {
      gp = popen(GNUPLOT_COMMAND,"w");
      if (gp==NULL)
      {
          printf("\n Could not open a pipe to gnuplot, check the definition of GNUPLOT_COMMAND in file global.h\n");
          printf("\n Edit the string to suit your system configuration and rerun the program\n");
          exit(1);
      }
      if (nsga2Params.nobj==2)
      {
          nsga2Params.obj1 = obj1;
          if (nsga2Params.obj1<1 || nsga2Params.obj1>nsga2Params.nobj)
          {
              printf("\n Wrong value of X objective entered, value entered was %d\n",nsga2Params.obj1);
              exit(1);
          }
          nsga2Params.obj2 = obj2;
          if (nsga2Params.obj2<1 || nsga2Params.obj2>nsga2Params.nobj)
          {
              printf("\n Wrong value of Y objective entered, value entered was %d\n",nsga2Params.obj2);
              exit(1);
          }
          nsga2Params.obj3 = -1;
      }
      else
      {
          nsga2Params.choice = choice;
          if (nsga2Params.choice!=2 && nsga2Params.choice!=3)
          {
              printf("\n Entered the wrong nsga2Params.choice, hence exiting, nsga2Params.choice entered was %d\n",nsga2Params.choice);
              exit(1);
          }
          if (nsga2Params.choice==2)
          {
              nsga2Params.obj1 = obj1;
              if (nsga2Params.obj1<1 || nsga2Params.obj1>nsga2Params.nobj)
              {
                  printf("\n Wrong value of X objective entered, value entered was %d\n",nsga2Params.obj1);
                  exit(1);
              }

              nsga2Params.obj2 = obj2;
              if (nsga2Params.obj2<1 || nsga2Params.obj2>nsga2Params.nobj)
              {
                  printf("\n Wrong value of Y objective entered, value entered was %d\n",nsga2Params.obj2);
                  exit(1);
              }
              nsga2Params.obj3 = -1;
          }
          else
          {
              nsga2Params.obj1 = obj1;
              if (nsga2Params.obj1<1 || nsga2Params.obj1>nsga2Params.nobj)
              {
                  printf("\n Wrong value of X objective entered, value entered was %d\n",nsga2Params.obj1);
                  exit(1);
              }
              nsga2Params.obj2 = obj2;
              if (nsga2Params.obj2<1 || nsga2Params.obj2>nsga2Params.nobj)
              {
                  printf("\n Wrong value of Y objective entered, value entered was %d\n",nsga2Params.obj2);
                  exit(1);
              }
              nsga2Params.obj3;
              if (nsga2Params.obj3<1 || nsga2Params.obj3>nsga2Params.nobj)
              {
                  printf("\n Wrong value of Z objective entered, value entered was %d\n",nsga2Params.obj3);
                  exit(1);
              }
              nsga2Params.angle1 = angle1;
              if (nsga2Params.angle1<0 || nsga2Params.angle1>180)
              {
                  printf("\n Wrong value for first angle entered, hence exiting \n");
                  exit(1);
              }
              nsga2Params.angle2 = angle2;
              if (nsga2Params.angle2<0 || nsga2Params.angle2>360)
              {
                  printf("\n Wrong value for second angle entered, hence exiting \n");
                  exit(1);
              }
          }
      }
  }
  obj_empty = new double[nsga2Params.nobj];
  void *inp = NULL;
  void *out = NULL;
  InitNSGA2( inp, out);
}


nsga2_optimization::~nsga2_optimization()
{
  // Closing the files and freeing up memories...
  fflush(stdout);
  fflush(all_pop);
  fflush(params);

  fclose(params);
  fclose(all_pop);

  if (nsga2Params.choice!=0)
  {
      pclose(gp);
  }
  if (nsga2Params.nreal!=0)
  {
      free (nsga2Params.min_realvar);
      free (nsga2Params.max_realvar);
  }
  if (nsga2Params.nbin!=0)
  {
      free (nsga2Params.min_binvar);
      free (nsga2Params.max_binvar);
      free (nsga2Params.nbits);
  }
  deallocate_memory_pop (parent_pop, nsga2Params.popsize);
  deallocate_memory_pop (child_pop, nsga2Params.popsize);
  deallocate_memory_pop (mixed_pop, 2*nsga2Params.popsize);
  free (parent_pop);
  free (child_pop);
  free (mixed_pop);
}

void nsga2_optimization::ReadParameters(int argc, char **argv)
{
    int i;

    if (argc<2)
    {
        printf("\n Usage ./nsga2r random_seed \n");
        exit(1);
    }
    nsga2Params.seed = (double)atof(argv[1]);
    if (nsga2Params.seed<=0.0 || nsga2Params.seed>=1.0)
    {
        printf("\n Entered seed value is wrong, seed value must be in (0,1) \n");
        exit(1);
    }
    printf("\n Enter the problem relevant and algorithm relevant parameters ... ");
    printf("\n Enter the population size (a multiple of 4) : ");
    scanf("%d",&nsga2Params.popsize);
    if (nsga2Params.popsize<4 || (nsga2Params.popsize%4)!= 0)
    {
        printf("\n population size read is : %d",nsga2Params.popsize);
        printf("\n Wrong population size entered, hence exiting \n");
        exit (1);
    }
    printf("\n Enter the number of generations : ");
    scanf("%d",&nsga2Params.ngen);
    if (nsga2Params.ngen<1)
    {
        printf("\n number of generations read is : %d",nsga2Params.ngen);
        printf("\n Wrong nuber of generations entered, hence exiting \n");
        exit (1);
    }
    printf("\n Enter the number of objectives : ");
    scanf("%d",&nsga2Params.nobj);
    if (nsga2Params.nobj<1)
    {
        printf("\n number of objectives entered is : %d",nsga2Params.nobj);
        printf("\n Wrong number of objectives entered, hence exiting \n");
        exit (1);
    }
    printf("\n Enter the number of constraints : ");
    scanf("%d",&nsga2Params.ncon);
    if (nsga2Params.ncon<0)
    {
        printf("\n number of constraints entered is : %d",nsga2Params.ncon);
        printf("\n Wrong number of constraints enetered, hence exiting \n");
        exit (1);
    }
    printf("\n Enter the number of real variables : ");
    scanf("%d",&nsga2Params.nreal);
    if (nsga2Params.nreal<0)
    {
        printf("\n number of real variables entered is : %d",nsga2Params.nreal);
        printf("\n Wrong number of variables entered, hence exiting \n");
        exit (1);
    }
    if (nsga2Params.nreal != 0)
    {
        nsga2Params.min_realvar = (double *)malloc(nsga2Params.nreal*sizeof(double));
        nsga2Params.max_realvar = (double *)malloc(nsga2Params.nreal*sizeof(double));
        for (i=0; i<nsga2Params.nreal; i++)
        {
            printf ("\n Enter the lower limit of real variable %d : ",i+1);
            scanf ("%lf",&nsga2Params.min_realvar[i]);
            printf ("\n Enter the upper limit of real variable %d : ",i+1);
            scanf ("%lf",&nsga2Params.max_realvar[i]);
            if (nsga2Params.max_realvar[i] <= nsga2Params.min_realvar[i])
            {
                printf("\n Wrong limits entered for the min and max bounds of real variable, hence exiting \n");
                exit(1);
            }
        }
        printf ("\n Enter the probability of crossover of real variable (0.6-1.0) : ");
        scanf ("%lf",&nsga2Params.pcross_real);
        if (nsga2Params.pcross_real<0.0 || nsga2Params.pcross_real>1.0)
        {
            printf("\n Probability of crossover entered is : %e",nsga2Params.pcross_real);
            printf("\n Entered value of probability of crossover of real variables is out of bounds, hence exiting \n");
            exit (1);
        }
        printf ("\n Enter the probablity of mutation of real variables (1/nsga2Params.nreal) : ");
        scanf ("%lf",&nsga2Params.pmut_real);
        if (nsga2Params.pmut_real<0.0 || nsga2Params.pmut_real>1.0)
        {
            printf("\n Probability of mutation entered is : %e",nsga2Params.pmut_real);
            printf("\n Entered value of probability of mutation of real variables is out of bounds, hence exiting \n");
            exit (1);
        }
        printf ("\n Enter the value of distribution index for crossover (5-20): ");
        scanf ("%lf",&nsga2Params.eta_c);
        if (nsga2Params.eta_c<=0)
        {
            printf("\n The value entered is : %e",nsga2Params.eta_c);
            printf("\n Wrong value of distribution index for crossover entered, hence exiting \n");
            exit (1);
        }
        printf ("\n Enter the value of distribution index for mutation (5-50): ");
        scanf ("%lf",&nsga2Params.eta_m);
        if (nsga2Params.eta_m<=0)
        {
            printf("\n The value entered is : %e",nsga2Params.eta_m);
            printf("\n Wrong value of distribution index for mutation entered, hence exiting \n");
            exit (1);
        }
    }
    printf("\n Enter the number of binary variables : ");
    scanf("%d",&nsga2Params.nbin);
    if (nsga2Params.nbin<0)
    {
        printf ("\n number of binary variables entered is : %d",nsga2Params.nbin);
        printf ("\n Wrong number of binary variables entered, hence exiting \n");
        exit(1);
    }
    if (nsga2Params.nbin != 0)
    {
        nsga2Params.nbits = (int *)malloc(nsga2Params.nbin*sizeof(int));
        nsga2Params.min_binvar = (double *)malloc(nsga2Params.nbin*sizeof(double));
        nsga2Params.max_binvar = (double *)malloc(nsga2Params.nbin*sizeof(double));
        for (i=0; i<nsga2Params.nbin; i++)
        {
            printf ("\n Enter the number of bits for binary variable %d : ",i+1);
            scanf ("%d",&nsga2Params.nbits[i]);
            if (nsga2Params.nbits[i] < 1)
            {
                printf("\n Wrong number of bits for binary variable entered, hence exiting");
                exit(1);
            }
            printf ("\n Enter the lower limit of binary variable %d : ",i+1);
            scanf ("%lf",&nsga2Params.min_binvar[i]);
            printf ("\n Enter the upper limit of binary variable %d : ",i+1);
            scanf ("%lf",&nsga2Params.max_binvar[i]);
            if (nsga2Params.max_binvar[i] <= nsga2Params.min_binvar[i])
            {
                printf("\n Wrong limits entered for the min and max bounds of binary variable entered, hence exiting \n");
                exit(1);
            }
        }
        printf ("\n Enter the probability of crossover of binary variable (0.6-1.0): ");
        scanf ("%lf",&nsga2Params.pcross_bin);
        if (nsga2Params.pcross_bin<0.0 || nsga2Params.pcross_bin>1.0)
        {
            printf("\n Probability of crossover entered is : %e",nsga2Params.pcross_bin);
            printf("\n Entered value of probability of crossover of binary variables is out of bounds, hence exiting \n");
            exit (1);
        }
        printf ("\n Enter the probability of mutation of binary variables (1/nsga2Params.nbits): ");
        scanf ("%lf",&nsga2Params.pmut_bin);
        if (nsga2Params.pmut_bin<0.0 || nsga2Params.pmut_bin>1.0)
        {
            printf("\n Probability of mutation entered is : %e",nsga2Params.pmut_bin);
            printf("\n Entered value of probability  of mutation of binary variables is out of bounds, hence exiting \n");
            exit (1);
        }
    }
    if (nsga2Params.nreal==0 && nsga2Params.nbin==0)
    {
        printf("\n Number of real as well as binary variables, both are zero, hence exiting \n");
        exit(1);
    }
    nsga2Params.choice=0;
    printf("\n Do you want to use gnuplot to display the results realtime (0 for NO) (1 for yes) : ");
    scanf("%d",&nsga2Params.choice);
    if (nsga2Params.choice!=0 && nsga2Params.choice!=1)
    {
        printf("\n Entered the wrong nsga2Params.choice, hence exiting, nsga2Params.choice entered was %d\n",nsga2Params.choice);
        exit(1);
    }
    if (nsga2Params.choice==1)
    {
        gp = popen(GNUPLOT_COMMAND,"w");
        if (gp==NULL)
        {
            printf("\n Could not open a pipe to gnuplot, check the definition of GNUPLOT_COMMAND in file global.h\n");
            printf("\n Edit the string to suit your system configuration and rerun the program\n");
            exit(1);
        }
        if (nsga2Params.nobj==2)
        {
            printf("\n Enter the objective for X axis display : ");
            scanf("%d",&nsga2Params.obj1);
            if (nsga2Params.obj1<1 || nsga2Params.obj1>nsga2Params.nobj)
            {
                printf("\n Wrong value of X objective entered, value entered was %d\n",nsga2Params.obj1);
                exit(1);
            }
            printf("\n Enter the objective for Y axis display : ");
            scanf("%d",&nsga2Params.obj2);
            if (nsga2Params.obj2<1 || nsga2Params.obj2>nsga2Params.nobj)
            {
                printf("\n Wrong value of Y objective entered, value entered was %d\n",nsga2Params.obj2);
                exit(1);
            }
            nsga2Params.obj3 = -1;
        }
        else
        {
            printf("\n #obj > 2, 2D display or a 3D display ?, enter 2 for 2D and 3 for 3D :");
            scanf("%d",&nsga2Params.choice);
            if (nsga2Params.choice!=2 && nsga2Params.choice!=3)
            {
                printf("\n Entered the wrong nsga2Params.choice, hence exiting, nsga2Params.choice entered was %d\n",nsga2Params.choice);
                exit(1);
            }
            if (nsga2Params.choice==2)
            {
                printf("\n Enter the objective for X axis display : ");
                scanf("%d",&nsga2Params.obj1);
                if (nsga2Params.obj1<1 || nsga2Params.obj1>nsga2Params.nobj)
                {
                    printf("\n Wrong value of X objective entered, value entered was %d\n",nsga2Params.obj1);
                    exit(1);
                }
                printf("\n Enter the objective for Y axis display : ");
                scanf("%d",&nsga2Params.obj2);
                if (nsga2Params.obj2<1 || nsga2Params.obj2>nsga2Params.nobj)
                {
                    printf("\n Wrong value of Y objective entered, value entered was %d\n",nsga2Params.obj2);
                    exit(1);
                }
                nsga2Params.obj3 = -1;
            }
            else
            {
                printf("\n Enter the objective for X axis display : ");
                scanf("%d",&nsga2Params.obj1);
                if (nsga2Params.obj1<1 || nsga2Params.obj1>nsga2Params.nobj)
                {
                    printf("\n Wrong value of X objective entered, value entered was %d\n",nsga2Params.obj1);
                    exit(1);
                }
                printf("\n Enter the objective for Y axis display : ");
                scanf("%d",&nsga2Params.obj2);
                if (nsga2Params.obj2<1 || nsga2Params.obj2>nsga2Params.nobj)
                {
                    printf("\n Wrong value of Y objective entered, value entered was %d\n",nsga2Params.obj2);
                    exit(1);
                }
                printf("\n Enter the objective for Z axis display : ");
                scanf("%d",&nsga2Params.obj3);
                if (nsga2Params.obj3<1 || nsga2Params.obj3>nsga2Params.nobj)
                {
                    printf("\n Wrong value of Z objective entered, value entered was %d\n",nsga2Params.obj3);
                    exit(1);
                }
                printf("\n You have chosen 3D display, hence location of eye required \n");
                printf("\n Enter the first angle (an integer in the range 0-180) (if not known, enter 60) :");
                scanf("%d",&nsga2Params.angle1);
                if (nsga2Params.angle1<0 || nsga2Params.angle1>180)
                {
                    printf("\n Wrong value for first angle entered, hence exiting \n");
                    exit(1);
                }
                printf("\n Enter the second angle (an integer in the range 0-360) (if not known, enter 30) :");
                scanf("%d",&nsga2Params.angle2);
                if (nsga2Params.angle2<0 || nsga2Params.angle2>360)
                {
                    printf("\n Wrong value for second angle entered, hence exiting \n");
                    exit(1);
                }
            }
        }
    }
}

void nsga2_optimization::InitNSGA2( void *inp, void *out)
{
    int i;

    printf("\n == InitNSGA2 ==");

    // WRITING OUT PARAMETERS
    FILE* initial_pop;

    int load_flow_type = bin_params[2];
    int hasTapChanger = bin_params[3];

    string initial_pop_string = "./output/NSGA-II/original/initial_pop";
    string params_string = "./output/NSGA-II/original/params";
    string all_pop_string = "./output/NSGA-II/original/all_pop";
  	if( load_flow_type == 1 )
    {
  		initial_pop_string += "_constP";
      params_string += "_constP";
      all_pop_string += "_constP";
    }
  	else if (load_flow_type == 2)
    {
  		initial_pop_string += "_constI";
      params_string += "_constI";
      all_pop_string += "_constI";
    }

  	if (hasTapChanger == 0)
    {
  		initial_pop_string += "_noRTR_";
      params_string += "_noRTR_";
      all_pop_string += "_noRTR_";
    }
  	else
    {
  		initial_pop_string += "_wRTR_";
      params_string += "_wRTR_";
      all_pop_string += "_wRTR_";
    }

  	initial_pop_string += to_string(bin_params[0]);
    params_string += to_string(bin_params[0]);
    all_pop_string += to_string(bin_params[0]);

  	initial_pop_string += "DG_Type";
    params_string += "DG_Type";
    all_pop_string += "DG_Type";

  	initial_pop_string += to_string(bin_params[1]);
    params_string += to_string(bin_params[1]);
    all_pop_string += to_string(bin_params[1]);

  	initial_pop_string += ".out";
    params_string += ".out";
    all_pop_string += ".out";

    initial_pop = fopen(initial_pop_string.c_str(), "w");
    params = fopen(params_string.c_str(),"w");
    all_pop = fopen(all_pop_string.c_str(),"w");

    fprintf(all_pop,"# This file contains the data of all generations\n");
    fprintf(all_pop,"# of objectives = %d, # of constraints = %d, # of real_var = %d, # of bits of bin_var = %d, constr_violation, rank, crowding_distance\n",nsga2Params.nobj,nsga2Params.ncon,nsga2Params.nreal,nsga2Params.bitlength);
    fprintf(initial_pop,"# This file contains the data of initial population\n");
    fprintf(initial_pop,"# of objectives = %d, # of constraints = %d, # of real_var = %d, # of bits of bin_var = %d, constr_violation, rank, crowding_distance\n",nsga2Params.nobj,nsga2Params.ncon,nsga2Params.nreal,nsga2Params.bitlength);


    fprintf(params,"# This file contains information about inputs as read by the program\n");
    fprintf(params,"\n Population size = %d",nsga2Params.popsize);
    fprintf(params,"\n Number of generations = %d",nsga2Params.ngen);
    fprintf(params,"\n Number of objective functions = %d",nsga2Params.nobj);
    fprintf(params,"\n Number of constraints = %d",nsga2Params.ncon);
    fprintf(params,"\n Number of real variables = %d",nsga2Params.nreal);
    if (nsga2Params.nreal!=0)
    {
        for (i=0; i<nsga2Params.nreal; i++)
        {
            fprintf(params,"\n Lower limit of real variable %d = %e",i+1,nsga2Params.min_realvar[i]);
            fprintf(params,"\n Upper limit of real variable %d = %e",i+1,nsga2Params.max_realvar[i]);
        }
        fprintf(params,"\n Probability of crossover of real variable = %e",nsga2Params.pcross_real);
        fprintf(params,"\n Probability of mutation of real variable = %e",nsga2Params.pmut_real);
        fprintf(params,"\n Distribution index for crossover = %e",nsga2Params.eta_c);
        fprintf(params,"\n Distribution index for mutation = %e",nsga2Params.eta_m);
    }
    fprintf(params,"\n Number of binary variables = %d",nsga2Params.nbin);
    if (nsga2Params.nbin!=0)
    {
        for (i=0; i<nsga2Params.nbin; i++)
        {
            fprintf(params,"\n Number of bits for binary variable %d = %d",i+1,nsga2Params.nbits[i]);
            fprintf(params,"\n Lower limit of binary variable %d = %e",i+1,nsga2Params.min_binvar[i]);
            fprintf(params,"\n Upper limit of binary variable %d = %e",i+1,nsga2Params.max_binvar[i]);
        }
        fprintf(params,"\n Probability of crossover of binary variable = %e",nsga2Params.pcross_bin);
        fprintf(params,"\n Probability of mutation of binary variable = %e",nsga2Params.pmut_bin);
    }
    fprintf(params,"\n Seed for random number generator = %e",nsga2Params.seed);
    nsga2Params.bitlength = 0;
    if (nsga2Params.nbin!=0)
    {
        for (i=0; i<nsga2Params.nbin; i++)
        {
            nsga2Params.bitlength += nsga2Params.nbits[i];
        }
    }
    nsga2Params.nbinmut = 0;
    nsga2Params.nrealmut = 0;
    nsga2Params.nbincross = 0;
    nsga2Params.nrealcross = 0;

    // Initializing the populations
    parent_pop = (population *)malloc(sizeof(population));
    child_pop = (population *)malloc(sizeof(population));
    mixed_pop = (population *)malloc(sizeof(population));
    allocate_memory_pop (parent_pop, nsga2Params.popsize);
    allocate_memory_pop (child_pop, nsga2Params.popsize);
    allocate_memory_pop (mixed_pop, 2*nsga2Params.popsize);

    // Preparing first Population
    randomize(nsga2Params.seed);
    initialize_pop (parent_pop);
    printf("\n Initialization done, now performing first generation");
    decode_pop(parent_pop);
    evaluate_pop (parent_pop, inp, out);
    assign_rank_and_crowding_distance (parent_pop);

    report_pop (initial_pop);
    fprintf(initial_pop,"# gen = 1\n");
    report_pop(all_pop);
    printf("\n -- Generation 1 --");

    fflush(stdout);
    if (nsga2Params.choice!=0)    onthefly_display (parent_pop,gp,1);

    fflush(params);
    fflush(all_pop);

}

void nsga2_optimization::set_real_params(double* real_params, int real_params_num)
{
  this->real_params = new double[real_params_num];
  for(int i=0; i<real_params_num; i++)
    this->real_params[i] = real_params[i];
}

void nsga2_optimization::set_bin_params(int* bin_params, int bin_params_num)
{
  this->bin_params = new int[bin_params_num];
  for(int i=0; i<bin_params_num; i++)
    this->bin_params[i] = bin_params[i];
}

void nsga2_optimization::set_params(double* real_params, int real_params_num, int* bin_params, int bin_params_num)
{
  set_real_params(real_params, real_params_num);
  set_bin_params(bin_params, bin_params_num);
}


int nsga2_optimization::NSGA2( void *inp, void *out)
{
    int i;

    for (i=2; i<=nsga2Params.ngen; i++)
    {
        selection (parent_pop, child_pop);
        mutation_pop ( child_pop);
        decode_pop( child_pop);
        evaluate_pop( child_pop, inp, out);
        merge ( parent_pop, child_pop, mixed_pop);
        fill_nondominated_sort ( mixed_pop, parent_pop);
        /* Comment following four lines if information for all
         generations is not desired, it will speed up the execution */
        fprintf(all_pop,"# gen = %d\n",i);
        report_pop(all_pop);
        fflush(all_pop);
        if (nsga2Params.choice!=0)    onthefly_display (parent_pop,gp,i);
        printf("\n -- Generation %d --", i);
    }
    printf("\n Generations finished, now reporting solutions");

    if (nsga2Params.nreal!=0)
    {
        fprintf(params,"\n Number of crossover of real variable = %d",nsga2Params.nrealcross);
        fprintf(params,"\n Number of mutation of real variable = %d",nsga2Params.nrealmut);
    }
    if (nsga2Params.nbin!=0)
    {
        fprintf(params,"\n Number of crossover of binary variable = %d",nsga2Params.nbincross);
        fprintf(params,"\n Number of mutation of binary variable = %d",nsga2Params.nbinmut);
    }


    printf("\n Algorithm was succesfully executed \n");
    return (0);
}

void nsga2_optimization::print_nsga2Params()
{
    int i;

    printf("NSGA2 Parameters:\n");
    printf("\n seed number is : %lf", nsga2Params.seed);
    printf("\n population size read is : %d",nsga2Params.popsize);
    printf("\n number of generations read is : %d",nsga2Params.ngen);
    printf("\n number of objectives entered is : %d",nsga2Params.nobj);
    printf("\n number of constraints entered is : %d",nsga2Params.ncon);
    printf("\n number of real variables entered is : %d",nsga2Params.nreal);
    printf("\n variables bounds: ");
    for (i=0; i<nsga2Params.nreal; i++){
        printf("[%lf", nsga2Params.min_realvar[i]);
        printf(" %lf], ", nsga2Params.max_realvar[i]);
    }
    printf("\n Probability of crossover entered is : %e",nsga2Params.pcross_real);
    printf("\n Probability of mutation entered is : %e",nsga2Params.pmut_real);
    printf("\n The eta_c entered is : %e",nsga2Params.eta_c);
    printf("\n The eta_m entered is : %e",nsga2Params.eta_m);
    if (nsga2Params.nbin != 0){
        printf ("\n number of binary variables entered is : %d",nsga2Params.nbin);
        printf("\n Probability of crossover entered is : %e",nsga2Params.pcross_bin);
        printf("\n Probability of mutation entered is : %e",nsga2Params.pmut_bin);
    }
}
