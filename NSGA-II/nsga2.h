/* This file contains the variable and function declarations */


# ifndef NSGA_INCLUDED
# define NSGA_INCLUDED

# include <stdio.h>
# include <stdlib.h>
# include <math.h>
# include <unistd.h>
# include "../NetworkModel/network.h"

# define INF 1.0e14
# define EPS 1.0e-14
# define eul  2.71828182845905
# define pi 3.14159265358979
# define GNUPLOT_COMMAND "gnuplot -persist"


typedef struct
{
    int rank;
    double constr_violation;
    double *xreal;
    int **gene;
    double *xbin;
    double *obj;
    double *constr;
    double crowd_dist;
  }
  individual;

typedef struct
{
    individual *ind;
}
population;

typedef struct lists
{
    int index;
    struct lists *parent;
    struct lists *child;
} list;


class nsga2_optimization
{
  private:
    typedef struct NSGA2Type
    {
        double seed;
        int nreal;
        int nbin;
        int nobj;
        int ncon;
        int popsize;
        double pcross_real;
        double pcross_bin;
        double pmut_real;
        double pmut_bin;
        double eta_c;
        double eta_m;
        int ngen;
        int nbinmut;
        int nrealmut;
        int nbincross;
        int nrealcross;
        int *nbits;
        double *min_realvar;
        double *max_realvar;
        double *min_binvar;
        double *max_binvar;
        int bitlength;
        int choice;
        int obj1;
        int obj2;
        int obj3;
        int angle1;
        int angle2;
    } NSGA2Type;

    double* obj_empty;

    NSGA2Type nsga2Params;

    FILE *gp;
    FILE *params;
    FILE *all_pop;

    population *parent_pop;
    population *child_pop;
    population *mixed_pop;

    double* real_params;
    int* bin_params;

    int feasnum;
    double* f_min, *f_max;
    /**
     *  allocate.c
     */
    void allocate_memory_pop ( population *pop, int size);
    void allocate_memory_ind ( individual *ind);
    void deallocate_memory_pop ( population *pop, int size);
    void deallocate_memory_ind ( individual *ind);

    // Methods
    double maximum (double a, double b);
    double minimum (double a, double b);
    void swap_int (int*, int*);

    void crossover ( individual *parent1, individual *parent2, individual *child1, individual *child2);
    void realcross ( individual *parent1, individual *parent2, individual *child1, individual *child2);
    void bincross ( individual *parent1, individual *parent2, individual *child1, individual *child2);

    void assign_crowding_distance_list (  population *pop, list *lst, int front_size);
    void assign_crowding_distance_indices (  population *pop, int c1, int c2);
    void assign_crowding_distance ( population *pop, int *dist, int **obj_array, int front_size);

    void decode_pop ( population *pop);
    void decode_ind ( individual *ind);

    void onthefly_display ( population *pop, FILE *gp, int ii);

    int check_dominance ( individual *a, individual *b);

    void evaluate_pop ( population *pop, void *, void *);
    void evaluate_ind ( individual *ind, void *, void *);


    void fill_nondominated_sort ( population *mixed_pop, population *new_pop);
    void crowding_fill ( population *mixed_pop, population *new_pop, int count, int front_size, list *cur);

    void initialize_pop (population *pop);
    void initialize_ind (individual *ind);

    void insert (list *node, int x);
    list* del (list *node);

    void merge( population *pop1, population *pop2, population *pop3);
    void copy_ind (individual *ind1, individual *ind2);

    void mutation_pop (population *pop);
    void mutation_ind (individual *ind);
    void bin_mutate_ind ( individual *ind);
    void real_mutate_ind ( individual *ind);

    void test_problem (double *xreal, double *xbin, int **gene, double *obj, double *constr, int *bin_params = NULL, double *real_params = NULL);

    void assign_rank_and_crowding_distance ( population *new_pop);

    double calculate2Dhypervolume( population *pop);
    double calculateSpacing( population *pop);

    int fuzzy_choice(population *pop);
    int fuzzy_choice_within_range (population *pop);
    int circle_choice(population *pop);

    void quicksort_front_obj(population *pop, int objcount, int obj_array[], int obj_array_size);
    void q_sort_front_obj(population *pop, int objcount, int obj_array[], int left, int right);
    void quicksort_dist(population *pop, int *dist, int front_size);
    void q_sort_dist(population *pop, int *dist, int left, int right);

    void selection (population *old_pop, population *new_pop);
    individual* tournament (individual *ind1, individual *ind2);

    void ReadParameters(int argc, char **argv);
    void InitNSGA2( void *, void *);

  /**
   * nsga2.c
   */
  public:
    nsga2_optimization() {}
    nsga2_optimization(int argc, char **argv, double* real_params, int real_params_num, int* bin_params, int bin_params_num);

    nsga2_optimization(double seed, int popsize, int ngen, int nobj, int ncon, int nreal,
      double* min_realval, double* max_realval, double pcross_real, double pmut_real,
      double dist_ind_cross_real, double dist_ind_mut_real, int nbin, int* nbits,
      double* min_binvar, double* max_binvar, double pcross_bin, double pmut_bin,
      double* real_params, int real_params_num, int* bin_params, int bin_params_num,
      int choice = 0, int obj1 = 1, int obj2 = 2, int obj3 = 3, int angle1 = 90, int angle2 = 90);

    ~nsga2_optimization();

    void set_real_params(double* real_params, int real_params_num) ;
    void set_bin_params(int* bin_params, int bin_params_num) ;
    void set_params(double* real_params, int real_params_num, int* bin_params, int bin_params_num) ;

    double getFmin(int i) {return f_min[i];}
    double getFmax(int i) {return f_max[i];}

    int NSGA2( void *, void *);

    void report_pop (FILE *fpt);
    void report_feasible (FILE *fpt);
    void report_ind (FILE *fpt);

    void print_nsga2Params();

    int calculateFeasnum(population* pop);
    void calculateMinMax(population* pop);

    void createParetoOctaveGraph(FILE *fpt);
    void log_table_form(FILE *fpt);
    void log_master_form(ofstream& os);

    void log_pareto_size(ofstream& os) {os<<calculateFeasnum(parent_pop)<<endl;}
    void log_MinMax(ofstream& os, int i);
    void log_2Dhypervolume(ofstream& os) {os<<calculate2Dhypervolume(parent_pop)<<endl;}
    void log_Spacing(ofstream& os) {os<<calculateSpacing(parent_pop)<<endl;}
};
# endif
