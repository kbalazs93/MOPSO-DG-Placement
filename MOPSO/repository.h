#ifndef GBEST_REPOSITORY_HEADER
#define GBEST_REPOSITORY_HEADER

#include "particle.h"
#include <vector>

class PSO_particle;

class gbest_repository
{
    int NREP;
    int Ngrid;
    int fitness_vals;
    vector<PSO_particle> REP_particles;
    vector<int> hypercube_coordinates;
    double* grid_mins;
    double* grid_maxes;
    double* alpha;
    double* f_max;
    double* f_min;
    double alphamax;
    int alphamax_ind;
  public:
    gbest_repository();
    gbest_repository(int NREP, int Ngrid);
//    gbest_repository(gbest_repository&);
    void setFitnessVals(const int);

    int repositorySize() {return REP_particles.size();}
    int roulette_wheel(int mode=0);
    void archive_controler(PSO_particle&);
    void adaptive_grid(bool reconstruct = false);
    PSO_particle select_lead();
    PSO_particle fuzzy_choice();
    PSO_particle circle_choice();
    int fuzzy_choice_within_range();

    int rep_size() {return REP_particles.size();}
    void calculateMinMax();

    double calculate2Dhypervolume(int funcs);
    double calculateSpacing(int switchnum);

    void log_repository_particles(const char* fname);
    void log_table_form(const char* fname);
    void log_master_form(ofstream& os);
    void log_MinMax(ofstream& os, int i);

    void createParetoOctaveGraph(const char* fname, int DGnum, int DGType);

};


#endif
