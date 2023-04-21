#ifndef PARTICLE_SWARM_H
#define PARTICLE_SWARM_H

#include "../NetworkModel/tempel.h"
#include "../NetworkModel/network.h"
#include "repository.h"
#include <iostream>
#include <vector>

using namespace std;

class network;

class particle_swarm
{
      network* netw;
      const double alpha = 0.1;
      const double beta = 2;
      const double gamma = 2;
      const double mu = 0.1;
      int NOP;
      int MAXITER;


      int cap_num;
      int disp_num;
      int solar_num;
      int wind_num;
      int fitness_num;
      int dispGType;

      int switchnum;
      int scenarionum;
      double W;
      double Wdamp;
      double C1;
      double C2;

      PSO_particle **particles;
      gbest_repository REP;

    public:
      bool doMutation(int iter);
      particle_swarm();
      particle_swarm(network* netw, int NOP = 500, int MAXITER = 200, double C1 = 2.0, double C2 = 2.0, double W = 0.5, double Wdamp = 0.99, int NREP = 100, int Ngrid = 7, bool hasTapChanger = false);
      particle_swarm(particle_swarm&);

      void initialize_scenario1(int switchnum, int NOP = 500, int MAXITER = 200, double C1 = 2.0, double C2 = 2.0, double W=0.5, double Wdamp = 0.99, bool hasTapChanger = false);  // Study with mostly renewable generators and capacitors
      void initialize_scenario2(int numOG, int switchnum, int dispType = 3, int NOP = 200, int MAXITER = 200, double C1 = 1.0, double C2 = 2.0, double W=0.5, double Wdamp = 0.99, bool hasTapChanger = false);  // Study with the dispachable generators
      void clear_particles();
      void particle_swarm_optimalization();
      PSO_particle fuzzy_choice();
      PSO_particle circle_choice() {return REP.circle_choice();}
      double calculate2Dhypervolume(int funcs) {return REP.calculate2Dhypervolume(funcs); }
      double calculateSpacing(int switchnum) {return REP.calculateSpacing(switchnum);}

      void log_particle_positions(const char*);
      void log_particle_velocities(const char*);
      void log_particle_bests(const char*);
      void log_particle_fitness_vals(const char*);

      void log_table_form(const char* fname) { REP.log_table_form(fname); }
      void log_master_form(ofstream& os) {REP.log_master_form(os);}

      void log_rep_size(ofstream& os) {os<<REP.rep_size()<<endl;}
      void log_MinMax(ofstream& os, int i) { REP.log_MinMax(os, i);}
      void log_2DHypervolume(ofstream& os, int switchnum) {os<<REP.calculate2Dhypervolume(switchnum)<<endl;}
      void log_Spacing(ofstream& os, int switchnum) {os<<REP.calculateSpacing(switchnum)<<endl;}

      void createParetoOctaveGraph(const char* fname) {REP.createParetoOctaveGraph(fname, disp_num, dispGType); }
};

#endif
