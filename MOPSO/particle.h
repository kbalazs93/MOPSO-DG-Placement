#ifndef PSO_PARTICLE_H
#define PSO_PARTICLE_H

#include "../NetworkModel/tempel.h"
#include "../NetworkModel/network.h"
#include <vector>

using namespace std;

class network;

class PSO_particle
{
    double W;
    double Wdamp;
    double C1;
    double C2;

    int minNode;
    int maxNode;
    int dispMinP;
    int dispMaxP;
    int dispMinQ;
    int dispMaxQ;
    double solarP;
    double solarQ;
    double windP;
    double windQ;
    double capQ;
    int minCapUnitNum;
    int maxCapUnitNum;
    int minWindUnitNum;
    int maxWindUnitNum;
    int minSolarUnitNum;
    int maxSolarUnitNum;
    int minDispUnitNum;
    int maxDispUnitNum;
    int minTapChanger;
    int maxTapChanger;

    int disp_type;

    bool relaxed_constraints;
    double Vmin;
    double Vmax;
    double Vmin_lim;
    double Vmax_lim;


    network* netw;

    vector<temp_capacitor> capacitors;
    vector<temp_generator> disp_generators;
    vector<temp_generator> wind_generators;
    vector<temp_generator> solar_generators;
    bool hasTapChanger;
    int tapChanger;

    vector<temp_capacitor> cap_vel;
    vector<temp_generator> disp_gen_vel;
    vector<temp_generator> wind_gen_vel;
    vector<temp_generator> solar_gen_vel;
    int tapChanger_vel;

    vector<int> cap_guide;
    vector<int> disp_guide;
    vector<int> wind_guide;
    vector<int> solar_guide;

    vector<temp_capacitor> cap_PBEST;
    vector<temp_generator> disp_gen_PBEST;
    vector<temp_generator> wind_gen_PBEST;
    vector<temp_generator> solar_gen_PBEST;
    int tapChanger_PBEST;

    vector<double> base_fitnes_vals;
    vector<double> fitness_values;
    vector<double> PBEST_fitness_vals;
    bool hasValidFitness;
    bool hasValidPBEST;

    int scenarionum;

    int switchnum;
    bool considerActiveLosses;
    bool considerVSI;
    bool considerVD;

    int round(double d);
    void swap_int(int*, int*);


  public:

    PSO_particle();
    PSO_particle(network* netw, int cap_num, int disp_num, int wind_num, int solar_num, int disp_Type, int switchnum, int scenarionum, double C1=2, double C2=2, double W=0.5, double Wdamp = 0.99, bool hasTapChanger = false);
    PSO_particle(const PSO_particle&);

    PSO_particle& operator=(const PSO_particle&);
    void setParticle(int capNum, int* capNodes, int windNum, int* windNodes, int solarNum, int* solarNodes, int dispNum, int* dispNodes, double* dispPvalues, double* dispQvalues, int tapChangerPos = 0);

    bool getFitnessValidation() {return hasValidFitness; }
    void getFitnessValues( vector<double>* );
    void getBaseFitnessValues (vector<double>* );
    double getFitnessValue(int i) {return fitness_values[i]; }
    int getFitnessNum()const {return fitness_values.size();}

    string getPowerFlowType() ;

    int getSwitchNum()const {return switchnum;}
    bool activeLossesConsidered()const {return considerActiveLosses;}
    bool voltageStabilityIndexConsidered()const {return considerVSI;}
    bool voltageDeviationConsidered()const {return considerVD;}

    void inertiaDamping() {W *= Wdamp;}

    void reset_particle(int dispGType);
    void evaluate_fitness_values();
    void evaluate_PBEST();
    void move_particle_path_based(PSO_particle&);
    void move_particle_layer_index_based(PSO_particle&);
    void clearParticle();

    void createOctaveGraph(string filename, string season);

    void log_particle_positions(const char*);
    void log_particle_velocities(const char*);
    void log_particle_bests(const char*);
    void log_particle_fitness_vals(const char*);

    void log_table_form(ofstream& os);
    void log_master_form(ofstream& os);

    friend bool operator< (const PSO_particle&, const PSO_particle&);
};

#endif
