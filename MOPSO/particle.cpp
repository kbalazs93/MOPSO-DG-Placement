#include "particle.h"
#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <cmath>

using namespace std;

// private:

int PSO_particle::round(double d)
{
  int a = (int)d;
  if(d-a < 0.5)
    return a;
  else
    return a+1;
}

void PSO_particle::swap_int(int* i, int* j)
{
  int helper = *i;
  *i = *j;
  *j = helper;
}

// public:
PSO_particle::PSO_particle()
{
  relaxed_constraints = false;

  minNode = 1;
  maxNode = 32;
  dispMinP = 0;
  dispMaxP = 2e6;
  dispMinQ = 0;
  dispMaxQ = 2e6;
  solarP = .25e6;
  solarQ = 0;
  windP = .25e6;
  windQ = 0;
  capQ = .125e6;
  minTapChanger = 0;
  maxTapChanger = 0;

  disp_type = 1;

  int cap_num = 8;
  int disp_num = 1;
  int solar_num = 4;
  int wind_num = 4;
  int fitness_num = 3;
  int disp_Type = 3;

  scenarionum = 1;

  switchnum = 7;
  considerActiveLosses = true;
  considerVSI = true;
  considerVD = true;

  temp_capacitor new_capacitor;
  temp_generator new_generator;
  tapChanger = 0;

  for(int i=0; i<cap_num; i++)
  {
    capacitors.push_back(new_capacitor);
    cap_PBEST.push_back(new_capacitor);
    cap_vel.push_back(new_capacitor);
    cap_guide.push_back(0);
  }
  for(int i=0; i<disp_num; i++)
  {
    disp_generators.push_back(new_generator);
    disp_gen_PBEST.push_back(new_generator);
    disp_gen_vel.push_back(new_generator);
    disp_guide.push_back(0);
  }
  for(int i=0; i<wind_num; i++)
  {
    wind_generators.push_back(new_generator);
    wind_gen_PBEST.push_back(new_generator);
    wind_gen_vel.push_back(new_generator);
    wind_guide.push_back(0);
  }
  for(int i=0; i<solar_num; i++)
  {
    solar_generators.push_back(new_generator);
    solar_gen_PBEST.push_back(new_generator);
    solar_gen_vel.push_back(new_generator);
    solar_guide.push_back(0);
  }
  for(int i=0; i<fitness_num; i++)
  {
    fitness_values.push_back(0);
    PBEST_fitness_vals.push_back(0);
  }
  hasValidPBEST = false;
  C1 = 1;
  C2 = 1;
  W = 0.5;
  Wdamp = 0.99;
}

PSO_particle::PSO_particle(network* netw, int cap_num, int disp_num, int wind_num, int solar_num, int disp_type, int switchnum, int scenarionum, double C1, double C2, double W, double Wdamp, bool hasTapChanger)
{
    this->netw = netw;
    this->scenarionum = scenarionum;

    this->disp_type = disp_type;
    this->hasTapChanger = hasTapChanger;

    relaxed_constraints = false;

    minNode = 1;
    maxNode = netw->getNumONodes()-1;
    dispMinP = 0;
    dispMaxP = 2e6;
    dispMinQ = 0;
    dispMaxQ = 2e6;
    solarP = .25e6;
    solarQ = 0;
    windP = .25e6;
    windQ = 0;
    capQ = .125e6;

    minCapUnitNum = 1;
    maxCapUnitNum = 8;
    minSolarUnitNum = 1;
    maxSolarUnitNum = 4;
    minWindUnitNum = 1;
    maxWindUnitNum = 4;
    minDispUnitNum = 1;
    maxDispUnitNum = 1;

    if (hasTapChanger)
    {
      minTapChanger = -12;
      maxTapChanger = 12;
    }
    else
    {
      minTapChanger = 0;
      maxTapChanger = 0;
    }

    this->C1 = C1;
    this->C2 = C2;
    this->W = W;
    this->Wdamp = Wdamp;


    temp_capacitor new_capacitor;
    temp_generator new_generator;

    this->switchnum = switchnum;
    int fitness_num = 3;


    considerVD = false;
    considerVSI = false;
    considerActiveLosses = false;

    int remainder = switchnum % 2;
    switchnum /= 2;
    if (remainder == 1)
      considerVD = true;
    remainder = switchnum % 2;
    switchnum /= 2;
    if (remainder == 1)
      considerVSI = true;
    remainder = switchnum % 2;
    if (remainder == 1)
      considerActiveLosses = true;

    for(int i=0; i<cap_num; i++)
    {
      capacitors.push_back(new_capacitor);
      cap_PBEST.push_back(new_capacitor);
      cap_vel.push_back(new_capacitor);
      cap_guide.push_back(0);
    }
    for(int i=0; i<disp_num; i++)
    {
      disp_generators.push_back(new_generator);
      disp_gen_PBEST.push_back(new_generator);
      disp_gen_vel.push_back(new_generator);
      disp_guide.push_back(0);
    }
    for(int i=0; i<wind_num; i++)
    {
      wind_generators.push_back(new_generator);
      wind_gen_PBEST.push_back(new_generator);
      wind_gen_vel.push_back(new_generator);
      wind_guide.push_back(0);
    }
    for(int i=0; i<solar_num; i++)
    {
      solar_generators.push_back(new_generator);
      solar_gen_PBEST.push_back(new_generator);
      solar_gen_vel.push_back(new_generator);
      solar_guide.push_back(0);
    }
    for(int i=0; i<fitness_num; i++)
    {
      fitness_values.push_back(0);
      PBEST_fitness_vals.push_back(0);
    }

    double activeLosses = 0;
    double voltageStabilityIndex = 0;
    double voltageDeviation = 0;

    if (scenarionum == 1)
    {
      for (int j=0; j<4; j++)
      {
        string season;
        switch (j) {
          case 0: season = "summer"; break;
          case 1: season = "autumn"; break;
          case 2: season = "winter"; break;
          case 3: season = "spring"; break;
        }
        for(int i=1; i<=24; i++)
        {
          netw->solvePowerFlow(season, i);
          activeLosses += netw->getActiveLosses();
          voltageStabilityIndex += netw->getVoltageStabilityIndexSum();
          voltageDeviation += netw->getVoltageDeviationSum();
        }
      }
    }
    else if (scenarionum == 2)
    {
      netw->solvePowerFlow("summer", 12);
      activeLosses = netw->getActiveLosses();
      voltageStabilityIndex = netw->getVoltageStabilityIndexSum();
      voltageDeviation = netw->getVoltageDeviationSum();
      if( !netw->isWithinVoltageLimits() )
        hasValidFitness = false;
    }
    else
    {
      cerr<<"Error: Wrong scenario number: "<<scenarionum<<endl;
      exit(EXIT_SUCCESS);
    }
    base_fitnes_vals.push_back(activeLosses);
    base_fitnes_vals.push_back(1/voltageStabilityIndex);
    base_fitnes_vals.push_back(voltageDeviation);


    reset_particle(disp_type);
    hasValidPBEST = false;
}

PSO_particle::PSO_particle(const PSO_particle& p2)
{
  this->netw = p2.netw;

  this->Vmin = p2.Vmin;
  this->Vmax = p2.Vmax;
  this->Vmin_lim = p2.Vmin_lim;
  this->Vmax_lim = p2.Vmax_lim;

  this->disp_type = p2.disp_type;
  this->hasTapChanger = p2.hasTapChanger;


  this->relaxed_constraints = p2.relaxed_constraints;

  minNode = p2.minNode;
  maxNode = p2.maxNode;
  dispMinP = p2.dispMinP;
  dispMaxP = p2.dispMaxP;
  dispMinQ = p2.dispMinQ;
  dispMaxQ = p2.dispMaxQ;
  solarP = p2.solarP;
  solarQ = p2.solarQ;
  windP = p2.windP;
  windQ = p2.windQ;
  capQ = p2.capQ;

  minCapUnitNum = p2.minCapUnitNum;
  maxCapUnitNum = p2.maxCapUnitNum;
  minSolarUnitNum = p2.minSolarUnitNum;
  maxSolarUnitNum = p2.maxSolarUnitNum;
  minWindUnitNum = p2.minWindUnitNum;
  maxWindUnitNum = p2.maxWindUnitNum;
  minDispUnitNum = p2.minDispUnitNum;
  maxDispUnitNum = p2.maxWindUnitNum;

  minTapChanger = p2.minTapChanger;
  maxTapChanger = p2.maxTapChanger;

  tapChanger = p2.tapChanger;


  C1 = p2.C1;
  C2 = p2.C2;
  W = p2.W;
  Wdamp = p2.Wdamp;

  capacitors = p2.capacitors;
  disp_generators = p2.disp_generators;
  wind_generators = p2.wind_generators;
  solar_generators = p2.solar_generators;

  cap_vel = p2.cap_vel;
  disp_gen_vel = p2.disp_gen_vel;
  wind_gen_vel = p2.wind_gen_vel;
  solar_gen_vel = p2.solar_gen_vel;

  cap_PBEST = p2.cap_PBEST;
  disp_gen_PBEST = p2.disp_gen_PBEST;
  wind_gen_PBEST = p2.wind_gen_PBEST;
  solar_gen_PBEST = p2.solar_gen_PBEST;

  cap_guide = p2.cap_guide;
  disp_guide = p2.disp_guide;
  solar_guide = p2.solar_guide;
  wind_guide = p2.wind_guide;

  base_fitnes_vals = p2.base_fitnes_vals;
  fitness_values = p2.fitness_values;
  PBEST_fitness_vals = p2.PBEST_fitness_vals;

  hasValidPBEST = p2.hasValidPBEST;
  hasValidFitness = p2.hasValidFitness;

  considerActiveLosses = p2.considerActiveLosses;
  considerVSI = p2.considerVSI;
  considerVD = p2.considerVD;
  switchnum = p2.switchnum;

  scenarionum = p2.scenarionum;
}

PSO_particle& PSO_particle::operator=(const PSO_particle& p2)
{
  this->netw = p2.netw;

  this->Vmin = p2.Vmin;
  this->Vmax = p2.Vmax;
  this->Vmin_lim = p2.Vmin_lim;
  this->Vmax_lim = p2.Vmax_lim;

  this->disp_type = p2.disp_type;
  this->hasTapChanger = p2.hasTapChanger;

  this->relaxed_constraints = p2.relaxed_constraints;

  minNode = p2.minNode;
  maxNode = p2.maxNode;
  dispMinP = p2.dispMinP;
  dispMaxP = p2.dispMaxP;
  dispMinQ = p2.dispMinQ;
  dispMaxQ = p2.dispMaxQ;
  solarP = p2.solarP;
  solarQ = p2.solarQ;
  windP = p2.windP;
  windQ = p2.windQ;
  capQ = p2.capQ;


  minCapUnitNum = p2.minCapUnitNum;
  maxCapUnitNum = p2.maxCapUnitNum;
  minSolarUnitNum = p2.minSolarUnitNum;
  maxSolarUnitNum = p2.maxSolarUnitNum;
  minWindUnitNum = p2.minWindUnitNum;
  maxWindUnitNum = p2.maxWindUnitNum;
  minDispUnitNum = p2.minDispUnitNum;
  maxDispUnitNum = p2.maxWindUnitNum;

  minTapChanger = p2.minTapChanger;
  maxTapChanger = p2.maxTapChanger;

  tapChanger = p2.tapChanger;

  capacitors = p2.capacitors;
  disp_generators = p2.disp_generators;
  wind_generators = p2.wind_generators;
  solar_generators = p2.solar_generators;

  C1 = p2.C1;
  C2 = p2.C2;
  W = p2.W;
  Wdamp = p2.Wdamp;

  cap_vel = p2.cap_vel;
  disp_gen_vel = p2.disp_gen_vel;
  wind_gen_vel = p2.wind_gen_vel;
  solar_gen_vel = p2.solar_gen_vel;

  cap_PBEST = p2.cap_PBEST;
  disp_gen_PBEST = p2.disp_gen_PBEST;
  wind_gen_PBEST = p2.wind_gen_PBEST;
  solar_gen_PBEST = p2.solar_gen_PBEST;

  cap_guide = p2.cap_guide;
  disp_guide = p2.disp_guide;
  solar_guide = p2.solar_guide;
  wind_guide = p2.wind_guide;

  base_fitnes_vals = p2.base_fitnes_vals;
  fitness_values = p2.fitness_values;

  hasValidPBEST = p2.hasValidPBEST;
  hasValidFitness = p2.hasValidFitness;

  considerActiveLosses = p2.considerActiveLosses;
  considerVSI = p2.considerVSI;
  considerVD = p2.considerVD;
  switchnum = p2.switchnum;

  scenarionum = p2.scenarionum;

  return *this;
}

void PSO_particle::getFitnessValues( vector<double>* v )
{
  for(int i=0; i<fitness_values.size(); ++i)
    v->push_back(fitness_values[i]);
}

void PSO_particle::getBaseFitnessValues( vector<double>* v )
{
  for(int i=0; i<base_fitnes_vals.size(); ++i)
    v->push_back(base_fitnes_vals[i]);
}

string PSO_particle::getPowerFlowType()
{
  return netw->getPowerFlowType();
}


void PSO_particle::setParticle(int capNum, int* capNodes, int windNum, int* windNodes, int solarNum, int* solarNodes, int dispNum, int* dispNodes, double* dispPvalues, double* dispQvalues, int tapChangerPos)
{

  clearParticle();
  for(int i=0; i<capNum; i++)
  {
    temp_capacitor new_capacitor(capNodes[i], 1, .125e6);
    capacitors.push_back(new_capacitor);
  }
  for(int i=0; i<windNum; i++)
  {
    temp_generator new_generator(windNodes[i], 1, .25e6, 0);
    wind_generators.push_back(new_generator);
  }
  for(int i=0; i<solarNum; i++)
  {
    temp_generator new_generator(solarNodes[i], 1, .25e6, 0);
    new_generator.setNode(solarNodes[i]);
    solar_generators.push_back(new_generator);
  }
  for(int i=0; i<dispNum; i++)
  {
    temp_generator new_generator;
    new_generator.setNode(dispNodes[i]);
    new_generator.setPvalue(dispPvalues[i]);
    new_generator.setQvalue(dispQvalues[i]);
    disp_generators.push_back(new_generator);
  }

  tapChanger = tapChangerPos;
}

void PSO_particle::clearParticle()
{
  capacitors.clear();
  wind_generators.clear();
  solar_generators.clear();
  disp_generators.clear();

  cap_vel.clear();
  wind_gen_vel.clear();
  solar_gen_vel.clear();
  disp_gen_vel.clear();

  cap_PBEST.clear();
  wind_gen_PBEST.clear();
  solar_gen_PBEST.clear();
  disp_gen_PBEST.clear();
}


void PSO_particle::evaluate_fitness_values()
{
  netw->removeTempElements();
  for(int j=0; j<capacitors.size(); j++)
    netw->setCapacitorQ(capacitors[j].getNode(), capacitors[j].getUnitNum(), capacitors[j].getUnitValue());
  for(int j=0; j<disp_generators.size(); j++)
    netw->setDispGenerator(disp_generators[j].getNode(), disp_generators[j].getUnitNum(), disp_generators[j].getUnitPvalue(), disp_generators[j].getUnitQvalue());
  for(int j=0; j<solar_generators.size(); j++)
    netw->setSolarGenerator(solar_generators[j].getNode(), solar_generators[j].getUnitNum(), solar_generators[j].getUnitPvalue(), solar_generators[j].getUnitQvalue());
  for(int j=0; j<wind_generators.size(); j++)
    netw->setWindGenerator(wind_generators[j].getNode(), wind_generators[j].getUnitNum(), wind_generators[j].getUnitPvalue(), wind_generators[j].getUnitQvalue());
  netw->setTapChanger(tapChanger);

  hasValidFitness = true;
  double activeLosses = 0;
  double voltageStabilityIndex = 0;
  double voltageDeviation = 0;

  if (scenarionum == 1)
  {
    for (int j=0; j<4; j++)
    {
      string season;
      switch (j) {
        case 0: season = "summer"; break;
        case 1: season = "autumn"; break;
        case 2: season = "winter"; break;
        case 3: season = "spring"; break;
      }
      for(int i=1; i<=24; i++)
      {
        netw->solvePowerFlow(season, i);
        activeLosses += netw->getActiveLosses();
        voltageStabilityIndex += netw->getVoltageStabilityIndexSum();
        voltageDeviation += netw->getVoltageDeviationSum();
        if( !netw->isWithinVoltageLimits() )
          hasValidFitness = false;
      }
    }
  }
  else if (scenarionum == 2)
  {
    netw->solvePowerFlow("summer", 12);
    activeLosses = netw->getActiveLosses();
    voltageStabilityIndex = netw->getVoltageStabilityIndexSum();
    voltageDeviation = netw->getVoltageDeviationSum();
    if( !netw->isWithinVoltageLimits() )
      hasValidFitness = false;
  }
  else
  {
    cerr<<"Error: Wrong scenario number: "<<scenarionum<<endl;
    exit(EXIT_SUCCESS);
  }
  fitness_values[0] = activeLosses;
  fitness_values[1] = 1/voltageStabilityIndex;
  fitness_values[2] = voltageDeviation;

  /*
  if (considerActiveLosses && (fitness_values[0] > base_fitnes_vals[0]))
    hasValidFitness = false;
  if (considerVSI && (fitness_values[1] > base_fitnes_vals[1]))
    hasValidFitness = false;
  if (considerVD && (fitness_values[2] > base_fitnes_vals[2]))
    hasValidFitness = false;
  //*/
  Vmin = netw->getMinimumVoltageMagnitude();
  Vmax = netw->getMaximumVoltageMagnitude();
  Vmin_lim = netw->getVmin();
  Vmax_lim = netw->getVmax();
}

void PSO_particle::evaluate_PBEST()
{
  bool current_dominates = true;
  bool PBEST_dominates = true;
  for( int i=0; i<fitness_values.size(); i++)
    if( fitness_values[i] >= PBEST_fitness_vals[i])
    {
      current_dominates = false;
      break;
    }
  for( int i=0; i<fitness_values.size(); i++)
    if( PBEST_fitness_vals[i] >= fitness_values[i])
    {
      PBEST_dominates = false;
      break;
    }
  if( !current_dominates && !PBEST_dominates)
  {
    double random_num = fmod(rand()*.1, 1);
    if( random_num >= 0.5)
      current_dominates = true;
  }
  if ( !hasValidFitness )
  {
    current_dominates = false;
  }
  else
  {
    hasValidPBEST = true;
  }
  if(current_dominates || PBEST_fitness_vals[0] == 0 || !hasValidPBEST)
  {
    cap_PBEST = capacitors;
    disp_gen_PBEST = disp_generators;
    solar_gen_PBEST = solar_generators;
    wind_gen_PBEST = wind_generators;
    PBEST_fitness_vals = fitness_values;
    tapChanger_PBEST = tapChanger;
  }
}


void PSO_particle::move_particle_layer_index_based(PSO_particle& GBEST)
{
  double R1 = fmod(rand()/1e6, 1);
  double R2 = fmod(rand()/1e6, 1);

  double C1 = 1.0;
  double C2 = 2.0;
  for(int i=0; i<disp_generators.size(); i++)
  {
    disp_gen_vel[i].setPvalue(W*disp_gen_vel[i].getPvalue() + C1*R1*(disp_gen_PBEST[i].getPvalue() -  disp_generators[i].getPvalue()) + C2*R2*(GBEST.disp_generators[i].getPvalue() - disp_generators[i].getPvalue()));
    disp_gen_vel[i].setQvalue(W*disp_gen_vel[i].getQvalue() + C1*R1*(disp_gen_PBEST[i].getQvalue() - disp_generators[i].getQvalue()) + C2*R2*(GBEST.disp_generators[i].getQvalue() - disp_generators[i].getQvalue()));
    disp_generators[i].setPvalue( disp_generators[i].getPvalue() + disp_gen_vel[i].getPvalue());
    disp_generators[i].setQvalue( disp_generators[i].getQvalue() + disp_gen_vel[i].getQvalue());

    if (disp_generators[i].getPvalue() > dispMaxP) disp_generators[i].setPvalue(dispMaxP);
    if (disp_generators[i].getPvalue() < dispMinP) disp_generators[i].setPvalue(dispMinP);
    if (disp_generators[i].getQvalue() > dispMaxQ) disp_generators[i].setQvalue(dispMaxQ);
    if (disp_generators[i].getQvalue() < dispMinQ) disp_generators[i].setQvalue(dispMinQ);

    disp_gen_vel[i].setNode( round( W*disp_gen_vel[i].getNode() + C1*R1*(disp_gen_PBEST[i].getNode() - disp_generators[i].getNode()) + C2*R2*(GBEST.disp_generators[i].getNode() - disp_generators[i].getNode())));
    disp_generators[i].setNode( disp_generators[i].getNode() + disp_gen_vel[i].getNode());
    if (disp_generators[i].getNode() > maxNode) disp_generators[i].setNode(maxNode);
    if (disp_generators[i].getNode() < minNode) disp_generators[i].setNode(minNode);
  }
  for(int i=0; i<capacitors.size(); i++)
  {
    cap_vel[i].setNode( round(W*cap_vel[i].getNode() + C1*R1*(cap_PBEST[i].getNode() - capacitors[i].getNode()) + C2*R2*(GBEST.capacitors[i].getNode() - capacitors[i].getNode())));
    capacitors[i].setNode( capacitors[i].getNode() + cap_vel[i].getNode());
    if (capacitors[i].getNode() > maxNode) capacitors[i].setNode(maxNode);
    if (capacitors[i].getNode() < minNode) capacitors[i].setNode(minNode);

    cap_vel[i].setUnitNum( round(W*cap_vel[i].getUnitNum() + C1*R1*(cap_PBEST[i].getUnitNum() - capacitors[i].getUnitNum()) + C2*R2*(GBEST.capacitors[i].getUnitNum() - capacitors[i].getUnitNum())));
    capacitors[i].setUnitNum(capacitors[i].getUnitNum() + cap_vel[i].getUnitNum());
    if (capacitors[i].getUnitNum() > maxCapUnitNum) capacitors[i].setUnitNum(maxCapUnitNum);
    if (capacitors[i].getUnitNum() < minCapUnitNum) capacitors[i].setUnitNum(minCapUnitNum);
  }
  for(int i=0; i<wind_generators.size(); i++)
  {
    wind_gen_vel[i].setNode( round (W*wind_gen_vel[i].getNode() + C1*R1*(wind_gen_PBEST[i].getNode() - wind_generators[i].getNode()) + C2*R2*(GBEST.wind_generators[i].getNode() - wind_generators[i].getNode())));
    wind_generators[i].setNode( wind_generators[i].getNode() + wind_gen_vel[i].getNode());
    if (wind_generators[i].getNode() > maxNode) wind_generators[i].setNode(maxNode);
    if (wind_generators[i].getNode() < minNode) wind_generators[i].setNode(minNode);

    wind_gen_vel[i].setUnitNum( round( W*wind_gen_vel[i].getUnitNum() + C1*R1*(wind_gen_PBEST[i].getUnitNum() - wind_generators[i].getUnitNum()) + C2*R2*(GBEST.wind_generators[i].getUnitNum() - wind_generators[i].getUnitNum())));
    wind_generators[i].setUnitNum(wind_generators[i].getUnitNum() + wind_gen_vel[i].getUnitNum());
    if (wind_generators[i].getUnitNum() > maxWindUnitNum) wind_generators[i].setUnitNum(maxWindUnitNum);
    if (wind_generators[i].getUnitNum() < minWindUnitNum) wind_generators[i].setUnitNum(minWindUnitNum);
  }
  for(int i=0; i<solar_generators.size(); i++)
  {
    solar_gen_vel[i].setNode(round (W*solar_gen_vel[i].getNode() + C1*R1*(solar_gen_PBEST[i].getNode() - solar_generators[i].getNode()) + C2*R2*(GBEST.solar_generators[i].getNode() - solar_generators[i].getNode())));
    solar_generators[i].setNode( solar_generators[i].getNode() + solar_gen_vel[i].getNode());
    if (solar_generators[i].getNode() > maxNode) solar_generators[i].setNode(maxNode);
    if (solar_generators[i].getNode() < minNode) solar_generators[i].setNode(minNode);

    solar_gen_vel[i].setUnitNum( round (W*solar_gen_vel[i].getUnitNum() + C1*R1*(solar_gen_PBEST[i].getUnitNum() - solar_generators[i].getUnitNum()) + C2*R2*(GBEST.solar_generators[i].getUnitNum() - solar_generators[i].getUnitNum())));
    solar_generators[i].setUnitNum(solar_generators[i].getUnitNum() + solar_gen_vel[i].getUnitNum());
    if (solar_generators[i].getUnitNum() > maxSolarUnitNum) solar_generators[i].setUnitNum(maxSolarUnitNum);
    if (solar_generators[i].getUnitNum() < minSolarUnitNum) solar_generators[i].setUnitNum(minSolarUnitNum);
  }
  tapChanger_vel = round(W*tapChanger + C1*R1*(tapChanger_PBEST - tapChanger) + C2*R2*(GBEST.tapChanger - tapChanger));
  tapChanger = tapChanger + tapChanger_vel;
  if (tapChanger > maxTapChanger) tapChanger = maxTapChanger;
  if (tapChanger < minTapChanger) tapChanger = minTapChanger;
}

void PSO_particle::reset_particle(int dispGType)
{
  minNode = 1;
  maxNode = netw->getNumONodes()-1;
  int nodePosibilities = maxNode - minNode;
  dispMinP = 0;
  dispMaxP = 2e6;
  dispMinQ = 0;
  dispMaxQ = 2e6;
  solarP = .25e6;
  solarQ = 0;
  windP = .25e6;
  windQ = 0;
  capQ = .125e6;

  for(int j=0; j<capacitors.size(); j++)
  {
    int random_node = rand() % nodePosibilities + minNode;
    int random_unit_num = rand() % maxCapUnitNum + minCapUnitNum;
    capacitors[j].setTempCapacitor(random_node, random_unit_num, capQ);
    cap_vel[j].setTempCapacitor(0, 0, 0);
    cap_guide[j] = random_node;
  }
  for(int j=0; j<wind_generators.size(); j++)
  {
    int random_node = rand() % nodePosibilities + minNode;
    int random_unit_num = rand() % maxWindUnitNum + minWindUnitNum;
    wind_generators[j].setTempGenerator(random_node, random_unit_num, windP, windQ);
    wind_gen_vel[j].setTempGenerator(0,0,0,0);
    wind_guide[j] = random_node;
  }
  for(int j=0; j<solar_generators.size(); j++)
  {
    int random_node = rand() % nodePosibilities + minNode;
    int random_unit_num = rand() % maxWindUnitNum + minWindUnitNum;
    solar_generators[j].setTempGenerator(random_node, random_unit_num, solarP, solarQ);
    solar_gen_vel[j].setTempGenerator(0,0,0,0);
    solar_guide[j] = random_node;
  }
  for(int j=0; j<disp_generators.size(); j++)
  {
    double disp_Pvalue, disp_Qvalue;
    if (dispGType == 2)
      disp_Pvalue = 0;
    else
      disp_Pvalue = rand() % (dispMaxP - dispMinP) + dispMinP;
    if (dispGType == 1)
      disp_Qvalue = 0;
    else
      disp_Qvalue = rand() % (dispMaxQ - dispMinQ) + dispMinQ;
    int random_node = rand() % nodePosibilities + minNode;
    int disp_unit_num = 1;

    disp_generators[j].setTempGenerator(random_node, disp_unit_num, disp_Pvalue, disp_Qvalue);
    disp_gen_vel[j].setTempGenerator(0,0,0,0);
    disp_guide[j] = random_node;
  }
  if (maxTapChanger - minTapChanger != 0)
    tapChanger = rand() % (maxTapChanger - minTapChanger) + minTapChanger;
  else
    tapChanger = 0;
}

void PSO_particle::createOctaveGraph(string graphname, string season)
{

  evaluate_fitness_values();

  bool filecreated = false;
  if (considerActiveLosses)
  {
    string appendable("_Ploss");
    netw->createOctaveGraph(graphname, season, &network::getActiveLosses, !filecreated, appendable);
    filecreated = true;
  }
  if (considerVSI)
  {
    string appendable("_VSI");
    netw->createOctaveGraph(graphname, season, &network::getMinimumVSI, !filecreated, appendable);
    filecreated = true;
  }
  if (considerVD)
  {
    string appendable("_VD");
    netw->createOctaveGraph(graphname, season, &network::getMaximumVoltageDeviation, !filecreated, appendable);
    filecreated = true;
  }
}


void PSO_particle::log_particle_positions(const char* fname)
{
  ofstream os(fname, fstream::app);
  os<<"\tCAPACITORS:"<<endl;
  for (int j=0; j<capacitors.size(); j++)
  {
    os<<"\tNO "<<j+1<<"  NODE: "<<capacitors[j].getNode();
    os<<"  Q: "<<capacitors[j].getValue()/1e3<<" kVAr"<<endl;
  }
  os<<endl;
  os<<"\tSOLAR GENERATORS:"<<endl;
  for (int j=0; j<solar_generators.size(); j++)
  {
    os<<"\tNO "<<j+1<<"  NODE: "<<solar_generators[j].getNode();
    os<<"  P: "<<solar_generators[j].getPvalue()/1e3 <<" kW";
    os<<"  Q: "<<solar_generators[j].getQvalue()/1e3 <<" kVAr"<<endl;
  }
  os<<endl;
  os<<"\tWIND GENERATORS:"<<endl;
  for (int j=0; j<wind_generators.size(); j++)
  {
    os<<"\tNO "<<j+1<<"  NODE: "<<wind_generators[j].getNode();
    os<<"  P: "<<wind_generators[j].getPvalue()/1e3 <<" kW";
    os<<"  Q: "<<wind_generators[j].getQvalue()/1e3 <<" kVAr"<<endl;
  }
  os<<endl;
  os<<"\tDISPACHABLE GENERATORS:"<<endl;
  for (int j=0; j<disp_generators.size(); j++)
  {
    os<<"\tNO "<<j+1<<"  NODE: "<<disp_generators[j].getNode();
    os<<"  P: "<<disp_generators[j].getPvalue()/1e6 <<" MW";
    os<<"  Q: "<<disp_generators[j].getQvalue()/1e6 <<" MVAr"<<endl;
  }
  os<<endl;
}

void PSO_particle::log_particle_velocities(const char* fname)
{
  ofstream os(fname, fstream::app);
  os<<"\tCAPACITOR VELOCITIES:"<<endl;
  for (int j=0; j<cap_vel.size(); j++)
  {
    os<<"\tNO "<<j+1<<"  NODE: "<<cap_vel[j].getNode();
    os<<" (with guide "<<cap_guide[j]<<" )";
    os<<"  Q: "<<cap_vel[j].getValue()/1e3<<" kVAr"<<endl;
  }
  os<<endl;
  os<<"\tSOLAR GENERATOR VELOCITIES:"<<endl;
  for (int j=0; j<solar_gen_vel.size(); j++)
  {
    os<<"\tNO "<<j+1<<"  NODE: "<<solar_gen_vel[j].getNode();
    os<<" (with guide "<<solar_guide[j]<<" )";
    os<<"  P: "<<solar_gen_vel[j].getPvalue()/1e3 <<" kW";
    os<<"  Q: "<<solar_gen_vel[j].getQvalue()/1e3 <<" kVAr"<<endl;
  }
  os<<endl;
  os<<"\tWIND GENERATOR VELOCITIES:"<<endl;
  for (int j=0; j<wind_gen_vel.size(); j++)
  {
    os<<"\tNO "<<j+1<<"  NODE: "<<wind_gen_vel[j].getNode();
    os<<" (with guide "<<wind_guide[j]<<" )";
    os<<"  P: "<<wind_gen_vel[j].getPvalue()/1e3 <<" kW";
    os<<"  Q: "<<wind_gen_vel[j].getQvalue()/1e3 <<" kVAr"<<endl;
  }
  os<<endl;
  os<<"\tDISPACHABLE GENERATORS VELOCITIES:"<<endl;
  for (int j=0; j<disp_gen_vel.size(); j++)
  {
    os<<"\tNO "<<j+1<<"  NODE: "<<disp_gen_vel[j].getNode();
    os<<" (with guide "<<disp_guide[j]<<" )";
    os<<"  P: "<<disp_gen_vel[j].getPvalue()/1e6 <<" MW";
    os<<"  Q: "<<disp_gen_vel[j].getQvalue()/1e6 <<" MVAr"<<endl;
  }
  os<<endl;
}

void PSO_particle::log_particle_bests(const char *fname)
{
  ofstream os(fname, fstream::app);
  os<<"\tCAPACITOR:"<<endl;
  for (int j=0; j<cap_PBEST.size(); j++)
  {
    os<<"\tNO "<<j+1<<"  NODE: "<<cap_PBEST[j].getNode();
    os<<"  Q: "<<cap_PBEST[j].getValue()/1e3<<" kVAr"<<endl;
  }
  os<<endl;
  os<<"\tSOLAR GENERATOR VELOCITIES:"<<endl;
  for (int j=0; j<solar_gen_PBEST.size(); j++)
  {
    os<<"\tNO "<<j+1<<"  NODE: "<<solar_gen_PBEST[j].getNode();
    os<<"  P: "<<solar_gen_PBEST[j].getPvalue()/1e3 <<" kW";
    os<<"  Q: "<<solar_gen_PBEST[j].getQvalue()/1e3 <<" kVAr"<<endl;
  }
  os<<endl;
  os<<"\tWIND GENERATOR VELOCITIES:"<<endl;
  for (int j=0; j<wind_gen_PBEST.size(); j++)
  {
    os<<"\tNO "<<j+1<<"  NODE: "<<wind_gen_PBEST[j].getNode();
    os<<"  P: "<<wind_gen_PBEST[j].getPvalue()/1e3 <<" kW";
    os<<"  Q: "<<wind_gen_PBEST[j].getQvalue()/1e3 <<" kVAr"<<endl;
  }
  os<<endl;
  os<<"\tDISPACHABLE GENERATORS VELOCITIES:"<<endl;
  for (int j=0; j<disp_gen_PBEST.size(); j++)
  {
    os<<"\tNO "<<j+1<<"  NODE: "<<disp_gen_PBEST[j].getNode();
    os<<"  P: "<<disp_gen_PBEST[j].getPvalue()/1e6 <<" MW";
    os<<"  Q: "<<disp_gen_PBEST[j].getQvalue()/1e6 <<" MVAr"<<endl;
  }
  os<<endl;
  int fitn_pointer = 0;
  //if (considerActiveLosses)
    os<<"\tAVTIVE LOSSES: "<<fitness_values[fitn_pointer++]*1e3<<" kW"<<endl;
  //if (considerVSI)
    os<<"\tVOLTAGE STABILITY INDEX: "<<fitness_values[fitn_pointer++]<<endl;
  //if (considerVD)
    os<<"\tVOLTAGE DEVIATION: "<<fitness_values[fitn_pointer++]<<endl;
  os<<endl;
}

void PSO_particle::log_particle_fitness_vals(const char *fname)
{
  ofstream os(fname, fstream::app);
  if (true)
  {
    os<<"Vmin_lim = "<<Vmin_lim<<endl;
    os<<"Vmin = "<<Vmin<<endl;
    os<<"Vmax_lim = "<<Vmax_lim<<endl;
    os<<"Vmax = "<<Vmax<<endl;
  }
  os<<"FITNESS VALUES:"<<endl;
  int fitn_pointer = 0;
//  if (considerActiveLosses)
    os<<"\tAVTIVE LOSSES: "<<fitness_values[fitn_pointer++]*1e3<<" kW"<<endl;
//  if (considerVSI)
    os<<"\tVOLTAGE STABILITY INDEX: "<<1/fitness_values[fitn_pointer++]<<endl;
//  if (considerVD)
    os<<"\tVOLTAGE DEVIATION: "<<fitness_values[fitn_pointer++]<<endl;
  os<<endl;
}

void PSO_particle::log_table_form(ofstream& os)
{
    if (hasValidFitness)
      os<<"Particle satisfies the constraints of the problem"<<endl;
    else
      os<<"Particle does not satisfie the constraints of the problem"<<endl;
    os<<"Tap changer: "<<tapChanger<<", ";
    os<<"nodes: ";
    for(int i=0; i<disp_generators.size(); i++)
    {
      os<<netw->convertIndex_LayerB2Longitudinal(disp_generators[i].getNode());
      if (i < disp_generators.size() - 1)
        os<<", ";
    }
    os<<"   ";
    os<<"Inj P: ";
    for(int i=0; i<disp_generators.size(); i++ )
    {
      os<<disp_generators[i].getPvalue() / 1e6;
      if (i < disp_generators.size() - 1)
        os<<", ";
    }

    os<<"  ";
    os<<"Inj Q: ";
    for(int i=0; i<disp_generators.size(); i++ )
    {
      os<<disp_generators[i].getQvalue() / 1e6;
      if (i < disp_generators.size() - 1)
        os<<", ";
    }

    os<<"   ";
    os<<"Active Losses: "<<getFitnessValue(0)*1e3;

    os<<"   ";
    os<<"VSI: "<<1/getFitnessValue(1);
}

void PSO_particle::log_master_form(ofstream& os)
{
    os<<fixed;
    int sorted_indexes[disp_generators.size()];
    for(int i=0; i<disp_generators.size(); i++)
      sorted_indexes[i] = i;

    for(int i=0; i<disp_generators.size()-1; i++)
      for(int j=0; j<disp_generators.size()-i-1; j++)
      {
        if (disp_generators[sorted_indexes[j]].getNode() > disp_generators[sorted_indexes[j+1]].getNode())
          swap_int(sorted_indexes+j, sorted_indexes+j+1);
      }

    if(hasTapChanger)
      os<<-tapChanger<<" ";
    for(int i=0; i<disp_generators.size(); i++)
    {
      os<<disp_generators[sorted_indexes[i]].getNode();
      if(i<disp_generators.size()-1)
        os<<",";
    }
    os<<" ";

    if (disp_type == 3)
      os<<"P:";

    if(disp_type == 1 || disp_type == 3)
    {
      for(int i=0; i<disp_generators.size(); i++)
      {
        os<<setprecision(3)<<disp_generators[sorted_indexes[i]].getPvalue()*1e-6;
        if(i<disp_generators.size()-1)
          os<<",";
      }
      os<<" ";
    }

    if (disp_type == 3)
      os<<"Q:";
    if(disp_type == 2 || disp_type == 3)
    {
      for(int i=0; i<disp_generators.size(); i++)
      {
        os<<setprecision(3)<<disp_generators[sorted_indexes[i]].getQvalue()*1e-6;
        if(i<disp_generators.size()-1)
          os<<",";
      }
      os<<" ";
    }

    if (switchnum == 4 || switchnum == 5 || switchnum == 6 || switchnum == 7)
      os<<setprecision(2)<<fitness_values[0]*1e3<<" ";
    if (switchnum == 2 || switchnum == 3 || switchnum == 6 || switchnum == 7)
      os<<setprecision(2)<<1/fitness_values[1]<<" ";
    if (switchnum == 1 || switchnum == 3 || switchnum == 5 || switchnum == 7)
      os<<setprecision(2)<<fitness_values[2]<<" ";

    os<<endl;
}

bool operator< (const PSO_particle& p1, const PSO_particle& p2)
{
  if (!p1.hasValidFitness)
    return false;

  if (p1.considerActiveLosses && p1.fitness_values[0] > p2.fitness_values[0])
    return false;
  if (p1.considerVSI && p1.fitness_values[1] > p2.fitness_values[1])
    return false;
  if (p1.considerVD && p1.fitness_values[2] > p2.fitness_values[2])
    return false;

  return true;
}
