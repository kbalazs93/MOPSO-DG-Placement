#include "particleswarm.h"
#include <cstdlib>
#include <sstream>


particle_swarm::particle_swarm()
{
  cerr<<"Error: You need to add values to the particles"<<endl;
  exit(EXIT_SUCCESS);
}

particle_swarm::particle_swarm(network* netw, int NOP, int MAXITER, double C1, double C2, double W, double Wdamp, int NREP, int Ngrid, bool hasTapChanger): REP(NREP, Ngrid)
{
  this->netw = netw;
  this->NOP = NOP;
  this->MAXITER = MAXITER;
  this->C1 = C1;
  this->C2 = C2;
  this->W = W;
  this->Wdamp = Wdamp;

  particles = NULL;

  initialize_scenario1(7, NOP, MAXITER, C1, C2, W, Wdamp, hasTapChanger);
}

particle_swarm::particle_swarm(particle_swarm& ps)
{
  this->netw = ps.netw;
  cap_num = ps.cap_num;
  disp_num = ps.disp_num;
  solar_num = ps.solar_num;
  wind_num = ps.wind_num;
  fitness_num = ps.fitness_num;
  dispGType = ps.dispGType;
  this->NOP = ps.NOP;
  this->MAXITER = ps.MAXITER;

  C1 = ps.C1;
  C2 = ps.C2;
  W = ps.W;
  Wdamp = ps.Wdamp;


  for(int i=0; i<NOP; i++)
    *particles[i] = *ps.particles[i];
}

bool particle_swarm::doMutation(int iter)
{
  double rr = pow(1-((double)iter-1)/((double)MAXITER-1), (1/mu));
  double random_num = fmod(rand()/1e4, 1);
}

void particle_swarm::initialize_scenario1(int switchnum, int NOP, int MAXITER, double C1, double C2, double W, double Wdamp, bool hasTapChanger)
{

  clear_particles();
  cap_num = 1;
  disp_num = 1;
  solar_num = 1;
  wind_num = 1;
  dispGType = 1;
  this->switchnum = switchnum;
  scenarionum = 1;
  particles = new PSO_particle*[NOP];
  this->NOP = NOP;
  this->MAXITER = MAXITER;

  for (int i=0; i<NOP; i++)
  {
    particles[i] = new PSO_particle(netw, cap_num, disp_num, wind_num, solar_num, dispGType, switchnum, scenarionum, C1, C2, W, Wdamp, hasTapChanger);
    particles[i]->reset_particle(dispGType);
  }

  fitness_num = particles[0]->getFitnessNum();
  REP.setFitnessVals( fitness_num );
}

void particle_swarm::initialize_scenario2(int numOG, int switchnum, int dispType, int NOP, int MAXITER, double C1, double C2, double W, double Wdamp, bool hasTapChanger)
{
  clear_particles();
  cap_num = 0;
  disp_num = numOG;
  solar_num = 0;
  wind_num = 0;
  this->switchnum = switchnum;
  scenarionum = 2;
  dispGType = dispType;
  particles = new PSO_particle*[NOP];
  this->NOP = NOP;
  this->MAXITER = MAXITER;

  for (int i=0; i<NOP; i++)
  {
    particles[i] = new PSO_particle(netw, cap_num, disp_num, wind_num, solar_num, dispGType, switchnum, scenarionum, C1, C2, W, Wdamp, hasTapChanger);
    particles[i]->reset_particle(dispGType);
  }
  fitness_num = particles[0]->getFitnessNum();
  REP.setFitnessVals( fitness_num );
}

void particle_swarm::clear_particles()
{
  if (particles != NULL)
    delete[] particles;
}

PSO_particle particle_swarm::fuzzy_choice()
{
  return REP.fuzzy_choice();
}

void particle_swarm::particle_swarm_optimalization()
{
  for(int i=0; i<MAXITER; i++)
  {
    for(int j=0; j<NOP; ++j)
    {
      if ( i==0 || doMutation(i) )
      {
        particles[j]->reset_particle(dispGType);
      }
      else
      {
        PSO_particle leader = REP.select_lead();
        particles[j]->move_particle_layer_index_based(leader);
      }
      int num = 1;
      while (true)
      {
        particles[j]->evaluate_fitness_values();
        particles[j]->evaluate_PBEST();
        REP.archive_controler(*particles[j]);
        if (num % 1000 == 0)
        {
          netw->relax_constraints();
        }
        if ( REP.repositorySize() )
        {

          break;
        }
        else
        {
          particles[j]->reset_particle(dispGType);
          num++;
        }
      }
      particles[j]->inertiaDamping();
    }
    cout<<"\rStatus: "<<(i+1)*100/MAXITER<<" %"<<flush;
  }
  stringstream output;
  output<<"./output/MOPSO/original/repository"<<NOP<<MAXITER<<".txt";
  string out = output.str();
  REP.log_repository_particles(out.c_str());
}


void particle_swarm::log_particle_positions(const char* fname)
{
  ofstream os(fname);
  os.close();
  os.open(fname, ofstream::app);
  for(int i=0; i<NOP; i++)
  {
    os<<endl<<endl;
    os<<"PARTICLE "<<i+1<<":"<<endl<<endl;
    particles[i]->log_particle_positions(fname);
  }
}

void particle_swarm::log_particle_velocities(const char* fname)
{
  ofstream os(fname);
  os.close();
  os.open(fname, ofstream::app);
  for(int i=0; i<NOP; i++)
  {
    os<<endl<<endl;
    os<<"PARTICLE "<<i+1<<":"<<endl<<endl;
    particles[i]->log_particle_velocities(fname);
  }
}

void particle_swarm::log_particle_bests(const char* fname)
{
  ofstream os(fname);
  os.close();
  os.open(fname, ofstream::app);
  for(int i=0; i<NOP; i++)
  {
    os<<endl<<endl;
    os<<"PARTICLE "<<i+1<<":"<<endl<<endl;
    particles[i]->log_particle_bests(fname);
  }
}

void particle_swarm::log_particle_fitness_vals(const char* fname)
{
  ofstream os(fname);
  os.close();
  os.open(fname, ofstream::app);
  for(int i=0; i<NOP; i++)
  {
    os<<endl<<endl;
    os<<"PARTICLE "<<i+1<<":"<<endl<<endl;
    particles[i]->log_particle_fitness_vals(fname);
  }
}
