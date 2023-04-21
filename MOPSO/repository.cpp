#include "repository.h"
#include "../NetworkModel/network.h"
#include <iostream>
#include <iomanip>
#include <cmath>

using namespace std;

gbest_repository::gbest_repository()
{
  NREP = 100;
  Ngrid = 7;
  setFitnessVals(3);
  alpha = new double[NREP];
}

gbest_repository::gbest_repository(int NREP, int Ngrid)
{
  this->NREP = NREP;
  this->Ngrid = Ngrid;
  setFitnessVals(3);
  alpha = new double[NREP];
}

void gbest_repository::setFitnessVals(const int n)
{
    delete[] grid_mins;
    delete[] grid_maxes;
    delete[] f_min;
    delete[] f_max;
    fitness_vals = n;
    grid_mins = new double[n];
    grid_maxes = new double[n];
    for(int i=0; i<fitness_vals; i++)
        grid_mins[i]=0;
    for(int i=0; i<fitness_vals; i++)
        grid_maxes[i]=0;

    f_min = new double[fitness_vals];
    f_max = new double[fitness_vals];
  }


void gbest_repository::adaptive_grid(bool reconstruct)
{
  if( REP_particles.size() < 2 )
    return;
  vector<double> last_fitness_vals;
  REP_particles[REP_particles.size()-1].getFitnessValues(&last_fitness_vals);

  bool resizeNeeded = false;
  for(int i=0; i<fitness_vals; ++i)
  {
    if (grid_mins[i] > last_fitness_vals[i])
    {
      resizeNeeded = true;
      break;
    }
    if (grid_maxes[i] < last_fitness_vals[i])
    {
      resizeNeeded = true;
      break;
    }
  }

  if (reconstruct || resizeNeeded)
  {
    for(int i=0; i<fitness_vals; i++)
    {
      grid_mins[i] = REP_particles[0].getFitnessValue(i);
      grid_maxes[i] = REP_particles[0].getFitnessValue(i);
    }
    for(int i=0; i<REP_particles.size(); ++i)
    {
      for(int j=0; j<fitness_vals; ++j)
      {
        if(REP_particles[i].getFitnessValue(j) < grid_mins[j])
          grid_mins[j] = REP_particles[i].getFitnessValue(j);
        if(REP_particles[i].getFitnessValue(j) > grid_maxes[j])
          grid_maxes[j] = REP_particles[i].getFitnessValue(j);
      }
    }
    for(int i=0; i<fitness_vals; ++i)
    {
      grid_mins[i] -= (grid_maxes[i] - grid_mins[i]) / 50;
      grid_maxes[i] += (grid_maxes[i] - grid_mins[i]) / 50;
    }

    for(int i=0; i<REP_particles.size(); i++)
    {
      int cube_num = 0;
      for (int j=0; j<fitness_vals; j++)
      {
        cube_num += (REP_particles[i].getFitnessValue(j) - grid_mins[j] ) / (grid_maxes[j] - grid_mins[j]) * Ngrid;
        cube_num *= Ngrid;
      }
      cube_num /= Ngrid;
      hypercube_coordinates[i] = cube_num;
    }
  }
  int cube_num = 0;
  for (int j=0; j<fitness_vals; j++)
  {
    cube_num += (REP_particles[REP_particles.size()-1].getFitnessValue(j) - grid_mins[j] ) / (grid_maxes[j] - grid_mins[j]) * Ngrid;
    cube_num *= Ngrid;
  }
  cube_num /= Ngrid;
  hypercube_coordinates[REP_particles.size()-1] = cube_num;
}


int gbest_repository::roulette_wheel(int mode)
{
  const int cube_num = Ngrid*Ngrid*Ngrid;
  double hypercube_partnum[cube_num];
  for(int i=0; i<cube_num; i++)
  {
    hypercube_partnum[i] = 0;
  }
  for(int i=0; i<hypercube_coordinates.size(); i++)
  {
    hypercube_partnum[hypercube_coordinates[i]]++;
  }

  double roulette_numbers[cube_num];
  for(int i=0; i<cube_num; i++)
  {
    if(mode != 0 && hypercube_partnum[i] != 0)
        roulette_numbers[i] = 1.0 / hypercube_partnum[i];
    else
        roulette_numbers[i] = hypercube_partnum[i];
  }

  double divisor = 0;
  for(int i=0; i<cube_num; i++)
    divisor += roulette_numbers[i];

  double roulette = fmod(rand(), divisor);


  int selected;
  for (int i=0; i<cube_num; i++)
    if (hypercube_partnum[i] != 0)
    {
      selected = i;
      break;
    }
  for (int i=0; roulette > 0 && i<cube_num; i++)
  {
    roulette -= roulette_numbers[i];
    selected = i;
  }

  int offset = rand() % (int)hypercube_partnum[selected];
  for (int i=0; offset>=0; i++)
  {
    if (hypercube_coordinates[i] == selected)
    {
      if (offset == 0)
        return i;
      offset--;
    }
  }
}

void gbest_repository::archive_controler(PSO_particle& part)
{
  if ( !part.getFitnessValidation() )
    return;
  if ( REP_particles.size() == 0 )
  {
    REP_particles.push_back(part);
    hypercube_coordinates.push_back(0);
    return;
  }
  for( int i=0; i<REP_particles.size(); i++)
  {
    if (REP_particles[i]<part)
      return;
  }
  for( int i=0; i<REP_particles.size(); i++)
  {
    if (part<REP_particles[i])
    {
      REP_particles.erase(REP_particles.begin() + i);
      hypercube_coordinates.erase(hypercube_coordinates.begin() + i);
    }
  }
  adaptive_grid(true);

  REP_particles.push_back(part);
  hypercube_coordinates.push_back(0);
  adaptive_grid();

  if (REP_particles.size() > NREP)
  {
    int removable_index = roulette_wheel();
    REP_particles.erase(REP_particles.begin() + removable_index);
    hypercube_coordinates.erase(hypercube_coordinates.begin() + removable_index);
    adaptive_grid(true);
  }
}

PSO_particle gbest_repository::select_lead()
{
  int lead_index = roulette_wheel(1);
  return REP_particles[lead_index];
}

int gbest_repository::fuzzy_choice_within_range()
{
  vector<double> ref_point;
  REP_particles[0].getBaseFitnessValues(&ref_point);

  double first_ref = ref_point[0];
  double second_ref = ref_point[1];

  double f_min_wrange[fitness_vals], f_max_wrange[fitness_vals];

  int feasnum = 0;
  for(int i=0; i<REP_particles.size(); i++)
    if ( REP_particles[i].getFitnessValue(0) < first_ref &&
          REP_particles[i].getFitnessValue(1) < second_ref)
      {
        feasnum++;
      }

  double obj1[feasnum], obj2[feasnum];
  int original_ind[feasnum];

  int j = 0;
  for (int i=0; i<REP_particles.size(); i++)
  {
      if (REP_particles[i].getFitnessValue(0) < first_ref && REP_particles[i].getFitnessValue(1) < second_ref)
        {
            obj1[j] = REP_particles[i].getFitnessValue(0);
            obj2[j] = REP_particles[i].getFitnessValue(1);
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

    f_min_wrange[i] = cur_obj[0];
    f_max_wrange[i] = cur_obj[0];
    for(int j=0; j<REP_particles.size(); j++)
    {
      if (cur_obj[i] < f_min_wrange[i])
        f_min_wrange[i] = cur_obj[i];
      if (cur_obj[i] > f_max_wrange[i])
        f_max_wrange[i] = cur_obj[i];
    }
  }

  double alphamax_wrange=0;
  int alphamax_wrange_ind=0;
  double alpha_wrange[feasnum];

  for(int i=0; i<feasnum; i++)
  {
    alpha_wrange[i] = 0;
    for(int j=0; j<fitness_vals; j++)
    {
      double* cur_obj;
      if (j==0)
        cur_obj = obj1;
      else if (j==1)
        cur_obj = obj2;

      if (cur_obj[i] == f_min_wrange[j])
        alpha_wrange[i] += 1.0;
      else if (cur_obj[i] == f_max_wrange[j])
        continue;
      else
        alpha_wrange[i] += ( f_max_wrange[j] - cur_obj[i] ) / (f_max_wrange[j] - f_min_wrange[j]) ;
      if (alpha_wrange[i] > alphamax_wrange)
      {
        alphamax_wrange = alpha_wrange[i];
        alphamax_wrange_ind = i;
      }
    }
  }
  return original_ind[alphamax_wrange_ind];
}

PSO_particle gbest_repository::fuzzy_choice()
{
  for(int i=0; i<fitness_vals; i++)
  {
    f_min[i] = REP_particles[0].getFitnessValue(i);
    f_max[i] = REP_particles[0].getFitnessValue(i);
    for(int j=0; j<REP_particles.size(); j++)
    {
      if (REP_particles[j].getFitnessValue(i) < f_min[i])
        f_min[i] = REP_particles[j].getFitnessValue(i);
      if (REP_particles[j].getFitnessValue(i) > f_max[i])
        f_max[i] = REP_particles[j].getFitnessValue(i);
    }
  }

  alphamax=0;
  alphamax_ind=0;
  for(int i=0; i<REP_particles.size(); i++)
  {
    alpha[i] = 0;
    for(int j=0; j<fitness_vals; j++)
    {
      if (j==0 && !REP_particles[i].activeLossesConsidered())
        continue;
      if (j==1 && !REP_particles[i].voltageStabilityIndexConsidered())
        continue;
      if (j==2 && !REP_particles[i].voltageDeviationConsidered())
        continue;

      int swnum = REP_particles[i].getSwitchNum();

      if (REP_particles[i].getFitnessValue(j) == f_min[j])
        alpha[i] += 1.0;
      else if (REP_particles[i].getFitnessValue(j) == f_max[j])
        continue;
      else
        alpha[i] += ( f_max[j] - REP_particles[i].getFitnessValue(j) ) / (f_max[j] - f_min[j]) ;
      if (alpha[i] > alphamax)
      {
        alphamax = alpha[i];
        alphamax_ind = i;
      }
    }
  }
  return REP_particles[alphamax_ind];
}

PSO_particle gbest_repository::circle_choice()
{
  double ref_point[fitness_vals];
  for (int i=0; i<fitness_vals; i++)
  {
    ref_point[i] = REP_particles[0].getFitnessValue(i);
    for(int j=0; j<REP_particles.size(); j++)
    {
      if( REP_particles[j].getFitnessValue(i) < ref_point[i] )
        ref_point[i] = REP_particles[j].getFitnessValue(i);
    }
  }

  int min_ind = 0;
  double minimum_radius_square;
  vector<double> base_point;
  REP_particles[0].getBaseFitnessValues(&base_point);

  for (int i=0; i<REP_particles.size(); ++i)
  {
    double radius_squere = 0;
    for (int j=0; j<fitness_vals; j++)
      radius_squere += pow((REP_particles[i].getFitnessValue(j) - ref_point[j]) / base_point[j],2);
    if (i == 0)
      minimum_radius_square = radius_squere;
    else
      if (radius_squere < minimum_radius_square)
      {
        minimum_radius_square = radius_squere;
        min_ind = i;
      }
  }
  return REP_particles[min_ind];
}

void gbest_repository::calculateMinMax()
{
  if (f_min != NULL)
    delete[] f_min;
  if (f_max != NULL)
    delete[] f_max;
  f_min = new double[fitness_vals];
  f_max = new double[fitness_vals];

  f_min[0] = REP_particles[0].getFitnessValue(0);
  f_min[1] = REP_particles[0].getFitnessValue(1);
  f_max[0] = REP_particles[0].getFitnessValue(0);
  f_max[1] = REP_particles[0].getFitnessValue(1);

  for(int i=0; i<REP_particles.size(); i++)
  {
      if (REP_particles[i].getFitnessValue(0) < f_min[0])
        f_min[0] = REP_particles[i].getFitnessValue(0);
      if (REP_particles[i].getFitnessValue(0) > f_max[0])
        f_max[0] = REP_particles[i].getFitnessValue(0);

      if (REP_particles[i].getFitnessValue(1) < f_min[1])
        f_min[1] = REP_particles[i].getFitnessValue(1);
      if (REP_particles[i].getFitnessValue(1) > f_max[1])
        f_max[1] = REP_particles[i].getFitnessValue(1);
  }
}

double gbest_repository::calculate2Dhypervolume(int funcs) // MSB - activeLosses, middle - VSI, LSB - vD
{
    vector<double> first_func;
    vector<double> second_func;

    vector<double> ref_point;
    REP_particles[0].getBaseFitnessValues(&ref_point);
    double first_ref;
    double second_ref;

    int first_ind;
    int second_ind;
    if (funcs == 3)
    {
        first_ind = 1;
        second_ind = 2;
    }
    else if (funcs == 5)
    {
        first_ind = 0;
        second_ind = 2;
    }
    else if (funcs == 6)
    {
        first_ind = 0;
        second_ind = 1;
    }
    else
    {
    //  cerr<<"Error: The switch value is not correct in hypervolume calculation"<<endl;
      return 0.0;
    //  exit(EXIT_SUCCESS);
    }
    first_ref = ref_point[first_ind];
    second_ref = ref_point[second_ind];
    for(int i=0; i<REP_particles.size(); ++i)
    {
      if (REP_particles[i].getFitnessValue(first_ind) < first_ref
        && REP_particles[i].getFitnessValue(second_ind) < second_ref)
      {
        first_func.push_back(REP_particles[i].getFitnessValue(first_ind));
        second_func.push_back(REP_particles[i].getFitnessValue(second_ind));
      }
    }
    for(int i=0; i<first_func.size(); i++)
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
    cout<<"Particles of the hypervolume area: "<<first_func.size()<<endl;
    double hypervolume;
    hypervolume = ((first_ref - first_func[0])/first_ref) * ((second_ref - second_func[0])/second_ref) ;
    for (int i=1; i<first_func.size(); ++i)
    {
      hypervolume += ((first_ref - first_func[i])/first_ref) * ((second_func[i-1] - second_func[i])/second_ref) ;
    }

    return hypervolume;
}

double gbest_repository::calculateSpacing(int switchnum)
{
  double first_func[REP_particles.size()];
  double second_func[REP_particles.size()];
  vector<double> ref_point;
  REP_particles[0].getBaseFitnessValues(&ref_point);
  double first_ref;
  double second_ref;

  if (switchnum == 3)
  {
    first_ref = ref_point[1];
    second_ref = ref_point[2];
    for(int i=0; i<REP_particles.size(); ++i)
    {
      first_func[i] = REP_particles[i].getFitnessValue(1);
      second_func[i] = REP_particles[i].getFitnessValue(2);
    }
  }
  else if (switchnum == 5)
  {
    first_ref = ref_point[0];
    second_ref = ref_point[2];
    for(int i=0; i<REP_particles.size(); ++i)
    {
      first_func[i] = REP_particles[i].getFitnessValue(0);
      second_func[i] = REP_particles[i].getFitnessValue(2);
    }
  }
  else if (switchnum == 6)
  {
    first_ref = ref_point[0];
    second_ref = ref_point[1];
    for(int i=0; i<REP_particles.size(); ++i)
    {
      first_func[i] = REP_particles[i].getFitnessValue(0);
      second_func[i] = REP_particles[i].getFitnessValue(1);
    }
  }
  else
  {
      return 0.0;
  }

  for(int i=0; i<REP_particles.size(); i++)
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
  double d[REP_particles.size()];
  double d_below[REP_particles.size()];

  d[0] = abs((first_func[1] - first_func[0])/first_ref) + abs((second_func[1] - second_func[0])/second_ref);
  d_below[0] = d[0];
  for (int i=1; i<REP_particles.size()-1; i++)
  {
    d_below[i] = abs((first_func[i+1]  - first_func[i])/first_ref) + abs((second_func[i+1] - second_func[i])/second_ref);
    d[i] = fmin(d_below[i], d_below[i-1]);
  }
  d[REP_particles.size() - 1] = d_below[REP_particles.size() - 2];

  double d_mean = 0;
  for(int i=0; i<REP_particles.size(); i++)
  {
    d_mean += d[i];
  }
  d_mean /= REP_particles.size();

  double devSum = 0;
  for(int i=0; i<REP_particles.size(); i++)
  {
    devSum += pow(d_mean - d[i], 2);
  }


  return sqrt((1.0/(REP_particles.size()-1)) * devSum);
}


void gbest_repository::log_repository_particles(const char* fname)
{
  fuzzy_choice();
  ofstream os(fname);
  os.close();
  os.open(fname, fstream::app);
  os<<"Repository contains "<<REP_particles.size()<<" particles"<<endl;
  os<<"Fuzzy logic (!!) best value: "<<alphamax<<endl;
  os<<"Fuzzy logic (!!) best choice: "<<alphamax_ind+1<<endl;
  REP_particles[alphamax_ind].log_particle_fitness_vals(fname);
  os<<endl<<endl;
  for(int i=0; i<fitness_vals; i++)
  {
    if (i==0)
      os<<"Active losses min: ";
    if (i==1)
      os<<"VSI min: ";
    if (i==2)
      os<<"Voltage deviations min:";
    os<<f_min[i]<<endl;
    if (i==0)
      os<<"Active losses max: ";
    if (i==1)
      os<<"VSI max: ";
    if (i==2)
      os<<"Voltage deviations max:";
    os<<f_max[i]<<endl;
  }
  os<<endl;
  for(int i=0; i<REP_particles.size(); i++)
  {
    os<<"REPOSITORY PARTICLE "<<i+1<<endl;
    REP_particles[i].log_particle_positions(fname);
    REP_particles[i].log_particle_fitness_vals(fname);
    os<<"hypercube: "<<hypercube_coordinates[i]<<endl;
    os<<endl<<endl;
  }
  os.close();
}

void gbest_repository::log_table_form(const char* fname)
{
  ofstream os(fname);

  PSO_particle my_best_particle = circle_choice();
  my_best_particle.log_table_form(os);

  os<<endl<<endl<<endl;

  /*
  int circle_best_wrange = fuzzy_choice_within_range();

  os<<"Fuzzy best particle within range active Losses: "<<REP_particles[fuzzy_best_wrange].getFitnessValue(0)*1e3<<endl;
  os<<"Fuzzy best particle within range VSI: "<<1/REP_particles[fuzzy_best_wrange].getFitnessValue(1)<<endl;
  */

  os<<endl<<endl;

  os<<"Active Losses range: "<<f_min[0]*1e3<<" - "<<f_max[0]*1e3<<endl;
  os<<"VSI range: "<<1/f_max[1]<<" - "<<1/f_min[1]<<endl;

  os<<endl;

  vector<double> base_point;
  REP_particles[0].getBaseFitnessValues(&base_point);
  double first_ref = base_point[0];
  double second_ref = base_point[1];

  double feasnum = 0;
  for(int i=0; i<REP_particles.size(); ++i)
  {
    if (REP_particles[i].getFitnessValue(0) < first_ref
      && REP_particles[i].getFitnessValue(1) < second_ref)
    {
      feasnum++;
    }
  }
  os<<"reference active losses: "<<first_ref*1e3<<endl;
  os<<"reference active losses: "<<1/second_ref<<endl;
  os<<endl;

  os<<"hypervolume: "<<calculate2Dhypervolume(REP_particles[0].getSwitchNum())<<endl;
  os<<"feasable solutions: "<<feasnum<<"/"<<REP_particles.size()<<endl;
  os<<"spacing: "<<calculateSpacing(REP_particles[0].getSwitchNum())<<endl;
  os<<endl;
  os<<"Load Modelling in Power Flow: "<<REP_particles[0].getPowerFlowType();

  os.close();
}

void gbest_repository::log_master_form(ofstream& os)
{
  PSO_particle my_best_particle = circle_choice();
  my_best_particle.log_master_form(os);
}

void gbest_repository::log_MinMax(ofstream& os, int i)
{
  os<<fixed;
  calculateMinMax();
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


void gbest_repository::createParetoOctaveGraph(const char* fname, int DGnum, int DGType)
{

  ofstream os(fname);
  os<<"p_front_mopso"<<DGnum<<"xt"<<DGType<<" = [..."<<endl;
  for(int i=0; i<REP_particles.size(); i++)
  {
    os<<'\t'<<REP_particles[i].getFitnessValue(0)*1e3<<'\t';
    os<<1/REP_particles[i].getFitnessValue(1);
    if (i != REP_particles.size()-1)
      os<<";"<<endl;
    else
      os<<"];"<<endl;
  }
  os<<endl<<endl;
  fuzzy_choice();
  os<<"fuzzy_best_particle_mopso = [";
  os<<REP_particles[alphamax_ind].getFitnessValue(0)*1e3<<'\t';
  os<<1/REP_particles[alphamax_ind].getFitnessValue(1)<<"];"<<endl;
  os.close();
}
