# include "rand.h"
# include "nsga2.h"
#include <iostream>
#include <ctime>

using namespace std;

nsga2_optimization* setup_optimization(int argc, char **argv);
void run_DGPlacement(int argc, char** argv);

int main(int argc, char **argv)
{
	if (argc != 5 && argc != 6)
	{
		cerr<<"Too few or too many arguments are used"<<endl;
		cerr<<"\t1st - load_flow_type 1 - const power, 2 - const current"<<endl;
		cerr<<"\t2nd - is the first branch a reg transformator? (0-no, 1-yes)"<<endl;
		cerr<<"\t3rd - Num of DGs"<<endl;
		cerr<<"\t4th - Type of DGs (1 - only P, 2 - only Q, 3 P,Q)"<<endl;
		cerr<<"\t5th - [Optional] if you want to use gnu plot (0-no, 1-yes)"<<endl;
		exit(EXIT_SUCCESS);
	}


	run_DGPlacement(argc, argv);
	return 0;
}

void run_DGPlacement(int argc, char** argv)
{
	FILE* final_pop;
	FILE* best_pop;
	FILE* pareto;
	FILE* table_output;

	ofstream master_output;
	ofstream feasnum_output;
	ofstream activeLosses_range_output;
	ofstream VSI_range_output;
	ofstream hypervolume_output;
	ofstream spacing_output;
	ofstream time_output;

	int load_flow_type = atoi(argv[1]); // 1 - const power, 2 - const current, 3 - mixed model
	int hasTapChanger = atoi(argv[2]);

	string final_pop_string = "./output/NSGA-II/original/final_pop";
	string best_pop_string = "./output/NSGA-II/original/best_pop";
	string pareto_string = "./output/NSGA-II/pareto";
	string table_output_string = "./output/NSGA-II/original/table_output";
	string master_output_string = "./output/NSGA-II/table";
	string feasnum_output_string = "./output/NSGA-II/feasnum";
	string activeLosses_range_output_string = "./output/NSGA-II/activeLosses_range";
	string VSI_range_output_string = "./output/NSGA-II/VSI_range";
	string hypervolume_output_string = "./output/NSGA-II/hypervolume";
	string spacing_output_string = "./output/NSGA-II/spacing";
	string time_output_string = "./output/NSGA-II/time";


	if( load_flow_type == 1 )
	{
		final_pop_string += "_constP";
		best_pop_string += "_constP";
		pareto_string += "_constP";
		table_output_string += "_constP";
		master_output_string += "_constP";
		feasnum_output_string += "_constP";
		activeLosses_range_output_string += "_constP";
		VSI_range_output_string += "_constP";
		hypervolume_output_string += "_constP";
		spacing_output_string += "_constP";
		time_output_string += "_constP";
	}
	else if (load_flow_type == 2)
	{
		final_pop_string += "_constI";
		best_pop_string += "_constI";
		pareto_string += "_constI";
		table_output_string += "_constI";
		master_output_string += "_constI";
		feasnum_output_string += "_constI";
		activeLosses_range_output_string += "_constI";
		VSI_range_output_string += "_constI";
		hypervolume_output_string += "_constI";
		spacing_output_string += "_constI";
		time_output_string += "_constI";
	}

	if (hasTapChanger == 0)
	{
		final_pop_string += "_noRTR_";
		best_pop_string += "_noRTR_";
		pareto_string += "_noRTR_";
		table_output_string += "_noRTR_";
		master_output_string += "_noRTR_";
		feasnum_output_string += "_noRTR_";
		activeLosses_range_output_string += "_noRTR_";
		VSI_range_output_string += "_noRTR_";
		hypervolume_output_string += "_noRTR_";
		spacing_output_string += "_noRTR_";
		time_output_string += "_noRTR_";
	}
	else
	{
		final_pop_string += "_wRTR_";
		best_pop_string += "_wRTR_";
		pareto_string += "_wRTR_";
		table_output_string += "_wRTR_";
		master_output_string += "_wRTR_";
		feasnum_output_string += "_wRTR_";
		activeLosses_range_output_string += "_wRTR_";
		VSI_range_output_string += "_wRTR_";
		hypervolume_output_string += "_wRTR_";
		spacing_output_string += "_wRTR_";
		time_output_string += "_wRTR_";

	}

	final_pop_string += argv[3];
	best_pop_string += argv[3];
	pareto_string += argv[3];
	table_output_string += argv[3];
	time_output_string += argv[3];

	final_pop_string += "DG_Type";
	best_pop_string += "DG_Type";
	pareto_string += "DG_Type";
	table_output_string += "DG_Type";
	time_output_string += "DG_Type";

	final_pop_string += argv[4];
	best_pop_string += argv[4];
	pareto_string += argv[4];
	table_output_string += argv[4];
	time_output_string += argv[4];

	final_pop_string += ".out";
	best_pop_string += ".out";
	pareto_string += ".out";
	table_output_string += ".out";
	master_output_string += ".out";
	feasnum_output_string += ".out";
	activeLosses_range_output_string += ".out";
	VSI_range_output_string += ".out";
	hypervolume_output_string += ".out";
	spacing_output_string += ".out";
	time_output_string += ".out";


	final_pop = fopen(final_pop_string.c_str(),"w");
	best_pop = fopen(best_pop_string.c_str(),"w");
	pareto = fopen(pareto_string.c_str(), "w");
	table_output = fopen(table_output_string.c_str(), "w");


	master_output.open(master_output_string, fstream::app);
	feasnum_output.open(feasnum_output_string, fstream::app);
	activeLosses_range_output.open(activeLosses_range_output_string, fstream::app);
	VSI_range_output.open(VSI_range_output_string, fstream::app);
	hypervolume_output.open(hypervolume_output_string, fstream::app);
	spacing_output.open(spacing_output_string, fstream::app);
	time_output.open(time_output_string, fstream::app);

	fprintf(final_pop,"# This file contains the data of final population\n");
	fprintf(best_pop,"# This file contains the data of final feasible population (if found)\n");


	nsga2_optimization *nsgaDG = setup_optimization(argc, argv);

	void *inp = NULL;
	void *out = NULL;

	time_t t1 = time(0);
	nsgaDG->NSGA2(inp, out);
	time_t t2 = time(0);

	cout<<time_output_string<<endl;
	time_output<<t2-t1<<endl;

	nsgaDG->report_pop(final_pop);
	nsgaDG->report_feasible(best_pop);
	nsgaDG->log_pareto_size(feasnum_output);
	nsgaDG->log_MinMax(activeLosses_range_output, 0);
	nsgaDG->log_MinMax(VSI_range_output, 1);
	nsgaDG->log_2Dhypervolume(hypervolume_output);
	nsgaDG->log_Spacing(spacing_output);


	nsgaDG->createParetoOctaveGraph(pareto);


//	nsgaDG->log_table_form(table_output);

	nsgaDG->log_master_form(master_output);


	fflush(final_pop);
	fflush(best_pop);
	fflush(pareto);
	fflush(table_output);

	fclose(final_pop);
	fclose(best_pop);
	fclose(pareto);
	fclose(table_output);

	master_output.close();
	feasnum_output.close();
	time_output.close();

	delete nsgaDG;
}

nsga2_optimization* setup_optimization(int argc, char** argv)
{



	int load_flow_type = atoi(argv[1]); // 1 - const power, 2 - const current, 3 - mixed model
	int hasTapChanger = atoi(argv[2]);

	int	sc2_dgNum = atoi(argv[3]);
	int	sc2_dgType = atoi(argv[4]);

	int gnu_plot = 0;
	if (argc == 6)
		gnu_plot = atoi(argv[5]);

	int scenarionum = 2;
	int switchnum = 6; // MSB = activeLosses, middle = voltageStabilityIndex, LSB = voltageDeviation



	int params_num = 4;
	int sc2_params[params_num];
	sc2_params[0] = sc2_dgNum; // Num o DGs
	sc2_params[1] = sc2_dgType; // Type o DGs
	sc2_params[2] = load_flow_type;
	sc2_params[3] = hasTapChanger;
	//sc2_params[3] = 1; // Contains tap changer (0 - no, everything else - yes)


	double node_min = 1;
	double node_max = 32;
	double injection_min = 0;
	double injection_max = 2e6;

	int nreal;

	if (sc2_dgType == 1 || sc2_dgType == 2)
		nreal = 2*sc2_dgNum+1;
	else if (sc2_dgType == 3)
		nreal = 3*sc2_dgNum+1;
	else
	{
		printf("Error: Wrong DG type\n");
		exit(EXIT_SUCCESS);
	}

	double min_realvar[nreal];
	double max_realvar[nreal];

	for(int i=0; i<nreal; i++)
	{
		if (sc2_dgType == 1 || sc2_dgType == 2)
		{
			if (i%2 == 0)
			{
				min_realvar[i] = node_min;
				max_realvar[i] = node_max;
			}
			else
			{
				max_realvar[nreal-1] = 0.1;
				min_realvar[i] = injection_min;
				max_realvar[i] = injection_max;
			}
		}
		else if (sc2_dgType == 3)
		{
			if (i%3 == 0)
			{
				min_realvar[i] = node_min;
				max_realvar[i] = node_max;
			}
			else
			{
				min_realvar[i] = injection_min;
				max_realvar[i] = injection_max;
			}
		}

	}
	if(hasTapChanger)
	{
		min_realvar[nreal-1] = -12;
		max_realvar[nreal-1] = 12;
	}
	else
	{
		min_realvar[nreal-1] = -0.1;
		max_realvar[nreal-1] = 0.1;
	}

	double seed = 0.5;
	int popsize = 200;
	int ngen = 100;
	int nobj = 2;
 	int ncon = 2;
	double pcross_real = 0.1;
	double pmut_real = 0.9;
	double dist_ind_cross_real = 20;
	double dist_ind_mut_real = 50;

	int nbin = 0;
	int* nbits = NULL;
	double* min_binvar = NULL;
	double* max_binvar = NULL;
	double pcross_bin = 0;
	double pmut_bin = 0;

	nsga2_optimization* nsgaDG = new nsga2_optimization(seed, popsize, ngen, nobj, ncon, nreal,
		min_realvar, max_realvar, pcross_real, pmut_real,
		dist_ind_cross_real, dist_ind_mut_real, nbin, nbits,
		min_binvar, max_binvar, pcross_bin, pmut_bin,
		NULL, 0, sc2_params, params_num, gnu_plot);

	return nsgaDG;
}
