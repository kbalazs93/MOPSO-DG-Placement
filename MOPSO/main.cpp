#include <iostream>
#include <fstream>
#include "../NetworkModel/network.h"
#include "particle.h"
#include "particleswarm.h"
#include <ctime>
#include <unistd.h>

using namespace std;

void setParticle(PSO_particle* modified_particle)
{
	int capnum = 8;
	int windnum = 4;
	int solarnum = 4;
	int dispnum = 1;

	int capNodes[capnum];
	int windNodes[windnum];
	int solarNodes[solarnum];
	int dispNodes[dispnum];
	double dispPvalues[dispnum];
	double dispQvalues[dispnum];


	for(int i=0; i<capnum; i++)
		capNodes[i] = 22;
	for(int i=0; i<windnum; i++)
		windNodes[i] = 28;
	for(int i=0; i<solarnum; i++)
		solarNodes[i] = 28;
	for(int i=0; i<dispnum; i++)
	{
		dispNodes[i] = 19;
		dispPvalues[i] = 0.812e6;
		dispQvalues[i] = 0;
	}
	modified_particle->setParticle(capnum, capNodes, windnum, windNodes, solarnum, solarNodes, dispnum, dispNodes, dispPvalues, dispQvalues);
}

void setParticleTest(PSO_particle* modified_particle)
{
	int capnum = 0;
	int windnum = 0;
	int solarnum = 0;
	int dispnum = 2;
	int tapChangerPos;

	int capNodes[capnum];
	int windNodes[windnum];
	int solarNodes[solarnum];
	int dispNodes[dispnum];
	double dispPvalues[dispnum];
	double dispQvalues[dispnum];

	dispNodes[0] = 13; // 7
	dispPvalues[0] = 2e6;
	dispQvalues[0] = 0;
	dispNodes[1] = 24; // 12
	dispPvalues[1] = 0;
	dispQvalues[1] = 0;

	tapChangerPos = 0;

	modified_particle->setParticle(capnum, capNodes, windnum, windNodes, solarnum, solarNodes, dispnum, dispNodes, dispPvalues, dispQvalues, tapChangerPos);
}


int main(int argc, char **argv)
{

	/*{
	double Vmin = 0.95;
	double Vmax = 1.05;
	network netw("IEEE33.txt", Vmin, Vmax);
	netw.solvePowerFlowConstP("summer", 12);
	cout<<netw<<endl;
	netw.setTapChanger(12);
	netw.solvePowerFlowConstP("summer", 12);
	cout<<netw<<endl;
	cout<<"Active Losses: "<<netw.getActiveLosses()<<endl;
	cout<<"Losses: "<<netw.getLosses()<<endl;

	exit(EXIT_SUCCESS);
}//*/

	int NOP;
	int MAXITER;
	double Vmin;
	double Vmax;
	double C1;
	double C2;
	int NREP;
	int Ngrid;

	if (argc != 5 && argc != 7)
	{
		cout<<argc<<endl;
		cerr<<"Too few or too many arguments are used"<<endl;
		cerr<<"\t1st - scenarionum (1-Renewable, long, 2-Dispatchable, short)"<<endl;
		cerr<<"\t2nd - switchnum (MSB = activeLosses, middle = VSI, LSB = Voltage deviation)"<<endl;
		cerr<<"\t3rd - load_flow_type 1 - const power, 2 - const current"<<endl;
		cerr<<"\t4th - is the first branch a reg transformator? (0-no, 1-yes)"<<endl;
		cerr<<"Next arguments are only used, when scenarionum == 2"<<endl;
		cerr<<"\t5th - Num of DGs"<<endl;
		cerr<<"\t6th - Type of DGs (1 - only P, 2 - only Q, 3 P,Q)"<<endl;
		exit(EXIT_SUCCESS);
	}

	int scenarionum = atoi(argv[1]); // 1 - Renewable long interval scenario, 2 - Momentary just disp scenario
	if (scenarionum == 2 && argc != 7)
	{
		cerr<<"You have choosen the dispatchable scenario, but the two last arguments are missing"<<endl;
		cerr<<"\t5th - Num of DGs"<<endl;
		cerr<<"\t6th - Type of DGs (1 - only P, 2 - only Q, 3 P,Q)"<<endl;
		exit(EXIT_SUCCESS);
	}
	int switchnum = atoi(argv[2]); // MSB = activeLosses, middle = voltageStabilityIndex, LSB = voltageDeviation
	int load_flow_type = atoi(argv[3]); // 1 - const power, 2 - const current, 3 - mixed model

	int hasTapChanger = atoi(argv[4]);

	int sc2_dgNum;
	int sc2_dgType;
	if(scenarionum == 2)
	{
		sc2_dgNum = atoi(argv[5]);
		sc2_dgType = atoi(argv[6]);
	}

	double W = 0.5;
	double Wdamp = 0.99;


	if (scenarionum == 1)
	{
		NOP = 500;
	  MAXITER = 200;
	  Vmin = 0.95;
	  Vmax = 1.05;
		C1 = 2.0;
		C2 = 2.0;
		NREP = 100;
		Ngrid = 7;
	}
	else if (scenarionum == 2)
	{
		NOP = 200;
	  MAXITER = 200;
	  Vmin = 0.95;
	  Vmax = 1.05;
		C1 = 1.0;
		C2 = 2.0;
		NREP = 50;
		Ngrid = 7;
	}


	network netw("./NetworkModel/IEEE33.txt", Vmin, Vmax, load_flow_type);

	particle_swarm ps(&netw, NOP, MAXITER, C1, C2, W, Wdamp, NREP, Ngrid, hasTapChanger);


	/*
	{
		PSO_particle test_particle(&netw, 0, 2, 0, 0, 3, 6, 2, C1, C2, W, Wdamp, true);
		setParticleTest (&test_particle);
		test_particle.evaluate_fitness_values();
		cout<<netw<<endl;
		cout<<"Losses: "<<netw.getLosses()<<endl;
		cout<<"VSI: "<<netw.getVoltageStabilityIndexSum()<<endl;
		netw.powerFlowControlAlgorithm();

		ofstream teststream;
		teststream.open("test.txt");
		test_particle.log_table_form(teststream);
		teststream.close();
		exit(EXIT_SUCCESS);
	}//*/

	if (scenarionum == 1)
	{
		netw.createOctaveGraph(string("./output/MOPSO/Renewable_Scenario_Graphs/empty_summer"), string("summer"), &network::getActiveLosses, true, string("_Ploss"));
		netw.createOctaveGraph(string("./output/MOPSO/Renewable_Scenario_Graphs/empty_summer"), string("summer"), &network::getMinimumVSI, false, string("_VSI"));
		netw.createOctaveGraph(string("./output/MOPSO/Renewable_Scenario_Graphs/empty_summer"), string("summer"), &network::getMaximumVoltageDeviation, false, string("_VD"));

		netw.createOctaveGraph(string("./output/MOPSO/Renewable_Scenario_Graphs/empty_autumn"), string("autumn"), &network::getActiveLosses, true, string("_Ploss"));
		netw.createOctaveGraph(string("./output/MOPSO/Renewable_Scenario_Graphs/empty_autumn"), string("autumn"), &network::getMinimumVSI, false, string("_VSI"));
		netw.createOctaveGraph(string("./output/MOPSO/Renewable_Scenario_Graphs/empty_autumn"), string("autumn"), &network::getMaximumVoltageDeviation, false, string("_VD"));

		netw.createOctaveGraph(string("./output/MOPSO/Renewable_Scenario_Graphs/empty_winter"), string("winter"), &network::getActiveLosses, true, string("_Ploss"));
		netw.createOctaveGraph(string("./output/MOPSO/Renewable_Scenario_Graphs/empty_winter"), string("winter"), &network::getMinimumVSI, false, string("_VSI"));
		netw.createOctaveGraph(string("./output/MOPSO/Renewable_Scenario_Graphs/empty_winter"), string("winter"), &network::getMaximumVoltageDeviation, false, string("_VD"));

		netw.createOctaveGraph(string("./output/MOPSO/Renewable_Scenario_Graphs/empty_spring"), string("spring"), &network::getActiveLosses, true, string("_Ploss"));
		netw.createOctaveGraph(string("./output/MOPSO/Renewable_Scenario_Graphs/empty_spring"), string("spring"), &network::getMinimumVSI, false, string("_VSI"));
		netw.createOctaveGraph(string("./output/MOPSO/Renewable_Scenario_Graphs/empty_spring"), string("spring"), &network::getMaximumVoltageDeviation, false, string("_VD"));

		PSO_particle base_particle(&netw, 8, 1, 4, 4, 3, 7, 1);
		setParticle (&base_particle);

		base_particle.evaluate_fitness_values();


		//base_particle.log_particle_fitness_vals("test.txt");
		base_particle.createOctaveGraph("./output/MOPSO/Renewable_Scenario_Graphs/base_summer", "summer");
		base_particle.createOctaveGraph("./output/MOPSO/Renewable_Scenario_Graphs/base_autumn", "autumn");
		base_particle.createOctaveGraph("./output/MOPSO/Renewable_Scenario_Graphs/base_winter", "winter");
		base_particle.createOctaveGraph("./output/MOPSO/Renewable_Scenario_Graphs/base_spring", "spring");

	}
	/*
	else if (scenarionum == 2)
	{
		if (load_flow_type == 1)
			netw.solvePowerFlow("summer", 12);
		else if  (load_flow_type == 2)
			netw.solvePowerFlowConstI("summer", 12);


		cout<<"Empty network, active losses: "<<netw.getActiveLosses()*1e3<<" kW"<<endl;
		cout<<"Empty network, VSI: "<<netw.getVoltageStabilityIndexSum()<<endl;
		cout<<"Empty network, voltage deviation: "<<netw.getVoltageDeviationSum()<<endl;
	}*/

	if(scenarionum == 1)
		ps.initialize_scenario1(switchnum, NOP, MAXITER, C1, C2, W, Wdamp, hasTapChanger);
	else if(scenarionum == 2)
		ps.initialize_scenario2(sc2_dgNum, switchnum, sc2_dgType, NOP, MAXITER, C1, C2, W, Wdamp, hasTapChanger);

	time_t t1 = time(0);
	ps.particle_swarm_optimalization();
	PSO_particle my_best_particle = ps.circle_choice();
	my_best_particle.log_particle_positions("my_best_particle.txt");
	my_best_particle.evaluate_fitness_values();
//	my_best_particle.log_particle_fitness_vals("./output/MOPSO/test.txt");

	time_t t2 = time(0);
	int secs = t2 - t1;
	int mins = secs / 60;
	secs %= 60;
	int hours = mins / 60;
	mins %= 60;
	cout<<endl<<"Time needed: "<<hours<<" h  "<<mins<<" m  "<<secs<<" s"<<" ("<<mins*60+secs<<" s)"<<endl;
	cout<<"hypervolume: "<<ps.calculate2Dhypervolume(switchnum)<<endl;
	cout<<"spacing: "<<ps.calculateSpacing(switchnum)<<endl;

	if (scenarionum == 1)
	{
		my_best_particle.createOctaveGraph("./output/MOPSO/Renewable_Scenario_Graphs/my_summer", "summer");
		my_best_particle.createOctaveGraph("./output/MOPSO/Renewable_Scenario_Graphs/my_autumn", "autumn");
		my_best_particle.createOctaveGraph("./output/MOPSO/Renewable_Scenario_Graphs/my_winter", "winter");
		my_best_particle.createOctaveGraph("./output/MOPSO/Renewable_Scenario_Graphs/my_spring", "spring");
	}
	else if (scenarionum == 2 )
	{
		string table_string = "./output/MOPSO/table";
		string ind_num_string = "./output/MOPSO/ind_num";
		string activeLosses_range_output_string = "./output/MOPSO/activeLosses_range";
		string VSI_range_output_string = "./output/MOPSO/VSI_range";
		string hypervolume_output_string = "./output/MOPSO/hypervolume";
		string spacing_output_string = "./output/MOPSO/spacing";
		string time_output_string = "./output/MOPSO/times";

		if( load_flow_type == 1 )
		{
			table_string += "_constP";
			ind_num_string += "_constP";
			activeLosses_range_output_string += "_constP";
			VSI_range_output_string += "_constP";
			hypervolume_output_string += "_constP";
			spacing_output_string += "_constP";
			time_output_string += "_constP";
		}
		else if (load_flow_type == 2)
		{
			table_string += "_constI";
			ind_num_string += "_constI";
			activeLosses_range_output_string += "_constI";
			VSI_range_output_string += "_constI";
			hypervolume_output_string += "_constI";
			spacing_output_string += "_constI";
			time_output_string += "_constI";
		}

		if (hasTapChanger == 0)
		{
			table_string += "_noRTR_";
			ind_num_string += "_noRTR_";
			activeLosses_range_output_string += "_noRTR_";
			VSI_range_output_string += "_noRTR_";
			hypervolume_output_string += "_noRTR_";
			spacing_output_string += "_noRTR_";
			time_output_string += "_noRTR_";
		}
		else
		{
			table_string += "_wRTR_";
			ind_num_string += "_wRTR_";
			activeLosses_range_output_string += "_wRTR_";
			VSI_range_output_string += "_wRTR_";
			hypervolume_output_string += "_wRTR_";
			spacing_output_string += "_wRTR_";
			time_output_string += "_wRTR_";
		}

		table_string += "switchnum";
		ind_num_string += "switchnum";
		activeLosses_range_output_string += "switchnum";
		VSI_range_output_string += "switchnum";
		hypervolume_output_string += "switchnum";
		spacing_output_string += "switchnum";
		time_output_string += "switchnum";

		table_string += argv[2];
		ind_num_string += argv[2];
		activeLosses_range_output_string += argv[2];
		VSI_range_output_string += argv[2];
		hypervolume_output_string += argv[2];
		spacing_output_string += argv[2];
		time_output_string += argv[2];

		table_string += ".out";
		ind_num_string += ".out";
		activeLosses_range_output_string += ".out";
		VSI_range_output_string += ".out";
		hypervolume_output_string += ".out";
		spacing_output_string += ".out";
		time_output_string += ".out";


		if (switchnum != 1 && switchnum != 2 && switchnum != 4)
		{
			string pareto_string = "./output/MOPSO/pareto";
			if( load_flow_type == 1 )
				pareto_string += "_constP";
			else if (load_flow_type == 2)
				pareto_string += "_constI";

			if (hasTapChanger == 0)
				pareto_string += "_noRTR_";
			else
				pareto_string += "_wRTR_";

			pareto_string += argv[5];
			pareto_string += "DG_Type";
			pareto_string += argv[6];
			pareto_string += ".out";

			ofstream master_table_os(table_string, fstream::app);
			ofstream ind_num_os(ind_num_string, fstream::app);
			ofstream activeLosses_range_os(activeLosses_range_output_string, fstream::app);
			ofstream VSI_range_os(VSI_range_output_string, fstream::app);
			ofstream hypervolume_os(hypervolume_output_string, fstream::app);
			ofstream spacing_os(spacing_output_string, fstream::app);

			ps.log_master_form(master_table_os);
			ps.log_rep_size(ind_num_os);
			ps.log_MinMax(activeLosses_range_os, 0);
			ps.log_MinMax(VSI_range_os, 1);
			ps.log_2DHypervolume(hypervolume_os, switchnum);
			ps.log_Spacing(spacing_os, switchnum);
			ps.createParetoOctaveGraph(pareto_string.c_str());
		}
		ofstream time_os(time_output_string, fstream::app);
		time_os<<hours*3600 + mins*60 + secs<<endl;
	
		ps.log_table_form("./output/MOPSO/original/results_table_form.txt");
	}

	cout<<"Program ends sucesfully"<<endl;
	return 0;
}
