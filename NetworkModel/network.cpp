#include "network.h"
#include <fstream>
#include <iostream>
#include <cstring>
#include <cstdlib>
#include <cmath>
#include <sstream>

using namespace std;

// CONSTRUCTORS
network::network()
{
	numONodes = 2;
	numOBranches = 1;
	nodes = new node[numONodes];
	branches = new branch[numOBranches];

	powerFlowType = "";

	tapChanger = 0;
	tapChangerDelta = 1.25 / 100;

	node firstNode;
	branch Branch;
	node secondNode;

	nodes[0] = firstNode;
	nodes[1] = secondNode;
	branches[0] = Branch;

	Sb = 1.0;
	double Vb1 = 1.0;

	Vb.push_back(Vb1);
	Vb.push_back(Vb1);

	filename = NULL;
}

network::network(network& n)
{
	delete[] branches;
	delete[] nodes;

	this->Vmin = n.Vmin;
	this->Vmax = n.Vmax;

	tapChanger = n.tapChanger;
	tapChangerMin = n.tapChangerMin;
	tapChangerMax = n.tapChangerMax;
	tapChangerDelta = n.tapChangerDelta;
	firstBranchImpedance = n.firstBranchImpedance;

	numOBranches = n.numOBranches;
	numONodes = n.numONodes;

	powerFlowType = n.powerFlowType;

	branches = new branch[n.numOBranches];
	nodes = new node[n.numONodes];

	for(int i=0; i<numOBranches; i++)
		branches[i] = n.branches[i];
	for(int i=0; i<numONodes; i++)
		nodes[i] = n.nodes[i];

	Sb = n.Sb;
	Vb = n.Vb;

	emptyLossesSummer12 = n.emptyLossesSummer12;
	emptyVSISummer12 = n.emptyVSISummer12;

	cout<<capacitors.size()<<endl;

	filename = new char[strlen(n.filename)];
	strcpy (filename, n.filename);
}

network::~network()
{
	delete[] branches;
	delete[] nodes;
	Vb.clear();
	capacitors.clear();
	wind_generators.clear();
	solar_generators.clear();
	disp_generators.clear();
}

network::network(const char* fname, double Vmin, double Vmax, int powerFlowType) : inputstream(fname, ifstream::in)
{
	cout<<"Initializing..."<<endl;

	this->Vmin = Vmin;
	this->Vmax = Vmax;

	tapChanger = 0;
	tapChangerMin = -12;
	tapChangerMax = 12;
	tapChangerDelta = 1.25 / 100;

	if (powerFlowType == 1)
		this->powerFlowType = "Constant Power";
	else if (powerFlowType == 2)
		this->powerFlowType = "Constant Current";
	else if (powerFlowType == 3)
		this->powerFlowType = "Mixed Model";

	filename = new char[strlen(fname)];
	strcpy(filename, fname);

	char ch;

	word = new char[50] ;
	inputstream.getline(word, 50);
	bool got_Base = false ;
	bool got_Element = false;


	numOBranches = 0;
	while( (ch = inputstream.get()) != -1 )
	{
		inputstream.putback(ch);
		if (!strncmp(word, "TRxx", 2 ) || !strncmp(word, "VODxx", 3 ) )
			numOBranches++;
		inputstream.getline(word, 50);
	}
	numONodes = numOBranches + 1;
	inputstream.close();
	nodes = new node[numONodes];
	branches = new branch[numOBranches];

	node firstNode;
	nodes[0] = firstNode;
	branchPointer = 0;
	nodePointer = 1;

	inputstream.open(fname, fstream::in);
	while( (ch = inputstream.get()) != -1 )
	{
		inputstream.putback(ch);
		if ( !strncmp(word, "BASExx", 4 ) )
		{
			readBase();
			got_Base = true;
		}
		else if (!strncmp(word, "TRxx", 2 ) )
		{
			readTR();
			got_Element = true;
		}
		else if ( !strncmp(word, "VODxx", 3 ) )
		{
			readSection();
			got_Element = true;
		}
		else if ( !strncmp(word, "NODExx", 4 ) )
		{
			readNode();
		}
		else
		{
			inputstream.getline(word, 50);
		}
	}
	if ( !got_Base )
		error("Error: There are no base values added. Program stops...");


	bool consistent = true; // Consistent means, the node index in class is equal with the vector index


	for(int i=0; i<numONodes; i++ )
		consistent *= (nodes[i].getNodeNum() == i);
	for(int i=0; i<numOBranches; i++ )
		consistent *= (branches[i].getEndPoint() == i+1);

	if (consistent)
		cout<<"\tNodes are consistently initialized"<<endl;
	else error("Error: The network is not consistent. Program stops..");

	bool numedByLayers = true;
	int lastIndex = 0;
	for(int i=0; i<numOBranches; i++ )
	{
		 numedByLayers *= (branches[i].getBeginPoint() >= lastIndex);
		 lastIndex = branches[i].getBeginPoint() ;
	}
	if (numedByLayers)
		cout<<"\tNodes are numed by layers"<<endl;
	else error("Error: Nodes are not numed by layers. Program stops..");
	cout<<endl;
	inputstream.close();

	string diagram_name;
	for(int i=0; i<3; i++)
	{
		switch (i) {
			case 0:
				diagram_name = "./NetworkModel/LoadDiagram.m";
				break;
			case 1:
				diagram_name = "./NetworkModel/windPowerDiagram.m";
				break;
			case 2:
				diagram_name = "./NetworkModel/solarPowerDiagram.m";
				break;
		}

		inputstream.open(diagram_name.c_str(), fstream::in);
		while( (ch = inputstream.get())  != -1  )
		{
			inputstream.putback(ch);
			while ( !isalpha(ch = inputstream.get() ) && ch != -1 && !inputstream.eof() )
				;
			inputstream.putback(ch);
			inputstream>>word;
			double* arrey_pointer;

			if ( i == 0 && !strncmp(word, "summer", 6 ) )
			{
				arrey_pointer = load_summer;
			}
			else if (i == 0 && !strncmp(word, "autumn", 6 ) )
			{
				arrey_pointer = load_autumn;
			}
			else if (i == 0 && !strncmp(word, "winter", 6 ) )
			{
				arrey_pointer = load_winter;
			}
			else if (i == 0 && !strncmp(word, "spring", 6 ) )
			{
				arrey_pointer = load_spring;
			}

			else if ( i == 1 && !strncmp(word, "summer", 6 ) )
			{
				arrey_pointer = wind_summer;
			}
			else if (i == 1 && !strncmp(word, "autumn", 6 ) )
			{
				arrey_pointer = wind_autumn;
			}
			else if (i == 1 && !strncmp(word, "winter", 6 ) )
			{
				arrey_pointer = wind_winter;
			}
			else if (i == 1 && !strncmp(word, "spring", 6 ) )
			{
				arrey_pointer = wind_spring;
			}

			else if ( i == 2 && !strncmp(word, "summer", 6 ) )
			{
				arrey_pointer = solar_summer;
			}
			else if (i == 2 && !strncmp(word, "autumn", 6 ) )
			{
				arrey_pointer = solar_autumn;
			}
			else if (i == 2 && !strncmp(word, "winter", 6 ) )
			{
				arrey_pointer = solar_winter;
			}
			else if (i == 2 && !strncmp(word, "spring", 6 ) )
			{
				arrey_pointer = solar_spring;
			}
			else
			{
				continue;
			}
			for (int j=0; ch != ']'; j++)
			{
				while ( !isdigit(ch = inputstream.get() ) && ch != -1 && ch != ']' && !inputstream.eof() )
					;
				if (ch == -1 || ch == ']' || inputstream.eof())
				{
					break;
				}
				inputstream.putback(ch);
				double val;
				inputstream>>val;
				if( j % 2 == 1)
					*(arrey_pointer + j/2) = val;
			}
		}
		inputstream.close();

		solvePowerFlow("summer", 12);
		emptyLossesSummer12 = getActiveLosses();
		emptyVSISummer12 = getVoltageStabilityIndexSum();
	}
	firstBranchImpedance = branches[0].getImpedance();
}

//*** READ-HELPER FUNCTIONS
void network::readBase()
{
	char ch;
	int num;

	double Vb;

	bool got_Sb = false;
	bool got_Vb = false;

	while (true)
	{
		while ( !isalpha(ch = inputstream.get() ) && ch != -1 && !inputstream.eof() )
			;
		inputstream.putback(ch);
		inputstream>>word;

		double* pointer;
		bool* got_flag;

		if ( ch != -1 && !inputstream.eof() && !strcmp(word, "Sb") )
		{
			pointer = &Sb;
			got_flag = &got_Sb;
		}
		else if ( ch != -1 && !inputstream.eof() && !strncmp(word, "Vbx", 2) )
		{
			pointer = &Vb;
			got_flag = &got_Vb;
		}
		else if ( ch == -1 || inputstream.eof() || !strcmp(word, "BASE") || !strncmp(word, "VODxx", 3) || !strncmp(word, "TRxx", 2) || !strncmp(word, "NODExx", 4) )
		{
			break;
		}
		else
		{
			continue;
		}

		if ( !strcmp(word, "Sb") || !strncmp(word, "Vbx", 2))
		{
			*got_flag = true;
			while ( !isdigit(ch = inputstream.get() )  && ch != -1 && !inputstream.eof())
				;
			inputstream.putback(ch);
			inputstream>>*pointer;
			while ( isspace(ch = inputstream.get() ) &&  ch != -1 && !inputstream.eof() )
				;
			inputstream.putback(ch);
			inputstream>>word;
			if (!strncmp(word,"kxx",1) )
				*pointer *= 1e3;
			else if (!strncmp(word,"Mxx",1) )
				*pointer *= 1e6;
			if (pointer == &Vb)
				this->Vb.push_back(Vb);
		}
	}
	if (!got_Sb || !got_Vb)
		error("Error: Some base parameters are missing");
}

void network::readTR()
{
	char ch;

	double beginpoint;
	double Sn;
	double Vn1;
	double uk;
	double Pcu;
	double X;
	double R;
	double tmax;
	double du;

	bool got_begpoint = false;
	bool got_Sn = false;
	bool got_Vn1 = false;
	bool got_uk = false;
	bool got_Pcu = false;
	bool got_tmax = false;
	bool got_du = false;

	while (true)
	{
		while ( !isalpha(ch = inputstream.get() )  && ch != -1 && !inputstream.eof() )
			;
		inputstream.putback(ch);
		inputstream>>word;

		double* pointer;
		bool* got_flag;

		if ( ch != -1 && !inputstream.eof() && !strcmp(word, "beginpoint") )
		{
			pointer = &beginpoint;
			got_flag = &got_begpoint;
		}
		else if ( ch != -1 && !inputstream.eof() && !strcmp(word, "Sn") )
		{
			pointer = &Sn;
			got_flag = &got_Sn;
		}
		else if ( ch != -1 && !inputstream.eof() && !strcmp(word, "Vn1") )
		{
			pointer = &Vn1;
			got_flag = &got_Vn1;
		}

		else if (  ch != -1 && !inputstream.eof() && !strcmp(word, "uk") )
		{
			pointer = &uk;
			got_flag = &got_uk;
		}
		else if (  ch != -1 && !inputstream.eof() && !strcmp(word, "Pcu") )
		{
			pointer = &Pcu;
			got_flag = &got_Pcu;
		}
		else if (  ch != -1 && !inputstream.eof() && !strcmp(word, "tmax") )
		{
			pointer = &tmax;
			got_flag = &got_tmax;
		}
		else if (  ch != -1 && !inputstream.eof() && !strcmp(word, "du") )
		{
			pointer = &du;
			got_flag = &got_du;
		}
		else if (  ch == -1 && inputstream.eof() || !strcmp(word, "BASE") || !strncmp(word, "VODxx", 3) || !strncmp(word, "TRxx", 2) || !strncmp(word, "NODExx", 4) )
		{
			break;
		}
		else
			continue;

		if ( *got_flag )
			error("Error: Transformer parameter cannot be initialized twice. Program stops..");
		*got_flag = true;
		while ( !isdigit(ch = inputstream.get() ) &&  ch != -1 && !inputstream.eof() )
			;
		inputstream.putback(ch);
		inputstream>>*pointer;
	}

	if ( !got_begpoint || !got_Sn || !got_Vn1 || !got_uk )
		error("Error: Transformer: a vital parameter is missing");
	if ( (got_tmax && !got_du) || (!got_tmax && got_du) )
		error("Error: some transformer regulation parameter is missing" );

	X = uk*pow(Vn1,2) / (100*Sn);
	R = 0;
	if (got_Pcu)
	{
		R = Pcu*pow(Vn1,2) / (3*pow(Sn,2)) ;
		X = sqrt( pow(X,2) - pow(R,2) );
	}
	if (got_tmax && got_du)
	{
		branch trbranch(beginpoint, R, X, tmax, du);
		branches[branchPointer++] = trbranch;
	}
	else
	{
		branch trbranch(beginpoint, R, X);
		branches[branchPointer++] = trbranch;
	}
	node endNode;
	nodes[nodePointer++] = endNode;
}

void network::readSection()
{
	char ch;

	double beginpoint;
	double endpoint;
	double r=0;
	double x=0;
	double g=0;
	double b=0;
	double l=1;

	bool got_begpoint = false;
	bool got_endpoint = false;
	bool got_r = false;
	bool got_x = false;
	bool got_g = false;
	bool got_b = false;
	bool got_l = false;

	while (true)
	{
		while ( !isalpha(ch = inputstream.get() ) &&  ch != -1 && !inputstream.eof() )
			;
		inputstream.putback(ch);
		inputstream>>word;

		double* pointer;
		bool* got_flag;

		if (  ch != -1 && !inputstream.eof() && !strcmp(word, "beginpoint") )
		{
			pointer = &beginpoint;
			got_flag = &got_begpoint;
		}
		else if (  ch != -1 && !inputstream.eof() && !strcmp(word, "endpoint") )
		{
			pointer = &endpoint;
			got_flag = &got_endpoint;
		}
		else if (  ch != -1 && !inputstream.eof() && !strcmp(word, "x") )
		{
			pointer = &x;
			got_flag = &got_x;
		}
		else if (  ch != -1 && !inputstream.eof() && !strcmp(word, "r") )
		{
			pointer = &r;
			got_flag = &got_r;
		}
		else if (  ch != -1 && !inputstream.eof() && !strcmp(word, "b") )
		{
			pointer = &b;
			got_flag = &got_b;
		}
		else if (  ch != -1 && !inputstream.eof() && !strcmp(word, "g") )
		{
			pointer = &g;
			got_flag = &got_g;
		}
		else if (  ch != -1 && !inputstream.eof() && !strcmp(word, "l") )
		{
			pointer = &l;
			got_flag = &got_l;
		}
		else if (  ch == -1 && inputstream.eof() || !strcmp(word, "BASE") || !strncmp(word, "VODxx", 3) || !strncmp(word, "TRxx", 2) || !strncmp(word, "NODExx", 4) )
		{
			break;
		}
		else
			continue;

		if ( *got_flag )
			error("Error: A section cannot be initialized twice. Program stops..");
		*got_flag = true;
		while ( !isdigit(ch = inputstream.get() ) &&  ch != -1 && !inputstream.eof() )
			;
		inputstream.putback(ch);
		inputstream>>*pointer;
	}
	if (  !got_x)
		error("Error: Section: a vital parameter is missing");

	nodes[(int)beginpoint].setShuntAdmitance(nodes[(int)beginpoint].getG() + g*l/2, nodes[(int)beginpoint].getB() + b*l/2);

	branch lineBranch(beginpoint, endpoint, r*l, x*l);
	branches[branchPointer++] = lineBranch;

	node endNode(g*l/2, b*l/2);
	nodes[nodePointer++] = endNode;
}

void network::readNode()
{
	char ch;

	double num = 0;
	double Pp = 0;
	double Qp = 0;
	double Pmax = 0;
	double Qmax = 0;
	double loadType = 0;

	bool got_Num = false;
	bool got_Pp = false;
	bool got_Qp = false;
	bool got_Pmax = false;
	bool got_Qmax = false;
	bool got_LoadType = false;

	while (true)
	{
		while ( !isalpha(ch = inputstream.get() ) && !inputstream.eof() )
			;
		inputstream.putback(ch);
		inputstream>>word;
		double* pointer;
		bool* got_flag;
		if ( ch != -1 && !inputstream.eof() && !strcmp(word, "nodeNum") )
		{
			pointer = &num;
			got_flag = &got_Num;
		}
		else if ( ch != -1 && !inputstream.eof() && !strcmp(word, "Pp"))
		{
			pointer = &Pp;
			got_flag = &got_Pp;
		}
		else if ( ch != -1 && !inputstream.eof() && !strcmp(word, "Qp") )
		{
			pointer = &Qp;
			got_flag = &got_Qp;
		}
		else if ( ch != -1 && !inputstream.eof() && !strcmp(word, "Pmax") )
		{
			pointer = &Pmax;
			got_flag = &got_Pmax;
		}
		else if ( ch != -1 && !inputstream.eof() && !strcmp(word, "Qmax") )
		{
			pointer = &Qmax;
			got_flag = &got_Qmax;
		}
		else if ( ch != -1 && !inputstream.eof() && !strcmp(word, "loadtype") )
		{
			pointer = &loadType;
			got_flag = &got_LoadType;
		}
		else if ( !strcmp(word, "BASE") || !strncmp(word, "VODxx", 3) || !strncmp(word, "TRxx", 2) || !strncmp(word, "NODExx", 4) || inputstream.eof() || ch == -1)
		{
			break;
		}
		else
		{
			continue;
		}

		if ( *got_flag )
			error("Error: A node parameter cannot be initialized twice. Program stops..");
		*got_flag = true;

		while ( !isdigit(ch = inputstream.get() ) && ch != -1 && !inputstream.eof()  && ch != '-')
			;

		inputstream.putback(ch);
		inputstream>>*pointer;
	}
	if (!got_Num )
		error("Error: Node: a vital parameter is missing");

  nodes[(int)num].setNodeNum((int)num);
	nodes[(int)num].setPmax(Pmax);
	nodes[(int)num].setQmax(Qmax);
	nodes[(int)num].setLoadPower(Pp, Qp);
	nodes[(int)num].setLoadType(loadType);
}

int network::convertIndex_LayerB2Longitudinal(int i)
{
		int indexes[33];
		indexes[0] = 1;
		indexes[1] = 2;
		indexes[2] = 19;
		indexes[3] = 3;
		indexes[4] = 20;
		indexes[5] = 4;
		indexes[6] = 23;
		indexes[7] = 21;
		indexes[8] = 5;
		indexes[9] = 24;
		indexes[10] = 22;
		indexes[11] = 6;
		indexes[12] = 25;
		indexes[13] = 7;
		indexes[14] = 26;
		indexes[15] = 8;
		indexes[16] = 27;
		indexes[17] = 9;
		indexes[18] = 28;
		indexes[19] = 10;
		indexes[20] = 29;
		indexes[21] = 11;
		indexes[22] = 30;
		indexes[23] = 12;
		indexes[24] = 31;
		indexes[25] = 13;
		indexes[26] = 32;
		indexes[27] = 14;
		indexes[28] = 33;
		indexes[29] = 15;
		indexes[30] = 16;
		indexes[31] = 17;
		indexes[32] = 18;

		return indexes[i];
}


void network::error(const char* errorstring)const
{
	cerr<<errorstring<<endl;
	exit(EXIT_SUCCESS);
}


double network::getVb(int i)const
{
	if (i>=0 && i<Vb.size() )
	{
		return Vb[i];
	}
	error("Error: Wrong indexparameter in getVb()");
	return 0;
}

double network::getActiveConsumption()
{
	complex<double> S = getConsumption();
	return S.real();
}

double network::getReactiveConsumption()
{
	complex<double> S = getConsumption();
	return S.imag();
}


double network::getMinimumVoltageMagnitude()const
{
	double minMagn = getVoltageMagnitude(0);
	for(int i=1; i<numONodes; i++)
		if (minMagn > getVoltageMagnitude(i))
			minMagn = getVoltageMagnitude(i);
	return minMagn;
}

double network::getMaximumVoltageMagnitude()const
{
	double maxMagn = getVoltageMagnitude(0);
	for(int i=1; i<numONodes; i++)
		if (maxMagn < getVoltageMagnitude(i))
			maxMagn = getVoltageMagnitude(i);
	return maxMagn;
}

complex<double> network::getLosses()
{
	complex<double> losses = complex<double>(0,0);
	for(int i=0; i<numONodes; i++)
		losses += pow(nodes[i].getVoltageMagnitude(),2) * conj(nodes[i].getShuntAdmitanse());
	for(int i=0; i<numOBranches; i++)
		losses += pow(branches[i].getCurrentMagnitude(), 2) * branches[i].getImpedance();

	return losses;
}

double network::getActiveLosses()
{
	double Pp = 0;
	for (int i=0; i<numONodes; i++)
		Pp += nodes[i].getActivePowerInjection();
	return getActiveConsumption() - Pp;
}

double network::getReactiveLosses()
{
	double Qp = 0;
	for (int i=0; i<numONodes; i++)
		Qp += nodes[i].getReactivePowerInjection();
	return getReactiveConsumption() - Qp;
}

double network::getVoltageStabilityIndex(int i)
{
	double VSI;
	int j;
	j = branches[i].getBeginPoint();
	//*
	VSI = pow(nodes[j].getVoltageMagnitude(), 4)
		- 4.0*pow(getActivePowerNode(i+1)*branches[i].getReaktance()-getReactivePowerNode(i+1)*branches[i].getResistance(), 2)
		- 4.0*(getActivePowerNode(i+1)*branches[i].getResistance()+getReactivePowerNode(i+1)*branches[i].getReaktance())*pow(nodes[j].getVoltageMagnitude(),2);
	/*/
	VSI = pow(nodes[j].getVoltageMagnitude(), 4)
		- 4.0*pow(nodes[i+1].getActivePowerInjection()*branches[i].getReaktance()-nodes[i+1].getReactivePowerInjection()*branches[i].getResistance(), 2)
		- 4.0*(nodes[i+1].getActivePowerInjection()*branches[i].getResistance()+nodes[i+1].getReactivePowerInjection()*branches[i].getReaktance())*pow(nodes[j].getVoltageMagnitude(),2);
		//*/
	return VSI;
}

double network::getVoltageDeviation(int i)
{
	return abs(1 - real(nodes[i].getVoltage()));
}

double network::getAverageActiveLosses()
{
	double avgLosses = 0;
	for(int i=0; i<getNumOBranches(); i++)
		avgLosses += getLineActiveLosses(i);
	avgLosses /= getNumOBranches();
	return avgLosses*1e3;
}

double network::getMinimumVSI()
{
	double minVSI = getVoltageStabilityIndex(0);
	for(int i=0; i < getNumOBranches(); i++)
		if (minVSI > getVoltageStabilityIndex(i) )
			minVSI = getVoltageStabilityIndex(i);
	return minVSI;
}

double network::getMaximumVoltageDeviation()
{
	double maxVD = getVoltageDeviation(0);
	for(int i=0; i < getNumONodes(); i++)
	{
		if (maxVD < getVoltageDeviation(i) )
			maxVD = getVoltageDeviation(i);
	}
	return 1-maxVD;
}

bool network::isWithinVoltageLimits()
{
	for(int i=0; i<getNumONodes(); i++)
	{
		if (nodes[i].getVoltageMagnitude() > Vmax || nodes[i].getVoltageMagnitude() < Vmin)
			return false;
	}
	return true;
}

double network::getVoltageStabilityIndexSum()
{
	double sum = 0;
	for(int i=0; i<getNumOBranches(); i++)
		sum += getVoltageStabilityIndex(i);
	return sum;
}

double network::getVoltageDeviationSum()
{
	double sum = 0;
	for(int i=0; i<getNumONodes(); i++)
		sum += getVoltageDeviation(i);
	return sum;
}



void network::getVoltagePNode(ostream& os)
{
	for (int i=0; i<numONodes; i++)
		os<<nodes[i].getVoltageMagnitude()<<"\t"<<nodes[i].getVoltagePhase()<<endl;
}

void network::getLoadPowerPNode(ostream& os)
{
	for (int i=0; i<numONodes; i++)
		os<<nodes[i].getActivePowerInjection()<<"\t"<<nodes[i].getReactivePowerInjection()<<endl;
}


string network::getPowerFlowType()
{
	if (powerFlowType == "")
		return "Power flow was not calculated so far";
	return powerFlowType;
}



void network::createOctaveGraph(string graphname, string season, double (network::*f)(), bool cleanfile, string appendable)
{
	string fname = graphname;
	graphname.append(appendable);
	fname.append(".m");
	if (cleanfile)
	{
		ofstream os(fname.c_str());
		os.close();
	}
	ofstream os(fname.c_str(), fstream::app);
	os<<endl;
	os<<graphname<<" = [..."<<endl;
	for(int i=1; i<=24; i++)
	{
		solvePowerFlow(season, i);
		os<<"\t"<<i<<"  "<<(this->*f)();
		if (i<24)
			os<<";"<<endl;
		else
			os<<"];"<<endl;
	}
	os<<endl;
}

void network::relax_constraints()
{
	Vmin -= 0.005;
	Vmax += 0.005;
}

void network::printFileName ()
{
	if (strlen(filename) != 0)
		cout<<filename<<endl;
	else
		cout<<"A network has no initializer file"<<endl;
}

void network::printNode (int i)
{
	if (i>=0 && i<numONodes )
		cout<<nodes[i]<<endl;
	else
		error("Error: printNode(): Wrong indexing");
}

void network::printBranch (int i)
{
	if (i>=0 && i<numONodes )
		cout<<nodes[i]<<endl;
	else
		error("Error: printBranch(): Wrong indexing");
}

void network::setTapChanger(int tPos)
{
	// Subtracting old value
	complex<double> unitComlpex(1.0, 0);

	double trRatio = 1 + tapChanger * tapChangerDelta;
	branches[0].setZ( branches[0].getImpedance() * trRatio ) ;
	nodes[0].setShuntAdmitance( nodes[0].getShuntAdmitanse() - (1-trRatio)*(unitComlpex/branches[0].getImpedance()) );
	nodes[1].setShuntAdmitance( nodes[1].getShuntAdmitanse() - trRatio*(trRatio-1)*(unitComlpex/branches[0].getImpedance() ) );

	// Adding new value
	tapChanger = tPos;
	trRatio = 1 + tapChanger * tapChangerDelta;
	//cout<<"trRatio: "<<trRatio<<endl;
	nodes[0].setShuntAdmitance( nodes[0].getShuntAdmitanse() + (1-trRatio)*(unitComlpex/branches[0].getImpedance()) );
	nodes[1].setShuntAdmitance( nodes[1].getShuntAdmitanse() + trRatio*(trRatio-1)*(unitComlpex/branches[0].getImpedance() ) );
	branches[0].setZ( branches[0].getImpedance() / trRatio ) ;
}


void network::setCapacitorQ(int i, int unitNum, double value)
{
	value /= pow(Vb[0],2);
	value /= Sb/pow(Vb[0],2);
	temp_capacitor new_C(i, unitNum, value);
	capacitors.push_back(new_C);
	double old_cap_value = nodes[i].getB();
	nodes[i].setCapacitorB(old_cap_value + unitNum*value);
}

void network::setDispGenerator(int i, int unitNum, double Pvalue, double Qvalue)
{
	Pvalue /= Sb;
	Qvalue /= Sb;
	temp_generator new_gen(i, unitNum, Pvalue, Qvalue);
	disp_generators.push_back(new_gen);
}

void network::setWindGenerator(int i, int unitNum, double Pvalue, double Qvalue)
{
	Pvalue /= Sb;
	Qvalue /= Sb;
	temp_generator new_gen(i, unitNum, Pvalue, Qvalue);
	wind_generators.push_back(new_gen);
}

void network::setSolarGenerator(int i, int unitNum, double Pvalue, double Qvalue)
{
	Pvalue /= Sb;
	Qvalue /= Sb;
	temp_generator new_gen(i, unitNum, Pvalue, Qvalue);
	solar_generators.push_back(new_gen);
}

void network::setAllGenerators(int capNum, int* capNodes, int windNum, int* windNodes, int solarNum, int* solarNodes, int dispNum, int* dispNodes, double* dispPavlues, double* dispQvalues)
{
	removeTempElements();
	for(int i=0; i<capNum; i++)
	{
		setCapacitorQ(capNodes[i], 1);
	}
	for(int i=0; i<windNum; i++)
	{
		setWindGenerator(windNodes[i], 1);
	}
	for(int i=0; i<solarNum; i++)
	{
		setSolarGenerator(solarNodes[i], 1);
	}
	for(int i=0; i<dispNum; i++)
	{
		setDispGenerator(dispNodes[i], 1, dispPavlues[i], dispQvalues[i]);
	}
}

void network::setLoadFlowType(int powerFlowType)
{
	if (powerFlowType == 1)
		this->powerFlowType = "Constant Power";
	else if (powerFlowType == 2)
		this->powerFlowType = "Constant Current";
	else if (powerFlowType == 3)
		this->powerFlowType = "Mixed Model";
}

void network::setVoltageLimits(double Vmin, double Vmax)
{
	this->Vmin = Vmin;
	this->Vmax = Vmax;
}


void network::removeTempElements()
{

	for(int i=0; i<capacitors.size(); i++)
	{
		double minus_val = capacitors[i].getValue();
		double old_val = nodes[capacitors[i].getNode()].getB();
		nodes[capacitors[i].getNode()].setCapacitorB( old_val - minus_val );
	}
	capacitors.clear();
	disp_generators.clear();
	wind_generators.clear();
	solar_generators.clear();
}

void network::solvePowerFlow(string season, int hour)
{
	if ( !strcmp(powerFlowType.c_str(), "Constant Power"))
		solvePowerFlowConstP(season, hour);
	else if ( !strcmp(powerFlowType.c_str(), "Constant Current"))
		solvePowerFlowConstI(season, hour);
	else if ( !strcmp(powerFlowType.c_str(), "Mixed Model"))
		solvePowerFlowMixed(season, hour);
	else
		error("Error: Unknown power flow argument");
}

void network::solvePowerFlowConstP(string season, int hour)
{
	powerFlowType = "Constant Power";

	// Set up loads
	double* load_arrey;
	double* wind_arrey;
	double* solar_arrey;
	if (!season.compare("summer"))
	{
		load_arrey = load_summer;
		wind_arrey = wind_summer;
		solar_arrey = solar_summer;
	}
	else if (!season.compare("autumn"))
	{
		load_arrey = load_autumn;
		wind_arrey = wind_autumn;
		solar_arrey = solar_autumn;
	}
	else if (!season.compare("winter"))
	{
		load_arrey = load_winter;
		wind_arrey = wind_winter;
		solar_arrey = solar_winter;
	}
	else if (!season.compare("spring"))
	{
		load_arrey = load_spring;
		wind_arrey = wind_spring;
		solar_arrey = solar_spring;
	}
	else
		error("Error: Load flow wrong argument");


	for (int i=0; i<numONodes; i++)
	{
		nodes[i].setVoltage(1,0);
		nodes[i].setLoadPower(load_arrey[hour-1]*nodes[i].getPmax(), load_arrey[hour-1]*nodes[i].getQmax());
	}

	for (int i=0; i<disp_generators.size(); i++)
	{
		double old_Pvalue = nodes[ disp_generators[i].getNode() ].getActivePowerInjection();
		double old_Qvalue = nodes[ disp_generators[i].getNode() ].getReactivePowerInjection();
		nodes[ disp_generators[i].getNode() ].setLoadPower( old_Pvalue-disp_generators[i].getPvalue(), old_Qvalue-disp_generators[i].getQvalue());
	}

	for (int i=0; i<wind_generators.size(); i++)
	{
		double old_Pvalue = nodes[ wind_generators[i].getNode() ].getActivePowerInjection();
		double old_Qvalue = nodes[ wind_generators[i].getNode() ].getReactivePowerInjection();
		double wind_active = wind_arrey[hour-1] * wind_generators[i].getPvalue();
		double wind_reactive = wind_arrey[hour-1] * wind_generators[i].getQvalue();
		nodes[ wind_generators[i].getNode() ].setLoadPower( old_Pvalue - wind_active, old_Qvalue - wind_reactive);
	}

	for (int i=0; i<solar_generators.size(); i++)
	{
		double old_Pvalue = nodes[ solar_generators[i].getNode() ].getActivePowerInjection();
		double old_Qvalue = nodes[ solar_generators[i].getNode() ].getReactivePowerInjection();
		double solar_active = solar_arrey[hour-1] * solar_generators[i].getPvalue();
		double solar_reactive = solar_arrey[hour-1] * solar_generators[i].getQvalue();
		nodes[ solar_generators[i].getNode() ].setLoadPower( old_Pvalue - solar_active, old_Qvalue - solar_reactive);
	}

	// Shirmohammadi algorithm
	int f = 0;
	double epsilon = 1e-5;

	while(true)
	{
		// Calculate injection currents
		for (int i=0; i<numONodes; i++)
		{
			nodes[i].setLoadCurrent (
				conj( nodes[i].getLoadPower() / nodes[i].getVoltage() ) + nodes[i].getVoltage()*nodes[i].getShuntAdmitanse() );
		}

		// Calculate branch currents
		for (int i=numOBranches; i>0; i--)
		{
			branches[i-1].setCurrent( nodes[i].getLoadCurrent() );
			for (int j=i; j<numOBranches; j++)
				if (branches[j].getBeginPoint() == branches[i-1].getEndPoint() )
				{
					complex<double> surplus = branches[i-1].getCurrent();
					surplus += branches[j].getCurrent();
					branches[i-1].setCurrent(surplus);
				}
		}

		// Calculate voltage drops
		for (int i=1; i<numONodes; i++)
		{
			nodes[i].setVoltage(
				nodes[ branches[i-1].getBeginPoint()].getVoltage() - branches[i-1].getImpedance()* branches[i-1].getCurrent() );
			if ( branches[i-1].isRegulative() )
			{
				double a = 1 + branches[i-1].getTapChangerPos() * branches[i-1].getDeltaU()/100 ;
				nodes[i].setVoltage( nodes[i].getVoltage() / a) ;
			}
		}

		//Testing delta
		bool pass = true;
		for (int i=0; i<numONodes; i++)
		{
			complex<double> S = nodes[i].getVoltage()*conj(nodes[i].getLoadCurrent()) - conj(nodes[i].getShuntAdmitanse())*pow(
				nodes[i].getVoltageMagnitude(), 2) ;
			complex<double> dS = (S - nodes[i].getLoadPower()) ;

			if ( abs(dS.real()) > epsilon || abs(dS.imag()) > epsilon )
			{
				pass = false;
				break;
			}
		}
		f++;
		if (pass )
			break;
		if (f >= 1000)
		{
			cout<<*this<<endl;
			error("Error: The power flow failed to convergate");
		}
	}
}

void network::solvePowerFlowConstI(string season, int hour)
{
	// Set up loads
	powerFlowType = "Constant Current";

	double* load_arrey;
	double* wind_arrey;
	double* solar_arrey;
	if (!season.compare("summer"))
	{
		load_arrey = load_summer;
		wind_arrey = wind_summer;
		solar_arrey = solar_summer;
	}
	else if (!season.compare("autumn"))
	{
		load_arrey = load_autumn;
		wind_arrey = wind_autumn;
		solar_arrey = solar_autumn;
	}
	else if (!season.compare("winter"))
	{
		load_arrey = load_winter;
		wind_arrey = wind_winter;
		solar_arrey = solar_winter;
	}
	else if (!season.compare("spring"))
	{
		load_arrey = load_spring;
		wind_arrey = wind_spring;
		solar_arrey = solar_spring;
	}
	else
		error("Error: Load flow wrong argument");

	//cout<<"Fingilingi"<<endl;
	for (int i=0; i<numONodes; i++)
	{
		nodes[i].setVoltage(1,0);
		nodes[i].setLoadPower(load_arrey[hour-1]*nodes[i].getPmax(), load_arrey[hour-1]*nodes[i].getQmax());
	//	cout<<"Load Power"<<nodes[i].getLoadPower()<<endl;
	}


	for (int i=0; i<disp_generators.size(); i++)
	{
		double old_Pvalue = nodes[ disp_generators[i].getNode() ].getActivePowerInjection();
		double old_Qvalue = nodes[ disp_generators[i].getNode() ].getReactivePowerInjection();
		nodes[ disp_generators[i].getNode() ].setLoadPower( old_Pvalue-disp_generators[i].getPvalue(), old_Qvalue-disp_generators[i].getQvalue());
	}

	for (int i=0; i<wind_generators.size(); i++)
	{
		double old_Pvalue = nodes[ wind_generators[i].getNode() ].getActivePowerInjection();
		double old_Qvalue = nodes[ wind_generators[i].getNode() ].getReactivePowerInjection();
		double wind_active = wind_arrey[hour-1] * wind_generators[i].getPvalue();
		double wind_reactive = wind_arrey[hour-1] * wind_generators[i].getQvalue();
		nodes[ wind_generators[i].getNode() ].setLoadPower( old_Pvalue - wind_active, old_Qvalue - wind_reactive);
	}

	for (int i=0; i<solar_generators.size(); i++)
	{
		double old_Pvalue = nodes[ solar_generators[i].getNode() ].getActivePowerInjection();
		double old_Qvalue = nodes[ solar_generators[i].getNode() ].getReactivePowerInjection();
		double solar_active = solar_arrey[hour-1] * solar_generators[i].getPvalue();
		double solar_reactive = solar_arrey[hour-1] * solar_generators[i].getQvalue();
		nodes[ solar_generators[i].getNode() ].setLoadPower( old_Pvalue - solar_active, old_Qvalue - solar_reactive);
	}


	// Shirmohammadi algorithm
	int f = 0;
	double epsilon = 1e-5; // 31.5 kW
	double current_mods[numONodes];
	double current_args[numONodes];
	double fi[numONodes]; // From power factor

	for (int i=0; i<numONodes; i++)
	{
		nodes[i].setVoltage(1.0,0);
		if (nodes[i].getActivePowerInjection() == 0 && nodes[i].getReactivePowerInjection() == 0)
			fi[i] = 0;
		else if (nodes[i].getActivePowerInjection() == 0)
			fi[i] = M_PI/2;
		else
			fi[i] = atan(nodes[i].getReactivePowerInjection() / nodes[i].getActivePowerInjection());
		current_mods[i] = load_arrey[hour-1]*(sqrt(pow(nodes[i].getActivePowerInjection(),2) + pow(nodes[i].getReactivePowerInjection(),2)) / nodes[i].getVoltageMagnitude());
		if ( nodes[i].getActivePowerInjection() < 0)
			fi[i] = M_PI + fi[i];
	}


	while(true)
	{
		for (int i=0; i<numONodes; i++)
		{
			current_args[i] = nodes[i].getVoltagePhase() - fi[i];
			nodes[i].setLoadCurrentModArg(current_mods[i], current_args[i]);
		}


		// Calculate injection currents
		for (int i=0; i<numONodes; i++)
		{
			nodes[i].setLoadCurrent(nodes[i].getLoadCurrent() + nodes[i].getVoltage()*nodes[i].getShuntAdmitanse() );
		}

		// Calculate branch currents
		for (int i=numOBranches; i>0; i--)
		{
			branches[i-1].setCurrent( nodes[i].getLoadCurrent() );
			for (int j=i; j<numOBranches; j++)
				if (branches[j].getBeginPoint() == branches[i-1].getEndPoint() )
				{
					complex<double> surplus = branches[i-1].getCurrent();
					surplus += branches[j].getCurrent();
					branches[i-1].setCurrent(surplus);
				}
		}

		// Calculate voltage drops
		for (int i=1; i<numONodes; i++)
		{
			nodes[i].setVoltage(
				nodes[ branches[i-1].getBeginPoint()].getVoltage() - branches[i-1].getImpedance()* branches[i-1].getCurrent() );
			if ( branches[i-1].isRegulative() )
			{
				double a = 1 + branches[i-1].getTapChangerPos() * branches[i-1].getDeltaU()/100 ;
				nodes[i].setVoltage( nodes[i].getVoltage() / a) ;
			}
		}

		//Testing delta
		bool pass = true;

		for (int i=0; i<numONodes; i++)
		{
			nodes[i].setLoadCurrent(nodes[i].getLoadCurrent() - nodes[i].getVoltage()*nodes[i].getShuntAdmitanse());
			nodes[i].setLoadPower( nodes[i].getVoltage()*conj(nodes[i].getLoadCurrent()) );
		}

		for(int i=0; i<numONodes; i++)
		{
			double fi_calc;
			if (nodes[i].getLoadPower().imag() == 0 && nodes[i].getLoadPower().real() == 0)
				fi_calc = 0;
			else if (nodes[i].getLoadPower().real() == 0)
				fi_calc = M_PI/2;
			else
				fi_calc = atan( nodes[i].getReactivePowerInjection()/nodes[i].getActivePowerInjection());
			if (nodes[i].getActivePowerInjection() < 0)
				fi_calc = M_PI + fi_calc;
			double dI = nodes[i].getCurrentMagnitude() - current_mods[i];
			double dcosfi = cos( fi_calc ) - cos(fi[i]);

			if ( abs(dI) > epsilon || abs(dcosfi) > epsilon )
			{
				pass = false;
				break;
			}
		}

		f++;
		if (pass )
			break;
		if (f >= 1000)
		{
			cout<<*this<<endl;
			error("Error: The power flow failed to convergate");
		}
	}
}

void network::solvePowerFlowMixed(string season, int hour)
{
	// Set up loads
	powerFlowType = "Mixed Model";

	double* load_arrey;
	double* wind_arrey;
	double* solar_arrey;
	if (!season.compare("summer"))
	{
		load_arrey = load_summer;
		wind_arrey = wind_summer;
		solar_arrey = solar_summer;
	}
	else if (!season.compare("autumn"))
	{
		load_arrey = load_autumn;
		wind_arrey = wind_autumn;
		solar_arrey = solar_autumn;
	}
	else if (!season.compare("winter"))
	{
		load_arrey = load_winter;
		wind_arrey = wind_winter;
		solar_arrey = solar_winter;
	}
	else if (!season.compare("spring"))
	{
		load_arrey = load_spring;
		wind_arrey = wind_spring;
		solar_arrey = solar_spring;
	}
	else
		error("Error: Load flow wrong argument");

	// Shirmohammadi algorithm
	int f = 0;
	double epsilon = 1e-5; // 31.5 kW
	double current_mods[numONodes];
	double current_args[numONodes];
	double fi[numONodes]; // From power factor


	for (int i=0; i<numONodes; i++)
	{
		nodes[i].setVoltage(1.0,0);
		if (nodes[i].getQmax() == 0 && nodes[i].getPmax() == 0)
			fi[i] = 0;
		else if (nodes[i].getPmax() == 0)
			fi[i] = M_PI/2;
		else
			fi[i] = atan(nodes[i].getQmax() / nodes[i].getPmax());
		if (nodes[i].getPmax() < 0)
			fi[i] = M_PI + fi[i];
		current_mods[i] = load_arrey[hour-1]*(sqrt(pow(nodes[i].getPmax(),2) + pow(nodes[i].getQmax(),2)) / nodes[i].getVoltageMagnitude());
	}


	while(true)
	{
		for (int i=0; i<numONodes; i++)
		{
			current_args[i] = nodes[i].getVoltagePhase() - fi[i];
			nodes[i].setLoadCurrentModArg(current_mods[i], current_args[i]);
		}

		for (int i=0; i<disp_generators.size(); i++)
		{
			complex<double> old_Current = nodes[ disp_generators[i].getNode()].getLoadCurrent();
			complex<double> genPower = complex<double> (disp_generators[i].getPvalue(), disp_generators[i].getQvalue());
			nodes[ disp_generators[i].getNode() ].setLoadCurrent( old_Current - conj(genPower / nodes[ disp_generators[i].getNode() ].getVoltage() ));

	//		cout<<"current aft: "<<nodes[disp_generators[i].getNode()].getLoadCurrent()<<endl;
	//		cout<<*this<<endl;
	//		exit(EXIT_SUCCESS);
		}

		for (int i=0; i<wind_generators.size(); i++)
		{
			complex<double> old_Current = nodes[ wind_generators[i].getNode()].getLoadCurrent();
			complex<double> genPower = complex<double>(wind_arrey[hour-1]*wind_generators[i].getPvalue(), wind_arrey[hour-1]*wind_generators[i].getQvalue());
			nodes[ wind_generators[i].getNode() ].setLoadCurrent( old_Current - conj(genPower / nodes[wind_generators[i].getNode()].getVoltage() ));
		}

		for (int i=0; i<solar_generators.size(); i++)
		{
			complex<double> old_Current = nodes[ solar_generators[i].getNode()].getLoadCurrent();
			complex<double> genPower = complex<double>(solar_arrey[hour-1]*solar_generators[i].getPvalue(), solar_arrey[hour-1]*solar_generators[i].getQvalue());
			nodes[ solar_generators[i].getNode() ].setLoadCurrent( old_Current - conj(genPower / nodes[ solar_generators[i].getNode()].getVoltage() ));
		}

		// Calculate injection currents
		for (int i=0; i<numONodes; i++)
		{
			nodes[i].setLoadCurrent(nodes[i].getLoadCurrent() + nodes[i].getVoltage()*nodes[i].getShuntAdmitanse() );
		}

  	// Calculate branch currents
		for (int i=numOBranches; i>0; i--)
		{
			branches[i-1].setCurrent( nodes[i].getLoadCurrent() );
			for (int j=i; j<numOBranches; j++)
				if (branches[j].getBeginPoint() == branches[i-1].getEndPoint() )
				{
					complex<double> surplus = branches[i-1].getCurrent();
					surplus += branches[j].getCurrent();
					branches[i-1].setCurrent(surplus);
				}
		}

		// Calculate voltage drops
		for (int i=1; i<numONodes; i++)
		{
			nodes[i].setVoltage(
				nodes[ branches[i-1].getBeginPoint()].getVoltage() - branches[i-1].getImpedance()* branches[i-1].getCurrent() );
			if ( branches[i-1].isRegulative() )
			{
				double a = 1 + branches[i-1].getTapChangerPos() * branches[i-1].getDeltaU()/100 ;
				nodes[i].setVoltage( nodes[i].getVoltage() / a) ;
			}
		}

		//Testing delta
		bool pass = true;
		complex<double> I[numONodes];
		complex<double> Sload[numONodes];
		for (int i=0; i<numONodes; i++)
		{
			I[i] = nodes[i].getLoadCurrent() - nodes[i].getVoltage()*nodes[i].getShuntAdmitanse();
			nodes[i].setLoadPower( nodes[i].getVoltage()*conj(I[i]) );
			Sload[i] = nodes[i].getLoadPower();
		}

			// Kivonni a generatorok teljesit az I bol
		for (int i=0; i<disp_generators.size(); i++)
		{
			complex<double> old_Current = nodes[ disp_generators[i].getNode()].getLoadCurrent();
			complex<double> old_Power = nodes[ disp_generators[i].getNode()].getLoadPower();
			complex<double> gen_Power = complex<double> (disp_generators[i].getPvalue(), disp_generators[i].getQvalue());
			I[disp_generators[i].getNode()] = old_Current + conj(gen_Power / nodes[disp_generators[i].getNode()].getVoltage() );
			Sload[disp_generators[i].getNode()] = old_Power + gen_Power;
		}
		for (int i=0; i<wind_generators.size(); i++)
		{
			complex<double> old_Current = nodes[ wind_generators[i].getNode()].getLoadCurrent();
			complex<double> old_Power = nodes[ wind_generators[i].getNode()].getLoadPower();
			complex<double> gen_Power = complex<double>(wind_arrey[hour-1]*wind_generators[i].getPvalue(), wind_arrey[hour-1]*wind_generators[i].getQvalue());
			I[wind_generators[i].getNode()] = old_Current + conj(gen_Power / nodes[wind_generators[i].getNode()].getVoltage() );
			Sload[wind_generators[i].getNode()] = old_Power + gen_Power;
		}
		for (int i=0; i<solar_generators.size(); i++)
		{
			complex<double> old_Current = nodes[ solar_generators[i].getNode()].getLoadCurrent();
			complex<double> old_Power = nodes[ solar_generators[i].getNode()].getLoadPower();
			complex<double> gen_Power = complex<double>(solar_arrey[hour-1]*solar_generators[i].getPvalue(), solar_arrey[hour-1]*solar_generators[i].getQvalue());
			I[solar_generators[i].getNode()] = old_Current + conj(gen_Power / nodes[solar_generators[i].getNode()].getVoltage() );
			Sload[solar_generators[i].getNode()] = old_Power + gen_Power;
		}

		for(int i=0; i<numONodes; i++)
		{
			double fi_calc;
			if (Sload[i].imag() == 0 && Sload[i].real() == 0)
				fi_calc = 0;
			else if (Sload[i].real() == 0)
				fi_calc = M_PI/2;
			else
				fi_calc = atan( Sload[i].imag()/Sload[i].real());
			if (Sload[i].real() < 0)
				fi_calc = M_PI + fi_calc;
			double dI = abs(I[i]) - current_mods[i];
			double dcosfi = cos( fi_calc ) - cos(fi[i]);
		//	cout<<i<<endl;
		//	cout<<"dI: "<<dI<<endl;
		//	cout<<"dcosfi: "<<dcosfi<<endl;

			if ( abs(dI) > epsilon || abs(dcosfi) > epsilon )
			{
				pass = false;
				break;
			}
		}

		f++;
		if (pass )
			break;
		if (f >= 1000)
		{
			cout<<*this<<endl;
			error("Error: The power flow failed to convergate");
		}

	}
	/*for(int i=0; i<numONodes; i++)
	{
		nodes[i].setLoadPower( conj(original_currents[i])*nodes[i].getVoltage() );
	}
	for (int i=0; i<disp_generators.size(); i++)
	{
		complex<double> old_Power = nodes[ disp_generators[i].getNode()].getLoadPower();
		complex<double> genPower = complex<double> (disp_generators[i].getPvalue(), disp_generators[i].getQvalue());
		nodes[ disp_generators[i].getNode() ].setLoadPower(old_Power - genPower);
	}

	for (int i=0; i<wind_generators.size(); i++)
	{
		complex<double> old_Power = nodes[ disp_generators[i].getNode()].getLoadPower();
		complex<double> genPower = complex<double>(wind_arrey[hour-1]*wind_generators[i].getPvalue(), wind_arrey[hour-1]*wind_generators[i].getQvalue());
		nodes[ wind_generators[i].getNode() ].setLoadPower(old_Power - genPower);
	}

	for (int i=0; i<solar_generators.size(); i++)
	{
		complex<double> old_Power = nodes[ solar_generators[i].getNode()].getLoadCurrent();
		complex<double> genPower = complex<double>(solar_arrey[hour-1]*solar_generators[i].getPvalue(), solar_arrey[hour-1]*solar_generators[i].getQvalue());
		nodes[ solar_generators[i].getNode() ].setLoadPower(old_Power - genPower);
	}*/
	ofstream additiona_info("results_table_form.txt", ofstream::app);

	additiona_info<<"Load Model for Power Flow: "<<powerFlowType<<endl;
	additiona_info.close();
}


void network::powerFlowControlAlgorithm()
{
		for(int i=0; i<getNumOBranches(); i++)
		{
			complex<double> sum = complex<double>(0,0);
			for(int j=i+1; j<getNumOBranches(); j++)
				if(branches[j].getBeginPoint() == i+1)
					sum += branches[j].getCurrent();
			sum += nodes[i+1].getLoadCurrent();
			if (abs(branches[i].getCurrent() - sum) > 1e-5)
			{
				cerr<<"Error section: "<<i<<endl;
				error("Error, the algorithm is not consistent: Current");
			}
		}
		cout<<"Currents are OK"<<endl;

		for(int i=0; i<getNumOBranches(); i++)
		{
			int j = branches[i].getBeginPoint();
			if( abs(nodes[j].getVoltage() - nodes[i+1].getVoltage() - branches[i].getCurrent()*branches[i].getImpedance()) > 1e-5)
			{
				cerr<<"Error section: "<<i<<endl;
				error("Error, the algorithm is not consistent: Voltage");
			}
		}
		cout<<"Voltages are OK"<<endl;
		cout<<"Algorithm is OK"<<endl;
}

ostream& operator<< (ostream& os, network& n)
{
	os<<"Network from file: "<<n.filename<<endl<<endl;

	os<<"Sb = "<<n.Sb<<" W"<<endl;
	for (int i=0; i<n.Vb.size(); i++)
		os<<"Vb["<<i<<"] = "<<n.Vb[i]<<" V"<<endl;
	os<<"Power flow type: "<<n.powerFlowType<<endl;
	os<<endl;

	os<<"NODES\t[nodeNum, shuntAdmitanse, voltage, loadPower, loadCurrent, loadType]"<<endl;
	for (int i=0; i<n.numONodes; i++)
		os<<"\t"<<n.nodes[i]<<endl;
	os<<endl;
	os<<"BRANCHES\t[beginpoint, endpoint, impedance, current, tapChanger, du]"<<endl;
	for (int i=0; i<n.numOBranches; i++)
		os<<"\t"<<n.branches[i]<<endl;
	os<<endl;
	os<<"CAPACITORS\t[node, unitNum, unitValue, value]"<<endl;
	for (int i=0; i<n.capacitors.size(); i++)
		os<<"\t"<<n.capacitors[i].getNode()<<"  "<<n.capacitors[i].getUnitNum()<<"  "<<n.capacitors[i].getUnitValue()<<"  "<<n.capacitors[i].getValue()<<endl;
	os<<endl;

	os<<"DISP_GENERATORS\t[node, unitNum, Pvalue, Qvalue]"<<endl;
	for (int i=0; i<n.disp_generators.size(); i++)
		os<<"\t"<<n.disp_generators[i].getNode()<<"  "<<n.disp_generators[i].getUnitNum()<<"  "<<n.disp_generators[i].getPvalue()<<"  "<<n.disp_generators[i].getQvalue()<<endl;
	os<<endl;

	os<<"WIND_GENERATORS\t[node, unitNum, unitPvalue, Pvalue, unitQvalue, Qvalue]"<<endl;
	for (int i=0; i<n.wind_generators.size(); i++)
		os<<"\t"<<n.wind_generators[i].getNode()<<"  "<<n.wind_generators[i].getUnitNum()<<"  "<<n.wind_generators[i].getUnitPvalue()<<"  "<<n.wind_generators[i].getPvalue()<<"  "<<n.wind_generators[i].getUnitQvalue()<<"  "<<n.wind_generators[i].getQvalue()<<endl;

	os<<"SOLAR_GENERATORS\t[node, unitNum, unitPvalue, Pvalue, unitQvalue, Qvalue]"<<endl;
	for (int i=0; i<n.solar_generators.size(); i++)
		os<<"\t"<<n.solar_generators[i].getNode()<<"  "<<n.solar_generators[i].getUnitNum()<<"  "<<n.solar_generators[i].getUnitPvalue()<<"  "<<n.solar_generators[i].getPvalue()<<"  "<<n.solar_generators[i].getUnitQvalue()<<"  "<<n.solar_generators[i].getQvalue()<<endl;

	return os;
}


void network::optimalElementPlacement()
{

}
