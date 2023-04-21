#ifndef NETWORK_H_INCLUDED
#define NETWORK_H_INCLUDED

#include "node.h"
#include "branch.h"
#include "tempel.h"
#include <vector>
#include <cmath>
#include <fstream>

using namespace std;

class particle_swarm;
class PSO_particle;

class network
{
		node* nodes;
		branch* branches;

		int numONodes;
		int numOBranches;

		int branchPointer;
		int nodePointer;

		int tapChanger;
		int tapChangerMin;
		int tapChangerMax;
		double tapChangerDelta;
		complex<double> firstBranchImpedance;

		double Sb;
		vector<double> Vb;
		double Vmin;
		double Vmax;

		char* filename;
		ifstream inputstream;

		char* word;
		void readTR();
		void readSection();
		void readNode();
		void readBase();
		void error(const char*)const;

		string powerFlowType;

		double load_summer[24];
		double load_autumn[24];
		double load_winter[24];
		double load_spring[24];

		double wind_summer[24];
		double wind_autumn[24];
		double wind_winter[24];
		double wind_spring[24];

		double solar_summer[24];
		double solar_autumn[24];
		double solar_winter[24];
		double solar_spring[24];

		double emptyLossesSummer12;
    double emptyVSISummer12;

		vector<temp_generator> disp_generators;
		vector<temp_generator> wind_generators;
		vector<temp_generator> solar_generators;
		vector<temp_capacitor> capacitors;

	public:
		network();
		network(network&);
		network(const char* filename, double Vmin = 0.95, double Vmax = 1.05, int powerFlowType = 1);
		~network();

		double getSb()const {return Sb; }
		double getVb(int i)const;
		double getNumONodes()const {return numONodes; }
		double getNumOBranches()const {return numOBranches; }
		void printFileName();
		void printNode(int i);
		void printBranch(int i);

		int convertIndex_LayerB2Longitudinal(int i);

		int getTapChangerPos()const {return tapChanger; }
		double getCapacitorC(int i)const {return nodes[i].getCapacitor(); }
		double getCapacitorB(int i)const {return nodes[i].getCapacitor()*100*M_PI; }
		complex<double> getVoltage(int i)const {return nodes[i].getVoltage(); }
		double getVoltageMagnitude(int i)const {return nodes[i].getVoltageMagnitude(); }
		double getVoltagePhase (int i)const {return nodes[i].getVoltagePhase(); }
		double getLoadActivePower( int i)const {return nodes[i].getActivePowerInjection(); }
		double getLoadReactivePower( int i)const {return nodes[i].getReactivePowerInjection(); }
		double getActivePowerNode(int i)const {return real(nodes[i].getVoltage()*conj(branches[i-1].getCurrent())); }
		double getReactivePowerNode(int i)const {return imag(nodes[i].getVoltage()*conj(branches[i-1].getCurrent())); }
		double getVmin()const {return Vmin; }
		double getVmax()const {return Vmax; }
		double getMinimumVoltageMagnitude()const ;
		double getMaximumVoltageMagnitude()const ;

		double getPmax (int i)const {return nodes[i].getPmax(); }
		int getLoadType (int i)const {return nodes[i].getLoadType(); }
		double getRootCurrent()const {return branches[0].getCurrentMagnitude(); }

		complex<double> getConsumption() {return nodes[0].getVoltage()*conj(branches[0].getCurrent() + nodes[0].getVoltage()*nodes[0].getShuntAdmitanse() ); }
		double getActiveConsumption() ;
		double getReactiveConsumption() ;

		complex<double> getLosses();
		double getActiveLosses() ;
		double getReactiveLosses();
		double getVoltageStabilityIndex(int);
		double getVoltageDeviation(int);
		double getLineActiveLosses(int i) { return real(nodes[branches[i-1].getBeginPoint()].getVoltage()*conj(branches[i-1].getCurrent())) - real(nodes[i].getVoltage()*conj(branches[i-1].getCurrent()));}
		bool isWithinVoltageLimits();
		double getAverageActiveLosses();
		double getMinimumVSI();
		double getMaximumVoltageDeviation();

		double getVoltageStabilityIndexSum();
		double getVoltageDeviationSum();

		complex<double> getComplexVoltage (int i) {return nodes[i].getVoltage(); }
		double getVoltageMagnitude(int i) {return nodes[i].getVoltageMagnitude(); }
		void getVoltagePNode(ostream& os);
		void getLoadPowerPNode(ostream& os);
		void createOctaveGraph(string graphname, string season, double (network::*f)(), bool cleanfile=false, string appendable="");

		double getEmptyLossesSummer12() {return emptyLossesSummer12;}
		double getEmptyVSISummer12() {return emptyVSISummer12;}

		string getPowerFlowType();



		void setTapChanger(int tPos) ;
		void setTapChangerOld (int tPos) {branches[0].setTapChanger(tPos); }

		void setCapacitorQ(int i, int unitNum, double value = 125000);
		void setDispGenerator(int i, int unitNum, double Pvalue, double Qvalue);
		void setWindGenerator(int i, int unitNum, double Pvalue = 250000, double Qvalue = 0);
		void setSolarGenerator(int i, int unitNum, double Pvalue = 250000, double Qvalue = 0);
		void setAllGenerators(int capNum, int* capNodes, int windNum, int* windNodes, int solarNum, int* solarNodes, int dispNum, int* dispNodes, double* dispPavlues, double* dispQvalues);

		void removeTempElements();
		void relax_constraints();

		void setLoadPower (int i, complex<double> value) {  nodes[i].setLoadPower(value); }
		void setLoadPower (int i, double P, double Q) {nodes[i].setLoadPower(P, Q); }
		void setLoadActivePower (int i, double P) { nodes[i].setLoadActivePower(P); }
		void setLoadReactivePower (int i, double Q) { nodes[i].setLoadReactivePower(Q); }

		void setLoadFlowType (int i) ;
		void setVoltageLimits (double Vmin, double Vmax);

		void solvePowerFlow(string, int) ;
		void solvePowerFlowConstP(string, int) ;
		void solvePowerFlowConstI(string, int) ;
		void solvePowerFlowMixed(string, int) ;

		void powerFlowControlAlgorithm();

		void optimalElementPlacement();

		friend ostream& operator<< (ostream& os, network& n);
};

#endif // NETWORK_H_INCLUDED
