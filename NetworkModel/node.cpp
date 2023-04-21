#include "node.h"
#include <cmath>
#include <iostream>

using namespace std;

int node::numONodes = 0;


node::node()
{
	num = 0;
    shuntAdmitanse = complex<double> (0,0);
    voltage = complex<double> (1,0);
    loadPower = complex<double> (0,0);
    loadCurrent = complex<double> (0,0);
    capacitor = 0;
	capacitorSwitchPos = false;
	Pmax = 0;
	loadType = 0;
}

node::node(int i, complex<double> adm)
{
    num = i;
    shuntAdmitanse = adm;
    voltage = complex<double>(1,0);
    loadPower = complex<double>(0,0);
    loadCurrent = complex<double> (0,0);
	Pmax = 0;

	loadType = 0;
    capacitor = 0;
	capacitorSwitchPos = false;
}

node::node(int i, double G, double B)
{
    num = i;
    shuntAdmitanse = complex<double>(G, B);
    voltage = polar(1,0);
    loadPower = complex<double> (0,0);
    loadCurrent = complex<double> (0,0);
	Pmax = 0;

    capacitor = 0;
	capacitorSwitchPos = false;
	loadType = 0;
}

node::node(const node &n)
{
	num = this->num;
    shuntAdmitanse = n.shuntAdmitanse;
    voltage = n.voltage;
    loadPower = n.loadPower;
    loadCurrent = n.loadCurrent;
    capacitor = n.capacitor;
	capacitorSwitchPos = n.capacitorSwitchPos;
	Pmax = n.Pmax;
	loadType = n.loadType;
}

node& node::operator= (const node &n)
{
	num = this->num;
  shuntAdmitanse = n.shuntAdmitanse;
  voltage = n.voltage;
  loadPower = n.loadPower;
  loadCurrent = n.loadCurrent;
  capacitor = n.capacitor;
	capacitorSwitchPos = n.capacitorSwitchPos;
	Pmax = n.Pmax;
	loadType = n.loadType;
}

void node::setCapacitorC(double C)
{
  shuntAdmitanse -= complex<double> (0, 100*M_PI*capacitor);
  shuntAdmitanse += complex<double> (0, 100*M_PI*C);
  capacitor = C;
	capacitorSwitchPos = true;
}

void node::setCapacitorB(double B)
{
	double C = B / 100/M_PI ;
  shuntAdmitanse -= complex<double> (0, 100*M_PI*capacitor);
  shuntAdmitanse += complex<double> (0, B);
  capacitor = C;
	capacitorSwitchPos = true;
}

void node::setLoadCurrentModArg(const double mod, const double arg)
{
	double ReI = mod * cos(arg);
	double ImI = mod * sin(arg);
	setLoadCurrent(ReI, ImI);
}

void node::removeCapacitor()
{
  shuntAdmitanse -= complex<double> (0, 100*M_PI*capacitor);
  capacitor = 0;
	capacitorSwitchPos = false;
}

ostream& operator<< (ostream& os, node& n)
{
	os<<"[ "<<n.num<<", "<<n.shuntAdmitanse<<", "<<n.voltage<<", "<<n.loadPower<<", "<<n.loadCurrent<<", "<<n.loadType<<" ] ";
	return os;
}
