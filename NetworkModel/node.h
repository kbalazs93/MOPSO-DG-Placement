#ifndef NODE_H_INCLUDED
#define NODE_H_INCLUDED

#include <complex>

using namespace std;

class node
{
        static int numONodes;
        int num;
        complex<double> shuntAdmitanse;
        complex<double> voltage;
        complex<double> loadPower;
        complex<double> loadCurrent;
        double capacitor;
		    double Pmax;
        double Qmax;
		    bool capacitorSwitchPos;
		    int loadType;
    public:
        node();
        node(int i, complex<double> adm);
        node(int i, double G, double B);
        node(const node&);
        node& operator= (const node&);


        int getNumONodes()const {return numONodes; }
        int getNodeNum()const {return num; }

        complex<double> getShuntAdmitanse()const {return shuntAdmitanse;}
        complex<double> getVoltage()const {return voltage;}
        complex<double> getLoadPower()const {return loadPower;}
        complex<double> getLoadCurrent()const {return loadCurrent;}
        double getCapacitor()const {return capacitor;}
		    bool getCapacitorSwitchPos()const {return capacitorSwitchPos;}

        double getG()const {return shuntAdmitanse.real();}
        double getB()const {return shuntAdmitanse.imag();}
        double getVoltageMagnitude()const {return abs(voltage);}
        double getVoltagePhase()const {return arg(voltage);}
        double getActivePowerInjection()const {return loadPower.real();}
        double getReactivePowerInjection()const {return loadPower.imag();}
        double getCurrentMagnitude()const {return abs(loadCurrent);}
        double getCurrentAngle()const {return arg(loadCurrent);}
        double getCurrentReal()const {return loadCurrent.real(); }
        double getCurrentImag()const {return loadCurrent.imag(); }
        double getPowerFactor()const {return loadPower.real() / abs(loadPower);}
		    double getPmax()const {return Pmax; }
        double getQmax()const {return Qmax; }
		    int getLoadType()const {return loadType; }

		    void setShuntAdmitance(const complex<double>& Y) { shuntAdmitanse = Y; }
		    void setShuntAdmitance(double G, double B) { shuntAdmitanse = complex<double>(G,B); }
        void setVoltage(const complex<double>& V) {voltage = V;}
        void setVoltage(double V, double fi) {voltage = polar(V,fi); }
        void setLoadPower(const complex<double>& S) {loadPower = S;}
        void setLoadPower(double P, double Q) {loadPower = complex<double>(P,Q); }
		    void setLoadActivePower (double P) { loadPower = complex<double> (P, loadPower.imag() ); }
		    void setLoadReactivePower (double Q) { loadPower = complex<double> (loadPower.real(), Q ); }
        void setLoadCurrent(const complex<double>& I) {loadCurrent = I;}
        void setLoadCurrent(const double ReI, const double ImI) {loadCurrent = complex<double>(ReI,ImI);}
        void setLoadCurrentModArg(const double mod, const double arg);
		    void setPmax(double value) {Pmax = value; }
        void setQmax(double value) {Qmax = value; }
		    void setLoadType(int lt) {loadType = lt; }
        void setNodeNum(int n) {num = n;}

        void setCapacitorC(double C);
        void setCapacitorB(double B);
        void removeCapacitor();
		    void turnCapacitorON();
		    void turnCapacitorOFF();

        friend ostream& operator<< (ostream&, node&);
};



#endif // NODE_H_INCLUDED
