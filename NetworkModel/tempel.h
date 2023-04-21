#ifndef TEMP_ELEMENT_H
#define TEMP_ELEMENT_H

class temp_capacitor
{
    int node;
    double value;
    int unitNum;
  public:
    temp_capacitor();
    temp_capacitor(int node, int unitNum, double value);
    temp_capacitor(const temp_capacitor& );
    int getNode()const;
    void setNode(const int);
    double getValue()const; // Get: unitNum*value
    double getUnitValue()const {return value;}
    void setValue(const double); // Set: single value
    int getUnitNum() {return unitNum;}
    void setUnitNum(int num) {unitNum = num;}
    void setTempCapacitor(const int, const int, const double);
    temp_capacitor& operator= (const temp_capacitor&);
};

class temp_generator
{
    int node;
    double Pvalue;
    double Qvalue;
    int unitNum;
  public:
    temp_generator();
    temp_generator(int node, int unitNum, double Pvalue, double Qvalue);
    temp_generator(const temp_generator& );
    int getNode()const;
    void setNode(const int);
    double getPvalue()const;    // Get: unitNum*value
    double getUnitPvalue()const {return Pvalue;}
    void setPvalue(const double); // Set: single value
    double getQvalue()const;
    double getUnitQvalue()const {return Qvalue;}
    void setQvalue(const double);
    void setPQ(const double, const double);
    int getUnitNum() {return unitNum;}
    void setUnitNum(int num) {unitNum = num;}
    void setTempGenerator(const int, const int, const double, const double);
    temp_generator& operator= (const temp_generator&);
};

#endif
