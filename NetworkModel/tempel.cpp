#include "tempel.h"

using namespace std;

temp_capacitor::temp_capacitor()
{
  node = 1;
  value = 0;
  unitNum = 1;
}

temp_capacitor::temp_capacitor(const temp_capacitor& tc)
{
  node = tc.node;
  value = tc.value;
  unitNum = tc.unitNum;
}

temp_capacitor::temp_capacitor(int node, int unitNum, double value)
{
  this->node = node;
  this->value = value;
  this->unitNum = unitNum;
}

int temp_capacitor::getNode()const
{
  return node;
}

void temp_capacitor::setNode(const int node)
{
  this->node = node;
}

double temp_capacitor::getValue()const
{
  return unitNum*value;
}

void temp_capacitor::setValue(const double value)
{
  this->value = value;
}

void temp_capacitor::setTempCapacitor(const int node, const int unitNum, const double value)
{
  this->node = node;
  this->value = value;
  this->unitNum = unitNum;
}

temp_capacitor& temp_capacitor::operator=(const temp_capacitor& tc)
{
  this->node = tc.node;
  this->value = tc.value;
  this->unitNum = tc.unitNum;
  return *this;
}



temp_generator::temp_generator()
{
  node = 1;
  Pvalue = 0;
  Qvalue = 0;
  unitNum = 1;
}

temp_generator::temp_generator(const temp_generator& tg)
{
  this->node = tg.node;
  this->Pvalue = tg.Pvalue;
  this->Qvalue = tg.Qvalue;
  this->unitNum = tg.unitNum;
}

temp_generator::temp_generator(int node, int unitNum, double Pvalue, double Qvalue)
{
  this->node = node;
  this->Pvalue = Pvalue;
  this->Qvalue = Qvalue;
  this->unitNum = unitNum;
}

int temp_generator::getNode()const
{
  return node;
}

void temp_generator::setNode(const int node)
{
  this->node = node;
}

double temp_generator::getPvalue()const
{
  return unitNum*Pvalue;
}

void temp_generator::setPvalue(const double Pvalue)
{
  this->Pvalue = Pvalue;
}

double temp_generator::getQvalue()const
{
  return unitNum*Qvalue;
}

void temp_generator::setQvalue(const double Qvalue)
{
  this->Qvalue = Qvalue;
}

void temp_generator::setPQ(const double Pvalue, const double Qvalue)
{
  this->Pvalue = Pvalue;
  this->Qvalue = Qvalue;
}


void temp_generator::setTempGenerator(const int node, int unitNum, const double Pvalue, const double Qvalue)
{
  this->node = node;
  this->Pvalue = Pvalue;
  this->Qvalue = Qvalue;
  this->unitNum = unitNum;
}

temp_generator& temp_generator::operator= (const temp_generator& tc)
{
  this->node = tc.node;
  this->Pvalue = tc.Pvalue;
  this->Qvalue = tc.Qvalue;
  this->unitNum = tc.unitNum;

  return *this;
}
