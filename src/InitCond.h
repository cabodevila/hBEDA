#include <iostream>
#include <fstream>
#include <sstream>
#include <string>

#ifndef MANAGE_DATA_H
#include "ManageData.h"
#endif

using namespace std;

#define INIT_H



class InitCond{
 private:
  //parameters
  bool _defaultQ;//is the default initial conditions?
  string _type;
  unsigned int _nf;
  
  //thermal quantities
  double _T, _mu;//mu is for all parton number with only 2<->2; mu = mu_B otherwise
  double _eg, _eq, _eqb;//energy densities
  double _ng, _nq, _nqb;//number densities
  double _sg, _sq, _sqb;//entropy densities

  //default, only for F=Fb
  double _fzero, _Qs;
  double _F0, _Fb0;//for quark and antiquark distributions

  //self defined functions
  double (*_fg)(double), (*_fq)(double), (*_fqb)(double);

  // Distiguish when reading from file
  bool _read = false;
  ManageData *_md;

 public:
  //constructor
  InitCond();
  InitCond(double, double=1.0, unsigned int=0);//inputs: fzero, Qs and nf
  InitCond(ManageData *, unsigned int=0); //inputs: ManageData of input data and nf

  //initialize
  void initDefault(double, double, unsigned int);
  void initSelfDef(double (*)(double), double (*)(double)=NULL, double (*)(double)=NULL);
  void setThermal(double, double, double, double, double, double, double, double, double, double);
  void initData(ManageData *, unsigned int);

  double getDefaultT();

  //set parameters
  void setfzero(double F0){_fzero = F0;};
  void setF0(double F0){_F0 = F0;};
  void setFb0(double F0){_Fb0 = F0;};
  //return parameters
  double Qs(){return _Qs;};
  double T(){return _T;};
  double fzero(){return _fzero;};
  double F0(){return _F0;};
  double Fb0(){return _Fb0;};

  //default initial distributions
  double step(double p);
  
  //thermal distributions
  double thermalg(double, double);
  double thermalq(double, double);

  //distributions for g, q, qb
  double pf0(double);
  double F0(double p);//initial distribution for q
  double Fb0(double p);//initial distribution for qb
  
  void outputT(ostream &out);//output thermal quantities
  void outputT();

  //deal with types
  void setType(const char*);
  string type();//return the type of initial conditions.

  //overload ()
  double operator()(int, double);
  double operator()(int, int, double);

  //summary
  void summary(ostream &out);//summarize the parameters 

  // Distiguish when reading from file
  bool get_read(){return _read;};
};
