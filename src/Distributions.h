#ifndef KTDATA_H
#include "KTData.h"
#endif

#ifndef INIT_H
#define INIT_H
#include "InitCond.h"
#endif

#define DIST_H

class Distributions: public KTData{
 private:  
  //\hat{q} = 4aS^2N_c Ia/pi, m_D^2 = 2 aS Ib/pi
  double _Ia,  _Ib, _Ic, _Ibc, _L2;//L2 = ln(pt^2/mD^2)
  //separate contributions from gluon and (anti)quarks
  double _Iag, _Iaq, _Iaqb, _Ibg, _Ibq, _Ibqb;
  double _Icg, _Icq, _Ibcqb;
  
 public:
  Distributions(Lattice *, unsigned int=0, bool=false);//latt, nf and BasymQ

  //initialize
  void initialize(InitCond *);//initialize the distribution functions

  //return the values of Ia and Ib
  double Ia(){return _Ia;};
  double Ib(){return _Ib;};
  double Ic(){return _Ic;};
  double Ibc(){return _Ibc;};
  double L2(){return _L2;};

  //interpolation functions
  double interp(double*, double);//interpolation of an array f in p
  double atp(int, double);//f at p for int = 0 (pf), 1 (F) and 2 (Fb)
  
  void calcCoefs();//calculate Is
  void showCoefs();

  Lattice *latt(){return _latt;};

  //output macroscopic quantities
  void infoMacs(ostream &out);
  void outputMacs(ostream &out);
  void outputChem(ostream &out);
  void outputMacs(const char*);

  //deal with time
  double t(){return _t;};//return _t
  void setTime(double t){_t =t;};

};
