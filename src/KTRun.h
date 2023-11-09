#include<iostream>
#include <cstring>
#include<cstdlib>
#include<cmath>
using namespace std;

#include "Distributions.h"
#include "KinTran.h"
#include "Elastic.h"
#include "InelasticLPM.h"

#ifndef MANAGE_DATA_H
#include "ManageData.h"
#endif

class KTRun{
private:
  //overall parameters
  unsigned int _nf;
  bool _BasymQ;
  double _t0, _dt, _tMax, _tSave;	//setup time variables
  double _aS;
  bool _elQ, _inelQ;

  string _version;//keep version number
  string _dir;//directory to save
  //initial distributions
  double _Qs, _fzero, _F0, _Fb0;//Initial
  InitCond *_init;
  string _fNInit;//record initial file name

  //grids
  int 	_n, _nmin, _nx; 	//nmin of refined grids
  double _pMax;//the maximum
  Lattice *_latt;

  //fs
  Distributions* _fs;
  
  //kernels
  int _nKerns;
  Kernel **_kern; Elastic *_el; Inelastic *_inel;

  //time integration
  KinTran *_kt;

  //save data
  ostringstream _ostr;
  bool _saveFsQ; //false by default, determines whether one wants to keep pf, F and Fb.
  bool _saveCsQ;//save Kernels
  float _tSaveCoef;//determine when to save the date
  
  //number of threads
  int _numThreads;
public:
  //constructor & destructor
  KTRun();
  KTRun(int argc, char **argv);
  ~KTRun();

  //set parameters
  void setDefaultParas();//default values of the parameters
  //There are two ways to modify the values of the parameters:
  //1. command line
  void takeParas(int argc, char **argv);//command line
  void setDir(const char* dir);
  //2. set parameters using
  void setaS(double aS) {_aS = aS;};
  void setfzero(double fzero, double Qs=1.0){_fzero=fzero; _Qs=Qs;};
  void setF0(double f0){_F0=f0;};
  void setFb0(double f0){_Fb0=f0;};
  
  void setns(unsigned int n, unsigned int nmin, unsigned int nx){_n=n; _nmin=nmin; _nx=nx;};
  void setnx(unsigned int nx){_nx=nx;};
  void setNfBQ(unsigned int nf, bool Ba){_nf=nf; _BasymQ=Ba;};
  void set_t0(double t){_t0=t;};
  void set_dt(double dt){_dt=dt;};
  void set_tMax(double tMax){_tMax=tMax;};
  void setKerns(bool elQ, bool inelQ){_elQ=elQ; _inelQ=inelQ;};
  void set_tSave(double tSave){_tSave=tSave;};
  void set_tSaveCoef(double tSaveCoef){_tSaveCoef=tSaveCoef;};
  void set_nmin(int nmin){_nmin = nmin;};
  void set_pMax(double pmax){_pMax = pmax;};
  //Then, run initialize
  void initialize();
  void initialize(string);
  void initialize_kernels();
  void outputParas(ostream&);

  //time integration
  void run();
  
  //save data
  void createFiles();
  void save();
  void saveFs(bool fsQ){_saveFsQ=fsQ;};
  void saveCs(bool CsQ){_saveCsQ=CsQ;};
  void saveCompTime();

  void calcKerns();
  void setNumThreads(int n){_numThreads=n;};
};
