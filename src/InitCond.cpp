#include <cmath>
#include <cstdlib>
using namespace std;

#include "InitCond.h"

//constructors
InitCond::InitCond(){
  _fzero = 0.0; _F0 = 0.0; _Fb0=0.0; _Qs = 1.0; _defaultQ = NULL;
  _fg = NULL; _fq = NULL; _fqb = NULL;
  setThermal(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0);
}

InitCond::InitCond(double fzero, double Qs, unsigned int nf){
  initDefault(fzero, Qs, nf);
}

InitCond::InitCond(ManageData *md, unsigned int nf){
  _fzero = 0.0; _F0 = 0.0; _Fb0=0.0; _Qs = 1.0; _defaultQ = NULL;
  _fg = NULL; _fq = NULL; _fqb = NULL;
  setThermal(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0);
  _type = "file";

  _read = true;
  _md = md;
  _nf = nf;
}

//initialize
double InitCond::getDefaultT(){
  return 0.80050698377213744*_Qs*pow(_fzero/(10.666666666666666 + 7.0*_nf),0.25);
}

void InitCond::initDefault(double fzero, double Qs, unsigned int nf){
  _fzero = fzero; _Qs = Qs; _nf = nf;
  _defaultQ = true; _type = "step";
  _fg = NULL; _fq = NULL; _fqb = NULL; 

  //        including _eg, _eq, _eqb, _ng, _nq, _nqb, _sg, _sq, _sqb
  _T = getDefaultT(); _mu = 0.0;

  double Tn = pow(_T, 4.0);
  _eg = 5.2637890139143246*Tn;
  _eq = 3.4543615403812755*Tn*_nf/2.0; _eqb = _eq;

  Tn = pow(_T,3.0);
  _ng = 1.9487012517371693*Tn;
  _nq = 1.0961444541021579*_nf*Tn/2.0; _nqb = _nq;

  _sg = 7.018385351014531*Tn;
  _sq = 4.60581538723398*_nf*Tn/2.0; _sqb = _sq;
}

void InitCond::initSelfDef(double (*fg)(double), double (*fq)(double), double (*fqb)(double)){//self-defined functions
  _fzero = 0.0; _Qs = 1.0; _defaultQ = NULL;
  _fg = fg; _fq = fq; _fqb = fqb; _defaultQ = false;

  setThermal(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0);
  cout << "InitCond: the self-defined functions are used.\n";
  cout << "          Please evalue the thermal quantites using setThermal." << endl;
}

void InitCond::setType(const char* type){
  cout << "Warning: The default distributions are used.\n";
  cout << "          _type is to be changed from 'step' to '";
  cout << type << "'." << endl; 
  _type = string(type);
}

string InitCond::type(){
  return _type;
}

void InitCond::setThermal(double T, double eg, double eq, double eqb, double ng, double nq, double nqb, double sg, double sq, double sqb){
  _T=T; _eg=eg; _eq=eq; _eqb=eqb;
  _ng=ng; _nq=nq; _nqb=nqb; _sg=sg; _sq=sq; _sqb=sqb;
}

//some implemented functions
double InitCond::step(double p){
  double f0;
  
  if(p < _Qs)
      f0=_fzero;
    else
      f0 = 0.0;
  
  return f0;
}

double InitCond::thermalg(double p, double mu){
  return 1.0/(exp((p-mu)/_T)-1.0);
}

double InitCond::thermalq(double p, double mu){
  return 1.0/(exp((p-mu)/_T)+1.0);
}

//initial functions
double InitCond::pf0(double p){	
  double f0i;
  if(_defaultQ){
    f0i = p*step(p);
  }
  else{
    if(_fg==NULL){
      cout << "Error in InitCond: please initialize self-defined fg." << endl;
      exit(EXIT_FAILURE);
    }else{
      f0i = p*(*_fg)(p);
    }
  }
  return f0i;
}

double InitCond::F0(double p){	
  double f0i;
  if(_defaultQ){
    f0i = 0.0;
    if(p < _Qs) f0i=_F0;
  }else{
    if(_fq==NULL){
      cout << "Error in InitCond: please initialize self-defined fq." << endl;
      exit(EXIT_FAILURE);
    }else{
      f0i = (*_fq)(p);
    }
  }
  return f0i;
}

double InitCond::Fb0(double p){	
  double f0i;
  if(_defaultQ){
    f0i = 0.0;
    if(p < _Qs) f0i=_Fb0;
  }else{
    if(_fqb==NULL){
      cout << "Error in InitCond: please initialize self-defined fqb." << endl;
      exit(EXIT_FAILURE);
    }else{
      f0i = (*_fqb)(p);
    }
  }
  return f0i;
}

void InitCond::outputT(ostream &out){
  out << "# type = " << _type;
  if(_defaultQ){
    out << ": f0 = " << _fzero << ", Qs = " << _Qs << ", nf = " << _nf <<endl;
  }
  out << "# T = " << _T << ", mu = " << _mu << "\n";
  out << "# ng = " << _ng  << ", nq = " << _nq   << ", nqb = " << _nqb << "\n";
  out << "# eg = " << _eg  << ", eq = " << _eq   << ", eqb = " << _eqb << "\n";
  out << "# sg = " << _sg  << ", sq = " << _sq   << ", sqb = " << _sqb << endl;
}

void InitCond::outputT(){
  ostringstream out;
  outputT(out);
  cout << out.str() << endl;
}

double InitCond::operator()(int flv, double p){
  double res = 0.0;
  switch(flv){
  case 0:
    res = pf0(p);
    break;
  case 1:
    res = F0(p);
    break;
  case 2:
    res = Fb0(p);
    break;
  }
  return res;
}

double InitCond::operator()(int flv, int ip, double p){
  return (*_md)(ip, flv+1);
}

void InitCond::summary(ostream& out){
  out << "# Summary of InitCond:\n";
  if(_defaultQ){
    outputT(out);
  }else{
    out << "# type = " << _type << ", Nf = " << _nf;
  }
  cout << "\n";
}
