#include <iostream>
#include <cmath>
#include <cstdlib>
#include <iomanip>
using namespace std;

#include "KTData.h"

//constructor & destructor
/*constructors always with _nf = 0 and _Basym = false*/
KTData::KTData(){
  setLattice(NULL); setNIndep(0); initData();
}

KTData::KTData(Lattice *latt){
  setLattice(latt); setNIndep(0); initData();
}

KTData::KTData(Lattice *latt, unsigned int nIndep){
  setLattice(latt); setNIndep(nIndep); initData();
}

//constructor for f and C: _nIndep = 1, 2, 3 so far
KTData::KTData(Lattice *latt, unsigned int nf, bool BQ){
  setLattice(latt); setNfBQ(nf, BQ);
  initData(); allocate();
  cout << "# KTData has been initialized with nf = " << _nf << ", Baryon asymetry = ";
  if(_BasymQ) cout << "true";
  else cout << "false";
  cout << endl;
  //cout << "# Hint: please then run allocate to allocate memory for _data." << endl;
}

KTData::~KTData(){
  delete [] _data;
}

//initialize
void KTData::setLattice(Lattice *latt){
  if(latt==NULL){
    _latt = NULL; _np = 0;
  }else{
    _latt = latt; _np = _latt->n()+1;
  }
}

void KTData::setNIndep(unsigned int nIndep){
  _nIndep = nIndep; _nf = 0; _BasymQ = false;
}

void KTData::setNfBQ(unsigned int nf, bool BQ){
  _nf = nf; _BasymQ = BQ;
  _nIndep = 1;//g
  if(nf>0){
    if(_BasymQ){
      _nIndep = 3;//g, q, qb
    }else{
      _nIndep = 2;//g, q
    }
  }
}

void KTData::initData(){
  _data = NULL; _n = 0;
}

void KTData::allocate(){
  bool status;
  if((_latt!=NULL)&&(_nIndep>0)){
    if(_data!=NULL){
      delete [] _data;
    }
    _n = _nIndep*_np;
    _data = new double[_n];
    cout << "# KTData: " << _n*sizeof(double)<< " bytes has been successfully allocated for _data." << endl;
  }else{
    cout << "# Please evaluate _latt and _nIndep in KTData before memory allocation." << endl;
    status = false;
  }
}

//access to _data
unsigned int KTData::nf(){
  return _nf;
}

bool KTData::Basym(){
  return _BasymQ;
}

unsigned int KTData::nIndep(){
  return _nIndep;
}

unsigned int KTData::n(){
  return _n;
}

unsigned int KTData::np(){
  return _np;
}

int KTData::index(int flv, int ip){
  return flv*_np + ip;
}

double& KTData::operator[](int i){
  if(i < 0 || i>=_n){
    cout << "Error in KTData: out of the range of _data\n";
    cout << i << " is not in the range [0, " << _n << "]\n" << endl;
    exit(EXIT_FAILURE);
  }
  return _data[i];
}

double& KTData::operator()(int flv, int ip){
  
  if(flv < 0 || flv>=_nIndep){
    cout << "Error in KTData: out of the range of _data(int,int)\n";
    cout << "The nIndep index =" <<  flv;
    cout << " is not in the range [0, " << _nIndep << "]\n" << endl;
    exit(EXIT_FAILURE);
  }
  
  if(ip < 0 || ip>=_np){
    cout << "Error in KTData: out of the range of _data(int,int)\n";
    cout << "The p index =" <<  ip;
    cout << " is not in the range [0, " << _np-1 << "]\n" << endl;
    exit(EXIT_FAILURE);
  }

  return _data[index(flv, ip)];
}

void KTData::output(ostream &out){
  out << "# t = " << _t << "\n";
  for(int ip = 0; ip < _np; ip++){
    out << _latt->p(ip);
    for(int flv = 0; flv < _nIndep; flv++){
      out << " " << operator()(flv, ip);
    }
    out << "\n";
  }
  out << "\n" << endl;
}

void KTData::output(ostream &out, int flv){
  out << "# t = " << _t << "\n";
  for(int ip = 0; ip < _np; ip++){
    out << _latt->p(ip);
    out << " " << operator()(flv, ip);
    out << endl;
  }
}

void KTData::output(){
  if(_n>0 && _np>0){
    output(cout);
  }else{
    cout << "There is not data to show.\n" << endl;
  }
}

void KTData::output(int flv){
  if(_n>0 && _np>0){
    output(cout, flv);
  }else{
    cout << "There is not data to show.\n" << endl;
  }
}

void KTData::output(const char* fn){
  if(_n>0 && _np>0){
    ofstream out;
    out.precision(15);
    out.open(fn,ios::app);
    output(out);
    out.close();
  }else{
    cout << "Empty _data: no output to " << fn << endl;
  }
}

void KTData::output(const char* fn, int flv){
  if(_n>0 && _np>0){
    ofstream out;
    out.precision(15);
    out.open(fn,ios::app);
    output(out, flv);
    out.close();
  }else{
    cout << "Empty _data: no output to " << fn << endl;
  }
}

//calculate macroscopic quantites

double KTData::get_ng(){
  double ng = 0.0;
  for(int i=0;i<_np;i++){
    ng += _data[i]/_latt->areaInv(i);
  }
  return 8.0*ng/(M_PI*M_PI);
}

double KTData::get_nq(int flv){
  double n = 0.0;
  if(_nf>0){
    double *F;
    if(_BasymQ) F = _data + flv*_np;
    else F = _data + _q*_np;
    for(int i=0;i< _np;i++){
      n += _latt->vol(i)*F[i];
    }
    n = 3.0*_nf*n/(M_PI*M_PI);
  }
  return n;  
}

double KTData::get_eg(){
  double eg = 0.0, *_pf = _data;
  for(int i=0;i<=_latt->n();i++){
    eg += _latt->vol(i)*_pf[i];
  }
  return 8.0*eg/(M_PI*M_PI);
}

double KTData::get_eq(int flv){
  double eq = 0.0;
  if(_nf>0){
    double *F;
    if(_BasymQ) F = _data + flv*_np;
    else F = _data + _q*_np;
    for(int i=0;i< _np;i++){
      eq += _latt->vol(i)*_latt->p(i)*F[i];
    }
    eq = 3.0*_nf*eq/(M_PI*M_PI);
  }
  return eq;  
}

double KTData::get_sg(){
  double sg = 0.0, f, f1, *_pf = _data;
  for(int i=0;i<=_latt->n();i++){
    f = _pf[i]/_latt->p(i); f1 = 1.0 + f;
    if(f>0){
      sg += _latt->vol(i)*(f1*log(f1) - f*log(f));
    }else{
      sg += _latt->vol(i)*f1*log(f1);
    }
  }
  return 8.0*sg/(M_PI*M_PI);
}

double KTData::get_sq(int flv){
  double s = 0.0;
  if(_nf>0){
    double F1, *F;
    if(_BasymQ) F = _data + flv*_np;
    else F = _data + _q*_np;
    for(int i=0;i<_np;i++){
      F1 = 1.0 - F[i];
      if(F[i]>0&&F1>0){
	s -= _latt->vol(i)*(F[i]*log(F[i]) + F1*log(F1));
      }
    }
    s = 3.0*_nf*s/(M_PI*M_PI);
  }
  return s;
}

