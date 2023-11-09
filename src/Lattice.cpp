#include <iostream>
#include <cmath>
#include <cstdlib>
using namespace std;

#include "Lattice.h"

Lattice::Lattice(int n, int nmin, double pMax){
  _n = n; _nmin = nmin;

  _dp = pMax/((double)_n); _dpr = 0.1*_dp;	//refined dp

  _prmin = _nmin*_dpr;

  _p = new double[_n+1];
  _ps = new double[_n+1];
  _vol = new double[_n+1];
  _areaInv = new double[_n+1];
  _dps = new double[_n+1];
  _dpInv = new double[_n+1];
  _psInterp = new double[_n+1];
  _psdpInv = new double[_n+1];

  generateGrids();
  generateElements();
}

Lattice::Lattice(ManageData *md){
  _n = md->number_rows()-1;
  //_dpr = 0.1*_dp; _prmin = _nmin*_dpr; //refined dp, replaced by
  _nmin = 0; _dpr = 0.0; _prmin = _nmin*_dpr;//For now, for the refined grid to be zero

  _p = new double[_n+1];
  _ps = new double[_n+1];
  _vol = new double[_n+1];
  _areaInv = new double[_n+1];
  _dps = new double[_n+1];
  _dpInv = new double[_n+1];
  _psInterp = new double[_n+1];
  _psdpInv = new double[_n+1];

  generateGrids(md);
  generateElements();
}

Lattice::~Lattice(){
  delete [] _p; delete [] _ps;
  delete [] _vol; delete [] _areaInv; delete [] _dps;
  delete [] _dpInv; delete [] _psInterp; delete [] _psdpInv;
}

void Lattice::generateGrids(){
  int i;
  for(i=0;i<=_n;i++){
    if(i<_nmin)
      _ps[i]=(i+1)*_dpr;
    else
      _ps[i]=_ps[i-1] + _dp;
    if(i==0) _p[0] = 0.5*_ps[0];
    else _p[i]=0.5*(_ps[i] + _ps[i-1]);
  }
}

void Lattice::generateGrids(ManageData *md){
  // Generate grid readed in md
  int i;
  for(i=0;i<=_n;i++){
      _p[i]=(*md)(i, 0);
    if(i>=1) _ps[i-1] = 0.5*(_p[i] + _p[i-1]);
    if(i==_n) _ps[i]=2.*_p[i] - _ps[i-1];
  }
  //_dp = pMax/((double)_n); //replaced by
  _dp = _ps[_nmin+1] - _ps[_nmin];//Again, this is only for homogeneous grids
}

void Lattice::generateElements(){
  for(int i=0;i<=_n;i++){
    double vol;
    if(i==0){
      _vol[i]=pow(_ps[0],3.0)/3.0;
      _areaInv[0]=2.0/pow(_ps[0],2.0);
      _dps[0]=_ps[0];
    }
    else{
      _vol[i]=(pow(_ps[i],3.0)-pow(_ps[i-1],3.0))/3.0;
      _areaInv[i]=2.0/( pow(_ps[i],2.0)-pow(_ps[i-1],2.0) );
      _dps[i]=_ps[i]-_ps[i-1];
    }
    if(i<_n){
      _dpInv[i]=1.0/(_p[i+1]-_p[i]);
    }else{
      _dpInv[_n] = _dpInv[_n-1];
    }
    _psInterp[i] = (_ps[i]-_p[i])*_dpInv[i];
    _psdpInv[i]=_ps[i]*_dpInv[i];	//Multiplication is cheaper than division.
  } 
}

int Lattice::n(){
  return _n;
}

int Lattice::nmin(){
  return _nmin;
}

double Lattice::dp(){
  return _dp;
}

double Lattice::pMax(){
  return _ps[_n];
}

double Lattice::dpr(){
  return _dpr;
}

int Lattice::getIdx(double p){
  int idx;
  if(p < _prmin){
    idx = static_cast<int>(floor(p/_dpr - 0.5));
  }else{
    idx = static_cast<int>(floor((p - _prmin)/_dp - 0.5 + _nmin));
  }

  //there might be some rounding error
  if(idx == _n){//it must be rounding error since p is less than p[n] > p[0] in Kernel
    idx = _n - 1;
    //cout << "\n#Warning: at p = " << p << ", idx = _n with p - _p[_n] = " << p - _p[_n] << endl;
  }else if(idx == -1){
    idx = 0;
    //cout << "\n#Warning:  at p = " << p << ", idx = -1 with p - _p[0] = " << p - _p[0] << endl;
  }else{
    if(p < _p[idx]){
      //cout << "\n#Warning:  at p = " << p << ", p < _p[idx] with p - _p[idx] = " << p - _p[idx] << endl;
      idx -= 1;
    }
  }

  return idx;
}

void Lattice::output(double *p){
  for(int i=0;i<=_n;i++){
    cout << p[i];
    if(i < _n)
      cout << ", ";
  }
  cout << endl;
}

void Lattice::output_p(){
  cout << "p = "; output(_p);
}

void Lattice::output_ps(){
  cout << "ps = "; output(_ps);
}

double Lattice::prmin(){
  return _prmin;
}

double Lattice::p(int i){
  if(i>=0&&i<=_n){
    return _p[i];
  }else{
    cout << "Index " << i << " is out of bound of _p.\n";
    exit(1);
  }
}

double Lattice::ps(int i){
  if(i>=0&&i<=_n){
    return _ps[i];
  }else{
    cout << "Index " << i << " is out of bound of _ps.\n";
    exit(1);
  }
}

double Lattice::vol(int i){
  if(i>=0&&i<=_n){
    return _vol[i];
  }else{
    cout << "Index " << i << "is out of bound of _vol.\n";
    exit(1);
  }
}

double Lattice::areaInv(int i){
  if(i>=0&&i<=_n){
    return _areaInv[i];
  }else{
    cout << "Index " << i << " is out of bound of _areaInv.\n";
    exit(1);
  }
}

double Lattice::dps(int i){
  if(i>=0&&i<=_n){
    return _dps[i];
  }else{
    cout << "Index " << i << " is out of bound of _dps.\n";
    exit(1);
  }
}

double Lattice::dpInv(int i){
  if(i>=0&&i<=_n){
    return _dpInv[i];
  }else{
    cout << "Index " << i << " is out of bound of _dpInv.\n";
    exit(1);
  }
}

double Lattice::psInterp(int i){
  if(i>=0&&i<=_n){
    return _psInterp[i];
  }else{
    cout << "Index " << i << " is out of bound of _psInterp.\n";
    exit(1);
  }
}


double Lattice::psdpInv(int i){
  if(i>=0&&i<=_n){
    return _psdpInv[i];
  }else{
    cout << "Index " << i << " is out of bound of _psdpInv.\n";
    exit(1);
  }
}
