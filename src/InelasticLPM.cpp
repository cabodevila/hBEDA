#include <iostream>
#include <cmath>
#include <cstdio>
#include <cstdlib>
using namespace std;

#include "InelasticLPM.h"
#include "omp.h"

InelasticLPM::InelasticLPM(Distributions *fs, int nx): Inelastic(fs, nx){
  //cout << "_nIndep in InelasticLPM = " << _nIndep << endl;
  for(int i=0;i<=_nx;i++){
    _Pg2gg[i] = Pg2gg(_x[i]); _xm2p5[i] = pow(_x[i], -2.5);
    if(_nf>0){
      _Pg2qqb[i] = Pg2qqb(_x[i]); _Pq2qg[i] = Pq2qg(_x[i]); _Pq2gq[i] = Pq2gq(_x[i]);
    }
  }
}

//description
string InelasticLPM::subtype(){
  return string("LPM");
}

//splitting functions
double InelasticLPM::Pg2gg(double x){
  double sp = 0.0, num;
  if(x>0.0&&x<1.0){
    num = 1.0 - x + x*x;
    sp=num*pow(num/(x - x*x), 1.5);
  }
  return 6.0*sp;//2Nc
}

double InelasticLPM::Pg2qqb(double x){
  double x1mx = x*(1.0-x);
  return 0.5*(x*x + (1.0-x)*(1.0-x))*sqrt((4.0/9.0-x1mx)/x1mx);
}

double InelasticLPM::Pq2gq(double x){
  double x1mx = x*(1.0-x), CF = 4.0/3.0;
  return CF*(1.0+(1.0-x)*(1.0-x))*sqrt((1.0-x+CF*x*x/3.0)/x1mx)/x;
}

double InelasticLPM::Pq2qg(double x){
  return Pq2gq(1.0-x);
}

void InelasticLPM::calc(int pIdx){
  //cout << "On " << omp_get_thread_num() << ": " << pIdx << endl;
  for(int flv=0; flv<_nIndep; flv++) _data[index(flv, pIdx)] = 0.0;
  //for(int flv=0; flv<_nIndep; flv++) cout << "Clear: " << _data[index(flv, pIdx)] << endl;
  FsInterp fsI;
  for(int i=0;i<=_nx;i++){
    fsInterp(i, pIdx, &fsI);
    for(int flv=0; flv<_nIndep; flv++){
      _data[index(flv, pIdx)] += Cin(i, flv, &fsI)*_dx[i];
    }
  }
  for(int flv=0; flv<_nIndep; flv++){
    int idata = index(flv, pIdx);
    if(flv==_g) _data[idata] *= (_rate*sqrt(_latt->p(pIdx)));
    else _data[idata] *= (_rate/sqrt(_latt->p(pIdx)));
    if(!isfinite(_data[idata])){
      //cout << "At i = " << idata << ", LPM C[idata] = " << _data[idata] <<endl;
      _data[idata] = 0.0;exit(1);
    }
  }
  //for(int flv=0; flv<_nIndep; flv++) cout << _data[index(flv, pIdx)] << endl;
}

void InelasticLPM::calc(){
  //Please run _fs->calcCoefs first
  if(_fs->Ia()>0.0)
    _rate = sqrt(3.0*_fs->Ia()*_fs->L2()/M_PI);
  else{
    cout << "Warning: Ia is negative." << endl;
    _rate = 0.0;
  }
  //#pragma omp parallel for
  int nThreads = omp_get_max_threads();
  //cout << nThreads << endl;
  int nPer = _np/nThreads;
  if(nPer*nThreads < _np) nPer++;
#pragma omp parallel
  {
   int thread = omp_get_thread_num();
   int ns = thread*nPer, ne = (thread+1)*nPer;
   if(ne>_np) ne = _np;
   //cout << _np << ": " << ns << ", " << ne << endl;
   for(int i=ns;i<ne;i++){
     calc(i);
   }
  }
}

