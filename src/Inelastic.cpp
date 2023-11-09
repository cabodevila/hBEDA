#include <iostream>
#include <cmath>
#include <cstdio>
#include <cstdlib>
using namespace std;

#include "Inelastic.h"
#include "omp.h"

Inelastic::Inelastic(Distributions *fs, int nx): Kernel(fs){

  _nxMid = nx/2;_nx = 2*_nxMid;

  //inelastic kernel
  _x = new double[_nx + 1]; _dx = new double[_nx + 1];
  _Pg2gg = new double[_nx + 1]; _xm2p5 = new double[_nx + 1];
  if(_nf>0){
    _Pg2qqb = new double[_nx + 1];
    _Pq2qg = new double[_nx + 1];
    _Pq2gq = new double[_nx + 1];
  }

  double dlx = ( log(0.5)-log(_fs->latt()->p(0)/_fs->latt()->p(_fs->latt()->n()) ) )/((double)_nxMid);
  //cout << _latt->p(0)/_latt->p(_latt->n()) << endl;

  _x[_nxMid] = 0.5;
  for(int i=0;i<_nxMid;i++){
    _x[i] = 0.5*exp(-(_nxMid-i)*dlx);
    _x[_nx-i] = 1.0 - _x[i];
  }

  _dx[0] = 0.5*(_x[1] - _x[0]);
  _dx[_nx] = 0.5*(_x[_nx] - _x[_nx - 1]);
  for(int i=1;i<_nx;i++){
    _dx[i] = 0.5*(_x[i+1] - _x[i-1]);
  }

  //for(int i=0;i<=_nx;i++) _xm2p5[i] = pow(_x[i], -2.5);

}

Inelastic::~Inelastic(){
  delete [] _x; delete [] _dx; delete [] _Pg2gg; delete [] _xm2p5;
  if(_nf>0){
    delete [] _Pg2qqb; delete [] _Pq2qg; delete [] _Pq2gq;
  }
}

double Inelastic::x(int xidx){
  return _x[xidx];
}

double Inelastic::dx(int xidx){
  return _dx[xidx];
}

int Inelastic::nx(){
  return _nx;
}

string Inelastic::type(){
  return string("inel");
}

void Inelastic::setNumThreads(int num){
  omp_set_num_threads(num);
}

//The collinear approximation
void Inelastic::fsInterp(int xIdx, int pIdx, FsInterp* fsI){
  double x = _x[xIdx], p = _latt->p(pIdx);
  double px = x*p, p1mx = (1.0-x)*p, pox = p/x, p1mxox = (1.0-x)*pox;
  fsI->_fp = (*_fs)(_g, pIdx)/p; fsI->_fpox = _fs->atp(_g, pox)/pox;
  fsI->_fp1mxox =  _fs->atp(_g, p1mxox)/p1mxox;
  fsI->_fpx =  _fs->atp(_g, px)/px; fsI->_fp1mx = _fs->atp(_g, p1mx)/p1mx;
  if(_nf>0){
    fsI->_Fp = (*_fs)(_q, pIdx); fsI->_Fpox = _fs->atp(_q, pox);
    fsI->_Fp1mxox =  _fs->atp(_q, p1mxox);
    fsI->_Fpx = _fs->atp(_q, px); fsI->_Fp1mx = _fs->atp(_q, p1mx);
    if(_BasymQ){
      fsI->_Fbp = (*_fs)(_qb, pIdx); fsI->_Fbpox = _fs->atp(_qb, pox);
      fsI->_Fbp1mxox =  _fs->atp(_qb, p1mxox);
      fsI->_Fbpx = _fs->atp(_qb, px); fsI->_Fbp1mx = _fs->atp(_qb, p1mx);
    }else{
      fsI->_Fbp = fsI->_Fp; fsI->_Fbpox = fsI->_Fpox;
      fsI->_Fbp1mxox =  fsI->_Fp1mxox;
      fsI->_Fbpx = fsI->_Fpx; fsI->_Fbp1mx = fsI->_Fp1mx;
    }
  }
}

double Inelastic::Fabc(double fa, double ea, double fb, double eb, double fc, double ec){
  return fa*(1.0+eb*fb)*(1.0+ec*fc) - fb*fc*(1.0+ea*fa);
}

double Inelastic::Cin(int xIdx, int flv, FsInterp* fsI){
  double c = 0.0;
  if(flv==_g){
    if(_x[xIdx]>0.5)
      c = _Pg2gg[xIdx]*(_xm2p5[xIdx]*Fabc(fsI->_fpox, 1.0, fsI->_fp, 1.0, fsI->_fp1mxox, 1.0)\
			-  Fabc(fsI->_fp, 1.0, fsI->_fpx, 1.0, fsI->_fp1mx, 1.0));
    else
      c = _Pg2gg[xIdx]*_xm2p5[xIdx]*Fabc(fsI->_fpox, 1.0, fsI->_fp, 1.0, fsI->_fp1mxox, 1.0);

    if(_nf>0){
      double cq2gq, CF = 4.0/3.0;
      cq2gq = Fabc(fsI->_Fpox, -1.0, fsI->_fp, 1.0, fsI->_Fp1mxox, -1.0);
      if(_fs->Basym()) cq2gq += Fabc(fsI->_Fbpox, -1.0, fsI->_fp, 1.0, fsI->_Fbp1mxox, -1.0);
      else cq2gq *= 2.0;
      c += 0.5*_nf*_xm2p5[xIdx]*_Pq2gq[xIdx]*cq2gq/CF;
      c -= _nf*_Pg2qqb[xIdx]*Fabc(fsI->_fp, 1.0, fsI->_Fpx, -1.0, fsI->_Fbp1mx, -1.0);
    }
  }else if(flv==_q){
    double CF = 4.0/3.0;
    c = _Pq2qg[xIdx]*(_xm2p5[xIdx]*Fabc(fsI->_Fpox, -1.0, fsI->_Fp, -1.0, fsI->_fp1mxox, 1.0)\
		      - Fabc(fsI->_Fp, -1.0, fsI->_Fpx, -1.0, fsI->_fp1mx, 1.0));
    c += 2.0*CF*_xm2p5[xIdx]*_Pg2qqb[xIdx]*Fabc(fsI->_fpox, 1.0, fsI->_Fp, -1.0, fsI->_Fbp1mxox, -1.0);
  }else if(flv==_qb){
    double CF = 4.0/3.0;
    c = _Pq2qg[xIdx]*(_xm2p5[xIdx]*Fabc(fsI->_Fbpox, -1.0, fsI->_Fbp, -1.0, fsI->_fp1mxox, 1.0)\
		      - Fabc(fsI->_Fbp, -1.0, fsI->_Fbpx, -1.0, fsI->_fp1mx, 1.0));
    c+=2.0*CF*_xm2p5[xIdx]*_Pg2qqb[xIdx]*Fabc(fsI->_fpox, 1.0, fsI->_Fbp, -1.0, fsI->_Fp1mxox, -1.0);
  }else{
    cout << "Inelastic Warning: call C for unknown partons!" << endl;
  }

  //cout << c << endl;
  return c;
}
