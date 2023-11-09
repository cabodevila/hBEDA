#ifndef KERN_H
#define KERN_H
#include "Kernel.h"
#endif
struct FsInterp{
  double _fp, _fpox, _fp1mxox, _fpx, _fp1mx;
  double _Fp, _Fpox, _Fp1mxox, _Fpx, _Fp1mx;
  double _Fbp, _Fbpox, _Fbp1mxox, _Fbpx, _Fbp1mx;
};

class Inelastic: public Kernel{
 protected:
  double *_x;//arrays for x, x^{-2.5} and (1-x+x^2)^2.5/(x-x^2)^1.5
                                      //we only need consider x >= xmin
  //double *_xm2p5;
  int _nx, _nxMid;//the number of x grids (-1),
  double *_dx;//the measure for numerical integration over x
  double *_xm2p5;
  
  double* _Pg2gg, *_Pg2qqb, *_Pq2qg, *_Pq2gq;

public:
  Inelastic(Distributions *, int);//fs and nx
  ~Inelastic();

  //about x grids
  int getIdx(double);//getIdex of omega0/p
  int nx();//return nx
  double x(int);//read x grids
  double dx(int);//read x grids
  
  virtual string type();
  virtual string subtype()=0;

  //parallel
  void setNumThreads(int);//set the number of threads through opm

  //virtual methods
  virtual double Pg2gg(double)=0;//splitting function for g->gg
  virtual double Pg2qqb(double)=0;//splitting function for g->qqb
  virtual double Pq2gq(double)=0;//splitting function for q/qb->gq/qb
  virtual double Pq2qg(double)=0;//splitting function for q/qb->q/qbg

  //Collinear 1<->2 kernels
  double Fabc(double, double, double, double, double, double);
  double Cin(int xIdx, int flv, FsInterp*);//Here, flv = _g, _q or _qb;

  void fsInterp(int, int, FsInterp*);//calculate all fs for given _x and _p
};
