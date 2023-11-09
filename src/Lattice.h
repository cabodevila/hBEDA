#ifndef MANAGE_DATA_H
#include "ManageData.h"
#endif

#define LATT_H

class Lattice{
 private:
  int _n, _nmin;//nmin of refined grids
  double _dp,	//spacing in p grids
        _dpr,	//refined dp
        _prmin;//the critical point that seperates refined and course grids in ps

  double *_p, *_ps;//p grids, which can be shared with other classes.

  //lend->dpInv, len->dps, psp->psInterp, pslend->psdpInv
  double *_vol,		//Keeps the inverse the volume of each momentum shell
	 *_areaInv,		//Keeps the (inverse) of area of each momentum shell
	 *_dps,		//Keeps the length of each momentum shell
         *_dpInv,		//Keeps the length of each momentum shell for differentiation
         *_psInterp,             //used for the interpolation of ps: psp[i] = (ps[i]-p[i])/(p[i+1]-p[i])
         *_psdpInv;             //used for calcualte the flux: _pslend[i] = _ps[i]/(p[i+1]-p[i])

  void output(double *);

public:
  
  Lattice(int, int=20, double=10.0);
  Lattice(ManageData *);
  ~Lattice();

  void generateGrids();
  void generateGrids(ManageData *);
  void generateElements();
  
  int n();
  int nmin();
  double dp();
  double dpr();
  double prmin();
  double pMax();

  int getIdx(double p);//given p, return the index i such that p lies within p[i; and p[i+1]

  //show p grids
  void output_p();
  void output_ps();

  //share data
  double p(int);double ps(int); double psdpInv(int);
  double vol(int); double dps(int); double dpInv(int);
  double areaInv(int); double psInterp(int);
};
