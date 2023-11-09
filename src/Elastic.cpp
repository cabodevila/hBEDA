#include <iostream>
#include <cmath>
#include <cstdio>
#include <cstdlib>
using namespace std;

#include "Elastic.h"

//constructor & destructor
Elastic::Elastic(Distributions *fs): Kernel(fs){
  //elastic kernel
  _condensateQ = false; _isBECQ = false;
}

Elastic::~Elastic(){
}

//BEC
void Elastic::checkBEC(double f0){
  _isBECQ = (f0>f0c());
}

void Elastic::checkCondensate(){
  _condensateQ = (abs((*_fs)[0]*_fs->Ib()/_fs->Ia()-1.0)<1.0e-4);
}

//calculate kernel
void Elastic::calcFluxG(){
  /*
    Calculate gluon jf = p^2 Jg from Ia, Ib and pf.
   */
  double pfi;
  for(int i=0;i<_latt->n();i++){
    pfi=(*_fs)(_g, i)+ ((*_fs)(_g, i+1)-(*_fs)(_g, i))*_latt->psInterp(i);//get pf at ps by interpolation
    _data[i]=_fs->Ia()*( _latt->psdpInv(i)*((*_fs)(_g, i+1)-(*_fs)(_g, i)) - pfi );
    _data[i]+= _fs->Ib()*pfi*( _latt->ps(i) + pfi);
  }
  _data[_latt->n()]=0.0;	//Boundary condtion jf & jF = 0 at p = pmax
}

//calculate kernel
void Elastic::calcFluxQ(int q){
  /*
    Calculate jF = p^2 Jq from quarks (q=_q) or antiquarks (q=_qb).
   */
  double *jF;
  jF = _data + q*_np;
  double Fi;
  for(int i=0;i<_latt->n();i++){
    double dF = (*_fs)(q, i+1)-(*_fs)(q, i);
    Fi=(*_fs)(q,i)+ dF*_latt->psInterp(i);//get pf at ps by interpolation
    jF[i]= (_fs->Ia()*dF/(_latt->p(i+1)-_latt->p(i)) +  _fs->Ib()*Fi*(1.0 - Fi))*_latt->ps(i)*_latt->ps(i);
  }
  jF[_latt->n()]=0.0;	//Boundary condtion jf & jF = 0 at p = pmax
}

void Elastic::calcDiffusionDragG(){
  /*
    Diffusion and drag terms for the number densities of gluons
  */
  calcFluxG();//now save in _data
  
  if(_isBECQ&&(!_condensateQ)) checkCondensate();
  
  double jf0;
  for(int i=_latt->n(); i >= 0;i--){
    if(i==0){//i==0
      if(_condensateQ){//after onset of BEC
	jf0=-_fs->Ia()*(*_fs)(_g, 0) + _fs->Ib()*(*_fs)(_g, 0)*(*_fs)(_g, 0);
      }
      else{//before onset of BEC
	jf0=0.0;
      }
      //jf0=-_Ia*pf[0] + _Ib*pf[0]*pf[0];
      _data[i] = (_data[i]-jf0);		//Boundary condition for jf =0 at p = 0
    }//end of i==0
    else{//if i>0
      _data[i] = (_data[i]-_data[i-1]);
    }//end of i>0
    if(!isfinite(_data[i])){
      cout << "At i = " << i << ", _Cg[i] = " << _data[i] <<endl;
      _data[i]=0.0;exit(1);
    }
    _data[i]*=(3.0);//CA for g
  }

}

void Elastic::calcDiffusionDragQ(int q){
  /*
    Diffusion and drag terms for the number densities of quarks (q=_q) and antiquarks (q=_qb)
   */
  calcFluxQ(q);

  double* Cq; Cq = _data + q*_np;

  double CF = 4.0/3.0;
  for(int i=_latt->n(); i > 0;i--){
    Cq[i] = CF*(Cq[i]-Cq[i-1]);
    if(!isfinite(Cq[i])){
      cout << "At q = " << _latt->p(i) << ", Cq[i] = " << Cq[i] <<endl; exit(1);
    }
  }
  Cq[0] *= (CF);
}

void Elastic::addConversionTerms(){
  /*
    Include the the source terms for number densities
   */

  double *Cg, *Cq;
  double CF = 4.0/3.0;
  Cg = _data; Cq = _data + _np*_q;
  for(int i=0; i<_np; i++){
    double sq = _fs->Ic()*(*_fs)(_g, i)*(1.0 - (*_fs)(_q, i));
    sq-=_fs->Ibc()*(*_fs)(_q, i)*(_latt->p(i) + (*_fs)(_g, i));//always keep Ibc
    sq *= (_latt->dps(i)*CF);
    if(_fs->Basym()){
      double *Cqb; Cqb = _data + _np*_qb;
      double sqb = _fs->Ibc()*(*_fs)(_g, i)*(1.0 - (*_fs)(_qb, i));
      sqb-=_fs->Ic()*(*_fs)(_qb, i)*(_latt->p(i) + (*_fs)(_g, i));//always keep Ibc
      sqb *= (_latt->dps(i)*CF);
      Cg[i] -= (0.5*_nf*(sq+sqb));
      Cq[i] += (CF*sq); Cqb[i] += (CF*sqb);
    }else{
      Cg[i] -= (_nf*sq);
      Cq[i] += (CF*sq);      
    }
  }
}

void Elastic::calc(){
  //calculate _jf=p^2Jg
  //please run _fs->calcCoefs first.
  calcDiffusionDragG();
  if(_fs->nf()>0){
    calcDiffusionDragQ(_q);
    if(_fs->Basym()) calcDiffusionDragQ(_qb);
    addConversionTerms();
  }

  //we need convert number densities to phase-space distributions (pf for g)
  for(int i=0; i<_n; i++){
    _data[i]*=_fs->L2();
    if(i<_np) _data[i] *= _latt->areaInv(i);
    else _data[i] /= _latt->vol(i%_np);
  }
}

string Elastic::type(){
  return string("el");
}

double Elastic::f0c(){
  return 0.2315750139367749*pow(5.333333333333333 + 3.0*_fs->nf(), 4.0)/pow(10.66666666666667 + 7.0*_fs->nf(),3.0);
}

