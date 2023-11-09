#include <iostream>
#include <cmath>
#include <cstdlib>
using namespace std;

#include "Distributions.h"

//constructor
Distributions::Distributions(Lattice *latt, unsigned int nf, bool BQ): KTData(latt, nf, BQ){
  _Ia = 0.0; _Ib = 0.0; _Ic = 0.0; _Ibc = 0.0; _L2 = 1.0;
  cout << "# Distrubtuions has been initialized!\n" << endl;
}

//initialization
void Distributions::initialize(InitCond *init){
  for(int flv=0; flv<_nIndep; flv++){
    for(int ip=0; ip <_np; ip++){
      //_data[index(flv, ip)] = (*init)(flv, _latt->p(ip));
      if (!(*init).get_read())
        _data[index(flv, ip)] = (*init)(flv, _latt->p(ip));
      else
        _data[index(flv, ip)] = (*init)(flv, ip, _latt->p(ip));
    }
  }
}

//calculate macroscropic quantities
void Distributions::showCoefs(){//show Is
  outputMacs(cout);
}
    
void Distributions::calcCoefs(){//calculate Is
  //univeral part for all
  
  //#at refined p 
  double m2Gref=0.0, qhatGref=0.0;
  double *pf; pf = _data;
  int i;
  for(i=0;i<_latt->nmin();i++){
    m2Gref+=pf[i];qhatGref+=pf[i]*(_latt->p(i)+pf[i]);
  }
  //#at regular p
  double m2G=0.0, qhatG=0.0;
  for(i=_latt->nmin();i<_np;i++){
    m2G+=pf[i];qhatG+=pf[i]*(_latt->p(i)+pf[i]);
  }
  m2G = m2G*_latt->dp()+m2Gref*_latt->dpr();//pf integrated over p

  _Ibg = 6.0*m2G;//2Nc
  _Ib = _Ibg;
  _Iag = 3.0*(qhatG*_latt->dp()+qhatGref*_latt->dpr());//Nc
  _Ia = _Iag;

  //all quanlitites proportional to _nf are set to be zero:
  _Iaq = 0.0; _Iaqb = 0.0; _Ibq = 0.0; _Ibqb = 0.0;
  _Ic = 0.0; _Ibc = 0.0; _Icq = 0.0; _Ibcqb = 0.0;

  if(_nf>0){//If there are quarks
    //integf=integf*dp+integfnm*dpr;integF=(plen*F).sum();
    double m2Q=0.0, qhatQ=0.0;//pF, p^2*F(1-F) integrated over p
    double *F;
    F = _data+_q*_np;
    for(i=0; i<_np; i++){
      qhatQ+=_latt->vol(i)*F[i]*(1.0-F[i]);//p^2F(1-F) integrated over p
      m2Q+=_latt->p(i)*_latt->dps(i)*F[i];//pF integrated over p
    }
    if(_BasymQ){//If the distributions of quarks and antiquarks are different
      double *Fb;
      double m2Qb=0.0, qhatQb=0.0;//pFb, p^2*Fb(1-Fb) integrated p
      Fb = _data+_qb*_np;
      for(i=0; i<_np; i++){
	qhatQb+=_latt->vol(i)*Fb[i]*(1.0-Fb[i]);
	m2Qb+=_latt->p(i)*_latt->dps(i)*Fb[i];
      }

      double pfdFref = 0.0, pfdF = 0.0;//pf*(F-Fb) integrated p
      for(i=0;i<_latt->nmin();i++){
	pfdFref+=pf[i]*(F[i]-Fb[i]);
      }
      for(i=_latt->nmin();i<_np;i++){
	pfdF+=pf[i]*(F[i]-Fb[i]);
      }
      pfdF = pfdF*_latt->dp()+pfdFref*_latt->dpr();

      _Iaq = 0.5*_nf*(qhatQ); _Iaqb = 0.5*_nf*(qhatQb); _Ia += _Iaq + _Iaqb;
      _Ibq = _nf*(m2Q); _Ibqb = _nf*(m2Qb); _Ib += _Ibq + _Ibqb;
      _Icg = m2G; _Icq = m2Q; _Ic = m2G + m2Q + pfdF;
      _Ibcqb = m2Qb; _Ibc = m2G + m2Qb - pfdF;
    }else{
      _Iaq = 0.5*_nf*qhatQ; _Iaqb = _Iaq; _Ia += _nf*qhatQ;
      _Ibq = _nf*m2Q; _Ibqb = _Ibq; _Ib += 2.0*_nf*m2Q;
      _Icg = m2G; _Icq = m2Q; _Ibcqb = _Icq; _Ic = m2G + m2Q; _Ibc = _Ic;
    }
  }

  //handle divergences
  if(!isfinite(_Ia)||!isfinite(_Ib)||!isfinite(_Ic)||!isfinite(_Ibc)){
    showCoefs(); exit(1);
  }
}

double Distributions::interp(double *f, double p){
  double pfp = 0.0;
  if(p <= _latt->p(0)){
    pfp = f[0];
  }else if( p >= _latt->p(_latt->n()) ){
    pfp = 0.0;
  }else{
    int idx = _latt->getIdx(p);
    if(idx==0||idx==_latt->n()-1||(idx==_latt->nmin()-1)){
      pfp = f[idx] + (p - _latt->p(idx))*(f[idx + 1] - f[idx])/(_latt->p(idx + 1) - _latt->p(idx));
    }else{//using 4 points
      double x1 = _latt->p(idx-1), x2 = _latt->p(idx), x3 = _latt->p(idx+1), x4 = _latt->p(idx+2);
      if(idx >=_latt->n()||idx <=0){
	cout << "Something wrong" << endl;
	cout << "idx = " << idx << endl;
	cout << x1 << ", " << x2 << ", " << x3 << ", " << x4 << endl;
      }
      pfp = f[idx-1]*(p-x2)*(p-x3)*(p-x4)/((x1-x2)*(x1-x3)*(x1-x4));
      pfp += f[idx]*(p-x1)*(p-x3)*(p-x4)/((x2-x1)*(x2-x3)*(x2-x4));
      pfp += f[idx+1]*(p-x1)*(p-x2)*(p-x4)/((x3-x1)*(x3-x2)*(x3-x4));
      pfp += f[idx+2]*(p-x1)*(p-x2)*(p-x3)/((x4-x1)*(x4-x2)*(x4-x3));
    }
    
  }
  return pfp;
}

double Distributions::atp(int flv, double p){//interpolation functions
  double fatp = 0.0;
  if((flv==_g)||(flv==_q)||(flv==_qb)){
    if(flv==_g){
      fatp = interp(_data + flv*_np ,p);
    }else{
      if(_nf>0){
	if(_BasymQ) fatp = interp(_data + flv*_np ,p);
	else fatp = interp(_data + _q*_np ,p);
      }
    }
  }
  return fatp;
}

//output macroscopic quantities
void Distributions::infoMacs(ostream& out){
  out << "#{t, Ia = qhatA/(alpha^2 N_c log/pi), ";
  out << "Ib = m_d^2/(2 alpha/pi), Ic, Ibc, ";
  out << "ng, nq, nqb, eg, eq, eqb, sg, sq, sqb}" << endl;		
}

void Distributions::outputMacs(ostream& out){
  //Format changed for V6.21, changed back to the formate for V6.1x
  out << _t << " " << _Ia << " " << _Ib << " " << _Ic  << " " << _Ibc;
  out << " " << get_ng() << " " << get_nq(_q) << " " << get_nq(_qb);
  out << " " << get_eg() << " " << get_eq(_q) << " " << get_eq(_qb);
  out << " " << get_sg() << " " << get_sq(_q) << " " << get_sq(_qb);
  out << endl;
}

void Distributions::outputChem(ostream& out){
  out << _t << " " << _Iag  << " " << _Iaq  << " " << _Iaqb  << " ";
  out << _Ibg  << " " << _Ibq  << " " << _Ibqb  << " ";
  out << _Icg  << " " << _Icq  << " " << _Ibcqb;
}

void Distributions::outputMacs(const char* fn){
  ofstream out; out.open(fn,ios::app);
  outputMacs(out); out.close();
}
