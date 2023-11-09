#include <cmath>
#include <iostream>
#include <limits>
#include <iomanip>

using namespace std;

#include "KinTran.h"
#include "omp.h"
#include "Inelastic.h"

KinTran::KinTran(Distributions *fs, Kernel **kern, int nKerns, double dt, double aS){
  _kern = kern; _nKerns = nKerns; _fs = fs;
  _dt = dt; _dt_ref = _dt; setaS(aS);
  
  resetStopWatch();
}

//deal with times

double KinTran::dt(){
  return _dt;
}

void KinTran::set_dt(double dt){
  _dt = dt;
}

//computation times
void KinTran::resetStopWatch(){
  _running_time=omp_get_wtime();//clock()/CLOCKS_PER_SEC;
}

double KinTran::getCompTime(){
  return omp_get_wtime() - _running_time;//clock()/CLOCKS_PER_SEC
}

void KinTran::outputCompTime(ostream& out){
  	//Compuation time
	_running_time=getCompTime();

	if(_running_time>=60.0){
		if(_running_time>=3600){
			out << "#Calculation is done and the running time is " << _running_time/3600.0 << " hours.\n\n\n";
		}
		else {
			out << "#Calculation is done and the running time is " << _running_time/60.0 << " minutes.\n\n\n";
		}
	}
	else {
		out << "#Calculation is done and the running time is " << _running_time << " seconds.\n\n\n";
	}
}

void KinTran::nextTime1st(){
  //time marching
  _fs->calcCoefs();
  for(int i=0; i<_nKerns; i++){
    _kern[i]->calc();
    for(int j=0; j< _fs->n(); j++){
      (*_fs)[j]+=_time_unit*_dt*(*_kern[i])[j];
    }
  }  
}

void KinTran::nextTime(){  
  nextTime1st(); //show_pf();
 
  _fs->setTime(_fs->t()+_dt);
}

//save data
void KinTran::infoChem(string fName){
  ofstream out;
  out.precision(15); out.open(string(fName+".chem.dat").c_str(),ios::app);
  out << "#{t,  Iag, Iaq, Iaqb, ";
  out << "Ibg, Ibq, Ibqb, Icg, Icq, Ibcqb";
  for(int i=0;i<_nKerns;i++){
    string nC = string("dn") + _kern[i]->type();
    out << ", " << nC + "g" <<  ", " << nC + "q" <<  ", " << nC + "qb" ;
    nC = string("de") + _kern[i]->type();
    out << ", " << nC + "g" <<  ", " << nC + "q" <<  ", " << nC + "qb" ;
    nC = string("ds") + _kern[i]->type();
    out << ", " << nC + "g" <<  ", " << nC + "q" <<  ", " << nC + "qb" ;
  }
  out << "}"; out.close();
}

void KinTran::save(string fName, bool saveFsQ, bool saveCsQ){
  _fs->outputMacs(string(fName+".mac.dat").c_str());

  //output contributions from each d.o.f.
  ofstream out;
  out.precision(15); out.open(string(fName+".chem.dat").c_str(),ios::app);
  _fs->outputChem(out);
  for(int i=0;i<_nKerns;i++) _kern[i]->outputChem(out);
  out << endl; out.close();

  if(saveFsQ) _fs->output(string(fName+".fs.dat").c_str());

  if(saveCsQ){
    for(int i=0;i<_nKerns;i++){
      _kern[i]->updateTime(); _kern[i]->output(string(fName+".C" + _kern[i]->type() +".dat").c_str());
    }
  }  
}
