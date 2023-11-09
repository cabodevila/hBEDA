//#include <time.h>
#include <fstream>
#include <ostream>
#include <sstream>

#ifndef DIST_H
#define DIST_H
#include "Distributions.h"
#endif

#ifndef KERN_H
#define KERN_H
#include "Kernel.h"
#endif

class KinTran{
 private:
  Distributions *_fs;
  Kernel **_kern; int _nKerns;//number of kernels
  
  double _dt, _dt_ref;	//setup time variables

  double _running_time;	//Start the stopwatch

  double _aS, _time_unit;//we take out an overall factor from kernels
 
public:
  
  KinTran(Distributions *fs, Kernel**, int=1, double=1.0e-6, double=0.1);
  void setaS(double aS){_aS=aS;_time_unit=_aS*_aS/M_PI;};

  //quantities at inital time
  void resetStopWatch();
  double getCompTime();

  //output
  void outputCompTime(ostream&);

  //calculate f at the next time step
  void nextTime1st();//1st order time integration, the same as condensate3.2.1.c
  void nextTime();

  //observables
  double dt();//get dt
  void set_dt(double);//set dt

  //save data
  void infoChem(string fNameBase);
  void save(string fNameBase, bool saveFsQ, bool saveCsQ);  

};
