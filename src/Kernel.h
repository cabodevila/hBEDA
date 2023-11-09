#ifndef DIST_H
#define DIST_H
#include "Distributions.h"
#endif
#include <fstream>

class Kernel: public KTData{
 protected:
  //pointer to the distributions
  Distributions *_fs;
  
public:
  Kernel();
  Kernel(Distributions *);
  
  void setFs(Distributions *);//set fs

  //update time from _fs
  void updateTime();//_t only used for saving data

  //access C
  virtual void calc()=0;//to be implemented by derived kernels

  //virtual double obsp(int) = 0;
  virtual string type()=0;

  //output
  void outputChem(ostream &out);
};
