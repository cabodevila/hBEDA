#ifndef KERN_H
#define KERN_H
#include "Kernel.h"
#endif

class Elastic: public Kernel{
 private:

  bool _isBECQ, _condensateQ;//is there BEC, is BEC onset
  
public:
  Elastic(Distributions *);
  ~Elastic();

  void calcFluxG();//calculate p^2 jg
  void calcDiffusionDragG();//calculate the total derivative term in Cg
  void calcFluxQ(int);//calculate p^2J_F with int = _q or _qb
  void calcDiffusionDragQ(int);//calculate the total derivative term with int = _q or _qb
  void addConversionTerms();//add the source terms for q/qb<->g
  virtual void calc();

  virtual string type();

  //BEC or not
  double f0c();//critical values for onset of BEC
  void checkBEC(double);//only for step function now
  void checkCondensate();
};
