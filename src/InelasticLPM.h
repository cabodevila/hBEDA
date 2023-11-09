#include "Inelastic.h"
/*
Here we use the LPM spectrum in the inelastic kernel.
BU stands for bottom-up
 */
class InelasticLPM: public Inelastic{
 private:
  double _rate;//1/lambda
  
 public:
  InelasticLPM(Distributions *, int=100);//both elastic and inelastic

  //description of this kernel
  virtual string subtype();
  
  //splitting functions
  virtual double Pg2gg(double);//splitting function for g->gg
  virtual double Pg2qqb(double);//splitting function for g->qqb
  virtual double Pq2gq(double);//splitting function for q/qb->gq/qb
  virtual double Pq2qg(double);//splitting function for q/qb->q/qbg

  void calc(int pIdx);
  virtual void calc();//calculate the inelastic kernel Cel from f

};
