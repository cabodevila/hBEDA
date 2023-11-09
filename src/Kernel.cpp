#include <iostream>
#include <cmath>
#include <cstdlib>
using namespace std;

#include "Kernel.h"

Kernel::Kernel(){
  _fs = NULL;
}

Kernel::Kernel(Distributions *fs){
  setFs(fs);
}

void Kernel::setFs(Distributions *fs){
  _fs = fs;
  setLattice(_fs->latt()); setNfBQ(_fs->nf(), _fs->Basym());
  initData(); allocate();
  cout << "# Kernel has been initialized!\n" << endl;
  
}

void Kernel::updateTime(){
  _t = _fs->t();
}

//output contributiosn from g, q, qb to dn/dt
void Kernel::outputChem(ostream &out){
  out << " " << get_ng() << " " << get_nq(_q) << " " << get_nq(_qb);
  out << " " << get_eg() << " " << get_eq(_q) << " " << get_eq(_qb);
  out << " " << get_sg() << " " << get_sq(_q) << " " << get_sq(_qb);
}
