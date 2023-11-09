#ifndef LATT_H
#define LATT_H
#include "Lattice.h"
#endif

#define KTDATA_H
#include <fstream>

enum { _g =0, _q =1, _qb = 2};

class KTData{
 protected:
  Lattice *_latt;

  double *_data;//point to data
  unsigned int _nIndep;//# of independent fs
  unsigned int _nf;//the number of quark flavors
  bool _BasymQ;//Is there baryon asymmetry?
  
  unsigned int _n;//the total number of data
  unsigned int _np;//number of p grids

  double _t;

 public:
  //constructor & destructor
  KTData();
  KTData(Lattice *);
  KTData(Lattice *, unsigned int);//set _latt and _flavor
  KTData(Lattice *, unsigned int, bool);//set _latt, _nf and _BasymQ  
  ~KTData();

  //initialize
  void setLattice(Lattice *);
  void setNIndep(unsigned int);
  void allocate();//allocate memory for _data;
  void setNfBQ(unsigned int, bool);//set _nf and _BasymQ to generate _flavor.
  void initData();//delete _data and set _n = 0
  
  //access to _data
  unsigned int n();//return _n;
  unsigned int np();//return _np;
  unsigned int nf();//return _nf
  unsigned int nIndep();//return _nIndep
  bool Basym();//return _BasymQ
  int index(int, int);//calculate index from flavor and ip

  double& operator[](int);//can be read and written
  double& operator()(int, int);//only read

  //calculate macroscopic quantities
  double get_ng();
  double get_eg();
  double get_sg();
  double get_nq(int);//int = _q or _qb
  double get_eq(int);
  double get_sq(int);

  //output _data
  void output(ostream &out);
  void output(ostream &out, int);//output one flavor only

  void output();//output to the screen (cout)
  void output(int);//output one flavor
  void output(const char*);//output to file
  void output(const char*, int i);//output one flavor to file
};
