#include "KTRun.h"
#include "omp.h"
#include <cmath>
#include <iostream>
#include <limits>
#include <iomanip>

using namespace std;

//constructor & destructor
KTRun::KTRun(){
    setDefaultParas();
    cout << "KTRun: All the parameters have default values.\n";
    cout << "       Please use setxxx methods to modify them if needed.\n";
    cout << "       Then, run initialize to generate all the objects." <<endl;
}

KTRun::KTRun(int argc, char **argv){
    setDefaultParas();
    cout << "# KTRun: Paras are set to default values!\n";
    takeParas(argc, argv);
    cout << "# KTRun: Paras are taken from the command line!\n";
}

KTRun::~KTRun(){
  delete _init, _latt, _kt, _fs, _el, _inel; delete [] _kern;
}

//parameters
void KTRun::outputParas(ostream& out){
  out << "\n# Parameters:\n";
  out << "# Init: Qs = " << _Qs << ", f0 = " << _fzero;
  if(_F0!=0.0) out <<  ", F0 = " << _F0;
  if(_Fb0!=0.0) out <<  ", Fb0 = " << _Fb0;
  out << ", aS = " << _aS << ", nf = " << _nf;
  if(_BasymQ) out << " with baryon asymmetry";
  else out << " with baryon symmetry";
  out << ", nThreads = " << _numThreads;
  out << "\n# Lattice: n_p = " <<  _n << ", n_pmin = " << _nmin << ", n_x = " << _nx;
  out << ", p_max = " << _pMax << "\n";
  out << "# Times: t0 = " << _t0 << ", dt = " << _dt << ", tMax = " << _tMax << "\n";
  out << "# tSave = " << _tSave << ", f's ";
  if(_saveFsQ) out << "to be saved";
  else out << "not to be saved";
  out << ", C's ";
  if(_saveCsQ) out << "to be saved, ";
  else out << "not to be saved, ";
  out << "\n" << endl;
}

void KTRun::setDefaultParas(){
  _version = "V7.0";_dir = ""; _numThreads = 32;
  _Qs = 1.0; _aS=0.1; _nf = 3; _BasymQ = false;
  _fzero = 0.1; _F0 = 0.0; _Fb0 = 0.0; _n=200; _nmin=20; _nx = 400;
  _t0 = 0.0; _dt = 1.0e-5; _tMax=1000.0;	//setup time variables
  _saveFsQ = true; _saveCsQ = true; _tSave = 0.1; _tSaveCoef = 0.25;
  _elQ = true; _inelQ = true; _pMax = 5.0;
}

void KTRun::takeParas(int argc, char **argv){
  //Take the parameters
  if((argc == 1)||((strcmp(argv[1],"e")!=0)&&(strcmp(argv[1],"eB")!=0)&&(strcmp(argv[1],"i")!=0) \
		    &&(strcmp(argv[1],"b")!=0))){
    cout << "Usage 2: ./beda type fzero (aS) (np) (dt) (nx) (tMax) (Nf) (Ba) (pMax)\n";
    cout << "type: e=elastic only, i=inelastic only, b=both.\n";
    cout << "Ba:   T (true) or F/anything else (false), Ba = F by default.\n \n";
    exit(0);
  }else{
    if(strcmp(argv[1],"e")==0) {_elQ = true; _inelQ = false;}
    else if(strcmp(argv[1],"i")==0) {_elQ = false; _inelQ = true;}
    else if(strcmp(argv[1],"b")==0) {_elQ = true; _inelQ = true;}
    if(argc>=3){
      _fzero = atof(argv[2]);
      if(argc>=4){
	_aS = atof(argv[3]);
	if(argc>=5){
	  _n = atoi(argv[4]);
	  if(argc>=6){
	    _dt = atof(argv[5]);
	    if(argc>=7){
	      _nx = atoi(argv[6]);
	      if(argc>=8){
		_tMax = atof(argv[7]);
		if(argc>=9){
		  _nf = atoi(argv[8]);
		  if(argc>=10){
		    if(strcmp(argv[9],"T")==0){
		      _BasymQ = true;
		    }
		    if(argc>=11){
		      _pMax = atof(argv[10]);
		      cout << "p_max = " <<  _pMax << endl;
		    }
		  }
		}
	      }
	    }
	  }
	}
      }
    }
  }
}

void KTRun::initialize(){
  outputParas(cout);
  cout << "# KTRun: initialzing..." << endl;
  _init = new InitCond(_fzero, _Qs, _nf); _init->setF0(_F0);; _init->setFb0(_Fb0);
  _latt = new Lattice(_n, _nmin, _pMax);
  _fs = new Distributions(_latt, _nf, _BasymQ);
  _fs->initialize(_init); _fs->setTime(_t0);//_t keeps the initial time
    
  initialize_kernels();
  
  _kt = new KinTran(_fs, _kern, _nKerns,  _dt, _aS);
  createFiles();
  cout << "# KTRun: all objects are initialized." << endl;
}

void KTRun::initialize(string name){
  _fNInit = name;

  ManageData md(name);

  // Adapt default parameters to data
  _n = md.number_rows()-1;
  _nmin = md.get_nmin();
  _t0 = md.get_time();
  _nf = md.get_nf();
  if (md.number_columns() == 2) { _nf = 0; _BasymQ = false;}
  else if (md.number_columns() == 3) _BasymQ = false;
  else if (md.number_columns() == 4) _BasymQ = true;

  _ostr << md.get_name();

  cout << "# KTRun: initializing using external data..." << endl;

  _init = new InitCond(&md, _nf);
  _latt = new Lattice(&md); _pMax = _latt->pMax();
  _fs = new Distributions(_latt, _nf, _BasymQ);
  _fs->initialize(_init); _fs->setTime(_t0);//_t keeps the initial time
    
  outputParas(cout);
  initialize_kernels();
  
  _kt = new KinTran(_fs, _kern, _nKerns,  _dt, _aS);

  cout << "# KTRun: all objects are initialized." << endl;

}

void KTRun::initialize_kernels(){
  if(_elQ&&!_inelQ){
    _nKerns = 1;
    _el = new Elastic(_fs); _el->checkBEC(_fzero);
    _kern = new Kernel*[_nKerns]; _kern[0] = _el;
    cout << "# Elastic collision kernel only." << endl;
  }else if(_inelQ&&!_elQ){
    _nKerns = 1;
    _inel = new InelasticLPM(_fs, _nx);
    _inel->setNumThreads(_numThreads);

    _kern = new Kernel*[_nKerns]; _kern[0] = _inel;
    cout << "# Inelastic collision kernel only." << endl;  
  }else{
    _nKerns = 2; 
    _el = new Elastic(_fs);
    _inel = new InelasticLPM(_fs, _nx);
    _inel->setNumThreads(_numThreads);

    _kern = new Kernel*[_nKerns]; _kern[0] = _el; _kern[1] = _inel;
    cout << "# Both collision kernels are included." << endl;
  }
}

//save data
void KTRun::save(){
  _kt->save(_ostr.str(), _saveFsQ, _saveCsQ);
}

void KTRun::createFiles(){
  ofstream out;
  cout << "# KTRun: Creating files...\n";
  _ostr << _dir;
  if(_init->get_read()){
    cout << "file name:" << _fNInit << endl; 
    _ostr << _fNInit;
  }else{
    _ostr << "fzero"  << _init->fzero();
    if(_F0!=0.0) _ostr <<  "F0" << _F0;
    if(_Fb0!=0.0) _ostr <<  "Fb0" << _Fb0;
  }
  _ostr<< "aS"  << _aS << "Qs"  << _init->Qs();
  if(_nf>0){
    _ostr << "Nf" << _nf;
    if(_BasymQ) _ostr << "Basym";
    else _ostr << "Bsym";
  }
  _ostr << "N"  << _latt->n();
  for(int i=0;i<_nKerns;i++){
    if(_kern[i]->type().compare("inel")==0)
      _ostr << "Nx" << ((Inelastic*)_kern[i])->nx();
  }
  
  _ostr << "dtm"  << setfill('0') << setw(1) << int(-log10(_dt)) << _version;
  _ostr << "nThread" << _numThreads;
  for(int i=0;i<_nKerns;i++){
    _ostr << "." << _kern[i]->type();
    if(_kern[i]->type().compare("inel")==0)
      _ostr << ((Inelastic*)_kern[i])->subtype();
  }

  _ostr << "." << _init->type();

  if(_saveFsQ){
    out.open(string(_ostr.str()+".fs.dat").c_str(),ios::app);
    out << "# {p, p*f";
    if(_nf>0){out << ", F"; if(_BasymQ) out << ", Fb";}
    out << "}\n";
    out.close();
  }

  if(_saveCsQ){
    for(int i=0; i<_nKerns; i++){
      out.open(string(_ostr.str()+".C" + _kern[i]->type() +".dat").c_str(),ios::app);
      out << "# {p, Cg";
      if(_nf>0){out << ", Cq"; if(_BasymQ) out << ", Cqb";}
      out << "}\n";
      out.close();
    }
  }

  out.open(string(_ostr.str()+".mac.dat").c_str(),ios::app);
  out << "# "; _init->outputT(out); out << "\n";//output thermal quantites
  _fs->infoMacs(out); out.close(); _kt->infoChem(_ostr.str());

  out.open(string(_ostr.str()+".time.dat").c_str(),ios::app);
  outputParas(out);
  out.close();
  
  cout << "# KTRun: Files have been just created." << endl;
}

//time integration
void KTRun::run(){
  unsigned int nSave = (unsigned int)(_tSave/_dt), ns = nSave;//time to save the data

  int logtMax = ceil(log10(_tMax));
  int logdt = log10(_dt);
  
  do{
    _kt->nextTime();
    //cout << "# t = " << _fs->t() << endl;
    
    if(ns==nSave){
      save(); ns=0;
    }
    
    for(int i=logdt;i<logtMax;i++){
      double dts;
      if(i==logdt) dts = pow(10.0,i);
      else dts = _tSaveCoef*pow(10.0,i);

      if((_fs->t()>pow(10.0,i))&&(_fs->t()<=pow(10.0,i+1))&&(_tSave!=dts)){
	_tSave = dts; nSave = (unsigned int)(_tSave/_dt);
	if(ns!=0){
	  save(); ns=0;
	}
      }
    }
    ns++;
  }while(_fs->t()<=_tMax);

  saveCompTime();
}

void KTRun::saveCompTime(){
  ofstream out;
  out.open(string(_ostr.str()+".time.dat").c_str(),ios::app);
  _kt->outputCompTime(out);
  out.close();
}

void KTRun::calcKerns(){
  for(int i=0;i<_nKerns;i++){
    _kern[i]->updateTime(); _fs->calcCoefs(); _kern[i]->calc();
  } 
}

void KTRun::setDir(const char* dir){
  _dir = dir;
  if(_dir.back()!='/') _dir.append("/");
  cout << _dir << endl;
}
