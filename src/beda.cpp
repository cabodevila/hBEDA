#include "KTRun.h"

void commandLine(int argc, char **argv){
  KTRun kt(argc, argv);

  kt.setNumThreads(8);//default num 32
  kt.set_nmin(20);
  
  //initialize
  kt.initialize(); kt.calcKerns();kt.save();

  //time integration
  kt.run();
}

void readFile(string data){
  KTRun kt;//(argc, argv);

  //set parameters using default InitCond
  //kt.setDir("/home/bwu/Documents/beda/");
  kt.set_nmin(0);kt.setnx(800);
  kt.set_dt(1e-3); kt.setKerns(true, false);
  kt.set_tMax(1e2);
  kt.setNumThreads(8);//default num 32
  kt.set_tSaveCoef(0.1);
  //kt.set_nmin(0);
  
  //initialize
  kt.initialize(data); 
  kt.calcKerns();

  kt.save();

  //time integration
  kt.run();
}

int main(int argc, char **argv){
  //take parameters from the command line
  
  if(argc==1){
    cout << "Usage 1: ./beda filename.dat can initialize with an arbitrary distribution function.\n";
    cout << "See more details in the README. \n" << endl;
    commandLine(argc, argv);
  }
  else if(strlen(argv[1])>4 && !strcmp(strrchr(argv[1], '\0') - 4, ".dat"))
    readFile(argv[1]);
  else
    commandLine(argc, argv);

  //take initial distribution from a file
  //string data = "fzero1em4nf3t8500.dat";
  //readFile(data);
  
  return 0;
}
