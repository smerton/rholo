// A simple mesh generator for DG problem input

// Author S. R. Merton

#include <iostream>
#include "input.h"

void Output(char*s,Input*I); // RhoLo/MFEM format
void Output_old(char*s,Input*I); // old output format

using namespace std;

// main function starts here

int main(int argc, char** argv){

  bool new_format(true);

  if(argc!=3){
    cout<<"main(): please specify an input file and an output file."<<endl;
    exit(1);
  }

  cout<<"main(): input file: "<<argv[1]<<endl; // note [0] is the name of the program
  cout<<"main(): output file: "<<argv[2]<<endl;

  cout<<"main(): Starting up main loop..."<<endl;

// initialise a new input from file

  Input I((char*) argv[1]);

// write the mesh file

  if(new_format){

// new format for RhoLo/MFEM

    Output((char*) argv[2],&I);

  }else{

// old format for backward compatibility with ToyALE

    Output_old((char*) argv[2],&I);

  }

  cout<<"Normal termination."<<endl;

  return 0;
}
