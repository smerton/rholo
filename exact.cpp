// This function codes for an exact solution where applicable

#define VVD vector<VD>    // vector of VD
#define TAYLOR 1              // Taylor-Green vortex problem
#define RAYLEIGH 2            // Rayleigh-Taylor instability problem
#define NOH 3                 // Noh stagnation shock problem
#define SEDOV 4               // Sedov expanding shock problem
#define TRIPLE 5              // Triple point problem
#define SOD 6                 // Sod's shock tube problem
#define R2R 7                 // 123 problem
#define SALTZMANN 8           // Saltzmann piston problem

#include <iostream>
#include <vector>
#include <iomanip>    // floating point precision
#include <cmath>      // for sqrt
#include <fstream>    // for file io
#include <algorithm>  // min_element, max_element
#include "mesh.h"
#include "shape.h"
#include "riemann.h"
#include "exact.h"

using namespace std;

void exact(VVD const &s,VVD const &x,int const &test_problem,double const &time){

// file handle for output

  ofstream f1;

// establish the mesh limits

  double xmin(*min_element(x.at(0).begin(),x.at(0).end()));
  double xmax(*max_element(x.at(0).begin(),x.at(0).end()));
  double ymin(*min_element(x.at(1).begin(),x.at(1).end()));
  double ymax(*max_element(x.at(1).begin(),x.at(1).end()));

// decalre the lineout structure

  struct lineout_type {
    double x1,y1; // start point of each line
    double x2,y2; // end point of each line
    string filename; //filename to output
    string filehead; // file header
    int nsamples; // number of sample points on each line
  } lineout;

  vector<lineout_type> Lineout;

// set up an exact solution for each problems that has one

  switch(test_problem){

    case(TAYLOR):

// Taylor Green vortex

      return;;

      break;

    case(RAYLEIGH):

// Rayleigh-Taylor instability

      return;

      break;

    case(NOH):

// Noh stagnation shock

      return;;

      break;

    case(SEDOV):

// Sedov expanding shock

      return;;

      break;

    case(SOD):

// Sod's shock tube

      lineout.x1=xmin;
      lineout.x2=xmax;
      lineout.y1=0.5*(ymin+ymax);
      lineout.y2=lineout.y1;
      lineout.filename="exact.dat";
      lineout.filehead="# Sod exact solution from (0.0,0.5) to (1.0,0.5) : Columns are x d p e u";
      lineout.nsamples=NSAMPLES;

      {

        double l[4]={s.at(0)[0],s.at(0)[1],s.at(0)[3],s.at(0)[4]}; // left flux state
        double r[4]={s.at(1)[0],s.at(1)[1],s.at(1)[3],s.at(1)[4]}; // right flux state
        Riemann R(Riemann::exact,l,r); // start Riemann solver from initial flux states
        vector<double> rx;vempty(rx); // sample point coordinates for the Riemann solver
//        for(int i=0;i<NSAMPLES;i++){rx.push_back(0.0+(i*(1.0-0.0)/double(NSAMPLES)));} // sample points for the Riemann solver
        for(int i=0;i<NSAMPLES;i++){rx.push_back(i*(xmax-xmin)/NSAMPLES);} // sample points for the Riemann solver
        R.profile(&rx,time); // Riemann solution at NSAMPLE sample points

// append to the data structure

        Lineout.push_back(lineout);

        cout<<"exact(): exact solution writing to file "<<Lineout.at(0).filename<<" ..."<<endl;

// open the output file for the exact solution and write the header part

        f1.open(Lineout.at(0).filename);
        f1<<Lineout.at(0).filehead<<endl;

// output the exact solution

        f1<<fixed<<setprecision(10);
        for(int j=0;j<rx.size();j++){
          f1<<rx.at(j)<<" "<<R.density(j)<<" "<<R.pressure(j)<<" "<<R.energy(j)<<" "<<R.velocity(j)<<endl;
        }

// close the output file

        f1.close();

      }

      break;


    case(R2R):

// 123 problem

      lineout.x1=xmin;
      lineout.x2=xmax;
      lineout.y1=0.5*(ymin+ymax);
      lineout.y2=lineout.y1;
      lineout.filename="exact.dat";
      lineout.filehead="# 123 problem solution from (0.0,0.5) to (1.0,0.5) : Columns are x d p e u";
      lineout.nsamples=NSAMPLES;

      {

        double l[4]={s.at(0)[0],s.at(0)[1],s.at(0)[3],s.at(0)[4]}; // left flux state
        double r[4]={s.at(1)[0],s.at(1)[1],s.at(1)[3],s.at(1)[4]}; // right flux state
        Riemann R(Riemann::exact,l,r); // start Riemann solver from initial flux states
        vector<double> rx;vempty(rx); // sample point coordinates for the Riemann solver
        for(int i=0;i<NSAMPLES;i++){rx.push_back(i*(xmax-xmin)/NSAMPLES);} // sample points for the Riemann solver
        R.profile(&rx,time); // Riemann solution at NSAMPLE sample points

// append to the data structure

        Lineout.push_back(lineout);

        cout<<"exact(): exact solution writing to file "<<Lineout.at(0).filename<<" ..."<<endl;

// open the output file for the exact solution and write the header part

        f1.open(Lineout.at(0).filename);
        f1<<Lineout.at(0).filehead<<endl;

// output the exact solution

        f1<<fixed<<setprecision(10);
        for(int j=0;j<rx.size();j++){
          f1<<rx.at(j)<<" "<<R.density(j)<<" "<<R.pressure(j)<<" "<<R.energy(j)<<" "<<R.velocity(j)<<endl;
        }

// close the output file

        f1.close();

      }

      break;

    case(SALTZMANN):

// Saltzmann piston

      return;

      break;

  }

  return;

}

// function to empty a vector

void vempty(vector<double>&v){vector<double> e;v.swap(e);return;}
