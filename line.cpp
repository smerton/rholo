// Function definitions for the line class

// Author S. R. Merton

#include <iostream>
#include <vector>
#include "line.h"
//#include <cmath> // for sqrt()
#include <math.h>

template <typename T> int sgn(T val); // return type safe sign of the argument

using namespace std;

// constructor to instantiate a new Line object and all its attributes

Line::Line(vector<double> r1,vector<double> r2){

// set the end points

  mstart=r1;
  mend=r2;

// set the gradient

  for(int i=0;i<r1.size();i++){
    mm.push_back(r2.at(i)-r1.at(i));
  }

// set the length, this is just the 2-norm of the gradient

  mlength=0.0;

  for(int i=0;i<mm.size();i++){
    mlength+=mm.at(i)*mm.at(i);
  }

  mlength=sqrt(mlength);

// set a default number of segments

  mnsegments=1;

// allocate the coordinates of the segment end points

  mcoord.resize(2);

// set the segmemnt end point coordinates

  mcoord.at(0).push_back(r2.at(0));
  mcoord.at(1).push_back(r2.at(1));

}


// function to return coordinate idim of the start

double Line::start(int idim) const {return mstart.at(idim);}

// function to return coordinate idim of the end

double Line::end(int idim) const {return mend.at(idim);}

// function to return coordinate idim of the gradient

double Line::m(int idim) const {return mm.at(idim);}

// function to return the gradient

double Line::m() const {return mm.at(1)/(sgn(mm.at(0))*max(1.0e-10,mm.at(0)));}

// function to return the length of the line

double Line::length() const {return mlength;}

// function to return segment i end point coordinate idim

double Line::coord(int idim,int i) const {return mcoord.at(idim).at(i);}

// function to divide the line into n segments

void Line::divide(int n){

// set the number of segment

  mnsegments=n;

// set the end point of each segment

  vector<vector<double> > vtmp(2);

  for(int i=1;i<=nsegments();i++){
    double ri(i*length()/nsegments()); // distance from line origin to segment end point
    double theta(atan(m()));
    vtmp.at(0).push_back(start(0)+ri*cos(theta)); // coordinate 0 of the segment end point
//    cout<<"Line::divide(): theta= "<<theta<<" m()= "<<m()<<endl;
    vtmp.at(1).push_back(start(1)+ri*sin(theta)); // coordinate 1 of the segment end point
  }

// replace coordinate vector

  mcoord=vtmp;

  return;

}

// function to return the number of segments

int Line::nsegments() const {return mnsegments;}

// type safe function to return the sign of the argument

template <typename T> int sgn(T val) {return( (val>=T(0))?T(1):T(-1));} // -1 or 1

// destructor function to release storage associated with a Line class object

Line::~Line(){}
