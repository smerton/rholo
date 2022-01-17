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

// set the segment end point coordinates

  mcoord.at(0).push_back(r2.at(0));
  mcoord.at(1).push_back(r2.at(1));

// intercept

  if(start(0)!=0.0){
    double theta(atan(m())),r(abs(cos(theta)/start(0)));
    mc=(m()>=0.0)?start(1)-(r*sin(theta)):start(1)+(r*sin(theta));
  }else{
    mc=start(1);
  }

}


// function to return coordinate idim of the start

double Line::start(int idim) const {return mstart.at(idim);}

// function to return coordinate idim of the end

double Line::end(int idim) const {return mend.at(idim);}

// function to return coordinate idim of the gradient

double Line::m(int idim) const {return mm.at(idim);}

// function to return the gradient

double Line::m() const {return mm.at(1)/(sgn(mm.at(0))*max(1.0e-10,abs(mm.at(0))));}

// function to return the length of the line

double Line::length() const {return mlength;}

// function to return segment i end point coordinate idim

double Line::coord(int idim,int i) const {return mcoord.at(idim).at(i);}

// function to divide the line into n segments

void Line::divide(int n){

// set the number of segment

  mnsegments=n;

// set segment number of intersection, only relevent if intersection exists between two lines in which case it will be updated

  msegint=-1;

// set the end point of each segment

  vector<vector<double> > vtmp(2);

  for(int i=1;i<=nsegments();i++){
    double ri(i*length()/nsegments()); // distance from line origin to segment end point
    double theta(atan(m()));
    vtmp.at(0).push_back(start(0)+ri*cos(theta)); // coordinate 0 of the segment end point
    vtmp.at(1).push_back(start(1)+ri*sin(theta)); // coordinate 1 of the segment end point
  }

// replace coordinate vector

  mcoord=vtmp;

  return;

}

// function to return the number of segments

int Line::nsegments() const {return mnsegments;}

// function to return the intercept

double Line::c() const {return mc;}

// function to check if the line intersects the line l1 passed in

bool Line::intersects(Line L1) {

  for(int iseg=0;iseg<this->nsegments();iseg++){

    vector<vector<double> > s0(2),s1(2);

// place end points of the two lines in these vectors

//    s0.at(0).push_back(this->start(0));
//    s0.at(0).push_back(this->start(1));
//    s0.at(1).push_back(this->end(0));
//    s0.at(1).push_back(this->end(1));

    s0.at(0).push_back(this->start(0)); // origin of line
    s0.at(0).push_back(this->start(1)); // origin of line
    s0.at(1).push_back(this->coord(0,iseg)); // end of current segment
    s0.at(1).push_back(this->coord(1,iseg)); // end of current segment

    s1.at(0).push_back(L1.start(0));
    s1.at(0).push_back(L1.start(1));
    s1.at(1).push_back(L1.end(0));
    s1.at(1).push_back(L1.end(1));

// take dot product of 2 cross-products

    double dx0 = s0[1][0]-s0[0][0];
    double dx1 = s1[1][0]-s1[0][0];
    double dy0 = s0[1][1]-s0[0][1];
    double dy1 = s1[1][1]-s1[0][1];
    double p0 = dy1*(s1[1][0]-s0[0][0]) - dx1*(s1[1][1]-s0[0][1]);
    double p1 = dy1*(s1[1][0]-s0[1][0]) - dx1*(s1[1][1]-s0[1][1]);
    double p2 = dy0*(s0[1][0]-s1[0][0]) - dx0*(s0[1][1]-s1[0][1]);
    double p3 = dy0*(s0[1][0]-s1[1][0]) - dx0*(s0[1][1]-s1[1][1]);

    bool z1((p0*p1<=0)&&(p2*p3<=0));

// if segment iseg crosses L1 store it and return true

    if(z1){
      msegint=iseg;
      return z1;
    }

  }

  return false;

}

// function to return the segment numebr at the intersection with another line

int Line::segint() const {return msegint;}

// type safe function to return the sign of the argument

template <typename T> int sgn(T val) {return( (val>=T(0))?T(1):T(-1));} // -1 or 1

// destructor function to release storage associated with a Line class object

Line::~Line(){}
