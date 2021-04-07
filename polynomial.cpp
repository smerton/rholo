// Function definitions for members of the polynomial class
//
// This generates a degree-n polynomial Pn(x) through nloc nodes along x.
// The polynomial has a value =1.0 at node i1 (the compact support node)
// and a value =0.0 at every other node. This ensures orthogonality
// allowing a finite element of arbitary order to be constructed from
// this polynomial.
//
// Usage: declare a polynomial at each node in the element and convolve
// to generate 1D, 2D or 3D finite element spaces of arbitrary polyhedral
// order

// Author S. R. Merton

#include "polynomial.h"
#include "matrix.h"
#include <iostream>
#include <fstream>
#include <math.h>
#include <iomanip>      // std::setprecision
#include <vector>

#define NPOINTS 1000

using namespace std;

Polynomial::Polynomial(int n, int i1, double x1, double x2){

// Constructor to instantiate new degree-n polynomial object Pn(x)
// on the range x1<=x<x2 with a value Pn(x)=1.0 at x=x_i1
// This works by dividing the range (x1,x2) into n partitions
// with n+1 evenly spaced nodes. Note it is convenient to address
// the n+1 nodes (0,n) in this code, rather than (1,n+1).
// Node 0 is at x=x1 and node n at x=x2.
// A system of n+1 simlutaneous equations for the polynomial are
// placed in a matrix equation Ac=b and solved for the vector of
// coeficients c. The load vector b in this matrix equation is the
// value of the Pn(x) at each of the n+1 nodes, which is zero at all
// n+1 nodes except the node i1 at which the value is 1.0. Node i1
// is the node selected as the node at which Pn(x) has compact
// support.
//
// Arguments:
// n is type int containing the degree of the polynomial, the choice is entirely arbitrary for n>=0
// i1 is type int containing a local node at which the polynomial has compact support (Pn(x_i1)=1.0)
//    note i1 must lie in the range 0<=i1<n
// x1 is type double and is the x coordinate of the left endpoint of the polynomial
// x2 is type double and is the x coordinate of the right endpoint of the polynomial

// set the number of nodes in the domain of the polynomial

  nloc=n+1;

// exception

  if(i1<0 || i1>n){
    cout<<"ERROR in Polynomial::Polynomial(): Node "<<i1+1<<" not in range (1,"<<nloc<<")"<<endl;
    exit(1);
  }

// set the range of the polynomial

  SetRange(x1,x2);

// set the x value of each node

  double x[nloc];
  for(int i=0;i<nloc;i++){
    x[i]=x1+i*(x2-x1)/n;
  }

// select node i1 to have compact support

  double b[nloc];
  for(int i=0;i<nloc;i++){
    b[i]=0.0;
  }
  b[i1]=1.0;


// place the simultaneous equations into a matrix and solve for their coefficients

  Matrix M(nloc);

// populate the matrix with the simultaneous equations

  for(int iloc=0;iloc<nloc;iloc++){
    for(int jloc=0;jloc<nloc;jloc++){
      M.write(iloc,jloc,pow(x[iloc],jloc));
    }
  }

// solve for the coefficients

  double soln[nloc];
  M.solve(soln,b);

// commit coefficients to the class address space

  for(int iloc=0;iloc<nloc;iloc++){
    coef.push_back(soln[iloc]);
  }

  return;

}

double Polynomial::x(double u, double umin, double umax){

// Member function to return the coordinate xmin<=x(u)<=xmax in the range
// of the polynomial object, given some other coordinate u on a range (umin,umax)
// that differs from the range (xmin,xmax) of the polynomial.

// Arguments:
// u is type double containing an arbitrary coordinate in the range (umin,umax)
// umin is type double containing the left endpoint of the domain passed in
// umax is type double containing the right endpoint of the domain passed in
//
// Returns:
// type double containing a coordinate in the range (xmin,xmax) of the polynomial

  return xmin+((u-umin)/(umax-umin))*(xmax-xmin);
}

bool Polynomial::check(double x,double xl,double xr){

// Member function to check x is within a given range (xl,xr)
// i.e. to check a given coordinate lies within the endpoints xl,xr
//
// Arguments:
// x is type double containing a coordinate
// xl is type double containing the left endpoint
// xr is type double containing the right endpoint
//
// Returns:
// type bool which is true if x is within range (xl,xr)

  return (x>=xl && x<=xr);
}

double Polynomial::value(double xval){

// Member function to return the value of P(xval), this simply
// evaluates P(x) at the coordinate x=xval
//
// Arguments:
// xval of type double containing a value of x at which to evaluate the polynomial
//
// Returns:
// type double containing the value of the polynomial at the given value of x or 0.0
// if xval is not on the range (xmin,xmax) of the polynomial.

  double pval=0.0;

  if(check(xval,xmin,xmax)){

    for(int i=0;i<nloc;i++){
      pval+=coef[i]*pow(xval,i);
    }

  }

  return pval;

}

double Polynomial::dvalue(double xval){

// Member function to evaluate the first derivative dPn(x)/dx at x=xval
//
// Arguments:
// xval of type double containing a value of x at which to evaluate the derivative dPn(x)/dx
//
// Returns:
// type double containing the value of the first derivative at x=xval

  double dpval=0.0;

  if(check(xval,xmin,xmax)){

    for(int i=1;i<nloc;i++){
      dpval+=double(i)*coef[i]*pow(xval,i-1);
    }

  }

  return dpval;

}

void Polynomial::draw(char*myfilename){

// Member function to draw the polynomial. This outputs the value
// of Pn(x) at a large number of points on the range (xmin,xmax)
// in a format convenient for plotting. This is useful for
// visualising Pn(x) and for debugging.
//
// Arguments:
// myfilename is a pointer to a char containing the output filename
//

  double x;

// open the file

  ofstream myfile(myfilename); 

// evaluate Pn(x) at NPOINTS points between xmin and xmax

  for(int j=0;j<=NPOINTS;j++){

// x value of point j

    x=xmin+j*(xmax-xmin)/NPOINTS;

// output Pn(x)

    myfile<<fixed<<setprecision(12)<<x<<" "<<value(x)<<" "<<xmin*0.0<<endl;

  }

  myfile.close();

  return;

}

inline void Polynomial::SetRange(double new_xmin, double new_xmax){

// Member function to set the range of the polynomial
//
// Arguments:
// new_xmin type int containing the x coordinate of the required left endpoint
// new_xmax type int containing the x coordinate of the required right endpoint

  xmin=new_xmin;
  xmax=new_xmax;

  return;
}

// Destructor function to release storage associated with a Polynomial class object

Polynomial::~Polynomial(){}
