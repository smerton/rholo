// Export the signature of the polynomial class
//
// Arguments
// n: degree of the polynomial
// i1: node at which the polynomial has compact support (see description)
//
// Description
// This class generates a polynomial Pn(x) of a given degree n
// across the element xmin<=x<=xmax
// The element domain (xmin,xmax) is divided into n+1 nodes, each one of
// these is a root of the polynomial except node i1 at which point
// Pn(x_i1)=1.0. This is the node at which Pn(x) has compact support
// in a finite element sense since this polynomial is intended
// for constructing finite element bases in 1D, 2D or 3D and it is
// a condition of orthogonality that the function has a value =1.0 at one
// of the nodes in the element and a value=0.0 at every other node.
//
// The polynomial is defined as Pn(x)=sum(C_i * x^i) i=0,1,2,...,n
// This is written in matrix form as a set of simultaneous equations
// See the source file for details of how this has been implemented

#include <vector>

using namespace std;

class Polynomial{
  public:
    Polynomial(int n,int i1,double x1,double x2); // constructor for a new polynomial of degree n with value Pn(x_i1)=1 on range (x1,x2)
    ~Polynomial();              // destructor function to release class storage
    double value(double xval);  // accessor function to evaluate Pn(x)
    double dvalue(double xval); // accessor function to evaluate d/dx(Pn(x))
    double x(double u,double umin,double umax); // find xmin<=x(u)<=xmax given umin<=u<=umax
    void draw(char*myfilename);// member function to draw Pn(x)
  private:
    vector<double> coef;       // coefficients of the polynomial
    int nloc;                  // number of nodes in the element (=n+1)
    double xmin=-1.0;          // minimum value of x
    double xmax=1.0;           // maximum value of x
    bool check(double x,double xl,double xr); // check x is with the range xl,xr
    void SetRange(double new_xmin,double new_xmax); // set the value of xmin and xmax
};
