// Export the signature of the quadrature rule used in numerical integration
//
// This was taken from the following reference:
// https://people.sc.fsu.edu/~jburkardt/cpp_src/legendre_rule/legendre_rule.cpp

#include <cstdlib>
#include <cmath>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <ctime>
#include <cstring>

using namespace std;

class QuadratureRule{
  public:
    QuadratureRule(int n);  // constructor for a new qadrature rule with n weights
    ~QuadratureRule(); // destructor to release class storage
    int npoints; // number of points
    double*x; // 1-D array to hold the abscissae
    double*w; // 1-D array to hold the weights
    string name; // name of the quadrature set
    double GetL(); // accessor function to return left endpoint
    double GetR(); // accessor function to return right endpoint
  private:
    double const a=-1.0; // left endpoint
    double const b=1.0; // right endpoint
    double alpha=0.0; // alpha parameter
    double beta=0.0; // beta parameter
    void cdgqf ( int nt, int kind, double alpha, double beta, double t[], double wts[] );
    void cgqf ( int nt, int kind, double alpha, double beta, double a, double b, double t[], double wts[] );
    double class_matrix ( int kind, int m, double alpha, double beta, double aj[], double bj[] );
    void imtqlx ( int n, double d[], double e[], double z[] );
    void parchk ( int kind, int m, double alpha, double beta );
    double r8_epsilon ( );
    double r8_sign ( double x );
    void r8mat_write ( string output_filename, int m, int n, double table[] );
    void rule_write ( int order, string filename, double x[], double w[], double r[] );
    void scqf ( int nt, double t[], int mlt[], double wts[], int nwts, int ndx[], double swts[], double st[], int kind, double alpha, double beta, double a, double b );
    void sgqf ( int nt, double aj[], double bj[], double zemu, double t[], double wts[] );
    void timestamp ( );
};
