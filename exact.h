// Signatures for the exact solution to test problems where applicable

// Author S. R. Merton

#define NSAMPLES 1000         // number of sample points for the exact solution

void vempty(vector<double>&v);                                                       // empty a vector
void exact(VVD const &s,VVD const &x,int const &test_problem,double const &time);    // exact solution at t=time
