// Signature of the Riemann class containing the Riemann solvers
// prefix m denotes member data

// Author S. R. Merton

#include <vector>

using namespace std;

class Riemann{

  public:

  enum solver{pvrs=0,exact,UNKOWN};  // available solver types
  char*solvername[3] ={"pvrs","exact","not implemented"}; // their names

  Riemann(solver s,double*l,double*r); // constructor function
  ~Riemann(); // destructor to release class storage

  double pstar; // pressure in the star region
  double ustar; // velocity in the star region

  void profile(vector<double>*x,double t); // profiles at (x,t)

  double density(long i); // density field at position i
  double pressure(long i); // pressure field at position i
  double energy(long i); // energy field at position i
  double velocity(long i); // energy field at position i
  double region(long i); // region field at position i

  void sample(double s,double &p,double &u,double &d,double &e,double &r); // sample function to construct profile

  private:

  double mPl; // left state pressure field
  double mDl; // left state density field
  double mul; // left state velocity field
  double mgl; // left state adiabatic constant

  double mPr; // right state pressure field
  double mDr; // right state density field
  double mur; // right state velocity field
  double mgr; // right state adiabatic constant

  double cl,cr; // sound speeds left and right

  vector<double> mSCoords; // sample coordinates

  vector<double> mdensity; // density field
  vector<double> mpressure; // pressure field
  vector<double> menergy; // energy field
  vector<double> mvelocity; // velocity field
  vector<double> mregion; // region field 

// solvers

  void pvrs_solver(); // primitive variable (acoustic wave) approximation (pp 280, Ch. 9 Toro)
  void exact_solver(); // exact Riemann solver for ideal gases (pp 152, Ch.4 Toro)
  double starte(); // returns initial guess for the pressure iteration
  void prefun(double*f,double*fd,double p,double dk,double pk,double ck); // computes pressure function

// gammas

  double g0;
  double g1;
  double g2;
  double g3;
  double g4;
  double g5;
  double g6;
  double g7;
  double g8;
  double g9;

};
