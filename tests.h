// Signatures of various test problems and the functions to set them up and produce lineouts and exact solutions

// Author S. R. Merton

#define TAYLOR 1              // Taylor-Green vortex problem
#define RAYLEIGH 2            // Rayleigh-Taylor instability problem
#define NOH 3                 // Noh stagnation shock problem
#define SEDOV 4               // Sedov expanding shock problem
#define TRIPLE 5              // Triple point problem
#define SOD 6                 // Sod's shock tube problem
#define R2R 7                 // 123 problem
#define SALTZMANN 8           // Saltzmann piston problem

#define NSAMPLES 1000         // number of sample points for Riemann solutions

// input overides for the Taylor-Green vortex

void init_TAYLOR(Mesh const &M,Shape const &S,Shape const &T,double const &dpi,VD &dinit,
                 VVD &u0,VVD &u1,VD &e0,VD &e1,VVD const &xk,VVD const &xt,VD const &gamma,
                 vector<int> const &mat,VD &detJ0,VVVD &detDJ0,VD &detJ,
                 VVVD &detDJ,VD const &m);

// input overides for the Rayleigh-Taylor instability

void init_RAYLEIGH(Mesh const &M,Shape const &S,Shape const &T,double const &dpi,VD &dinit,
                   VVD &u0,VVD &u1,VD &e0,VD &e1,VVD const &xk,VVD const &xt,VD const &gamma,
                   vector<int> const &mat,VD &detJ0,VVVD &detDJ0,VD &detJ,
                   VVVD &detDJ,VD const &m);

// input overides for the Noh stagnation shock

void init_NOH(Mesh const &M,Shape const &S,Shape const &T,double const &dpi,VD &dinit,
              VVD &u0,VVD &u1,VD &e0,VD &e1,VVD const &xk,VVD const &xt,VD const &gamma,
              vector<int> const &mat,VD &detJ0,VVVD &detDJ0,VD &detJ,
              VVVD &detDJ,VD const &m);

// input overides for the Sedov explosion

void init_SEDOV(Mesh const &M,Shape const &S,Shape const &T,double const &dpi,VD &dinit,
                VVD &u0,VVD &u1,VD &e0,VD &e1,VVD const &xk,VVD const &xt,VD const &gamma,
                vector<int> const &mat,VD &detJ0,VVVD &detDJ0,VD &detJ,
                VVVD &detDJ,VD const &m);

// input overides for the Saltzmann piston

void init_SALTZMANN(Mesh const &M,Shape const &S,Shape const &T,double const &dpi,VD &dinit,
                VVD &u0,VVD &u1,VD &e0,VD &e1,VVD &xk,VVD &xt,VD const &gamma,
                vector<int> const &mat,VD &detJ0,VVVD &detDJ0,VD &detJ,
                VVVD &detDJ,VD const &m);

// exact solution at t=time for various problems

void exact(VVD const &s,VVD const &x,int const &test_problem,double const &time);

// top level function to choose between 1D and 2D lines

void lineouts(Mesh const &M,Shape const &S,Shape const &T,VD const &dinit,VD const &e,VVD const &xinit,VVD const &x,
              VVD const &xt,VVD const &u,int const &test_problem,vector<int> const &mat,VD const &g);

// function to code for 1D lineouts

void lineouts_1d(Mesh const &M,Shape const &S,Shape const &T,VD const &dinit,VD const &e,VVD const &xinit,VVD const &x,
              VVD const &xt,VVD const &u,int const &test_problem,vector<int> const &mat,VD const &g);

// function to code for 2D lineouts

void lineouts_2d(Mesh const &M,Shape const &S,Shape const &T,VD const &dinit,VD const &e,VVD const &xinit,VVD const &x,
              VVD const &xt,VVD const &u,int const &test_problem,vector<int> const &mat,VD const &g);

// impose a manufactured solution on the internal energy field

void manufactured_soln_ie(Mesh const &M,Shape const &S,Shape const &T,VVD const &x,int const &test_problem,VD &detJ,VVVD &detDJ,VD &b);
