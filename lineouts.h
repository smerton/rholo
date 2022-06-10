// Signatures for lineout functions

// Author S. R. Merton

#define VD vector<double> // vector of doubles
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
#include <iomanip>   // floating point precision
#include <cmath>     // for sqrt
#include <fstream>   // for file io
#include <algorithm> // min_element, max_element
#include "eos.h"     // eos lookups
#include "line.h"
#include "mesh.h"
#include "shape.h"

// top level function to choose between 1D and 2D lines

void lineouts(Mesh const &M,Shape const &S,Shape const &T,VD const &dinit,VD const &e,VVD const &xinit,VVD const &x,
              VVD const &xt,VVD const &u,int const &test_problem,vector<int> const &mat,VD const &g);

// function to code for 1D lineouts

void lineouts_1d(Mesh const &M,Shape const &S,Shape const &T,VD const &dinit,VD const &e,VVD const &xinit,VVD const &x,
              VVD const &xt,VVD const &u,int const &test_problem,vector<int> const &mat,VD const &g);

// function to code for 2D lineouts

void lineouts_2d(Mesh const &M,Shape const &S,Shape const &T,VD const &dinit,VD const &e,VVD const &xinit,VVD const &x,
              VVD const &xt,VVD const &u,int const &test_problem,vector<int> const &mat,VD const &g);
