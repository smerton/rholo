// Variant of RhoLo (Really High Order Lagrangian Operator - RhoLo)
// RhoLo is an experimental high-order 2D pure Lagrangian hydrodynamics test code used for
// making assessments of high-order methods and developing an implementation strategy for
// high-order methods.
// This code solves the Euler equations in their non-conservative form in the fluid frame (the Lagrangian frame)
// using a high-order finite element method (discontinuous thermodynamic variables p,rho,e and node 
// centred kinematic variables x,u,a) and bulk viscosity q to increase entropy across element boundaries, initial 
// implementation is only first order in time
//
// Author S. R. Merton
//
// how to set up push and pull from multiple remotes (e.g. home NAS and github):
//
// git remote add NAS ssh://smerton@192.168.1.79/shares/git/rholo.git
// git remote set-url --add NAS git@github.com:smerton/rholo.git
// git remote add github git@github.com:smerton/rholo.git
// git push -vu NAS finite-element
//
// then check url's have been added with something like
//
// git remote -v
//
// and also check .git/config
//
// then git pull --all will pull from 1st url in location 1 (NAS)
// and 1st url in location 2 (github)
//
// for graphics: convert -density 300 filename.png filename.pdf
//

#define DTSTART 0.001     // insert a macro for the first time step
#define ENDTIME 0.15      // insert a macro for the end time
#define ECUT 1.0e-8       // cut-off on the energy field
#define NSAMPLES 1000     // number of sample points for the exact solution
//#define VISFREQ 200     // frequency of the graphics dump steps
//#define OUTFREQ 50      // frequency of the output print steps
#define VISFREQ 0.02      // frequency of the graphics dump times
#define OUTFREQ 0.01      // frequency of the output print times
#define VD vector<double> // vector of doubles
#define VVD vector<VD>    // vector of VD
#define VVVD vector<VVD>  // vector of VVD
#define VI vector<int>    // vector of ints
#define COURANT 0.333     // Courant number for CFL condition
#define DTSFACTOR 0.1     // safety factor on time-step control
#define NROWS nknodes     // number of rows in the global matrix
#define NCOLS nknodes     // number of columns in the global matrix
#define NGI S.ngi()*n     // number of integration points on the mesh
#define GPNT i*T.ngi()+gi // global address of integration point gi
#define ROW M.GlobalNode_CFEM(i,iloc)  // row address in global matrix
#define COL M.GlobalNode_CFEM(i,jloc)  // column address in global matrix

#define VACUUM 1              // vacuum boundary
#define REFLECTIVE 2          // reflective boundary
#define VELOCITY 3            // velocity v.n applied to boundary
#define FLUID 4               // fluid on the boundary
#define ACCELERATION 5        // velocity a.n applied to boundary

#define TAYLOR 1              // Taylor-Green vortex problem
#define RAYLEIGH 2            // Rayleigh-Taylor instability problem
#define NOH 3                 // Noh stagnation shock problem
#define SEDOV 4               // Sedov expanding shock problem
#define TRIPLE 5              // Triple point problem
#define SOD 6                 // Sod's shock tube problem
#define R2R 7                 // 123 problem

#include <iostream>
#include <vector>
#include <iomanip>
#include <cmath>
#include <fstream>   // for file io
#include "matrix.h"  // matrix operations
#include "shape.h"   // signature of the shape class
#include <bits/stdc++.h>
#include "eos.h"     // eos lookups
#include "mesh.h"    // signature of the mesh class
#include <ctime>     // date and time
#include <stdlib.h>  // getenv
#include <limits.h>
#include <stdio.h>
#include <unistd.h>
#include "line.h"
#include "riemann.h"

// function signatures

string date();                                                                                                       // function to return the date and time
void header();                                                                                                       // header part
void vempty(vector<double>&v);                                                                                       // empty a vector
void jacobian(int const &i,VVD const &x,Mesh const &M,Shape const &S,VD &detJ,VVVD &detDJ);                          // calculate a jacobian and determinant
void jacobian(int const &i,VVD const &x0,VVD const &x,Mesh const &M,Shape const &S,VD &detJ,VVVD &Js);               // calculate a jacobian for the Lagrangian motion
void jacobian(int const &i,VVD const &x,Mesh const &M,Shape const &S,VD &detJ);                                      // calculate a jacobian at the local nodes
void sum_ke(double &ke,VVD const &u,VD const &dinit,Mesh const &M,VVD const &xinit,                                  // sum the global kinetic energy field
            VVD const &x,Shape const &S,Shape const &T,VD &detJ0,VVVD &detDJ0,VD &detJ,VVVD &detDJ);
void sum_ie(double &ie,VD const &e,VD const &dinit,Mesh const &M,VVD const &xinit,                                   // sum the global internal energy field
            VVD const &x,Shape const &S,Shape const &T,VD &detJ0,VVVD &detDJ0,VD &detJ,VVVD &detDJ);
void initial_data(int const n, long const nknodes,long const ntnodes,Shape const S,                                  // echo some initial information
                  int const ndims, int const nmats,Mesh const &M);
void lineouts(Mesh const &M,Shape const &S,Shape const &T,VD const &dinit,VD const &e,VVD const &xinit,VVD const &x, // line-outs
              VVD const &xt,VVD const &u,int const &test_problem,vector<int> const &mat,VD const &g);
void exact(VVD const &s,VVD const &x,int const &test_problem);                                                       // exact solution where applicable
void silo(VVD const &x,VVD const &xt,VVD const &xinit,VD const &d,VD const &l,VD const &V,VD const &e,               // silo graphics output
          VVD const &u,VI const &mat,int s, double t,Mesh const &M,VD const &g,Shape const &S,Shape const &T);
void state_print(int const n,int const ndims, int const nmats, VI const &mat,                                        // output material states
                  VD const &d, VD const &V, VD const &m, VD const &e, VD const &p, 
                  VVD const &x, VVD const &u, int const &s, double const &t,VD const &gamma);
void bc_insert(Matrix &A,Mesh const &M,Shape const &S,VD const &d,VD const &detJ,                                    // insert boundary conditions in to the mass matrix via row elimination
               VVD &u0,VVD &u1,VD &b0,VD &b1,long const &nnodes);
void init_TAYLOR(Mesh const &M,Shape const &S,Shape const &T,double const &dpi,VD &dinit,                            // input overides for the Taylor-Green vortex
                 VVD &u0,VVD &u1,VD &e0,VD &e1,VVD const &x,VD const &gamma,
                 vector<int> const &mat,VD &detJ0,VVVD &detDJ0,VD &detJ,
                 VVVD &detDJ,VD const &m);
void init_RAYLEIGH(Mesh const &M,Shape const &S,Shape const &T,double const &dpi,VD &dinit,                          // input overides for the Rayleigh-Taylor instability
                   VVD &u0,VVD &u1,VD &e0,VD &e1,VVD const &x,VD const &gamma,
                   vector<int> const &mat,VD &detJ0,VVVD &detDJ0,VD &detJ,
                   VVVD &detDJ,VD const &m);
void init_NOH(Mesh const &M,Shape const &S,Shape const &T,double const &dpi,VD &dinit,                               // input overides for the Noh stagnation shock
              VVD &u0,VVD &u1,VD &e0,VD &e1,VVD const &x,VD const &gamma,
              vector<int> const &mat,VD &detJ0,VVVD &detDJ0,VD &detJ,
              VVVD &detDJ,VD const &m);
void init_SEDOV(Mesh const &M,Shape const &S,Shape const &T,double const &dpi,VD &dinit,                             // input overides for the Sedov explosion
                VVD &u0,VVD &u1,VD &e0,VD &e1,VVD const &x,VD const &gamma,
                vector<int> const &mat,VD &detJ0,VVVD &detDJ0,VD &detJ,
                VVVD &detDJ,VD const &m);
template <typename T> int sgn(T val);                                                                                // return type safe sign of the argument

using namespace std;

int main(){

// global data

  Mesh M("mesh/r2r-100x1.mesh");                                  // load a new mesh from file
  Shape S(2,3,CONTINUOUS);                                       // load a shape function for the kinematics
  Shape T(1,sqrt(S.ngi()),DISCONTINUOUS);                        // load a shape function for the thermodynamics
  ofstream f1,f2,f3;                                             // files for output
  int const n(M.NCells()),ndims(M.NDims());                      // no. ncells and no. dimensions
  long const nknodes(M.NNodes(S.order(),S.type()));              // insert shape functions in to the mesh
  long const ntnodes(M.NNodes(T.order(),T.type()));              // insert shape functions in to the mesh
  int const nmats(M.NMaterials());                               // number of materials
  double cl,cq;                                                  // linear & quadratic coefficients for bulk viscosity (for weak (cl) and strong (cq) shock control)
  Matrix KMASS(2*NROWS),KMASSI(2*NROWS);                         // mass matrix for kinematic field
  vector<double> dinit(n),V0(n),V1(n),m(n);                      // density, volume & mass
  vector<double> e0(ntnodes),e1(ntnodes);                        // internal energy field
  vector<double> l(S.ngi()),d(S.ngi()),c(S.ngi());               // length scale, density and sound speed as functions
  vector<double> p(S.ngi()),q(S.ngi());                          // pressure and artificial viscosity as functions
  vector<vector<double> > u0(ndims),u1(ndims);                   // node velocity
  vector<vector<double> > x0(ndims),x1(ndims);                   // kinematic node coordinates
  vector<vector<double> > xt0(ndims),xt1(ndims);                 // thermodynamic node coordinates
  vector<vector<double> > xinit(ndims);                          // initial node coordinates
  vector<double> detJ0(S.ngi()),detJ(S.ngi()),detJs(S.ngi());    // determinant of jacobian at each integration point at time 0 and time-t
  vector<vector<vector<double> > > Js(S.ngi());                  // jacobian of the lagrangian motion to map from time-t to time-0
  vector<vector<vector<double> > > detDJ0(ndims),detDJ(ndims);   // determinant of jacobian for each derivative at time 0 and time-t
  vector<double> dt_cfl(NGI);                                    // time-step at each integration point
  vector<double> dts(2);                                         // time-step for each condition (0=CFL, 1=graphics hits)
  vector<int> mat(n);                                            // element material numbers
  vector<double> b0(M.NDims()*nknodes);                          // boundary conditions on the acceleration field
  vector<double> b1(M.NDims()*nknodes);                          // rows eliminated from the mass matrix
  vector<double> l0(n);                                          // initial length of each element
  double ls(0.0);                                                // element length scale at the quadrature point
  double ke(0.0),ie(0.0);                                        // kinetic and internal energy for conservation checks
  double time(0.0),dt(DTSTART);                                  // start time and time step
  int step(0);                                                   // step number
  int test_problem(0);                                           // set later for specific test cases that may need some overides
  double dpi(4.0*atan(1.0));                                     // definition of pi to double precision
  vector<double> gamma(M.NMaterials());                          // ratio of specific heats, set from material definition
  long nzeroes(0);                                               // number of non-zeroes in the CSR force matrix
  vector<double> F(nzeroes);                                     // force matrix in CSR format
  vector<long> frow(nzeroes);                                    // row addresses in the force matrix
  vector<long> fcol(nzeroes);                                    // column addresses in the force matrix
  vector<int> fdim(nzeroes);                                     // dimension addresses in the force matrix
  int length_scale_type(LS_LOCAL);                               // length scale definition: LS_AVERAGE,LS_LOCAL or LS_DIRECTIONAL (LS_PSEUDO_1D for 1D tests on 2D meshes)

// initial flux state in each material is in the form (d,ux,uy,p,gamma)

//  test_problem=SOD;length_scale_type=LS_PSEUDO_1D;cl=0.5;cq=4.0/3.0;       // set overides needed to run this problem
//  vector<vector<double> > state={{1.000, 0.000,0.000, 1.000,5.0/3.0},      // initial flux state in each material for Sod's shock tube 
//                                 {0.125, 0.000,0.000, 0.100,5.0/3.0}};

//  test_problem=SODSOD;length_scale_type=LS_PSEUDO_1D;cl=0.5;cq=4.0/3.0;    // set overides needed to run this problem
//  vector<vector<double> > state={{1.000, 0.000,0.000, 1.000,5.0/3.0},      // initial flux state in each material for double shock problem 
//                                 {0.125, 0.000,0.000, 0.100,5.0/3.0},
//                                 {1.000, 0.000,0.000, 1.000,5.0/3.0}};

  test_problem=R2R;length_scale_type=LS_PSEUDO_1D;cl=0.25;cq=2.0/3.0;      // set overides needed to run this problem
  vector<vector<double> > state={{1.000,-2.000,0.000, 0.400,1.4},          // initial flux state in each material for the 123 problem 
                                 {1.000, 2.000,0.000, 0.400,1.4}};

//  test_problem=BLASTWAVE;length_scale_type=LS_PSEUDO_1D;cl=0.5;cq=4.0/3.0; // set overides needed to run this problem
//  vector<vector<double> > state={{1.000,0.000,0.000, 1000.0,1.4},          // initial flux state in each material for the blast wave
//                                 {1.000t,0.000,0.000, 0.0100,1.4}};

//  test_problem=TAYLOR;length_scale_type=LS_AVERAGE;cl=0.0;cq=0.0;          // set overides needed to run this problem
//  vector<vector<double> > state={{1.000, 0.000,0.000, 1.000,5.0/3.0}};     // initial flux state in each material for Taylor problem

//  test_problem=NOH;length_scale_type=LS_LOCAL;cl=0.5;cq=4.0/3.0;             // set overides needed to run this problem
//  vector<vector<double> > state={{1.000, 0.000,0.000, 0.000,5.0/3.0}};       // initial flux state in each material for Noh problem

//  test_problem=SEDOV;length_scale_type=LS_AVERAGE;cl=0.5;cq=4.0/3.0;       // set overides needed to run this problem
//  vector<vector<double> > state={{1.000, 0.000,0.000, 1.000,1.4}};         // initial flux state in each material for Sedov problem

//  test_problem=TRIPLE;length_scale_type=LS_AVERAGE;cl=0.5;cq=4.0/3.0;      // set overides needed to run this problem
//  vector<vector<double> > state={{1.000, 0.000,0.000, 1.000,1.5},          // initial flux state in each material for triple-point problem
//                                 {1.000, 0.000,0.000, 0.100,1.4},
//                                 {0.125, 0.000,0.000, 0.100,1.5}};

// acquire adiabatic constant for each material

  for(int imat=0;imat<M.NMaterials();imat++){gamma.at(imat)=state[imat][4];}

// set boundary conditions on the edges of the mesh in the form (side,type,v.n) where side 0,1,2,3 = bottom,right,top,left

//  M.bc_set(0,VELOCITY,0.0);  // set boundary condition on bottom edge of mesh
//  M.bc_set(1,VELOCITY,0.0);  // set boundary condition on right edge of mesh
//  M.bc_set(2,VELOCITY,0.0);  // set boundary condition on top edge of mesh
//  M.bc_set(3,VELOCITY,0.0);  // set boundary condition on left edge of mesh
//  M.bc_set(0,VACUUM);  // set boundary condition on bottom edge of mesh
//  M.bc_set(1,VACUUM);  // set boundary condition on right edge of mesh
//  M.bc_set(2,VACUUM);  // set boundary condition on top edge of mesh
//  M.bc_set(3,VACUUM);  // set boundary condition on left edge of mesh
//  M.bc_set(0,VELOCITY,0.0);  // set boundary condition on bottom edge of mesh
//  M.bc_set(1,VACUUM);  // set boundary condition on right edge of mesh
//  M.bc_set(2,VACUUM);  // set boundary condition on top edge of mesh
//  M.bc_set(3,VELOCITY,0.0);  // set boundary condition on left edge of mesh
  M.bc_set(0,VELOCITY,0.0);  // set boundary condition on bottom edge of mesh
  M.bc_set(1,VELOCITY,2.0);  // set boundary condition on right edge of mesh
  M.bc_set(2,VELOCITY,0.0);  // set boundary condition on top edge of mesh
  M.bc_set(3,VELOCITY,-2.0);  // set boundary condition on left edge of mesh

// initialise the problem

  M.InitCoords(x0,S.order(),S.type());  // set initial coordinates for kinematic nodes
  x1=x0;
  xinit=x0;

  M.InitCoords(xt0,T.order(),T.type()); // set initial coordinates for thermodynamic nodes
  xt1=xt0;

  for(int i=0;i<M.NCells();i++){mat.at(i)=M.Material(i);}                                                                                             // load material numbers from the mesh class
  for(int i=0;i<M.NCells();i++){V0.at(i)=M.Volume(i);}                                                                                                // set initial element volume for start of step
  for(int i=0;i<M.NCells();i++){V1.at(i)=M.Volume(i);}                                                                                                // set initial element volume for end of step
  for(int i=0;i<M.NCells();i++){dinit.at(i)=state[mat.at(i)-1][0];}                                                                                   // load density field from initial flux state
  for(int i=0;i<M.NCells();i++){m.at(i)=dinit.at(i)*V0[i];}                                                                                           // calculate the initial mass field
  for(int i=0;i<M.NCells();i++){for(int j=0;j<T.nloc();j++){e0.at(M.GlobalNode_DFEM(i,j))=E(dinit.at(i),state[mat.at(i)-1][3],gamma[mat.at(i)-1]);}}  // invert the eos to start the energy field
  for(int i=0;i<e0.size();i++){e1.at(i)=e0.at(i);}
  for(int i=0;i<dt_cfl.size();i++){dt_cfl.at(i)=DTSTART;}                                                                                             // initial time-step

// initial element length

  M.InitLength(l0,S.order(),V0,length_scale_type);                                                                                                    // initialise element length scale

// allocate a determinant for each derivative

  for(int idim=0;idim<detDJ.size();idim++){
    detDJ0.at(idim).resize(S.nloc());
    detDJ.at(idim).resize(S.nloc());
    for(int j=0;j<S.nloc();j++){
      detDJ0.at(idim).at(j).resize(S.ngi());
      detDJ.at(idim).at(j).resize(S.ngi());
    }
  }

// allocate jacobian of the lagrangian motion

  for(int gi=0;gi<S.ngi();gi++){
    Js.at(gi).resize(ndims);
    for(int idim=0;idim<ndims;idim++){
      Js.at(gi).at(idim).resize(2);
    }
  }

// load velocity fields from initial flux state

  for(int idim=0;idim<ndims;idim++){

    vector<double> vtmp(x0.at(idim).size());

    for(int i=0;i<mat.size();i++){
      for(int iloc=0;iloc<S.nloc();iloc++){
        long j(M.GlobalNode_CFEM(i,iloc));
        vtmp.at(j)=(vtmp[j]==0.0)?(state[mat[i]-1][1+idim]):((vtmp[j]!=(state[mat[i]-1][1+idim]))?0.0:(state[mat[i]-1][1+idim]));
      }
    }

// commit to velocty address spaces

    u0.at(idim)=vtmp;
    u1.at(idim)=vtmp;

  }

// input overides needed to initialise certain test problems

  switch(test_problem){

    case(TAYLOR):

// Taylor Green vortex

      init_TAYLOR(M,S,T,dpi,dinit,u0,u1,e0,e1,xt0,gamma,mat,detJ0,detDJ0,detJ,detDJ,m);

      break;

    case(RAYLEIGH):

// Rayleigh-Taylor instability

      init_RAYLEIGH(M,S,T,dpi,dinit,u0,u1,e0,e1,x0,gamma,mat,detJ0,detDJ0,detJ,detDJ,m);

      break;

    case(NOH):

// Noh stagnation shock

      init_NOH(M,S,T,dpi,dinit,u0,u1,e0,e1,x0,gamma,mat,detJ0,detDJ0,detJ,detDJ,m);

      break;

    case(SEDOV):

// Sedov expanding shock

      init_SEDOV(M,S,T,dpi,dinit,u0,u1,e0,e1,x0,gamma,mat,detJ0,detDJ0,detJ,detDJ,m);

      break;

  }

// display the header so we have a date and time stamp in the output

  header();

// echo some initial information

  initial_data(n,nknodes,ntnodes,S,ndims,nmats,M);

// assign storage to the force matrix - use a compressed coordinate format

  for(int i=0, k=0;i<n;i++){
    for(int idim=0;idim<M.NDims();idim++){
      for(int iloc=0;iloc<S.nloc();iloc++){
        for(int jloc=0;jloc<T.nloc();jloc++,k++){
//          frow.push_back(idim*nknodes+M.GlobalNode_CFEM(i,iloc));
          frow.push_back(M.GlobalNode_CFEM(i,iloc)); // more convenient for loop structures where F[][] is accessed ?
//          fcol.push_back(idim*ntnodes+M.GlobalNode_DFEM(i,jloc));
          fcol.push_back(M.GlobalNode_DFEM(i,jloc)); // more convenient for loop structures where F[][] is accessed ?
          fdim.push_back(idim); // store the dimension idim associated with the address in F[][]
        }
      }
    }
  }

// set number of non-zeroes

  nzeroes=frow.size();
  F.resize(nzeroes);

  cout<<"Force matrix assigned, "<<nzeroes<<" non-zeroes in its address space: "<<fixed<<setprecision(2)<<F.size()*8.0/1024.0/1024.0<<" Mb acquired."<<endl;

// assemble mass matrix for acceleration field

  cout<<"Mass matrix for acceleration field assembly: ";
  for(int idim=0;idim<M.NDims();idim++){
    for(int i=0;i<M.NCells();i++){
      jacobian(i,x0,M,S,detJ,detDJ);
      for(int iloc=0;iloc<S.nloc();iloc++){
        for(int jloc=0;jloc<S.nloc();jloc++){
          double nn(0.0),nnx(0.0),nny(0.0),nx(0.0),ny(0.0),nlx(0.0),nly(0.0); // mass matrix
          for(int gi=0;gi<S.ngi();gi++){
            nn+=dinit.at(i)*S.value(iloc,gi)*S.value(jloc,gi)*detJ[gi]*S.wgt(gi);   // mass matrix
            nnx+=S.value(iloc,gi)*detDJ[0][jloc][gi]*detJ[gi]*S.wgt(gi);      // not used: left in as a diagnostic
            nny+=S.value(iloc,gi)*detDJ[1][jloc][gi]*detJ[gi]*S.wgt(gi);      // not used: left in as a diagnostic
            nlx+=detDJ[0][jloc][gi]*detJ[gi]*S.wgt(gi);                       // not used: left in as a diagnostic
            nly+=detDJ[1][jloc][gi]*detJ[gi]*S.wgt(gi);                       // not used: left in as a diagnostic
//          nx+=S.dvalue(0,iloc,gi)*(dy/2.0)*S.wgt(gi);                       // not used: left in as a diagnostic
//          ny+=S.dvalue(1,iloc,gi)*(dx/2.0)*S.wgt(gi);                       // not used: left in as a diagnostic
          }
          KMASS.add(idim*NROWS+ROW,idim*NROWS+COL,nn);
        }
      }
    }
  }
  cout<<"Done."<<endl;

// impose boundary constraints via row elimination

  bc_insert(KMASS,M,S,dinit,detJ,u0,u1,b0,b1,nknodes);

// invert the mass matrix

  cout<<"Inverting the mass matrix: ";

  KMASSI.inverse2(&KMASS); // lapack drivers dgetrf_ and dgetri_

  cout<<"Done."<<endl;

// set output precision

  cout<<fixed<<setprecision(17);

  cout<<"main(): Starting up main loop..."<<endl;

// time integration

  while(time<=ENDTIME){

// calculate a new stable time step that will impose the CFL limit on each integration point

    for(int i=0;i<n;i++){
      jacobian(i,xinit,M,S,detJ0,detDJ0);
      jacobian(i,x0,M,S,detJ,detDJ);
      jacobian(i,xinit,x1,M,S,detJs,Js);

// evaluate energy, divergence, length, density, pressure, sound speed and q at each integration point

      for(int gi=0;gi<S.ngi();gi++){
        double egi(0.0);
        for(int iloc=0;iloc<T.nloc();iloc++){
          egi+=e1.at(M.GlobalNode_DFEM(i,iloc))*T.value(iloc,gi);
        }
        egi=max(ECUT,egi);
        double divu(0.0);
        for(int idim=0;idim<M.NDims();idim++){
          for(int jloc=0;jloc<S.nloc();jloc++){divu+=S.dvalue(idim,jloc,gi)*u1.at(idim).at(M.GlobalNode_CFEM(i,jloc))/detJ.at(gi);}
        }
        l.at(gi)=M.UpdateLength(S.order(),V1.at(i),l0.at(i),detJs.at(gi),length_scale_type); // update element length scale at the quadrature point
        d.at(gi)=dinit.at(i)*detJ0.at(gi)/detJ.at(gi);
        p.at(gi)=M.UpdatePressure(d.at(gi),egi,gamma.at(mat.at(i)-1));
        c.at(gi)=M.UpdateSoundSpeed(gamma.at(mat.at(i)-1),p.at(gi),d.at(gi));
        q.at(gi)=M.UpdateQ(l.at(gi),d.at(gi),c.at(gi),cq,cl,divu);
        dt_cfl.at(GPNT)=(step==0)?DTSTART:(COURANT*S.wgt(gi)*(l.at(gi)/sqrt(c.at(gi)*c.at(gi)+2.0*q.at(gi)/d.at(gi)))); // impose the CFL limit on each quadrature point
      }
    }

// reduce time step across element and apply a saftey factor

    dts.at(0)=DTSFACTOR*(*min_element(dt_cfl.begin(), dt_cfl.end()));

// calculate a time step so we hit I/O event times smoothly and accurately

    double dt_g(-1.0*remainder(time,VISFREQ)); //  time until next graphics event
    dts.at(1)=1.0e6; // next I/O event too far away to care

    if( dt_g/(4.0*dt)<1.0 && dt_g>1.0e-12){
      if(dt_g/dt<4.000000001) {dts.at(1)=dt_g/4.0;}
      if(dt_g/dt<3.000000001) {dts.at(1)=dt_g/3.0;}
      if(dt_g/dt<2.000000001) {dts.at(1)=dt_g/2.0;}
      if(dt_g/dt<1.000000001) {dts.at(1)=dt_g;}
      if(abs(1.0-(dts[1]/dt))>1.0e-5){
        cout<<"       time-step cut on I/O event approach, dt lowered by "<<(1.0-dts[1]/dt)*100.00<<"% to "<<dts[1]<<endl;
      }
    }

// choose the smallest time step

//    dt=(*min_element(dts.begin(), dts.end()));
    dt=DTSTART;cout<<"DT HARDWIRED !! "<<endl;

// sum kinetic/internal energy fields for conservation checks

    sum_ke(ke,u1,dinit,M,xinit,x1,S,T,detJ0,detDJ0,detJ,detDJ);
    sum_ie(ie,e1,dinit,M,xinit,x1,S,T,detJ0,detDJ0,detJ,detDJ);

    cout<<fixed<<setprecision(5)<<"  step "<<step<<" time= "<<time<<" dt= "<<dt<<fixed<<setprecision(5)<<" energy (i/k/tot)= "<<ie<<" "<<ke<<" "<<ie+ke<<endl;

// graphics output

    if(abs(remainder(time,VISFREQ))<1.0e-12){
//      state_print(n,ndims,nmats,mat,dinit,V0,m,e0,p,x0,u0,step,time,gamma);
      silo(x0,xt0,xinit,dinit,l0,V0,e0,u0,mat,step,time,M,gamma,S,T);
    }else{
      if(abs(remainder(step,VISFREQ))==0){
//      state_print(n,ndims,nmats,mat,dinit,V0,m,e0,p,x0,u0,step,time,gamma);
      silo(x0,xt0,xinit,dinit,l0,V0,e0,u0,mat,step,time,M,gamma,S,T);
      }
    }

// debug
//  lineouts(M,S,T,d1,p,e1,q,xinit,x1,xt1,u1,test_problem,mat,gamma);
//  exit(1);
// debug

// advect the nodes to the full-step position

    M.UpdateCoords(x1,u0,dt);  // kinematics
    M.MapCoords(x1,xt1,S.order(),T.order()); // thermodynamics

// update volume field at the full-step

    M.UpdateVolume(V1,x1,S.order());

// assemble finite element energy field

    {Matrix A(T.nloc());vector<double> b(ntnodes);double bloc[T.nloc()],edot[T.nloc()];for(long i=0;i<b.size();i++){b.at(i)=0.0;}

// assemble the rhs of the energy equation from the force matrix using F^T dot (ux,uy)^T

      for(long iz=0;iz<nzeroes;iz++){b.at(fcol.at(iz))+=F.at(iz)*u1.at(fdim.at(iz)).at(frow.at(iz));}

// solve the energy equation locally in each cell

      for(int i=0;i<n;i++){

// update jacobians for this cell

        jacobian(i,xinit,M,S,detJ0,detDJ0);
        jacobian(i,x1,M,S,detJ,detDJ);

// update density field at each quadrature point in the cell

        for(int gi=0;gi<S.ngi();gi++){d.at(gi)=dinit.at(i)*detJ0.at(gi)/detJ.at(gi);}

// assemble a local mass matrix for the energy equation

        for(int iloc=0;iloc<T.nloc();iloc++){
          for(int jloc=0;jloc<T.nloc();jloc++){
            double nn(0.0); // mass matrix
            for(int gi=0;gi<S.ngi();gi++){
              nn+=d.at(gi)*T.value(iloc,gi)*T.value(jloc,gi)*detJ.at(gi)*S.wgt(gi);
            }
            A.write(iloc,jloc,nn);
          }
          bloc[iloc]=b.at(M.GlobalNode_DFEM(i,iloc)); // local sourcing of the energy field
        }

// solve the system for edot=de/dt

        A.solve(edot,bloc);

// advance the solution and commit to the global address space in the energy field

        for(int iloc=0;iloc<T.nloc();iloc++){e1.at(M.GlobalNode_DFEM(i,iloc))=e0.at(M.GlobalNode_DFEM(i,iloc))-edot[iloc]*dt;}

      }

    }

// assemble force matrix to connect thermodynamic/kinematic spaces, this can be used as rhs of both energy and momentum equations

    for(long k=0;k<nzeroes;k++){F.at(k)=0.0;}
    for(int i=0, k=0;i<n;i++){
      jacobian(i,xinit,x1,M,S,detJs,Js);
      jacobian(i,xinit,M,S,detJ0,detDJ0);
      jacobian(i,x1,M,S,detJ,detDJ);

// evaluate energy, divergence, length, density, pressure, sound speed and q at each integration point

      for(int gi=0;gi<S.ngi();gi++){
        double egi(0.0);
        for(int jloc=0;jloc<T.nloc();jloc++){
          egi+=e1.at(M.GlobalNode_DFEM(i,jloc))*T.value(jloc,gi);
        }
        egi=max(ECUT,egi);
        double divu(0.0);
        for(int idim=0;idim<M.NDims();idim++){
          for(int jloc=0;jloc<S.nloc();jloc++){
            divu+=detDJ[idim][jloc][gi]*u1.at(idim).at(M.GlobalNode_CFEM(i,jloc));
          }
        }


        l.at(gi)=M.UpdateLength(S.order(),V1.at(i),l0.at(i),detJs.at(gi),length_scale_type); // update element length scale at the quadrature point
        d.at(gi)=dinit.at(i)*detJ0.at(gi)/detJ.at(gi);
        p.at(gi)=P(d.at(gi),egi,gamma.at(mat.at(i)-1));
        c.at(gi)=M.UpdateSoundSpeed(gamma.at(mat.at(i)-1),p.at(gi),d.at(gi));
        q.at(gi)=M.UpdateQ(l.at(gi),d.at(gi),c.at(gi),cq,cl,divu);
      }

// construct force terms

      for(int idim=0;idim<M.NDims();idim++){
        for(int iloc=0;iloc<S.nloc();iloc++){
          for(int jloc=0;jloc<T.nloc();jloc++,k++){
            for(int gi=0;gi<S.ngi();gi++){
              F.at(k)+=(p.at(gi)+q.at(gi))*detDJ[idim][iloc][gi]*T.value(jloc,gi)*detJ[gi]*S.wgt(gi);
            }
          }
        }
      }

    }

// assemble acceleration field

    vector<double> b(M.NDims()*nknodes);
    {for(long i=0;i<b.size();i++){b.at(i)=-b1.at(i);}

// assemble the rhs of the momentum equation from the force matrix using F dot (unit vector)^T

      for(long iz=0;iz<nzeroes;iz++){b.at(fdim.at(iz)*nknodes+frow.at(iz))+=F.at(iz);}

// impose boundary constraints

      for(long i=0;i<b.size();i++){if(b0.at(i)!=0.0){b.at(i)=b0.at(i);}}

// solve global system

      for(int idim=0;idim<M.NDims();idim++){
        for(long i=0;i<nknodes;i++){
          double udot(0.0);long ladd(idim*nknodes);
          for(long j=0;j<nknodes;j++){
            udot+=KMASSI.read(ladd+i,ladd+j)*b.at(ladd+j);
          }

// advance the solution and commit to the global address space

          u1.at(idim).at(i)=u0.at(idim).at(i)+udot*dt;

        }
      }

    }

// advance the time step

    time+=dt;
    step++;

// advance the flux states for the new time step

    u0=u1;   // velocity field
    x0=x1;   // kinematic node coordinates
    xt0=xt1; // thermodynamic node coordinates
    e0=e1;   // internal energy
    V0=V1;   // cell volumes

// debug
//  if(step==5){
//    cout<<"debug stop."<<endl;
//    output(); // might want this ??
//    exit(1);
//  }
// debug

  }

// some output

  lineouts(M,S,T,dinit,e1,xinit,x1,xt1,u1,test_problem,mat,gamma);
  exact(state,x1,test_problem);

// estimate convergence rate in the L1/L2 norms using a Riemann solution as the exact solution

//  vector<double> rx;vempty(rx);  // sample points
//  double xstart(0.0),xstop(1.0); // sample range - should be whole mesh though bc.s may artificially reduce convergence
//  for(long i=0;i<ng+1;i++){if(x1[i]>=xstart&&x1[i]<=xstop){rx.push_back(x1[i]);}}
//  R0.profile(&rx,ENDTIME); // evolve a Riemann problem on the sample range to final time level
//
//  double l1(0.0),l2(0.0),l1r(0.0),l2r(0.0),l1d(0.0),l2d(0.0),l1n(0.0),l2n(0.0);
//  int ii(0);

// numerator and denominator for each norm

//  for(int i=0;i<ng+1;i++){
//    if(x1[i]>=xstart&&x1[i]<=xstop){
//      double err(abs(R0.velocity(ii)-u1[i])); // absolute error
//      l1d+=abs(R0.velocity(ii)); // l1 denominator
//      l2d+=R0.velocity(ii)*R0.velocity(ii); // l2 denominator
//      l1n+=err*V1[i]; // l1 numerator
//      l2n+=err*err*V1[i]; // l2 numerator
//      ii++;
//    }
//  }

// construct norms and relative errors

//  l1=l1n/ii; // L1 error norm
//  l1r=l1n/l1d; // L1 relative error norm
//  l2=sqrt(l2n/ii); // L2 error norm
//  l2r=sqrt(l2n/l2d); // L2 relative error norm

//  cout<<endl<<fixed<<setprecision(10)<<"  Error Estimates (grid spacing h= "<<(xstop-xstart)/ii<<"):"<<endl;

//  cout<<"  L1 norm= "<<l1<<" (relative error= "<<l1r<<")"<<endl;
//  cout<<"  L2 norm= "<<l2<<" (relative error= "<<l2r<<")"<<endl;
//  cout<<"  No. points sampled= "<<ii<<endl;
//  cout<<"  Range sampled= "<<xstart<<","<<xstop<<endl;

  cout<<"Normal termination."<<endl;

  return 0;
}

// this function codes for some line-outs

void lineouts(Mesh const &M,Shape const &S,Shape const &T,VD const &dinit,VD const &e,VVD const &xinit,VVD const &x,VVD const &xt,VVD const &u, int const &test_problem,vector<int> const &mat,VD const &g){

// file handle for output

  ofstream f1;

// establish the mesh limits

  double xmin(*min_element(x.at(0).begin(),x.at(0).end()));
  double xmax(*max_element(x.at(0).begin(),x.at(0).end()));
  double ymin(*min_element(x.at(1).begin(),x.at(1).end()));
  double ymax(*max_element(x.at(1).begin(),x.at(1).end()));

// decalre the line-out structure

  struct lineout_type {
    double x1,y1; // start point of each line
    double x2,y2; // end point of each line
    string filename; //filename to output
    string filehead; // file header
    int nsamples; // number of sample points on each line
  } lineout;

  vector<lineout_type> Lineout;

// set up different line-outs for different problems

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

      for(int iline=0;iline<4;iline++){

        switch(iline){

          case(0):

            lineout.x1=0.5*(xmin+xmax);
            lineout.x2=xmax;
            lineout.y1=0.5*(ymin+ymax);
            lineout.y2=ymax;
            lineout.filename="lineout_1.dat";
            lineout.filehead="# Sedov lineout from (0.0,0.0) to (1.2,1.2) : Columns are x d p e q u";
            lineout.nsamples=100;

            break;

          case(1):

            lineout.x1=0.5*(xmin+xmax);
            lineout.x2=xmin;
            lineout.y1=0.5*(ymin+ymax);
            lineout.y2=ymax;
            lineout.filename="lineout_2.dat";
            lineout.filehead="# Sedov lineout from (0.0,0.0) to (-1.2,1.2) : Columns are x d p e q u";
            lineout.nsamples=100;

            break;

          case(2):

            lineout.x1=0.5*(xmin+xmax);
            lineout.x2=xmin;
            lineout.y1=0.5*(ymin+ymax);
            lineout.y2=ymin;
            lineout.filename="lineout_3.dat";
            lineout.filehead="# Sedov lineout from (0.0,0.0) to (-1.2,-1.2) : Columns are x d p e q u";
            lineout.nsamples=100;

            break;

          case(3):

            lineout.x1=0.5*(xmin+xmax);
            lineout.x2=xmax;
            lineout.y1=0.5*(ymin+ymax);
            lineout.y2=ymin;
            lineout.filename="lineout_4.dat";
            lineout.filehead="# Sedov lineout from (0.0,0.0) to (1.2,-1.2) : Columns are x d p e q u";
            lineout.nsamples=100;

            break;

        }

      }

      break;

    case(SOD):

// Sod's shock tube

      lineout.x1=-0.001;
      lineout.x2=1.0001;
      lineout.y1=0.5*(ymin+ymax);
      lineout.y2=lineout.y1;
      lineout.filename="lineout_1.dat";
      lineout.filehead="# Sod lineout from (0.0,0.5) to (1.0,0.5) : Columns are x d e vx vy p";
      lineout.nsamples=100;

      break;

    case(R2R):

// 123 problem

      lineout.x1=xmin;
      lineout.x2=xmax;
      lineout.y1=0.5*(ymin+ymax);
      lineout.y2=lineout.y1;
      lineout.filename="lineout_1.dat";
      lineout.filehead="# 123 problem lineout from (xmin,0.5) to (xmax,0.5) : Columns are x d e vx vy p";
      lineout.nsamples=20;

      break;

  }

// append to the data structure

  Lineout.push_back(lineout);

// loop over the lineouts and produce the output

  for(int iline=0;iline<Lineout.size();iline++){

    cout<<"lineouts(): Lineout "<<iline<<" writing to file "<<Lineout.at(iline).filename<<" ..."<<endl;

    double p(0.0),q(0.0),d(0.0);

// open the output file for the lineout and write the header part

    f1.open(Lineout.at(iline).filename);
    f1<<Lineout.at(iline).filehead<<endl;

// set up line AB to sample along

    vector<double> A(2),B(2); // A and B are the two end points of the line
    vector<double> E(2),F(2); // E and F are the two end points of a segment EF on the line AB

    A.at(0)=Lineout.at(iline).x1;A.at(1)=Lineout.at(iline).y1;
    B.at(0)=Lineout.at(iline).x2;B.at(1)=Lineout.at(iline).y2;

    Line AB(A,B);
    AB.divide(Lineout.at(iline).nsamples);

// search for cells that AB intersects

    E.at(0)=A.at(0);E.at(1)=A.at(1); // start of segment 1 is origin of AB
    int celllist[AB.nsegments()];fill_n(celllist,AB.nsegments(),-1); // cell transit list

    for(int iseg=0;iseg<AB.nsegments();iseg++){

// create a line to represent this segment

      F.at(0)=AB.coord(0,iseg);F.at(1)=AB.coord(1,iseg); // segment end point
      Line EF(E,F); // this line is the segment

      for(int i=0;i<M.NCells();i++){

// set up line CD to represent each cell side in turn

        vector<double> C(2),D(2); // C and D are the two end points of a cell side
        int endloc[5]={0,S.sloc()-1,S.nloc()-1,S.nloc()-S.sloc(),0}; // node at end of face
        int nsides(0); // number of sides of i that are crossed by AB

        for(int j=0;j<S.nfaces();j++){

          C.at(0)=x[0].at(M.GlobalNode_CFEM(i,endloc[j])),C.at(1)=x[1].at(M.GlobalNode_CFEM(i,endloc[j]));
          D.at(0)=x[0].at(M.GlobalNode_CFEM(i,endloc[j+1])),D.at(1)=x[1].at(M.GlobalNode_CFEM(i,endloc[j+1]));

// create a line CD that is coincident with the cell side and check if it is intersected by the segment EF on the line AB

          Line CD(C,D);

          if(EF.intersects(CD)) {nsides++;}

        }

// if at least 1 side was crossed store the cell number against this segment

        if(nsides!=0){celllist[iseg]=i;}

      }

// update start of next segment

      E.at(0)=F.at(0);E.at(1)=F.at(1);

    }

// traverse the mesh along the line AB and interpolate data onto each segment end point
// Note we ignore last segment to avoid an extra point at end of line-out

    int i(0);
    for(int iseg=0;iseg<AB.nsegments()-1;iseg++){

      if(celllist[iseg]>=0){i=celllist[iseg];} // update cell number, <0 means iseg is fully contained within the cell

      double interpolated_value[5];

// cell vertices

      vector<vector<double> > rk(2);
      for(int j=0;j<S.nloc();j++){
        rk.at(0).push_back(x.at(0).at(M.GlobalNode_CFEM(i,j)));
        rk.at(1).push_back(x.at(1).at(M.GlobalNode_CFEM(i,j)));
      }

      vector<vector<double> > rt(2);
      for(int j=0;j<T.nloc();j++){
        rt.at(0).push_back(xt.at(0).at(M.GlobalNode_DFEM(i,j)));
        rt.at(1).push_back(xt.at(1).at(M.GlobalNode_DFEM(i,j)));
      }

// instantiate a shape function in global coordinates for each space

      Shape Gk(S.order(),rk);
      Shape Gt(T.order(),rt);

// global coordinates of interpolation

      vector<double> ri;
      ri.push_back(AB.coord(0,iseg));
      ri.push_back(AB.coord(1,iseg));

// local node positions

      double xpos[S.nloc()],ypos[S.nloc()],dx(2.0/S.order()),dy(2.0/S.order());

      for(int isuby=0,k=0;isuby<S.order()+1;isuby++){
        for(int isubx=0;isubx<S.order()+1;isubx++,k++){
          xpos[k]=-1.0+isubx*dx;
          ypos[k]=-1.0+isuby*dy;
        }
      }

// map global coordinates to local coordinates

      double xloc(0.0),yloc(0.0);
      for(int iloc=0;iloc<S.nloc();iloc++){
        xloc+=Gk.value(iloc,ri)*xpos[iloc];
        yloc+=Gk.value(iloc,ri)*ypos[iloc];
      }

// jacobian and determinant at local coordinates

      double dxdu(0.0),dxdv(0.0),dydu(0.0),dydv(0.0),dxdu0(0.0),dxdv0(0.0),dydu0(0.0),dydv0(0.0);

      for(int iloc=0;iloc<S.nloc();iloc++){
        dxdu0+=S.dvalue(0,iloc,xloc,yloc)*xinit.at(0).at(M.GlobalNode_CFEM(i,iloc));
        dxdv0+=S.dvalue(1,iloc,xloc,yloc)*xinit.at(0).at(M.GlobalNode_CFEM(i,iloc));
        dydu0+=S.dvalue(0,iloc,xloc,yloc)*xinit.at(1).at(M.GlobalNode_CFEM(i,iloc));
        dydv0+=S.dvalue(1,iloc,xloc,yloc)*xinit.at(1).at(M.GlobalNode_CFEM(i,iloc));
        dxdu+=S.dvalue(0,iloc,xloc,yloc)*x.at(0).at(M.GlobalNode_CFEM(i,iloc));
        dxdv+=S.dvalue(1,iloc,xloc,yloc)*x.at(0).at(M.GlobalNode_CFEM(i,iloc));
        dydu+=S.dvalue(0,iloc,xloc,yloc)*x.at(1).at(M.GlobalNode_CFEM(i,iloc));
        dydv+=S.dvalue(1,iloc,xloc,yloc)*x.at(1).at(M.GlobalNode_CFEM(i,iloc));
      }

      double detJ0(dxdu0*dydv0-dxdv0*dydu0),detJ(dxdu*dydv-dxdv*dydu);

// density at interpolation point in local coordinates

      interpolated_value[0]=dinit.at(i)*detJ0/detJ;

// evaluate energy field at global coordinate ri(x,y)

      for(int j=0;j<Gt.nloc();j++){interpolated_value[1]+=Gt.value(j,ri)*e.at(M.GlobalNode_DFEM(i,j));}

// evaluate velocity field x-component at global coordinate ri(x,y)

      for(int j=0;j<Gk.nloc();j++){interpolated_value[2]+=Gk.value(j,ri)*u.at(0).at(M.GlobalNode_CFEM(i,j));}

// evaluate velocity field y-component at global coordinate ri(x,y)

      for(int j=0;j<Gk.nloc();j++){interpolated_value[3]+=Gk.value(j,ri)*u.at(1).at(M.GlobalNode_CFEM(i,j));}

// evaluate pressure field at global coordinate ri(x,y)

      interpolated_value[4]=P(interpolated_value[0],interpolated_value[1],g.at(mat.at(i)-1));

// output the interpolated values along the line AB

      double xx(sqrt(pow((AB.coord(0,0)-ri.at(0)),2)+pow((AB.coord(1,0)-ri.at(1)),2))+(AB.coord(0,1)-AB.coord(0,0))); // distance along line-out from end of the first segment
      f1<<fixed<<setprecision(10)<<xx<<" ";
      for(int j=0;j<5;j++){
        f1<<interpolated_value[j]<<" ";
      }
      f1<<endl;

    }

// close the output file

    f1.close();

  }

  return;

}

// generate an exact solution where applicable

void exact(VVD const &s,VVD const &x,int const &test_problem){

// file handle for output

  ofstream f1;

// establish the mesh limits

  double xmin(*min_element(x.at(0).begin(),x.at(0).end()));
  double xmax(*max_element(x.at(0).begin(),x.at(0).end()));
  double ymin(*min_element(x.at(1).begin(),x.at(1).end()));
  double ymax(*max_element(x.at(1).begin(),x.at(1).end()));

// decalre the line-out structure

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
        R.profile(&rx,ENDTIME); // Riemann solution at NSAMPLE sample points

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
        R.profile(&rx,ENDTIME); // Riemann solution at NSAMPLE sample points

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
  }

  return;

}

// empty a vector

void vempty(vector<double>&v){vector<double> e;v.swap(e);return;}

// return today's date

string date(){time_t now = time(0);return ctime(&now);}

// return user name

string gethost(){

  char*hostname(getenv("HOSTNAME"));
  string computername;

  if(hostname!=NULL) {
    computername=hostname;
    hostname=NULL;
  }else{
    hostname=new char[512];
    if(gethostname(hostname,512)==0){ // success = 0, failure = -1
      computername=hostname;
    }
    delete[]hostname;
    hostname=NULL;
  }

  return computername;
}

// output the header with a date and time stamp

void header(){

  cout<<endl;
  cout<<"             RhoLo"<<endl;
  cout<<endl;
  cout<<"Date: "<<date();
  cout<<"Host: "<<gethost()<<endl;
  cout<<"User: "<<getenv("USER")<<endl;
  cout<<"Home: "<<getenv("HOME")<<endl;
  cout<<"pwd:  "<<getenv("PWD")<<endl;

  cout<<endl;

}

// output some initial data

void initial_data(int const n,long const nknodes,long const ntnodes,Shape const S,int const ndims, int const nmats, Mesh const &M){

  cout<<"Initial data for the simulation"<<endl;

  cout<<"Number of dimensions:             "<<ndims<<endl;
  cout<<"Number of cells:                  "<<n<<endl;
  cout<<"Number of kinematic nodes:        "<<nknodes<<endl;
  cout<<"Number of thermodynamic nodes:    "<<ntnodes<<endl;
  cout<<"Number of integration points:     "<<S.ngi()<<endl;
  cout<<"Number of materials:              "<<nmats<<endl;
  cout<<"Boundary conditions have been set on "<<M.nbcs()<<" edges :"<<endl;

  for(int i=0;i<M.nbcs();i++){

    string bcname;

    switch(M.bc_edge(i)){
      case(VACUUM):
        bcname="vacuum";
        break;
      case(VELOCITY):
        bcname="velocity";
        break;
      case(REFLECTIVE):
        bcname="reflective";
        break;
      case(FLUID):
        bcname="fluid";
        break;
      case(ACCELERATION):
        bcname="acceleration";
        break;
    }

    cout<<"Edge "<<i<<" boundary type : "<<bcname<<endl;
  }

  return;

}

// material state - args are passed by ref to prevent a copy, and are const to prevent changes

void state_print(int const n,int const ndims,int const nmats,VI const &mat,VD const &d,VD const &V,
                  VD const &m,VD const &e,VD const &p,VVD const &x, VVD const &u, int const &s, double const &t,VD const &gamma){

// material average data

  vector<double> vmat(nmats); // volume of material
  vector<double> pmat(nmats); // average pressure
  vector<double> mmat(nmats); // mass of material
  vector<double> dmat(nmats); // density of material  
  vector<double> iemat(nmats); // internal energy
  vector<double> kemat(nmats); // kinetic energy
  vector<double> temat(nmats); // total energy

// totals

  double vtot(0.0),ptot(0.0),mtot(0.0),dtot(0.0),ietot(0.0),ketot(0.0),tetot(0.0);

// state of each material

  for(int i=0;i<n;i++){
    int ii(mat[i]-1);
    vmat.at(ii)+=V[i];
    mmat.at(ii)+=m[i];
    iemat.at(ii)+=e[i]*m[i];
    kemat.at(ii)+=0.0;
    temat.at(ii)+=e[i]*m[i]+0.0;
  }

  for(int i=0;i<nmats;i++){
    dmat.at(i)=mmat[i]/vmat[i];
    pmat.at(i)=P(dmat[i],iemat[i]/mmat[i],gamma[mat[i]-1]);
  }

// sum across material to get domain totals

  for(int i=0;i<nmats;i++){
    vtot+=vmat[i];
    ptot+=pmat[i]*vmat[i];
    mtot+=mmat[i];
    ietot+=iemat[i];
    ketot+=kemat[i];
  }

  ptot=ptot/vtot;
  dtot=mtot/vtot;
  tetot=ietot+ketot;

// print out the material averages in a table

  cout<<endl;
  cout<<fixed<<setprecision(6)<<"Material state print at time "<<t<<" on step "<<s<<endl;
  cout<<"~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"<<endl;
  cout<<"Mat.  "<<"Density  "<<"Mass     "<<"Volume   "<<"Pressure "<<"Energy (Specific IE/IE/KE/Total)"<<endl;
  cout<<"~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"<<endl;
  for(int i=0;i<nmats;i++){
    cout<<fixed<<setprecision(6);
    cout<<i+1<<"     "<<dmat[i]<<" "<<mmat[i]<<" "<<vmat[i]<<" "<<pmat[i]<<" "<<iemat[i]/mmat[i]<<" "<<iemat[i]<<" "<<kemat[i]<<" "<<temat[i]<<endl;
  }

  cout<<"total "<<dtot<<" "<<mtot<<" "<<vtot<<" "<<P(dtot,ietot,1.4)<<" "<<ietot/mtot<<" "<<ietot<<" "<<ketot<<" "<<tetot<<endl;

  cout<<"~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"<<endl;

  cout<<endl;

  return;

}

// calculate a jacobian and the determinant

void jacobian(int const &i,VVD const &x,Mesh const &M,Shape const &S,VD &detJ,VVVD &detDJ){

// loop over quadrature points and calculate the jacobian

  for(int gi=0;gi<S.ngi();gi++){

    double dxdu(0.0),dydu(0.0),dxdv(0.0),dydv(0.0);

// derivatives of the physical coordinates at the quadrature points

    for(int iloc=0;iloc<S.nloc();iloc++){
      long gloc=(S.type()==CONTINUOUS)?M.GlobalNode_CFEM(i,iloc):M.GlobalNode_DFEM(i,iloc);
      dxdu+=x.at(0).at(gloc)*S.dvalue(0,iloc,gi); // dx/du
      dxdv+=x.at(0).at(gloc)*S.dvalue(1,iloc,gi); // dx/dv
      dydu+=x.at(1).at(gloc)*S.dvalue(0,iloc,gi); // dy/du
      dydv+=x.at(1).at(gloc)*S.dvalue(1,iloc,gi); // dy/dv
    }

// calculate the determinant at the quadrature point and commit to the vector

    detJ.at(gi)=dxdu*dydv-dxdv*dydu;

// determinants for the deriavtives at the quadrature points

    for(int iloc=0;iloc<S.nloc();iloc++){
      detDJ.at(0).at(iloc).at(gi)=(dydv*S.dvalue(0,iloc,gi)-dydu*S.dvalue(1,iloc,gi))/detJ[gi]; // original, Taylor runs better ??
//      detDJ.at(0).at(iloc).at(gi)=(dydv*S.dvalue(0,iloc,gi)-dxdv*S.dvalue(1,iloc,gi))/detJ[gi]; // Noh/triple run better ??
      detDJ.at(1).at(iloc).at(gi)=(-dxdv*S.dvalue(0,iloc,gi)+dxdu*S.dvalue(1,iloc,gi))/detJ[gi]; // original, Taylor runs better ??
//      detDJ.at(1).at(iloc).at(gi)=(-dydu*S.dvalue(0,iloc,gi)+dxdu*S.dvalue(1,iloc,gi))/detJ[gi];// Noh/triple run better ??
    }

  }

  return;

}

// calculate a jacobian for the Lagrangian motion and return the determinant

void jacobian(int const &i,VVD const &x0,VVD const &x,Mesh const &M,Shape const &S,VD &detJs,VVVD &Js){

// loop over quadrature points and calculate jacobians J0 and J

  for(int gi=0;gi<S.ngi();gi++){

    double dx0du(0.0),dy0du(0.0),dx0dv(0.0),dy0dv(0.0);
    double dxdu(0.0),dydu(0.0),dxdv(0.0),dydv(0.0);

// derivatives of the physical coordinates at the quadrature points

    for(int iloc=0;iloc<S.nloc();iloc++){

      long gloc=(S.type()==CONTINUOUS)?M.GlobalNode_CFEM(i,iloc):M.GlobalNode_DFEM(i,iloc);

// at time 0

//      dx0du+=x0.at(0).at(gloc)*S.dvalue(0,iloc,gi); // dx0/du
//      dx0dv+=x0.at(0).at(gloc)*S.dvalue(1,iloc,gi); // dx0/dv
//      dy0du+=x0.at(1).at(gloc)*S.dvalue(0,iloc,gi); // dy0/du
//      dy0dv+=x0.at(1).at(gloc)*S.dvalue(1,iloc,gi); // dy0/dv

//      dx0du=1.0; // mod for direction s
//      dx0dv=1.0; // mod for direction s
//      dy0du+=x0.at(1).at(gloc)*S.dvalue(0,iloc,gi); // dy0/du
//      dy0dv+=x0.at(1).at(gloc)*S.dvalue(1,iloc,gi); // dy0/dv

//      dx0du+=x0.at(0).at(gloc)*S.dvalue(0,iloc,gi); // dx0/du
//      dx0dv=1.0; // mod for direction s
//      dy0du=1.0; // mod for direction s
//      dy0dv+=x0.at(1).at(gloc)*S.dvalue(1,iloc,gi); // dy0/dv

//      dx0du=1.0; // mod for direction s
//      dx0dv=1.0; // mod for direction s
//      dy0du=1.0;
//      dy0dv+=x0.at(1).at(gloc)*S.dvalue(1,iloc,gi); // dy0/dv

      dx0du=1.0; // mod for direction s
      dx0dv+=x0.at(0).at(gloc)*S.dvalue(1,iloc,gi); // dx0/dv
      dy0du+=x0.at(1).at(gloc)*S.dvalue(0,iloc,gi); // dy0/du
      dy0dv=1.0; // mod for direction s

// at time t

//      dxdu+=x.at(0).at(gloc)*S.dvalue(0,iloc,gi); // dx/du
//      dxdv+=x.at(0).at(gloc)*S.dvalue(1,iloc,gi); // dx/dv
//      dydu+=x.at(1).at(gloc)*S.dvalue(0,iloc,gi); // dy/du
//      dydv+=x.at(1).at(gloc)*S.dvalue(1,iloc,gi); // dy/dv

//      dxdu=1.0; // mod for direction s
//      dxdv=1.0; // mod for direction s
//      dydu+=x.at(1).at(gloc)*S.dvalue(0,iloc,gi); // dy/du
//      dydv+=x.at(1).at(gloc)*S.dvalue(1,iloc,gi); // dy/dv

//      dxdu+=x.at(0).at(gloc)*S.dvalue(0,iloc,gi); // dx/du
//      dxdv=1.0; // mod for direction s
//      dydu=1.0; // mod for direction s
//      dydv+=x.at(1).at(gloc)*S.dvalue(1,iloc,gi); // dy/dv

//      dxdu=1.0; // mod for direction s
//      dxdv=1.0; // mod for direction s
//      dydu=1.0;
//      dydv+=x.at(1).at(gloc)*S.dvalue(1,iloc,gi); // dy/dv

      dxdu=1.0; // mod for direction s
      dxdv+=x.at(0).at(gloc)*S.dvalue(1,iloc,gi); // dx/dv
      dydu+=x.at(1).at(gloc)*S.dvalue(0,iloc,gi); // dy/du
      dydv=1.0; // mod for direction s

    }

// define the jacobian for time-0, this maps to the isoparametric element from the time-0 element

    vector<vector<double> > J0{{dx0du,dx0dv},
                               {dy0du,dy0dv}};

// determinant of J0

    double detJ0(dx0du*dy0dv-dx0dv*dy0du);

// the inverse of J0 is simply the adjugate divided by the determinant

    vector<vector<double> > invJ0{{ dy0dv/detJ0,-dx0dv/detJ0},
                                  {-dy0du/detJ0, dx0du/detJ0}};

// define the jacobian for time-t, this maps to the isoparametric element from the time-t element

    vector<vector<double> > J{{dxdu,dxdv},
                              {dydu,dydv}};

// determinant of J

    double detJ(dxdu*dydv-dxdv*dydu);

// define a jacobian of the lagrangian motion, this maps to the time-0 element from the time-t element

//    vector<vector<double> > Js{{invJ0[0][0]*J[0][0]+invJ0[0][1]*J[1][0],invJ0[0][0]*J[0][1]+invJ0[0][1]*J[1][1]},
//                               {invJ0[1][0]*J[0][0]+invJ0[1][1]*J[1][0],invJ0[1][0]*J[0][1]+invJ0[1][1]*J[1][1]}};

    Js.at(gi)[0][0]=invJ0[0][0]*J[0][0]+invJ0[0][1]*J[1][0];
    Js.at(gi)[0][1]=invJ0[0][0]*J[0][1]+invJ0[0][1]*J[1][1];
    Js.at(gi)[1][0]=invJ0[1][0]*J[0][0]+invJ0[1][1]*J[1][0];
    Js.at(gi)[1][1]=invJ0[1][0]*J[0][1]+invJ0[1][1]*J[1][1];

// determinant of Js

    detJs.at(gi)=(Js.at(gi)[0][0]*Js.at(gi)[1][1]-Js.at(gi)[0][1]*Js.at(gi)[1][0]);

  }

  return;

}

// calculate a jacobian at the local nodes

void jacobian(int const &i,VVD const &x,Mesh const &M,Shape const &S,VD &detJ){

// node positions and displacement in local coordinates

  double xpos[S.nloc()],ypos[S.nloc()],disp(2.0/S.order());

  for(int jsloc=0,kloc=0;jsloc<S.sloc();jsloc++){
    for(int isloc=0;isloc<S.sloc();isloc++,kloc++){
      xpos[kloc]=-1.0+isloc*disp;
      ypos[kloc]=-1.0+jsloc*disp;
    }
  }

// loop over local nodes and calculate the jacobian

  for(int iloc=0;iloc<S.nloc();iloc++){

    double dxdu(0.0),dydu(0.0),dxdv(0.0),dydv(0.0);

// derivatives of the physical coordinates at the nodes

    for(int jloc=0;jloc<S.nloc();jloc++){
      long gloc=(S.type()==CONTINUOUS)?M.GlobalNode_CFEM(i,jloc):M.GlobalNode_DFEM(i,jloc);
      dxdu+=x.at(0).at(gloc)*S.dvalue(0,jloc,xpos[iloc],ypos[iloc]); // dx/du
      dxdv+=x.at(0).at(gloc)*S.dvalue(1,jloc,xpos[iloc],ypos[iloc]); // dx/dv
      dydu+=x.at(1).at(gloc)*S.dvalue(0,jloc,xpos[iloc],ypos[iloc]); // dy/du
      dydv+=x.at(1).at(gloc)*S.dvalue(1,jloc,xpos[iloc],ypos[iloc]); // dy/dv
    }

// calculate the determinant at the quadrature point and commit to the vector

    detJ.at(iloc)=dxdu*dydv-dxdv*dydu;

  }

  return;

}

// input overides for the Taylor-Green vortex

void init_TAYLOR(Mesh const &M,Shape const &S,Shape const &T,double const &dpi,VD &dinit,VVD &u0,VVD &u1,VD &e0,VD &e1,VVD const &x,VD const &gamma,vector<int> const &mat,VD &detJ0,VVVD &detDJ0,VD &detJ,VVVD &detDJ,VD const &m){

  cout<<"init_TAYLOR(): Input overides for Taylor-Green..."<<endl;

// start the energy field

  for(int i=0;i<M.NCells();i++){
    for(int iloc=0;iloc<T.nloc();iloc++){
      double xval(x.at(0).at(M.GlobalNode_DFEM(i,iloc)));
      double yval(x.at(0).at(M.GlobalNode_DFEM(i,iloc)));
      e0.at(M.GlobalNode_DFEM(i,iloc))=(3.0*dpi/8.0)*cos(3.0*dpi*xval)*cos(dpi*yval)-cos(dpi*xval)*cos(3.0*dpi*yval);
      e1.at(M.GlobalNode_DFEM(i,iloc))=e0.at(M.GlobalNode_DFEM(i,iloc));
    }
  }

  return;

}

// input overides for the Rayleigh-Taylor instability

void init_RAYLEIGH(Mesh const &M,Shape const &S,Shape const &T,double const &dpi,VD &dinit,VVD &u0,VVD &u1,VD &e0,VD &e1,VVD const &x,VD const &gamma,vector<int> const &mat,VD &detJ0,VVVD &detDJ0,VD &detJ,VVVD &detDJ,VD const &m){

  cout<<"init_RAYLEIGH(): Input overides for the Rayleigh-Taylor instability test not coded yet."<<endl;

  exit(1);

  return;

}

// input overides for the Noh stagnation shock

void init_NOH(Mesh const &M,Shape const &S,Shape const &T,double const &dpi,VD &dinit,VVD &u0,VVD &u1,VD &e0,VD &e1,VVD const &x,VD const &gamma,vector<int> const &mat,VD &detJ0,VVVD &detDJ0,VD &detJ,VVVD &detDJ,VD const &m){

  cout<<"init_NOH(): Input overides for Noh..."<<endl;

// set origin

  double xorig(0.5*(M.Min(0)+M.Max(0))),yorig(0.5*(M.Min(1)+M.Max(1)));

// correct origin for reflection

  if((M.bc_edge(0)==VELOCITY)){yorig=x.at(1).at(0);} // ymin forced reflective
  if((M.bc_edge(3)==VELOCITY)){xorig=x.at(0).at(0);} // xmin forced reflective

  double origin[2]={xorig,yorig}; // origin coordinates

// start velocity field

  for(long i=0;i<u0.at(0).size();i++){

    double rx(x.at(0).at(i)-origin[0]);     // radial vector component from domain origin to node
    double ry(x.at(1).at(i)-origin[1]);     // radial vector component from domain origin to node

// length of radial vector from domain origin to node

    double rnorm(sqrt(rx*rx+ry*ry));

// velocity is a radial vector from degree of freedom towards the domain origin

    u0.at(0).at(i)=-rx/max(rnorm,1.0e-12);
    u1.at(0).at(i)=u0.at(0).at(i);

    u0.at(1).at(i)=-ry/max(rnorm,1.0e-12);
    u1.at(1).at(i)=u0.at(1).at(i);

  }

  return;

}

// input overides for the Sedov explosion

void init_SEDOV(Mesh const &M,Shape const &S,Shape const &T,double const &dpi,VD &dinit,VVD &u0,VVD &u1,VD &e0,VD &e1,VVD const &x,VD const &gamma,vector<int> const &mat,VD &detJ0,VVVD &detDJ0,VD &detJ,VVVD &detDJ,VD const &m){

  cout<<"init_SEDOV(): Input overides for Sedov..."<<endl;

// load energy field

  for(long i=0;i<M.NCells();i++){

    for(int iloc=0;iloc<T.nloc();iloc++){

    e0.at(M.GlobalNode_DFEM(i,iloc))=0.0;
    e1.at(M.GlobalNode_DFEM(i,iloc))=0.0;

// delta function at domain origin

      if( abs(x.at(0).at(M.GlobalNode_DFEM(i,iloc)))<1.0e-7 && abs(x.at(1).at(M.GlobalNode_DFEM(i,iloc)))<1.0e-7  ){
//        e0.at(i)=0.3014676/0.025; // drive ~ 12.058704, from another code see ref. paper for numbers
        e0.at(M.GlobalNode_DFEM(i,iloc))=0.25/m[i]; // place 1/4 of the drive in each of the 4 cells at the origin per unit mass
        e1.at(M.GlobalNode_DFEM(i,iloc))=e0.at(M.GlobalNode_DFEM(i,iloc));
      }

    }

  }

  return;

}

// sum the global kinetic energy field

void sum_ke(double &ke,VVD const &u,VD const &dinit,Mesh const &M,VVD const &xinit,VVD const &x,Shape const &S,Shape const &T,VD &detJ0,VVVD &detDJ0,VD &detJ,VVVD &detDJ){

// density as a function

  vector<double> d(S.ngi());

// sum global kinetic field field

  ke=0.0;

  for(int i=0;i<M.NCells();i++){

// update jacobians

    jacobian(i,xinit,M,S,detJ0,detDJ0);
    jacobian(i,x,M,S,detJ,detDJ);

// density at the Gauss points

    for(int gi=0;gi<T.ngi();gi++){d.at(gi)=(dinit.at(i)*detJ0.at(gi)/detJ.at(gi));}

// loop over element velocity field

    for(int iloc=0;iloc<S.nloc();iloc++){

// compute update nodal mass

      double nodmass(0.0);
      for(int gi=0;gi<T.ngi();gi++){nodmass+=S.value(iloc,gi)*d.at(gi)*S.wgt(gi);}

      for(int idim=0;idim<M.NDims();idim++){
        ke+=nodmass*u.at(idim).at(M.GlobalNode_CFEM(i,iloc))*u.at(idim).at(M.GlobalNode_CFEM(i,iloc));
      }

    }

  }

  ke=0.5*ke;

  return;

}

// sum the global internal energy field

void sum_ie(double &ie,VD const &e,VD const &dinit,Mesh const &M,VVD const &xinit,VVD const &x,Shape const &S,Shape const &T,VD &detJ0,VVVD &detDJ0,VD &detJ,VVVD &detDJ){

// density as a function

  vector<double> d(T.ngi());

// sum global kinetic field field

  ie=0.0;

  for(int i=0;i<M.NCells();i++){

// update jacobians

    jacobian(i,xinit,M,S,detJ0,detDJ0);
    jacobian(i,x,M,S,detJ,detDJ);

// density at the Gauss points

    for(int gi=0;gi<T.ngi();gi++){d.at(gi)=(dinit.at(i)*detJ0.at(gi)/detJ.at(gi));}

// loop over element energy field

    for(int iloc=0;iloc<T.nloc();iloc++){

// compute nodal mass - this should not change during the calculation

      double nodmass(0.0);
      for(int gi=0;gi<T.ngi();gi++){nodmass+=T.value(iloc,gi)*d.at(gi)*T.wgt(gi);}

      ie+=nodmass*e.at(M.GlobalNode_DFEM(i,iloc));

    }


  }

  return;

}

// function to insert boundary conditions into the mass matrix via row elimination

void bc_insert(Matrix &A,Mesh const &M,Shape const &S,VD const &d,VD const &detJ,VVD &u0,VVD &u1,VD &b0,VD &b1,long const &nknodes){

// loop over boundary elements and choose what type of boundary needs to be applied

  cout<<"There are "<<M.NSides()<<" element sides on the mesh boundary:"<<endl;

// initialise boundary vectors

  fill(b0.begin(),b0.end(),0.0); // value on the boundary
  fill(b1.begin(),b1.end(),0.0); // eliminated row

// find cell sides coincident with the edges of the mesh

  for(long i=0;i<M.NCells();i++){
    for(int iside=0;iside<M.NVertices(i);iside++){

      if(M.E2E(i,iside)<M.NCells()){continue;}

// side iside is on a mesh edge, impose boundary conditions

      string bcname;

      int idim=(iside==0||iside==2)?1:0; // direction perpendicular to mesh boundary

// acquire local node numbers on side iside

      int local_node[S.order()+1];

      switch(iside){

        case(0): // bottom edge of mesh

          for(int iloc=0;iloc<S.order()+1;iloc++){local_node[iloc]=iloc;}

          break;

        case(1): // right edge of mesh

          for(int iloc=0;iloc<S.order()+1;iloc++){local_node[iloc]=iloc*(S.order()+1)+S.order();}

          break;

        case(2): // top edge of mesh

          for(int iloc=0;iloc<S.order()+1;iloc++){local_node[iloc]=S.nloc()-S.order()-1+iloc;}

          break;

        case(3): // left edge of mesh

          for(int iloc=0;iloc<S.order()+1;iloc++){local_node[iloc]=iloc*(S.order()+1);}

          break;

      }

// select boundary condition type on side iside

      switch(M.bc_edge(iside)){

        case(VACUUM):

// do nothing so mesh expands into the void

          bcname="vacuum";

          break;

        case(REFLECTIVE):

// set v.n=0 on domain boundary and impose a constraint on the acceleration field

          bcname="reflective";

          break;

        case(VELOCITY):

// set v.n=<value> on domain boundary and impose a constraint on the acceleration field

          bcname="velocity";

// set v.n=<value> on domain boundary

          for(int isloc=0;isloc<S.order()+1;isloc++){
            long k(M.GlobalNode_CFEM(i,local_node[isloc]));
            u0.at(idim).at(k)=M.bc_value(iside);
            u1.at(idim).at(k)=M.bc_value(iside);
          }

// eliminate k'th solution as we are imposing a condition on it

          for(int isloc=0;isloc<S.order()+1;isloc++){

            long k(M.GlobalNode_CFEM(i,local_node[isloc])); // boundary node
            double bval(1.0e-200); // boundary value

// collect known information - is this double counting ??

            for(long irow=0;irow<nknodes;irow++){b1.at(idim*nknodes+irow)+=A.read(idim*nknodes+irow,idim*nknodes+k)*bval;}

// store boundary value

            b0.at(idim*NROWS+k)=bval;

          }

// modify mass matrix to restore symmetry, this is to try and avoid costly changes to the solution strategy

          for(int isloc=0;isloc<S.order()+1;isloc++){
            long k(M.GlobalNode_CFEM(i,local_node[isloc])); // boundary node
            for(long l=0;l<M.NDims()*nknodes;l++){
              if(l!=k){
                A.write(l,idim*nknodes+k,0.0);
                A.write(idim*nknodes+k,l,0.0);
              }
            }
            A.write(idim*nknodes+k,idim*nknodes+k,1.0);
          }

          break;

        case(FLUID):

// add in additional boundary fluid masses to give zero pressure gradient across the edges of the domain

          bcname="fluid";

          break;

        case(ACCELERATION):

// impose a.n=<value> on domain boundary via row elimination of the mass matrix

          bcname="acceleration";

          break;

        default:

          bcname="undefined";

          cout<<"bc_insert(): "<<bcname<<" boundary conditions not coded, stopping."<<endl;

          exit(1);

      }

      cout<<"Edge "<<iside<<" boundary type : "<<bcname<<endl;
    }
  }

  return;

}

// type safe function to return the sign of the argument

//template <typename T> int sgn(T val) {return(T(0)<val)-(val<T(0));} // -1,0 or 1
template <typename T> int sgn(T val) {return( (val>=T(0))?T(1):T(-1));} // -1 or 1
