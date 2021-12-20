// Finite element variant of RhoLo (Really High Order Lagrangian Operator - RhoLo)
// RhoLo is an ultra simple finite element (DG) hydrodynamics test code
// This finite element variant solves the Euler equations in their non-conservative form in the fluid frame (the Lagrangian frame)
// using a mixed continuous finite element method (cell-centred thermodynamic variables p,rho,e with node 
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
#define ENDTIME 0.5       // insert a macro for the end time
//#define ECUT 1.0e-8       // cut-off on the energy field
#define NSAMPLES 1000     // number of sample points for the exact solution
//#define VISFREQ 200     // frequency of the graphics dump steps
//#define OUTFREQ 50      // frequency of the output print steps
#define VISFREQ 0.01       // frequency of the graphics dump times
#define OUTFREQ 0.01       // frequency of the output print times
#define VD vector<double> // vector of doubles
#define VVD vector<VD>    // vector of VD
#define VVVD vector<VVD>  // vector of VVD
#define VI vector<int>    // vector of ints
#define COURANT 0.333     // Courant number for CFL condition
#define DTSFACTOR 0.5     // safety factor on time-step control
#define NROWS nnodes      // number of rows in the global matrix
#define NCOLS nnodes      // number of columns in the global matrix
#define ROW M.Vertex(i,iloc)  // row address in global matrix
#define COL M.Vertex(i,jloc)  // column address in global matrix

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

// function signatures

string date();                                                                              // function to return the date and time
void header();                                                                              // header part
void vempty(vector<double>&v);                                                              // empty a vector
double length(Mesh const &M,Shape const &S,int const i);                                    // element i length scale
void jacobian(int const &i,VVD const &x,Mesh const &M,Shape const &S,VD &detJ,VVVD &detDJ); // calculate a jacobian and determinant
void initial_data(int const n, int const nnodes, int const ndims, int const nmats,          // echo some initial information
                  Mesh const &M);
void lineouts(Mesh const &M, Shape const &S, VD const &d,VD const &p,VD const &e,           // line-outs
              VD const &q, VVD const &x, VVD const &u);
void silo(VVD const &x, VD const &d,VD const &p,VD const &e,VD const &q,VD const &c,        // silo graphics output
          VVD const &u,VI const &mat,int s, double t,Mesh const &M,VD const &g);
void state_print(int const n,int const ndims, int const nmats, VI const &mat,               // output material states
                  VD const &d, VD const &V, VD const &m, VD const &e, VD const &p, 
                  VVD const &x, VVD const &u, int const &s, double const &t,VD const &gamma);
void bc_insert(Matrix &A,Mesh const &M,Shape const &S,VD const &d,VD const &detJ,           // insert boundary conditions into the mass matrix
               VVD &u0,VVD &u1,VD &b0,VD &b1);
void bc_insert(Mesh const &M,Shape const &S,VD &b,VD const &b0,VD const &p,VD const &q,VVVD const &detDJ,VD const &detJ); // insert boundary conditions on acceleration field
void init_TAYLOR(Mesh const &M,Shape const &S,double const &dpi,VD &d0,VD &d1,VVD &u0,VVD &u1,VD &p,VD &e0,VD &e1,VD const &g,vector<int> const &mat);   // input overides for the Taylor-Green vortex
void init_RAYLEIGH(Mesh const &M,Shape const &S,double const &dpi,VD &d0,VD &d1,VVD &u0,VVD &u1,VD &p,VD &e0,VD &e1,VD const &g,vector<int> const &mat); // input overides for the Rayleigh-Taylor instability
void init_NOH(Mesh const &M,Shape const &S,double const &dpi,VD &d0,VD &d1,VVD &u0,VVD &u1,VD &p,VD &e0,VD &e1,VD const &g,vector<int> const &mat); // input overides for the Noh stagnation shock
void init_SEDOV(Mesh const &M,Shape const &S,double const &dpi,VD &d0,VD &d1,VVD &u0,VVD &u1,VD &p,VD &e0,VD &e1,VD const &g,vector<int> const &mat); // input overides for the Sedov explosion

using namespace std;

int main(){

  cout<<"main(): Starting up main loop..."<<endl;

// global data

  Mesh M("mesh/taylor-green-10x10.mesh");                               // load a new mesh from file
  Shape S(1);                                                    // load a p1 shape function
  ofstream f1,f2,f3;                                             // files for output
  int const n(M.NCells()),ndims(M.NDims());                      // no. ncells and no. dimensions
  int const nnodes(M.NNodes());                                  // no. nodes in the mesh
  int const nmats(M.NMaterials());                               // number of materials
  double const cl(0.3),cq(1.0);                                  // linear & quadratic coefficients for bulk viscosity
  Matrix KMASS(2*NROWS),KMASSI(2*NROWS);                         // mass matrix for kinematic field
  vector<double> d0(n),d1(n),V0(n),V1(n),m(n);                   // density, volume & mass
  vector<double> e0(n),e1(n);                                    // cell-centred energy field
  vector<double> c(n),p(n),q(n);                                 // element sound speed, pressure and bulk viscosity
  vector<vector<double> > u0(ndims),u1(ndims);                   // node velocity
  vector<vector<double> > x0(ndims),x1(ndims);                   // node coordinates
  vector<double> detJ(S.ngi());                                  // determinant of jacobian at each integration point
  vector<vector<vector<double> > > detDJ(ndims);                 // determinant of jacobian for each derivative
  vector<double> dt_cfl(n);                                      // element time-step
  vector<double> dts(2);                                         // time-step for each condition (0=CFL, 1=graphics hits)
  vector<int> mat(n);                                            // element material numbers
  vector<double> b0(2*NROWS),b1(2*NROWS);                        // vectors for boundary conditions (b0=value, b1=eliminated row)
  double ke(0.0),ie(0.0);                                        // kinetic and internal energy for conservation checks
  double time(0.0),dt(DTSTART);                                  // start time and time step
  int step(0);                                                   // step number
  int test_problem(0);                                           // set later for specific test cases that may need some overides
  double dpi(4.0*atan(1.0));                                     // definition of pi to double precision
  vector<double> gamma(M.NMaterials());                          // ratio of specific heats, set from material definition

// initial flux state in each material is in the form (d,ux,uy,p,gamma)

//  test_problem=SOD;                                                    // set overides needed to run this problem
//  vector<vector<double> > state={{1.000, 0.000,0.000, 1.000,5.0/3.0},  // initial flux state in each material for Sod's shock tube 
//                                 {0.125, 0.000,0.000, 0.100,5.0/3.0}};

//  test_problem=SODSOD;                                                 // set overides needed to run this problem
//  vector<vector<double> > state={{1.000, 0.000,0.000, 1.000,5.0/3.0},  // initial flux state in each material for double shock problem 
//                                 {0.125, 0.000,0.000, 0.100,5.0/3.0},
//                                 {1.000, 0.000,0.000, 1.000,5.0/3.0}};

//  test_problem=R2R;                                                    // set overides needed to run this problem
//  vector<vector<double> > state={{1.000,-2.000,0.000, 0.400,1.4},      // initial flux state in each material for the 123 problem 
//                                 {1.000, 2.000,0.000, 0.400,1.4}};

//  test_problem=BLASTWAVE;                                              // set overides needed to run this problem
//  vector<vector<double> > state={{1.000,0.000,0.000, 1000.0,1.4},      // initial flux state in each material for the blast wave
//                                 {1.000,0.000,0.000, 0.0100,1.4}};

//  test_problem=VACUUMBC;                                               // set overides needed to run this problem
//  vector<vector<double> > state={{1.000, 0.000,0.000, 1.000,1.4},      // initial flux state in each material for vacuum boundary test
//                                 {1.000, 0.000,0.000, 1.000,1.4}};

  test_problem=TAYLOR;                                                 // set overides needed to run this problem
  vector<vector<double> > state={{1.000, 0.000,0.000, 1.000,5.0/3.0}}; // initial flux state in each material for Taylor problem

//  test_problem=NOH;                                                    // set overides needed to run this problem
//  vector<vector<double> > state={{1.000, 0.000,0.000, 1.000,5.0/3.0}}; // initial flux state in each material for Noh problem

//  test_problem=SEDOV;                                                  // set overides needed to run this problem
//  vector<vector<double> > state={{1.000, 0.000,0.000, 1.000,1.4}};     // initial flux state in each material for Sedov problem

//  test_problem=TRIPLE;                                                 // set overides needed to run this problem
//  vector<vector<double> > state={{1.000, 0.000,0.000, 1.000,1.5},      // initial flux state in each material for triple-point problem
//                                 {1.000, 0.000,0.000, 0.100,1.4},
//                                 {0.125, 0.000,0.000, 0.100,1.5}};

// acquire adiabatic constant for each material

  for(int imat=0;imat<M.NMaterials();imat++){gamma.at(imat)=state[imat][4];}

// set boundary conditions on the edges of the mesh in the form (side,type,v.n) where side 0,1,2,3 = bottom,right,top,left

  M.bc_set(0,VELOCITY,0.0);  // set boundary condition on bottom edge of mesh
  M.bc_set(1,VELOCITY,0.0);  // set boundary condition on right edge of mesh
  M.bc_set(2,VELOCITY,0.0);  // set boundary condition on top edge of mesh
  M.bc_set(3,VELOCITY,0.0);  // set boundary condition on left edge of mesh

// initialise the problem

  M.InitCoords(x0,EXCLUDE_GHOSTS); // set initial coordinates
  M.InitCoords(x1,EXCLUDE_GHOSTS); // set initial coordinates

  for(int i=0;i<n;i++){mat.at(i)=M.Material(i);}
  for(int i=0;i<n;i++){V0.at(i)=M.Volume(i);}
  for(int i=0;i<n;i++){V1.at(i)=M.Volume(i);}
  for(int i=0;i<n;i++){int imat(mat[i]);p.at(i)=state[imat-1][3];}  // load pressure field from initial flux state
  for(int i=0;i<n;i++){int imat(mat[i]);d0.at(i)=state[imat-1][0];} // load density field from initial flux state
  for(int i=0;i<n;i++){int imat(mat[i]);d1.at(i)=state[imat-1][0];}
  for(int i=0;i<n;i++){m.at(i)=d0[i]*V0[i];}                        // calculate the initial mass field
  for(int i=0;i<n;i++){e0.at(i)=E(d0[i],p[i],gamma[mat[i]-1]);}     // invert the eos to start the energy field
  for(int i=0;i<n;i++){e1.at(i)=E(d1[i],p[i],gamma[mat[i]-1]);}
  for(int i=0;i<n;i++){q.at(i)=0.0;}
  for(int i=0;i<n;i++){c.at(i)=sqrt(gamma[mat[i]-1]*p[i]/d0[i]);}

// allocate a determinant for each derivative

  for(int idim=0;idim<ndims;idim++){
    detDJ.at(idim).resize(S.nloc());
    for(int j=0;j<S.nloc();j++){
      detDJ.at(idim).at(j).resize(S.ngi());
    }
  }

// load velocity fields from initial flux state

  for(int idim=0;idim<ndims;idim++){

    vector<double> vtmp(x0.at(idim).size());

    for(int i=0;i<n;i++){
      for(int iloc=0;iloc<S.nloc();iloc++){
        long j(M.Vertex(i,iloc));
        vtmp.at(j)=(vtmp[j]==0.0)?(state[mat[i]-1][1+idim]):((vtmp[j]!=(state[mat[i]-1][1+idim]))?0.0:(state[mat[i]-1][1+idim]));
      }
    }

// commit to velocty field address spaces

    u0.at(idim)=vtmp;
    u1.at(idim)=vtmp;

  }

// input overides needed to initialise certain test problems

  switch(test_problem){

    case(TAYLOR):

// Taylor Green vortex

      init_TAYLOR(M,S,dpi,d0,d1,u0,u1,p,e0,e1,gamma,mat);

      break;

    case(RAYLEIGH):

// Rayleigh-Taylor instability

      init_RAYLEIGH(M,S,dpi,d0,d1,u0,u1,p,e0,e1,gamma,mat);

      break;

    case(NOH):

// Noh stagnation shock

      init_NOH(M,S,dpi,d0,d1,u0,u1,p,e0,e1,gamma,mat);

    case(SEDOV):

// Sedov expanding shock

      init_SEDOV(M,S,dpi,d0,d1,u0,u1,p,e0,e1,gamma,mat);

      break;

  }

// display the header so we have a date and time stamp in the output

  header();

// echo some initial information

  initial_data(n,nnodes,ndims,nmats,M);

// assemble mass matrix for acceleration field

  for(int idim=0;idim<M.NDims();idim++){
    for(int i=0;i<n;i++){
      jacobian(i,x0,M,S,detJ,detDJ);
      for(int iloc=0;iloc<S.nloc();iloc++){
        for(int jloc=0;jloc<S.nloc();jloc++){
          double nn(0.0),nnx(0.0),nny(0.0),nx(0.0),ny(0.0),nlx(0.0),nly(0.0); // mass matrix
          for(int gi=0;gi<S.ngi();gi++){
            nn+=d0[i]*S.value(iloc,gi)*S.value(jloc,gi)*detJ[gi]*S.wgt(gi);   // mass matrix
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

// impose boundary constraints via row elimination

  bc_insert(KMASS,M,S,d0,detJ,u0,u1,b0,b1);

// invert the mass matrix

  cout<<"Inverting the mass matrix..."<<endl;

  KMASSI.inverse2(&KMASS); // lapack drivers dgetrf_ and dgetri_

  cout<<"Done."<<endl;

// set output precision

  cout<<fixed<<setprecision(17);

// time integration

  while(time<=ENDTIME){

// calculate a new stable time-step

    for(int i=0;i<n;i++){double l(length(M,S,i));dt_cfl.at(i)=(COURANT*l/sqrt((c[i]*c[i])+2.0*q[i]/d0[i]));} // impose the CFL limit on each element

// reduce across element and apply a saftey factor

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

    cout<<fixed<<setprecision(5)<<"  step "<<step<<" time= "<<time<<" dt= "<<dt<<fixed<<setprecision(5)<<" energy (i/k/tot)= "<<ie<<" "<<ke<<" "<<ie+ke<<endl;

// graphics output

    if(abs(remainder(time,VISFREQ))<1.0e-12){
//      state_print(n,ndims,nmats,mat,d0,V0,m,e0,p,x0,u0,step,time,gamma);
      silo(x0,d0,p,e0,q,c,u0,mat,step,time,M,gamma);
    }else{
      if(abs(remainder(step,VISFREQ))==0){
//        state_print(n,ndims,nmats,mat,d0,V0,m,e0,p,x0,u0,step,time,gamma);
        silo(x0,d0,p,e0,q,c,u0,mat,step,time,M,gamma);
      }
    }

// move the nodes to their full-step position

    M.UpdateCoords(x1,u0,dt);

// update cell volumes at the full-step

    M.UpdateVolume(V1,x1,S.order());

// update cell density at the full-step

    M.UpdateDensity(d1,V1,m);

// update cell energy at the full-step

    M.UpdateEnergy(e0,e1,p,q,V0,V1,m);

// update cell pressure at the full-step

    for(int i=0;i<n;i++){p.at(i)=P(d1[i],e1[i],gamma[mat[i]-1]);if(p[i]<0.0){cout<<"-'ve pressure detected in cell "<<i<<" e1= "<<e1[i]<<endl;exit(1);}}

// bulk q

    for(int i=0;i<n;i++){
      c.at(i)=sqrt(gamma[mat[i]-1]*p[i]/d1[i]);
      double l(length(M,S,i)),divu((d0[i]-d1[i])/(d1[i]*dt)); // element length and divergence field
      if(divu<0.0){
        q.at(i)=d1[i]*l*divu*((cq*l*divu)-cl*c[i]);
      }else{
        q.at(i)=0.0; // turn off q as cell divergence indicates expansion
      }
    }

//  for(int i=0;i<n;i++){
//
//    double l(sqrt(V1[i])),divu(0.0);                // length of element and divergence field
//    double dxdu(0.0),dydu(0.0),dxdv(0.0),dydv(0.0); // derivatives of the jacobian
//    vector<vector<double> > detJS(ndims,vector<double>(S.nloc()));  // jacobian for each component of the divergence term
//    c.at(i)=sqrt(GAMMA*p[i]/d1[i]);                 // sound speed

// compute a jacobian and determinant for the element

//    for(int j=0;j<S.nloc();j++){
//      dxdu+=x1.at(0).at(M.Vertex(i,j))*S.dvalue(0,j,0.0,0.0); // dx/du
//      dydu+=x1.at(1).at(M.Vertex(i,j))*S.dvalue(0,j,0.0,0.0); // dy/du
//      dxdv+=x1.at(0).at(M.Vertex(i,j))*S.dvalue(1,j,0.0,0.0); // dx/dv
//      dydv+=x1.at(1).at(M.Vertex(i,j))*S.dvalue(1,j,0.0,0.0); // dy/dv
//    }

//              a1   b3    b1   a3
//    double detj(dxdu*dydv-dydu*dxdv),area(4.0*detj);

//    double ui1(u1.at(0)[M.Vertex(i,0)]);
//    double ui2(u1.at(0)[M.Vertex(i,1)]);
//    double ui3(u1.at(0)[M.Vertex(i,2)]);
//    double ui4(u1.at(0)[M.Vertex(i,3)]);

//    double vi1(u1.at(1)[M.Vertex(i,0)]);
//    double vi2(u1.at(1)[M.Vertex(i,1)]);
//    double vi3(u1.at(1)[M.Vertex(i,2)]);
//    double vi4(u1.at(1)[M.Vertex(i,3)]);

//    double a1(dxdu);
//    double b3(dydv);
//    double b1(dydu);
//    double a3(dxdv);

// determinants for the deriavtives in the divergence term

//    for(int j=0;j<S.nloc();j++){
//      detJS.at(0).at(j)=(dydv*S.dvalue(0,j,0.0,0.0)-dydu*S.dvalue(1,j,0.0,0.0))/detj;
//      detJS.at(1).at(j)=(-dxdv*S.dvalue(0,j,0.0,0.0)+dxdu*S.dvalue(1,j,0.0,0.0))/detj;
//   }

// calculate element divergence field 

////    divu=(ui1*(b1-b3)+ui2*(b1+b3)+ui3*(-b1+b3)+ui4*(-b1-b3))/area;
////    divu+=(vi1*(-a1+a3)+vi2*(-a1-a3)+vi3*(a1-a3)+vi4*(a1+a3))/area;

//    for(int idim=0;idim<ndims;idim++){
//      for(int j=0;j<S.nloc();j++){
//        divu+=u1.at(idim)[M.Vertex(i,j)]*detJS[idim][j];
//      }
//    }

// artificial viscosity term

//    if(divu<0.0){
//      q.at(i)=d1[i]*l*divu*((cq*l*divu)-cl*c[i]);
////      q.at(i)=d1[i]*((cq*divu*divu*area) - cl*sqrt(detj*c[i]))*divu;
//    }else{
//      q.at(i)=0.0; // turn off q as cell divergence indicates expansion
//    }

//  }

// assemble acceleration field

  vector<double> b=vector<double> (M.NDims()*nnodes);for(int i=0;i<M.NDims()*nnodes;i++){b.at(i)=-b1[i];}
  for(int idim=0;idim<M.NDims();idim++){
    for(int i=0;i<n;i++){
      jacobian(i,x1,M,S,detJ,detDJ);
      for(int iloc=0;iloc<S.nloc();iloc++){
        for(int gi=0;gi<S.ngi();gi++){
          b.at(idim*NROWS+ROW)+=(p[i]+q[i])*detDJ[idim][iloc][gi]*detJ[gi]*S.wgt(gi);
        }
      }
    }
  }

// insert boundary conditions

  bc_insert(M,S,b,b0,p,q,detDJ,detJ);

// solve global system

  vector<double> udot=vector<double>(M.NDims()*nnodes);
  for(int i=0;i<M.NDims()*nnodes;i++){
    double x(0.0);
    for(int j=0;j<M.NDims()*nnodes;j++){
      x+=KMASSI.read(i,j)*b[j];
    }
    udot.at(i)=x;
  }

// advance the solution

//  for(int i=0;i<nnodes;i++){u1.at(idim).at(i)=u0.at(idim).at(i)+udot.at(i)*dt;}

  for(int idim=0;idim<M.NDims();idim++){
    for(int i=0;i<nnodes;i++){
      u1.at(idim).at(i)=u0.at(idim).at(i)+udot.at(idim*NROWS+i)*dt;
    }
  }

// advance the time step

    time+=dt;
    step++;

// advance the states for the new time step

    u0=u1;x0=x1;e0=e1;V0=V1;d0=d1;
    for(int i=0;i<n;i++){c.at(i)=sqrt(gamma[mat[i]-1]*(p[i]+q[i])/d1[i]);}

// debug
//  if(step==1){
//    cout<<"debug stop."<<endl;
//    output(); // might want this ??
//    exit(1);
//  }
// debug

  }

// some output

  lineouts(M,S,d1,p,e1,q,x1,u1);

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

// this function codes for some line-outs intended for use in the pseudo-1D tests

void lineouts(Mesh const &M,Shape const &S,VD const &d,VD const &p,VD const &e,VD const &q,VVD const &x,VVD const &u){

// some output

  ofstream f1,f2,f3; // file handles for output

  f1.open("dpe.dat");f2.open("q.dat");f3.open("u.dat");

// set up line-outs

  vector<int> plotdim={0}; // dimension to plot against (0 for x, 1 for y)
  vector<int> linedim={1}; // dimension in which the line-out is located (0 for x, 1 for y)
  vector<double> posline={0.5*(M.Min(1)+M.Max(1))}; // datum on dimension linedim
  vector<string> label={"YMID"}; // label for the line-out

// loop over line-outs

  for(int line=0;line<plotdim.size();line++){

// set the output precision

    f1<<fixed<<setprecision(10);
    f2<<fixed<<setprecision(10);
    f3<<fixed<<setprecision(10);

// write out a header so we know what each file contains

    f1<<"# "<<label[line]<<" line-out at "<<((linedim[line]==0)?"X=":"Y=")<<posline[line]<<": Coord Density Pressure Energy"<<endl;
    f2<<"# "<<label[line]<<" line-out at "<<((linedim[line]==0)?"X=":"Y=")<<posline[line]<<": Coord q"<<endl;
    f3<<"# "<<label[line]<<" line-out at "<<((linedim[line]==0)?"X=":"Y=")<<posline[line]<<": Coord velocity"<<endl;

// sweep mesh and use centroid to find cells along the line-out

    for(int i=0;i<M.NCells();i++){
    
// vertex coordinates of cell i

      double xpos[4],ypos[4],xc[2]={0.0,0.0};
      for(int j=0;j<S.nloc();j++){
        xpos[j]=x[0].at(M.Vertex(i,j));
        ypos[j]=x[1].at(M.Vertex(i,j));
      }

// get cell centroid

      for(int j=0;j<S.nloc();j++){
        xc[0]+=S.value(j,0.0,0.0)*xpos[j];
        xc[1]+=S.value(j,0.0,0.0)*ypos[j];
      }

// check cell coincides with the line-out

      if(abs(xc[linedim[line]]-posline[line])<1.0e-8){

        f1<<xc[plotdim[line]]<<" "<<d.at(i)<<" "<<p.at(i)<<" "<<e.at(i)<<endl;
        f2<<xc[plotdim[line]]<<" "<<q.at(i)<<endl;

// output velocity through the two nodes parallel to the line-out

        if(plotdim[line]==0){
          f3<<x.at(plotdim[line]).at(M.Vertex(i,0))<<" "<<u.at(plotdim[line]).at(M.Vertex(i,0))<<endl;
          f3<<x.at(plotdim[line]).at(M.Vertex(i,1))<<" "<<u.at(plotdim[line]).at(M.Vertex(i,1))<<endl;
        }else{
          f3<<x.at(plotdim[line]).at(M.Vertex(i,0))<<" "<<u.at(plotdim[line]).at(M.Vertex(i,0))<<endl;
          f3<<x.at(plotdim[line]).at(M.Vertex(i,2))<<" "<<u.at(plotdim[line]).at(M.Vertex(i,2))<<endl;
        }

      }

    }

  }

// close the output files

  f1.close();f2.close();f3.close();

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

void initial_data(int const n,int const nnodes,int const ndims, int const nmats, Mesh const &M){

  cout<<"Initial data for the simulation"<<endl;

  cout<<"Number of dimensions:             "<<ndims<<endl;
  cout<<"Number of cells:                  "<<n<<endl;
  cout<<"Number of nodes:                  "<<nnodes<<endl;
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

    for(int j=0;j<S.nloc();j++){
      dxdu+=x.at(0).at(M.Vertex(i,j))*S.dvalue(0,j,gi); // dx/du
      dxdv+=x.at(0).at(M.Vertex(i,j))*S.dvalue(1,j,gi); // dx/dv
      dydu+=x.at(1).at(M.Vertex(i,j))*S.dvalue(0,j,gi); // dy/du
      dydv+=x.at(1).at(M.Vertex(i,j))*S.dvalue(1,j,gi); // dy/dv
    }

// calculate the determinant at the quadrature point and commit to the vector

    detJ.at(gi)=dxdu*dydv-dxdv*dydu;

// determinants for the deriavtives at the quadrature points

    for(int j=0;j<S.nloc();j++){
      detDJ.at(0).at(j).at(gi)=(dydv*S.dvalue(0,j,gi)-dydu*S.dvalue(1,j,gi))/detJ[gi];
      detDJ.at(1).at(j).at(gi)=(-dxdv*S.dvalue(0,j,gi)+dxdu*S.dvalue(1,j,gi))/detJ[gi];
    }

  }

  return;

}

// compute the shortest distance a signal can take to cross element i

double length(Mesh const &M,Shape const &S,int const i){

  VD xydist; // distance between mid-points
  VVD xymid(M.NDims()); // midpoint coordinates for each side

// resize mid-point coordinate vector to hold information for 4 sides

  for(int idim=0;idim<M.NDims();idim++){xymid.at(idim).resize(4);}

// side mid-points using the finite element method

  for(int idim=0;idim<M.NDims();idim++){
    for(int j=0;j<S.nloc();j++){
      xymid[idim][0]+=S.value(j,0.0,-1.0)*M.Coord(idim,M.Vertex(i,j)); // bottom face
      xymid[idim][1]+=S.value(j,1.0,0.0)*M.Coord(idim,M.Vertex(i,j)); // right face
      xymid[idim][2]+=S.value(j,0.0,1.0)*M.Coord(idim,M.Vertex(i,j)); // top face
      xymid[idim][3]+=S.value(j,-1.0,0.0)*M.Coord(idim,M.Vertex(i,j)); // left face
    }
  }

// distances between the mid-points

  for(int iside=0;iside<3;iside++){
    for(int jside=iside+1;jside<4;jside++){
      double dx(abs(xymid[0][iside]-xymid[0][jside]));
      double dy(abs(xymid[1][iside]-xymid[1][jside]));
      xydist.push_back(sqrt(dx*dx+dy*dy));
    }
  }

// return minimum distance between the mid-points

  return (*min_element(xydist.begin(),xydist.end()));

}

// insert boundary conditions into the mass matrix

void bc_insert(Matrix &A,Mesh const &M,Shape const &S,VD const &d,VD const &detJ,VVD &u0,VVD &u1,VD &b0,VD &b1){

// loop over boundary elements and choose what type of boundary needs to be applied

  cout<<"There are "<<M.NSides()<<" elements on the mesh boundary:"<<endl;

  int nnodes(M.NNodes()); // needed for expansion of the NROWS macro

// initialise boundary vectors

  for(int i=0;i<M.NDims()*M.NNodes();i++){
    b0.at(i)=0.0; // value on the boundary
    b1.at(i)=0.0; // eliminated row
  }

// impose boundary constraints via row elimination of the global mass matrix

  for(int ib=0;ib<M.NSides();ib++){

    string bcname;

    int j(M.SideAttr(ib)); // element side coincident with mesh boundary
    int idim=(j==0||j==2)?1:0; // direction perpendicular to mesh boundary

    switch(M.bc_edge(M.SideAttr(ib))){
      case(VACUUM):

// do nothing so mesh expands into the void

        bcname="vacuum";

        break;

      case(REFLECTIVE):

// set v.n=0 on domain boundary and impose a constraint on the acceleration field

        bcname="reflective";

// set v.n=0 on domain boundary

        for(int iloc=0;iloc<M.NSideNodes(ib);iloc++){
          int k(M.SideNode(ib,iloc));
          u0.at(idim).at(k)=0.0;
          u1.at(idim).at(k)=0.0;
        }

// reflect the mass matrix

        for(int jdim=0;jdim<M.NDims();jdim++){
          for(int iloc=0;iloc<M.NSideNodes(ib);iloc++){
            for(int jloc=0;jloc<M.NSideNodes(ib);jloc++){
              double nn(0.0); // mass matrix
              for(int gi=0;gi<S.ngi();gi++){
                nn+=d[M.E2E(M.NCells()+ib,0)]*S.value(iloc,gi)*S.value(jloc,gi)*detJ[gi]*S.wgt(gi);
              }
              A.add(jdim*NROWS+M.SideNode(ib,iloc),jdim*NROWS+M.SideNode(ib,jloc),nn); // add boundary mass to the mass matrix
            }
          }
        }

// eliminate k'th solution as we are imposing a condition on it

        for(int iloc=0;iloc<M.NSideNodes(ib);iloc++){

// boundary node and boundary value

          int k(M.SideNode(ib,iloc));
          double rhs(0.0),bval(1.0e-200); // no net force acting on boundary node

// collect known information

          for(int i=0;i<M.NNodes();i++){rhs+=A.read(i,idim*NROWS+k)*bval;}

// store boundary value

          b0.at(idim*NROWS+k)=bval;

// move known information onto rhs

          b1.at(idim*NROWS+k)=rhs;

        }

// modify mass matrix and restore symmetry

        for(int iloc=0;iloc<M.NSideNodes(ib);iloc++){
          int k(M.SideNode(ib,iloc));
          for(int i=0;i<M.NDims()*M.NNodes();i++){
            if(i!=k){
              A.write(i,idim*NROWS+k,0.0);
              A.write(idim*NROWS+k,i,0.0);
            }
          }
          A.write(idim*NROWS+k,idim*NROWS+k,1.0);
        }

        break;

      case(VELOCITY):

// set v.n=<value> on domain boundary and impose a constraint on the acceleration field

        bcname="velocity";

// set v.n=<value> on domain boundary

        for(int iloc=0;iloc<M.NSideNodes(ib);iloc++){
          int k(M.SideNode(ib,iloc));
          u0.at(idim).at(k)=M.bc_value(M.SideAttr(ib));
          u1.at(idim).at(k)=M.bc_value(M.SideAttr(ib));
        }

// eliminate k'th solution as we are imposing a condition on it

        for(int iloc=0;iloc<M.NSideNodes(ib);iloc++){

// boundary node and boundary value

          int k(M.SideNode(ib,iloc));
          double rhs(0.0),bval(1.0e-200); // no net force acting on boundary node

// collect known information

          for(int i=0;i<M.NNodes();i++){rhs+=A.read(i,idim*NROWS+k)*bval;}

// store boundary value

          b0.at(idim*NROWS+k)=bval;

// move known information onto rhs

          b1.at(idim*NROWS+k)=rhs;

        }

// modify mass matrix and restore symmetry

        for(int iloc=0;iloc<M.NSideNodes(ib);iloc++){
          int k(M.SideNode(ib,iloc));
          for(int i=0;i<M.NDims()*M.NNodes();i++){
            if(i!=k){
              A.write(i,idim*NROWS+k,0.0);
              A.write(idim*NROWS+k,i,0.0);
            }
          }
          A.write(idim*NROWS+k,idim*NROWS+k,1.0);
        }

        break;

      case(FLUID):

// add in additional boundary fluid masses to give zero pressure gradient across the edges of the domain

        bcname="fluid";

        for(int jdim=0;jdim<M.NDims();jdim++){
          for(int iloc=0;iloc<M.NSideNodes(ib);iloc++){
            for(int jloc=0;jloc<M.NSideNodes(ib);jloc++){
              double nn(0.0); // mass matrix
              for(int gi=0;gi<S.ngi();gi++){
                nn+=d[M.E2E(M.NCells()+ib,0)]*S.value(iloc,gi)*S.value(jloc,gi)*detJ[gi]*S.wgt(gi);
              }
              A.add(jdim*NROWS+M.SideNode(ib,iloc),jdim*NROWS+M.SideNode(ib,jloc),nn); // add boundary mass to the mass matrix
            }
          }
        }

        break;

      case(ACCELERATION):

// impose a.n=<value> on domain boundary via row elimination of the mass matrix

        bcname="acceleration";

// eliminate k'th solution as we are imposing a condition on it

        for(int iloc=0;iloc<M.NSideNodes(ib);iloc++){

// boundary node and boundary value

          int k(M.SideNode(ib,iloc));
          double rhs(0.0),bval(M.bc_value(M.SideAttr(ib))); // net force acting on boundary node = m*(<value>).n

// collect known information

          for(int i=0;i<M.NNodes();i++){rhs+=A.read(i,idim*NROWS+k)*bval;}

// store boundary value

          b0.at(idim*NROWS+k)=bval;

// move known information onto rhs

          b1.at(idim*NROWS+k)=rhs;

        }

// modify mass matrix and restore symmetry

        for(int iloc=0;iloc<M.NSideNodes(ib);iloc++){
          int k(M.SideNode(ib,iloc));
          for(int i=0;i<M.NDims()*M.NNodes();i++){
            if(i!=k){
              A.write(i,idim*NROWS+k,0.0);
              A.write(idim*NROWS+k,i,0.0);
            }
          }
          A.write(idim*NROWS+k,idim*NROWS+k,1.0);
        }

        break;

      default:

        bcname="undefined";

        cout<<"bc_insert(): "<<bcname<<" boundary conditions not coded, stopping."<<endl;

        exit(1);

    }

    cout<<"Edge "<<ib<<" boundary type : "<<bcname<<endl;
  }

  return;

}

// insert boundary conditions on acceleration field

void bc_insert(Mesh const &M,Shape const &S,VD &b,VD const &b0,VD const &p,VD const &q,VVVD const &detDJ,VD const &detJ){

  int nnodes(M.NNodes()); // needed for expansion of the NROWS macro

// loop over boundary elements and choose what type of boundary needs to be applied

  for(int ib=0;ib<M.NSides();ib++){

    string bcname;

    int j(M.SideAttr(ib)); // element side coincident with mesh boundary
    int idim=(j==0||j==2)?1:0; // perpendicular direction

// modify source vector accordingly

    switch(M.bc_edge(M.SideAttr(ib))){

      case(VACUUM):

// do nothing so mesh expands into the void

        bcname="vacuum";

        break;

      case(REFLECTIVE):

// reflective boundary so a.n=0.0

        bcname="reflective";

        for(int i=0;i<M.NNodes();i++){if(b0[idim*NROWS+i]!=0.0){b.at(idim*NROWS+i)=-b0[idim*NROWS+i];}}

        break;

      case(VELOCITY):

// velocity has been imposed so prevent any perpendicular acceleration to retain the value

        bcname="velocity";

        for(int i=0;i<M.NNodes();i++){if(b0[idim*NROWS+i]!=0.0){b.at(idim*NROWS+i)=b0[idim*NROWS+i];}}

        break;

      case(ACCELERATION):

// acceleration has been imposed so enforce the value

        bcname="acceleration";

        for(int i=0;i<M.NNodes();i++){if(b0[idim*NROWS+i]!=0.0){b.at(idim*NROWS+i)=b0[idim*NROWS+i];}}

        cout<<"bc_insert(): "<<bcname<<" boundary conditions not coded, stopping."<<endl;

        exit(1);

        break;

      case(FLUID):

// add extra fluid mass on this boundary

        bcname="fluid";

        for(int jdim=0;jdim<M.NDims();jdim++){
          for(int iloc=0;iloc<M.NSideNodes(ib);iloc++){
            double nn(0.0);
            for(int gi=0;gi<S.ngi();gi++){
              nn+=detDJ[jdim][iloc][gi]*detJ[gi]*S.wgt(gi);
            }
//if(jdim==0){cout<<"Adding extra mass to side node "<<M.SideNode(ib,iloc)<<" : ele "<<M.E2E(M.NCells()+ib,0)<<" is copied across segment "<<ib<<endl;}
            b.at(jdim*NROWS+M.SideNode(ib,iloc))+=(p[M.E2E(M.NCells()+ib,0)]+q[M.E2E(M.NCells()+ib,0)])*nn;
          }
        }

        break;

      default:

        bcname="undefined";

        cout<<"bc_insert(): "<<bcname<<" boundary conditions not coded, stopping."<<endl;

        exit(1);
    }

  }

// debug
//  exit(1);
// debug
  return;

}

// input overides for the Taylor-Green vortex

void init_TAYLOR(Mesh const &M,Shape const &S,double const &dpi,VD &d0,VD &d1,VVD &u0,VVD &u1,VD &p,VD &e0,VD &e1,VD const &gamma,vector<int> const &mat){

  cout<<"init_TAYLOR(): Input overides for Taylor-Green..."<<endl;

// start x component of velocity field

  for(long i=0;i<M.NNodes();i++){
    u0.at(0).at(i)=sin(dpi*M.Coord(0,i))*cos(dpi*M.Coord(1,i));
    u1.at(0).at(i)=u0.at(0).at(i);
  }

// start y component of velocity field

  for(long i=0;i<M.NNodes();i++){
    u0.at(1).at(i)=-cos(dpi*M.Coord(0,i))*sin(dpi*M.Coord(1,i));
    u1.at(1).at(i)=u0.at(1).at(i);
  }

// load density field

  for(long i=0;i<M.NCells();i++){
    d0.at(i)=1.0;
    d1.at(i)=d0.at(i);
  }

// load pressure field

  for(long i=0;i<M.NCells();i++){
    double xc[2];xc[0]=0.0;xc[1]=0.0;
    for(int iloc=0;iloc<S.nloc();iloc++){
      xc[0]+=S.value(iloc,0.0,0.0)*M.Coord(0,M.Vertex(i,iloc));
      xc[1]+=S.value(iloc,0.0,0.0)*M.Coord(1,M.Vertex(i,iloc));
    }
    p.at(i)=0.25*d0.at(i)*(cos(2.0*dpi*xc[0])+cos(2.0*dpi*xc[1]))+1.0;
  }

// invert the eos to start the energy field
// energy should = (3.0*dpi/8.0)*cos(3.0*dpi*xc[0])*cos(dpi*xc[1])-cos(dpi*xc[0])*cos(3.0*dpi*xc[1])
// so we can use this as a check

  for(int i=0;i<M.NCells();i++){
    double xc[2];xc[0]=0.0;xc[1]=0.0;
    for(int iloc=0;iloc<S.nloc();iloc++){
      xc[0]+=S.value(iloc,0.0,0.0)*M.Coord(0,M.Vertex(i,iloc));
      xc[1]+=S.value(iloc,0.0,0.0)*M.Coord(1,M.Vertex(i,iloc));
    }
    double echeck((3.0*dpi/8.0)*(cos(3.0*dpi*xc[0])*cos(dpi*xc[1])-cos(dpi*xc[0])*cos(3.0*dpi*xc[1])));
    e0.at(i)=E(d0[i],p[i],gamma[mat[i]-1]);
    e1.at(i)=E(d1[i],p[i],gamma[mat[i]-1]);
  }

  return;

}

// input overides for the Rayleigh-Taylor instability

void init_RAYLEIGH(Mesh const &M,Shape const &S,double const &dpi,VD &d0,VD &d1,VVD &u0,VVD &u1,VD &p,VD const &gamma,vector<int> const &mat){

  cout<<"init_TAYLOR(): Input overides for the Rayleigh-Taylor instability test not coded yet."<<endl;

  exit(1);

  return;

}

// input overides for the Noh stagnation shock

void init_NOH(Mesh const &M,Shape const &S,double const &dpi,VD &d0,VD &d1,VVD &u0,VVD &u1,VD &p,VD &e0,VD &e1,VD const &gamma,vector<int> const &mat){

  cout<<"init_NOH(): Input overides for Noh..."<<endl;

// start x component of velocity field

  for(long i=0;i<M.NNodes();i++){

    double origin[2]={0.0,0.0};            // origin coordinates assuming we are on the domain [-1,1]
    double rx(M.Coord(0,i)-origin[0]);     // radial vector component from domain origin to node
    double ry(M.Coord(1,i)-origin[1]);     // radial vector component from domain origin to node
    double rnorm(sqrt(rx*rx+ry*ry));       // length of radial vector from domain origin to node

// velocity is a radial vector from degree of freedom towards the domain origin

    u0.at(0).at(i)=-rx/max(rnorm,1.0e-12);
    u1.at(0).at(i)=u0.at(0).at(i);

    u0.at(1).at(i)=-ry/max(rnorm,1.0e-12);
    u1.at(1).at(i)=u0.at(1).at(i);

  }

// load density field

  for(long i=0;i<M.NCells();i++){
    d0.at(i)=1.0;
    d1.at(i)=d0.at(i);
  }

// load pressure field

  for(long i=0;i<M.NCells();i++){
    p.at(i)=0.0;
  }

// start the energy field

  for(int i=0;i<M.NCells();i++){
    e0.at(i)=0.0;
    e1.at(i)=e0.at(i);
  }

  return;

}

// input overides for the Sedov explosion

void init_SEDOV(Mesh const &M,Shape const &S,double const &dpi,VD &d0,VD &d1,VVD &u0,VVD &u1,VD &p,VD &e0,VD &e1,VD const &gamma,vector<int> const &mat){

  cout<<"init_SEDOV(): Input overides for Sedov..."<<endl;

// load density field

  for(long i=0;i<M.NCells();i++){
    d0.at(i)=1.0;
    d1.at(i)=d0.at(i);
  }

// load energy field

  for(long i=0;i<M.NCells();i++){

    e0.at(i)=0.0;
    e1.at(i)=e0.at(i);

    for(int iloc=0;iloc<S.nloc();iloc++){

// delta function at domain origin

      if( abs(M.Coord(0,M.Vertex(i,iloc)))<1.0e-7 && abs(M.Coord(1,M.Vertex(i,iloc)))<1.0e-7  ){
        e0.at(i)=1.0;
        e1.at(i)=e0.at(i); 
      }

    }

  }

// load pressure field

  for(long i=0;i<M.NCells();i++){
    p.at(i)=P(d0[i],e0[i],gamma[mat[i]-1]);
  }

  return;

}
