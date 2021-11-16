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

#define DTSTART 0.0005    // insert a macro for the first time step
#define ENDTIME 0.2       // insert a macro for the end time
#define ECUT 1.0e-8       // cut-off on the energy field
#define NSAMPLES 1000     // number of sample points for the exact solution
//#define VISFREQ 200     // frequency of the graphics dump steps
//#define OUTFREQ 50      // frequency of the output print steps
#define VISFREQ 0.01      // frequency of the graphics dump times
#define OUTFREQ 0.04      // frequency of the output print times
#define VD vector<double> // vector of doubles
#define VVD vector<VD>    // vector of VD
#define VVVD vector<VVD>  // vector of VVD
#define VI vector<int>    // vector of ints
#define VTOL 1.0e-10      // threshold for volume errors
#define COURANT 0.333     // Courant number for CFL condition
#define DTSFACTOR 0.5     // safety factor on time-step control
#define NROWS nnodes      // number of rows in the global matrix
#define NCOLS nnodes      // number of columns in the global matrix
#define ROW M.Vertex(i,iloc)  // row address in global matrix
#define COL M.Vertex(i,jloc)  // column address in global matrix
#define VACUUM 1              // vacuum boundary
#define REFLECTIVE 2          // reflective boundary
#define FREE_SURFACE 3        // free surface (transmissive) boundary

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
double length(Mesh const &M,Shape const &S,int const i);                                          // element i length scale
void bc_insert(Matrix &A,Mesh const &M,Shape const &S,VD const &d,VD const &detJ);          // iinsert boundary conditions by modifying the mass matrix
void jacobian(int const &i,VVD const &x,Mesh const &M,Shape const &S,VD &detJ,VVVD &detDJ); // calculate a jacobian and determinant
void initial_data(int const n, int const nnodes, int const ndims, int const nmats,          // echo some initial information
                  Mesh const &M);
void lineouts(Mesh const &M, Shape const &S, VD const &d,VD const &p,VD const &e,           // line-outs
              VD const &q, VVD const &x, VVD const &u);
void silo(VVD const &x, VD const &d,VD const &p,VD const &e,VD const &q,VD const &c,        // silo graphics output
          VVD const &u,VI const &mat,int s, double t,Mesh const &M);
void state_print(int const n,int const ndims, int const nmats, VI const &mat,               // output material states
                  VD const &d, VD const &V, VD const &m, VD const &e, VD const &p, 
                  VVD const &x, VVD const &u, int const &s, double const &t);

using namespace std;

int main(){

  cout<<"main(): Starting up main loop..."<<endl;

// global data

  Mesh M("mesh/input-40cells.mesh");                             // load a new mesh from file
  Shape S(1);                                                    // load a p1 shape function
  ofstream f1,f2,f3;                                             // files for output
  int const n(M.NCells()),ndims(M.NDims());                      // no. ncells and no. dimensions
  int const nnodes(M.NNodes());                                  // no. nodes in the mesh
  int const nmats(M.NMaterials());                               // number of materials
  double const cl(0.3),cq(1.0);                                  // linear & quadratic coefficients for bulk viscosity
  Matrix KMASS(NROWS),KMASSI(NROWS);                             // mass matrix for kinematic field
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
  double ke(0.0),ie(0.0);                                        // kinetic and internal energy for conservation checks
  double time(0.0),dt(DTSTART);                                  // start time and time step
  int step(0);                                                   // step number

  vector<vector<double> > state={{1.000, 0.000,0.000, 1.000},    // initial flux state in each material for Sod's shock tube 
                                 {0.125, 0.000,0.000, 0.100}};   // where each flux state is in the form (d,ux,uy,p)

//  vector<vector<double> > state={{1.000,-2.000,0.000, 0.400},  // initial flux state in each material for the 123 problem 
//                                 {1.000, 2.000,0.000, 0.400}}; // where each flux state is in the form (d,ux,uy,p)

//  vector<vector<double> > state={{1.000,0.000,0.000, 1000.0},  // initial flux state in each material for the blast wave
//                                 {1.000,0.000,0.000, 0.0100}}; // where each flux state is in the form (d,ux,uy,p)

//  vector<vector<double> > state={{1.000, 0.000,0.000, 1.000},  // initial flux state in each material for vacuum boundary test
//                                 {1.000, 0.000,0.000, 1.000}}; // where each flux state is in the form (d,ux,uy,p)

  double l[3]={state[0][0],state[0][1],state[0][3]};             // left flux state for input to the Riemann solvers
  double r[3]={state[1][0],state[1][1],state[1][3]};             // right flux state for input to the Riemann solvers

// initialise the problem

  M.InitCoords(x0);                                                 // set initial coordinates
  M.InitCoords(x1);                                                 // set initial coordinates

  for(int i=0;i<n;i++){mat.at(i)=M.Material(i);}
  for(int i=0;i<n;i++){V0.at(i)=M.Volume(i);}
  for(int i=0;i<n;i++){V1.at(i)=M.Volume(i);}
  for(int i=0;i<n;i++){int imat(mat[i]);p.at(i)=state[imat-1][3];}  // load pressure field from initial flux state
  for(int i=0;i<n;i++){int imat(mat[i]);d0.at(i)=state[imat-1][0];} // load density field from initial flux state
  for(int i=0;i<n;i++){int imat(mat[i]);d1.at(i)=state[imat-1][0];}
  for(int i=0;i<n;i++){m.at(i)=d0[i]*V0[i];}                        // calculate the initial mass field
  for(int i=0;i<n;i++){e0.at(i)=E(d0[i],p[i]);}                     // invert the eos to start the energy field
  for(int i=0;i<n;i++){e1.at(i)=E(d1[i],p[i]);}
  for(int i=0;i<n;i++){q.at(i)=0.0;}
  for(int i=0;i<n;i++){c.at(i)=sqrt(GAMMA*p[i]/d0[i]);}

// set boundary condition on the edges of the mesh

  M.bc_set(REFLECTIVE);                                             // set boundary condition on bottom edge of mesh
  M.bc_set(VACUUM);                                                 // set boundary condition on right edge of mesh
  M.bc_set(REFLECTIVE);                                             // set boundary condition on top edge of mesh
  M.bc_set(VACUUM);                                                 // set boundary condition on left edge of mesh

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
    u0.at(idim)=vtmp;
    u1.at(idim)=vtmp;
  }

// display the header so we have a date and time stamp in the output

  header();

// echo some initial information

  initial_data(n,nnodes,ndims,nmats,M);

// assemble mass matrix for acceleration field

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
//          nx+=S.dvalue(0,iloc,gi)*(dy/2.0)*S.wgt(gi);                     // not used: left in as a diagnostic
//          ny+=S.dvalue(1,iloc,gi)*(dx/2.0)*S.wgt(gi);                     // not used: left in as a diagnostic
        }
        KMASS.add(ROW,COL,nn);
      }
    }
  }

// impose boundary constraints via row elimination

  bc_insert(KMASS,M,S,d0,detJ);

// invert the mass matrix

  KMASSI.inverse2(&KMASS); // lapack drivers dgetrf_ and dgetri_

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

    dt=(*min_element(dts.begin(), dts.end()));
//    dt=DTSTART;cout<<"DT HARDWIRED !! "<<endl;

    cout<<fixed<<setprecision(5)<<"  step "<<step<<" time= "<<time<<" dt= "<<dt<<fixed<<setprecision(5)<<" energy (i/k/tot)= "<<ie<<" "<<ke<<" "<<ie+ke<<endl;

// graphics output

    if(abs(remainder(time,VISFREQ))<1.0e-12){
      state_print(n,ndims,nmats,mat,d0,V0,m,e0,p,x0,u0,step,time);
      silo(x0,d0,p,e0,q,c,u0,mat,step,time,M);
    }else{
      if(abs(remainder(step,VISFREQ))==0){
        state_print(n,ndims,nmats,mat,d0,V0,m,e0,p,x0,u0,step,time);
        silo(x0,d0,p,e0,q,c,u0,mat,step,time,M);
      }
    }

// move the nodes to their full-step position

    for(int idim=0;idim<ndims;idim++){
      for(int i=0;i<nnodes;i++){
        x1.at(idim).at(i)=x0.at(idim)[i]+u0.at(idim)[i]*dt;
      }
    }

// update cell volumes at the full-step

    for(int i=0;i<n;i++){

// update the jacobian

      jacobian(i,x1,M,S,detJ,detDJ);

// reset the cell volume

      V1.at(i)=0.0;
      for(int gi=0;gi<S.ngi();gi++){
        V1.at(i)+=detJ[gi]*S.wgt(gi);
      }

      if(V1.at(i)<VTOL){cout<<"-'ve volume detected in cell "<<i<<endl;lineouts(M,S,d1,p,e1,q,x1,u1);exit(1);}

    }

// update cell density at the full-step

    for(int i=0;i<n;i++){d1.at(i)=m[i]/V1[i];}

// update cell energy at the full-step

    for(int i=0;i<n;i++){e1.at(i)=max(ECUT,e0[i]-((p[i]+q[i])*(V1[i]-V0[i]))/m[i]);}

// update cell pressure at the full-step

    for(int i=0;i<n;i++){p.at(i)=P(d1[i],e1[i]);if(p[i]<0.0){cout<<"-'ve pressure detected in cell "<<i<<" e1= "<<e1[i]<<endl;exit(1);}}

// bulk q

    for(int i=0;i<n;i++){
      c.at(i)=sqrt(GAMMA*p[i]/d1[i]);
      double l(length(M,S,i)),divu((d0[i]-d1[i])/(d1[i]*dt));
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

  for(int idim=0;idim<M.NDims();idim++){
    double b[nnodes];for(int i=0;i<nnodes;i++){b[i]=0.0;}
    for(int i=0;i<n;i++){
      for(int iloc=0;iloc<S.nloc();iloc++){
        for(int gi=0;gi<S.ngi();gi++){
          b[ROW]+=(p[i]+q[i])*detDJ[idim][iloc][gi]*detJ[gi]*S.wgt(gi);
if(idim==1){b[ROW]=0.0;} // reflective
        }
      }
    }

if(idim==0){
b[0]=0.0;b[40]=0.0;
b[41]=0.0;b[81]=0.0;
}


// solve global system

    for(int i=0;i<nnodes;i++){
      double x(0.0);
      for(int j=0;j<nnodes;j++){
        x+=KMASSI.read(i,j)*b[j];
      }

// advance the solution

      u1.at(idim).at(i)=u0.at(idim).at(i)+x*dt;

    }

  }






// old code:

//  Matrix A(n+3);double b[n+3],x[n+3];for(int i=0;i<n+3;i++){b[i]=0.0;x[i]=0.0;}
//  double m[2][2]={};m[0][0]=0.5;m[1][1]=0.5;// for mass lumping

// next block codes for matrix assembly with a continuous Galerkin finite element type

//  for(int iel=0;iel<ng;iel++){
//    for(int iloc=0;iloc<S.nloc();iloc++){
//      int i(iel+iloc); // column address in the global matrix
//      if((i>0&&i<ng)){for(int gi=0;gi<S.ngi();gi++){b[i-1]+=(p[iel]+q[iel])*S.dvalue(0,iloc,gi)*S.wgt(gi);}} // integrate the shape derivative for rhs
//      for(int jloc=0;jloc<S.nloc();jloc++){
//        double nn(0.0); // mass matrix
//        int j(iel+jloc); // row address in the global matrix
//        for(int gi=0;gi<S.ngi();gi++){
//          nn+=S.value(iloc,gi)*S.value(jloc,gi)*S.wgt(gi)*0.5*(x1[iel+1]-x1[iel]); // DG & notes use this - double check ??
//        }
//        if((i>0&&i<ng)&&(j>0&&j<ng)){
//          A.add(i-1,j-1,d1[iel]*nn);
////          A.add(i-1,j-1,d1[iel]*(x1[iel+1]-x1[iel])*m[iloc][jloc]); // use mass lumping
//        }
//      }
//    }
//
//  }

// solve global system

//  A.solve(x,b);

// advance the solution

//  for(int i=0;i<n+3;i++){u1.at(i+1)=u0[i+1]+x[i]*dt;}

// impose boundary constraints on the acceleration field

//  u1.at(1)=u1[2];u1.at(ng-1)=u1[ng-2];
//  u1.at(0)=u1[1];u1.at(ng)=u1[ng-1];






























// advance the time step

    time+=dt;
    step++;

// advance the states for the new time step

    u0=u1;x0=x1;e0=e1;V0=V1;d0=d1;
    for(int i=0;i<n;i++){c.at(i)=sqrt(GAMMA*(p[i]+q[i])/d1[i]);}

// debug
//  if(step==230){
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
  vector<double> posline={0.25}; // datum on dimension linedim
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
      case(FREE_SURFACE):
        bcname="free-surface";
        break;
      case(REFLECTIVE):
        bcname="forced-reflective";
        break;
    }

    cout<<"Edge "<<i<<" boundary type : "<<bcname<<endl;
  }

  return;

}

// material state - args are passed by ref to prevent a copy, and are const to prevent changes

void state_print(int const n,int const ndims,int const nmats,VI const &mat,VD const &d,VD const &V,
                  VD const &m,VD const &e,VD const &p,VVD const &x, VVD const &u, int const &s, double const &t){

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
    pmat.at(i)=P(dmat[i],iemat[i]/mmat[i]);
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

  cout<<"total "<<dtot<<" "<<mtot<<" "<<vtot<<" "<<P(dtot,ietot)<<" "<<ietot/mtot<<" "<<ietot<<" "<<ketot<<" "<<tetot<<endl;

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

// insert boundary masses into the global matrix if required

  void bc_insert(Matrix &A,Mesh const &M,Shape const &S,VD const &d,VD const &detJ){

// loop over boundary elements (note these are refered to here as "sides")
// and choose what type of boundary to apply

    for(int ib=0;ib<M.NSides();ib++){

      if(M.bc_edge(M.SideAttr(ib))!=VACUUM){

// material on this face - add in boundary masses

        int i(M.E2E(M.NCells()+ib,0)); // physical element connecting

        for(int iloc=0;iloc<M.NSideNodes(ib);iloc++){
          for(int jloc=0;jloc<M.NSideNodes(ib);jloc++){
            double nn(0.0); // mass matrix
            for(int gi=0;gi<S.ngi();gi++){
              nn+=d[i]*S.value(iloc,gi)*S.value(jloc,gi)*detJ[gi]*S.wgt(gi);
            }
            A.add(M.SideNode(ib,iloc),M.SideNode(ib,jloc),nn); // add boundary mass on to existing value in the mass matrix
          }
        }

      }
    }

    return;

  }
