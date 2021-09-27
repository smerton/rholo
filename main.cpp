// Finite element variant of RhoLo (Riemann-based Hydro in One-dimension at Low Order - RhoLo)
// RhoLo is an ultra simple 1-D discontinuous finite element (DG) hydrodynamics test code
// This finite element variant solves the Euler equations in their non-conservative form in the fluid frame (the Lagrangian frame)
// using a mixed continuous finite element method (cell-centred thermodynamic variable d,rho,e with node 
// centred kinematic variables u,a) and bulk viscosity q to increase entropy across element boundaries, initial 
// implementation is only first order in time

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

// for graphics: convert -density 300 filename.png filename.pdf

// Author S. R. Merton

#define DTSTART 0.0005    // insert a macro for the first time step
#define ENDTIME 0.20      // insert a macro for the end time
#define GAMMA 1.4         // ratio of specific heats for ideal gases
#define ECUT 1.0e-8       // cut-off on the energy field
#define NSAMPLES 500      // number of sample points for the exact solution
#define VISFREQ 10000     // frequency of the graphics dumps
#define VD vector<double> // vector of doubles
#define VTOL 1.0e-10      // threshold for volume errors
#define COURANT 0.333     // Courant number for CFL condition
#define DTSFACTOR 0.5     // safety factor on time-step control

#include <iostream>
#include <vector>
#include <iomanip>
#include <cmath>
#include <fstream>
#include "riemann.h"
#include "matrix.h"
#include "shape.h"
#include <bits/stdc++.h>

// sigantures for eos lookups

double P(double d,double e); // eos returns pressure as a function of energy
double E(double d,double p); // invert the eos to get energy if we only have pressure
void vempty(vector<double>&v); // signature for emptying a vector

using namespace std;

int main(){

  cout<<"main(): Starting up main loop..."<<endl;

// global data

  ofstream f1,f2,f3,f4,f5;                              // files for output
  int const n(320),ng(n+4);                             // no. ncells, no. ghosts
  Shape S(1);                                           // p1 shape function
  double const cl(0.3),cq(1.0);                         // linear & quadratic coefficients for bulk viscosity
  vector<double> d0(ng),d1(ng),V0(ng),V1(ng),m(ng);     // density, volume & mass
  vector<double> e0(ng),e1(ng);                         // cell-centred energy field
  vector<double> c(ng),p(ng),q(ng);                     // element sound speed, pressure and bulk viscosity
  vector<double> u0(ng+1),u1(ng+1),utmp(ng+1),u2(ng+1); // node velocity
  vector<double> x0(ng+1),x1(ng+1);                     // node coordinates
  vector<double> dt_cfl(ng);                            // element time-step
  double ke(0.0),ie(0.0);                               // kinetic and internal energy for conservation checks
  double time(0.0),dt(DTSTART);                         // start time and time step
  int step(0);                                          // step number
  double l[3]={1.0,0.0,1.0},r[3]={0.125,0.0,0.1};       // left/right flux states for the problem
//  double l[3]={1.0,-2.0,0.4},r[3]={1.0,2.0,0.4};      // left/right flux states for the problem
//  double l[3]={1.0,0.0,1000.0},r[3]={1.0,0.0,0.01};   // left/right flux states for the problem

// initialise the problem

  double dx(1.0/n);x0.at(0)=-2.0*dx;x1.at(0)=x0[0];
  for(int i=1;i<ng+1;i++){x0.at(i)=x0[i-1]+dx;x1.at(i)=x0[i];}
  for(int i=0;i<ng;i++){p.at(i)=(0.5*(x0[i]+x0[i+1])<=0.5)?l[2]:r[2];}
  for(int i=0;i<ng;i++){d0.at(i)=(0.5*(x0[i]+x0[i+1])<=0.5)?l[0]:r[0];}
  for(int i=0;i<ng;i++){d1.at(i)=(0.5*(x1[i]+x1[i+1])<=0.5)?l[0]:r[0];}
  for(int i=0;i<ng;i++){e0.at(i)=E(d0[i],p[i]);}
  for(int i=0;i<ng;i++){e1.at(i)=E(d1[i],p[i]);}
  for(int i=0;i<ng;i++){V0.at(i)=x0[i+1]-x0[i];}
  for(int i=0;i<ng;i++){V1.at(i)=x1[i+1]-x1[i];}
  for(int i=0;i<ng;i++){m.at(i)=d0[i]*V0[i];}
  for(int i=0;i<ng+1;i++){u0.at(i)=(x0[i]<=0.5)?l[1]:r[1];}
  for(int i=0;i<ng+1;i++){u1.at(i)=u0[i];}
  for(int i=0;i<ng+1;i++){utmp.at(i)=u0[i];}
  for(int i=0;i<ng;i++){q.at(i)=0.0;}
  for(int i=0;i<ng;i++){c.at(i)=sqrt(GAMMA*p[i]/d0[i]);}
  for(int i=1;i<ng;i++){ke+=0.25*(m[i-1]+m[i])*u0[i]*u0[i];};ke+=0.25*m[0]*u0[0]*u0[0];ke+=0.25*m[ng-1]*u0[ng]*u0[ng];
  for(int i=0;i<ng;i++){ie+=e1[i]*m[i];}

// start the Riemann solvers from initial flux states

  Riemann R0(Riemann::exact,l,r),R1(Riemann::exact,l,r);

// set output precision

  cout<<fixed<<setprecision(17);

// time integration

  while(time<=ENDTIME){

// calculate a new stable time-step

    for(int i=0;i<ng;i++){double l(x0[i+1]-x0[i]);dt_cfl.at(i)=(COURANT*l/sqrt((c[i]*c[i])+2.0*q[i]/d0[i]));} // impose the CFL limit on each element
    double dt=DTSFACTOR*(*min_element(dt_cfl.begin(), dt_cfl.end())); // reduce across element and apply a saftey factor

    cout<<fixed<<setprecision(5)<<"  step "<<step<<" time= "<<time<<" dt= "<<dt<<fixed<<setprecision(5)<<" energy (i/k/tot)= "<<ie<<" "<<ke<<" "<<ie+ke<<endl;

// move the nodes to their full-step position

    for(int i=0;i<ng+1;i++){x1.at(i)=x0[i]+u0[i]*dt;}

// update kinetic energy for conservation checks

    ke=0.0;for(int i=1;i<ng;i++){ke+=0.25*(m[i-1]+m[i])*u0[i]*u0[i];} // include level 1 halo
//    ke+=0.25*m[0]*u0[0]*u0[0];ke+=0.25*m[ng-1]*u0[ng]*u0[ng]; // update level 2 halo

// evolve the Riemann problems to the end of the time-step on the end of time-step meshes

    vector<double> r0x,rx;vempty(r0x); // sample point coordinates
//    for(long i=0;i<NSAMPLES;i++){r0x.push_back(x1[0]+(i*(x1[ng]-x1[0])/double(NSAMPLES)));} // sample points
    for(long i=0;i<NSAMPLES;i++){r0x.push_back(0.0+(i*(1.0-0.0)/double(NSAMPLES)));} // sample points
    R0.profile(&r0x,time+dt); // Riemann solution at the sample points along the mesh
    rx.clear();for(int i=0;i<ng;i++){rx.push_back(x1[i]);rx.push_back(0.5*(x1[i]+x1[i+1]));rx.push_back(x1[i+1]);}
    R1.profile(&rx,time+dt);

// update cell volumes at the full-step

    for(int i=0;i<ng;i++){V1.at(i)=x1[i+1]-x1[i];if(V1[i]<VTOL){cout<<"-'ve volume detected in cell "<<i<<endl;exit(1);}}

// update cell density at the full-step

    for(int i=0;i<ng;i++){d1.at(i)=m[i]/V1[i];}

// debug
//    for(int i=0;i<ng;i++){d1.at(i)=R1.density(3*i+1);} // 3*i+1 is cell-centre address
// debug

// update cell energy at the full-step

    for(int i=0;i<ng;i++){e1.at(i)=max(ECUT,e0[i]-((p[i]+q[i])*(V1[i]-V0[i]))/m[i]);}

// FDS

//    for(int i=1;i<ng;i++){
//      double de(((p[i-1]+q[i-1]-p[i]-q[i])/(0.5*(x1[i+1]-x1[i-1])))*u0[i]*dt);
//      double vol(0.5*(V1[i-1]+V1[i]));
//      utmp.at(i-1)-=de*vol;
//      utmp.at(i)+=de*vol;
//      e1.at(i-1)=max(ECUT,e0[i-1]-(de*vol)/m[i]);
//      e1.at(i)=max(ECUT,e0[i]+(de*vol)/m[i]);
//      cout<<fixed<<setprecision(17)<<0.5*(x1[i]+x1[i+1])<<" "<<(p[i]+q[i])*(V1[i]-V0[i])<<" "<<el*vl+er*vr<<endl;
//      cout<<fixed<<setprecision(17)<<0.5*(x1[i]+x1[i+1])<<" "<<(p[i]+q[i])*(V1[i]-V0[i])<<" "<<de*vol<<endl;
//      cout<<fixed<<setprecision(17)<<i<<" "<<el*vl<<" "<<er*vr<<endl;
//    }

// update internal energy for conservation checks

    ie=0.0;for(int i=0;i<ng;i++){ie+=e1[i]*m[i];}

// debug
//    for(int i=0;i<ng;i++){e1.at(i)=R1.energy(3*i+1);} // 3*i+1 is cell-centre address
// debug

// update cell pressure at the full-step

    for(int i=0;i<ng;i++){p.at(i)=P(d1[i],e1[i]);if(p[i]<0.0){cout<<"-'ve pressure detected in cell "<<i<<" e1= "<<e1[i]<<endl;exit(1);}}

// bulk q

    for(int i=0;i<ng;i++){
      c.at(i)=sqrt(GAMMA*p[i]/d1[i]);
      double l(x1[i+1]-x1[i]),divu((d0[i]-d1[i])/(d1[i]*dt));
      if(divu<0.0){
        q.at(i)=d0[i]*l*divu*((cq*l*divu)-cl*c[i]);
      }else{
        q.at(i)=0.0; // turn off q as cell divergence indicates expansion
      }
    }

// assemble acceleration field

  Matrix A(n+3);double b[n+3],x[n+3];for(int i=0;i<n+3;i++){b[i]=0.0;x[i]=0.0;}
  double m[2][2]={};m[0][0]=0.5;m[1][1]=0.5;// for mass lumping

// next block codes for matrix assembly with a continuous Galerkin finite element type

  for(int iel=0;iel<ng;iel++){
    for(int iloc=0;iloc<S.nloc();iloc++){
      int i(iel+iloc); // column address in the global matrix
      if((i>0&&i<ng)){for(int gi=0;gi<S.ngi();gi++){b[i-1]+=(p[iel]+q[iel])*S.dvalue(iloc,gi)*S.wgt(gi);}} // integrate the shape derivative for rhs
      for(int jloc=0;jloc<S.nloc();jloc++){
        double nn(0.0); // mass matrix
        int j(iel+jloc); // row address in the global matrix
        for(int gi=0;gi<S.ngi();gi++){
          nn+=S.value(iloc,gi)*S.value(jloc,gi)*S.wgt(gi)*0.5*(x1[iel+1]-x1[iel]); // DG & notes use this - double check ??
        }
        if((i>0&&i<ng)&&(j>0&&j<ng)){
          A.add(i-1,j-1,d1[iel]*nn);
//          A.add(i-1,j-1,d1[iel]*(x1[iel+1]-x1[iel])*m[iloc][jloc]); // use mass lumping
        }
      }
    }

  }

// solve global system

  A.solve(x,b);

// advance the solution

  for(int i=0;i<n+3;i++){u1.at(i+1)=u0[i+1]+x[i]*dt;}

// impose boundary constraints on the acceleration field

  u1.at(1)=u1[2];u1.at(ng-1)=u1[ng-2];
  u1.at(0)=u1[1];u1.at(ng)=u1[ng-1];

// some output

    f1.open("exact.dat");f2.open("dpe.dat");f3.open("q.dat");f4.open("u.dat");f5.open("r.dat");
    f1<<fixed<<setprecision(17);f2<<fixed<<setprecision(17);f3<<fixed<<setprecision(17);f4<<fixed<<setprecision(17);
    for(int i=0;i<NSAMPLES;i++){f1<<r0x[i]<<" "<<R0.density(i)<<" "<<R0.pressure(i)<<" "<<R0.velocity(i)<<" "<<R0.energy(i)<<endl;}
    for(int i=1;i<NSAMPLES;i++){if(R0.region(i)!=R0.region(i-1)){f5<<r0x[i]<<" -20000.0"<<endl<<r0x[i]<<" 20000.0"<<endl<<r0x[i]<<" -20000.0"<<endl;}} // sample region
    for(int i=2;i<n+2;i++){f2<<0.5*(x1[i]+x1[i+1])<<" "<<d1[i]<<" "<<p[i]<<" "<<e1[i]<<" "<<endl;}
    for(int i=2;i<n+2;i++){f3<<0.5*(x1[i]+x1[i+1])<<" "<<q[i]<<endl;}
    for(int i=2;i<n+3;i++){f4<<x1[i]<<" "<<u1[i]<<endl;}

    f1.close();f2.close();f3.close();f4.close();f5.close();

// advance the time step

    time+=dt;
    step++;

// debug
//    for(int i=0;i<ng;i++){u1.at(i)=R1.velocity(3*i);u1.at(i+1)=R1.velocity(3*i+2);} // 3*i,3*i+2 are nodal address
// debug

// advance the solution for the new time step

    for(int i=0;i<ng+1;i++){u0.at(i)=u1[i];}
    for(int i=0;i<ng+1;i++){x0.at(i)=x1[i];}
    for(int i=0;i<ng;i++){e0.at(i)=e1[i];}
    for(int i=0;i<ng;i++){V0.at(i)=V1[i];}
    for(int i=0;i<ng;i++){d0.at(i)=d1[i];}
//    for(int i=0;i<ng;i++){c.at(i)=sqrt(GAMMA*(p[i]+q[i])/d1[i]);}

  }

// estimate convergence rate in the L1/L2 norms using a Riemann solution as the exact solution

  vector<double> rx;vempty(rx);  // sample points
  double xstart(0.0),xstop(1.0); // sample range - should be whole mesh though bc.s may artificially reduce convergence
  for(long i=0;i<ng+1;i++){if(x1[i]>=xstart&&x1[i]<=xstop){rx.push_back(x1[i]);}}
  R0.profile(&rx,ENDTIME); // evolve a Riemann problem on the sample range to final time level

  double l1(0.0),l2(0.0),l1r(0.0),l2r(0.0),l1d(0.0),l2d(0.0),l1n(0.0),l2n(0.0);
  int ii(0);

// numerator and denominator for each norm

  for(int i=0;i<ng+1;i++){
    if(x1[i]>=xstart&&x1[i]<=xstop){
      double err(abs(R0.velocity(ii)-u1[i])); // absolute error
      l1d+=abs(R0.velocity(ii)); // l1 denominator
      l2d+=R0.velocity(ii)*R0.velocity(ii); // l2 denominator
      l1n+=err*V1[i]; // l1 numerator
      l2n+=err*err*V1[i]; // l2 numerator
      ii++;
    }
  }

// construct norms and relative errors

  l1=l1n/ii; // L1 error norm
  l1r=l1n/l1d; // L1 relative error norm
  l2=sqrt(l2n/ii); // L2 error norm
  l2r=sqrt(l2n/l2d); // L2 relative error norm

  cout<<endl<<fixed<<setprecision(10)<<"  Error Estimates (grid spacing h= "<<(xstop-xstart)/ii<<"):"<<endl;

  cout<<"  L1 norm= "<<l1<<" (relative error= "<<l1r<<")"<<endl;
  cout<<"  L2 norm= "<<l2<<" (relative error= "<<l2r<<")"<<endl;
  cout<<"  No. points sampled= "<<ii<<endl;
  cout<<"  Range sampled= "<<xstart<<","<<xstop<<endl;

  cout<<"Normal termination."<<endl;

  return 0;
}

// return pressure given the energy

double P(double d,double e){return (GAMMA-1.0)*d*e;}

// invert the eos to return energy given the pressure

double E(double d,double p){return p/((GAMMA-1.0)*d);}

// empty a vector

void vempty(vector<double>&v){vector<double> e;v.swap(e);return;}
