// Finite difference variant of RhoLo (Riemann-based Hydro in One-dimension at Low Order - RhoLo)
// RhoLo is an ultra simple 1-D discontinuous finite element (DG) hydrodynamics test code
// This finite difference variant solves the Euler equations in their non-conservative form in the fluid frame (the Lagrangian frame)
// using a staggered-grid finite difference method (cell-centred thermodynamic variable d,rho,e with node 
// centred kinematic variables u,a) and bulk viscosity q to increase entropy across element boundaries, initial 
// implementation is only first order in time
// for graphics: convert -density 300 filename.png filename.pdf

// Author S. R. Merton

#define DTSTART 0.0005    // insert a macro for the first time step
#define ENDTIME 0.25      // insert a macro for the end time
#define GAMMA 1.4         // ratio of specific heats for ideal gases
#define ECUT 1.0-8        // cut-off on the energy field
#define NSAMPLES 500     // number of sample points for the exact solution
#define VISFREQ 10000       // frequency of the graphics dumps
#define VD vector<double> // vector of doubles

#include <iostream>
#include <vector>
#include <iomanip>
#include <cmath>
#include "riemann.h"
#include "matrix.h"

// sigantures for eos lookups

double P(double d,double e); // eos returns pressure as a function of energy
double E(double d,double p); // invert the eos to get energy if we only have pressure
void vempty(vector<double>&v); // signature for emptying a vector

using namespace std;

int main(){

  cout<<"main(): Starting up main loop..."<<endl;

// global data

  int const n(100),ng(n+4);                               // no. ncells, no. ghosts
  double const cl(1.0),cq(0.1);                         // linear & quadratic coefficients for bulk viscosity
  vector<double> d(ng),p(ng),q(ng),V0(ng),V1(ng),m(ng); // pressure, bulk viscosity, density, volume & mass
  vector<double> e0(ng),e1(ng);                         // cell-centred energy field
  vector<double> u0(ng+1),u1(ng+1);                     // node velocity
  vector<double> x0(ng+1),x1(ng+1);                     // node coordinates
  double time(0.0),dt(DTSTART);                         // start time and time step
  int step(0);                                          // step number
  double l[3]={1.0,0.0,1.0},r[3]={0.125,0.0,0.1};       // left/right flux states for the problem
//  double l[3]={1.0,-2.0,0.4},r[3]={1.0,2.0,0.4};      // left/right flux states for the problem
//  double l[3]={1.0,0.0,1000.0},r[3]={1.0,0.0,0.01};   // left/right flux states for the problem

// initialise the problem

  double dx(1.0/n);x0.at(0)=-2.0*dx;
  for(int i=1;i<ng+1;i++){x0.at(i)=x0[i-1]+dx;}
  for(int i=0;i<ng;i++){p.at(i)=(0.5*(x0[i]+x0[i+1])<=0.5)?l[2]:r[2];}
  for(int i=0;i<ng;i++){d.at(i)=(0.5*(x0[i]+x0[i+1])<=0.5)?l[0]:r[0];}
  for(int i=0;i<ng;i++){e0.at(i)=E(d[i],p[i]);}
  for(int i=0;i<ng;i++){V0.at(i)=x0[i+1]-x0[i];}
  for(int i=0;i<ng;i++){m.at(i)=d[i]*V0[i];}
  for(int i=0;i<ng+1;i++){u0.at(i)=(x0[i]<=0.5)?l[1]:r[1];}
  for(int i=0;i<ng;i++){q.at(i)=0.0;} // bulk viscosity

// start the Riemann solvers from initial flux states

  Riemann R0(Riemann::exact,l,r),R1(Riemann::exact,l,r),R2(Riemann::exact,l,r);

// set output precision

  cout<<fixed<<setprecision(17);

// time integration

  while(time<ENDTIME+dt){

    cout<<"  step "<<step<<" time= "<<time<<" dt= "<<dt<<endl;

// move the nodes to their full-step position

    for(int i=0;i<ng+1;i++){x1.at(i)=x0[i]+u0[i]*dt;}

// evolve the Riemann problems to the end of the time-step on the end of time-step meshes

    vector<double> r0x,rx;vempty(r0x); // sample point coordinates
    for(long i=0;i<NSAMPLES;i++){r0x.push_back(x1[0]+(i*(x1[ng]-x1[0])/double(NSAMPLES)));} // sample points
    R0.profile(&r0x,time+dt); // Riemann solution at the sample points along the mesh
    rx.clear();for(int i=0;i<ng;i++){rx.push_back(0.5*(x1[i]+x1[i+1]));} // cell centres
    R1.profile(&rx,time+dt); // Riemann solution at the cell centres
    rx.clear();for(int i=0;i<ng+1;i++){rx.push_back(x1[i]);} // nodes
    R2.profile(&rx,time+dt); // Riemann solution at the nodes of the mesh

// update cell volumes at the full-step

    for(int i=0;i<ng;i++){V1.at(i)=x1[i+1]-x1[i];} 

// update cell density at the full-step

    for(int i=0;i<ng;i++){d.at(i)=m[i]/V1[i];} 

// debug
    for(int i=0;i<ng;i++){d.at(i)=R1.density(i);} // R1 is the Riemann solution at the cell centres
// debug

// update cell energy at the full-step

    for(int i=0;i<ng;i++){e1.at(i)=max(ECUT,e0[i]-(p[i]+q[i])*(V1[i]-V0[i])/m[i]);}

// debug
    for(int i=0;i<ng;i++){e1.at(i)=R1.energy(i);} // R1 is the Riemann solution at the cell centres
// debug

// update cell pressure at the full-step

    for(int i=0;i<ng;i++){p.at(i)=P(d[i],e1[i]);}

// update acceleration field

    for(int i=2;i<ng-1;i++){

      double udot(0.0);

// advance the solution

//      u1.at(i)=u0[i]+udot*dt;

    }

// debug
    for(int i=0;i<ng+1;i++){u1.at(i)=R2.velocity(i);} // R2 is the Riemann solution at the nodesd
// debug

// impose a constraint on the acceleration field at domain boundaries to stop the mesh taking off

    u1.at(1)=u1[2];u1.at(0)=u1[1];u1.at(ng-1)=u1[ng-2];u1.at(ng)=u1[ng-1];

// some output - toggle this to output either the exact solutions from the Riemann solver or the finite element solution generated by the code

//    for(int i=0;i<NSAMPLES;i++){cout<<r0x[i]<<" "<<R0.density(i)<<" "<<R0.pressure(i)<<" "<<R0.velocity(i)<<" "<<R0.energy(i)<<endl;} // exact solution from Riemann solver
    for(int i=2;i<n+2;i++){cout<<0.5*(x1[i]+x1[i+1])<<" "<<d[i]<<" "<<p[i]<<" "<<e1[i]<<" "<<0.5*(u1[i]+u1[i+1])<<endl;} // finite difference solution
//    for(int i=2;i<n+3;i++){cout<<x1[i]<<" "<<u1[i]<<endl;}


// advance the time step

    time+=dt;
    step++;

// advance the solution for the new time step

    for(int i=0;i<ng+1;i++){u0.at(i)=u1[i];}
    for(int i=0;i<ng+1;i++){x0.at(i)=x1[i];}
    for(int i=0;i<ng;i++){e0.at(i)=e1[i];}
    for(int i=0;i<ng;i++){V0.at(i)=V1[i];}

  }

  cout<<"Normal termination."<<endl;

  return 0;
}

// return pressure given the energy

double P(double d,double e){return (GAMMA-1.0)*d*e;}

// invert the eos to return energy given the pressure

double E(double d,double p){return p/((GAMMA-1.0)*d);}

// empty a vector

void vempty(vector<double>&v){vector<double> e;v.swap(e);return;}
