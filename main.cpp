// main program for rholo
// Riemann Hydro One-dimensional Low Order hydrodynamics test code

// Author S. R. Merton

#define DTSTART 0.002  // insert a macro for the first time step
#define ENDTIME 0.25 // insert a macro for the end time
#define GAMMA 1.4 // ration of specific heats for ideal gases

#include <iostream>
#include <vector>
#include <iomanip>
#include "riemann.h"
#include <cmath>

// sigantures for eos lookups

double P(double d,double e); // return pressure as a function of energy
double E(double d,double p); // return energy as a function of pressure

// signature for emptying a vector

void vempty(vector<double>&v);

using namespace std;

int main(){

  cout<<"main(): Starting up main loop..."<<endl;

// global data

  int const n(100),ng(n+2);            // no. ncells and ghosts
  vector<double> d(ng),p(ng),V0(ng),V1(ng),m(ng); // pressure, density, volume, mass
  vector<double> e0(2*ng),e1(2*ng);   // fe energy
  vector<double> ec0(ng),ec1(ng); // cell energy
  vector<double> u0(2*ng),u1(2*ng);   // velocity
  vector<double> x0(2*ng),x1(2*ng);   // coordinates
  double time(0.0),dt(DTSTART);       // start time and time step
  int step(0);                        // step number
  double l[3]={1.0,0.0,1.0},r[3]={0.125,0.0,0.1}; // left/right flux states for Riemann solver
  int const nloc(2);ngi(2);            // number of local nodes and integration points
  double N[nloc][ngi],NX[nloc][ngi];dw[ngi]; // fe shape function, derivatives and weight
  vector<double> normal(2*ng);         // face normal

// initialise the problem (Sod's shock tube)

  double dx(1.0/n);x0.at(0)=0.0-dx;x0.at(1)=0.0;
  for(int i=1;i<ng;i++){x0.at(2*i)=x0[2*(i-1)+1];x0.at(2*i+1)=x0[2*i]+dx;}
  for(int i=0;i<ng;i++){p.at(i)=(0.5*(x0[2*i]+x0[2*i+1])<=0.5)?1.0:0.1;}
  for(int i=0;i<ng;i++){d.at(i)=(0.5*(x0[2*i]+x0[2*i+1])<=0.5)?1.0:0.125;}
  for(int i=0;i<ng;i++){e0.at(2*i)=E(d[i],p[i]);e0.at(2*i+1)=E(d[i],p[i]);}
  for(int i=0;i<ng;i++){ec0.at(i)=E(d[i],p[i]);}
  for(int i=0;i<ng;i++){V0.at(i)=x0[2*i+1]-x0[2*i];}
  for(int i=0;i<ng;i++){m.at(i)=d[i]*V0[i];}
  for(int i=0;i<2*ng;i++){u0.at(i)=0.0;u1.at(i)=0.0;}
  for(int i=0;i<2*ng;i++){normal.at(i)=-1.0?i%2:1.0;}

// fe spaces

  N[0][0]=0.5*(1.0+1.0/sqrt(3.0));N[0][1]=0.5*(1.0-1.0/sqrt(3.0));
  N[1][0]=0.5*(1.0-1.0/sqrt(3.0));N[1][1]=0.5*(1.0+1.0/sqrt(3.0));
  NX[0][0]=-0.5;NX[0][1]=-0.5;NX[1][0]=0.5;NX[1][1]=0.5;

// start the Riemann solver

  Riemann R(Riemann::exact,l,r);

// set output precision

  cout<<fixed<<setprecision(17);

// time integration

  while(time<ENDTIME+dt){

    cout<<"  step "<<step<<" time= "<<time<<" dt= "<<dt<<endl;

// evolve the Riemann problem to start of the time step

    vector<double> rx;vempty(rx); // sample coordinates
    for(int i=1;i<=n;i++){rx.push_back(x0[2*i]);rx.push_back(x0[2*i+1]);}
    R.profile(&rx,time);

// move nodes to full step position

    for(long i=0;i<n;i++){x1.at(2*i+2)=x0[2*i+2]+R.velocity(2*i)*dt;x1.at(2*i+3)=x0[2*i+3]+R.velocity(2*i+1)*dt;}

//    l[0]=R.density(2*i);l[1]=R.velocity(2*i);l[2]=R.pressure(2*i);
//    r[0]=R.density(2*i+1);r[1]=R.velocity(2*i+1);r[2]=R.pressure(2*i+1);
//    Riemann R1(Riemann::exact,l,r);
//    for(long i=0;i<n;i++){x1.at(2*i+2)=x0[2*i+2]+R.velocity(2*i)*dt;x1.at(2*i+3)=x0[2*i+3]+R1.ustar*dt;}
//    exit(1);

// update cell volumes at the full step

    for(int i=0;i<ng;i++){V1.at(i)=x1[2*i+1]-x1[2*i];} 

// update cell density at the full step

    for(int i=0;i<ng;i++){d.at(i)=m[i]/V1[i];} 

// update cell energy at the full step

    for(int i=0;i<ng;i++){ec1.at(i)=ec0[i]-(p[i]*(V1[i]-V0[i]))/m[i];}

// update cell pressure at the full step

    for(int i=0;i<ng;i++){p.at(i)=P(d[i],ec1[i]);}

// update node velocity at the full step

// ...

// some output

    for(long i=0;i<2*n;i++){cout<<rx[i]<<" "<<R.density(i)<<" "<<R.pressure(i)<<" "<<R.velocity(i)<<" "<<R.energy(i)<<endl;}
//    for(long i=1;i<=n;i++){cout<<0.5*(x1[2*i]+x1[2*i+1])<<" "<<d[i]<<" "<<ec1[i]<<endl;}

// advance the time step

    time+=dt;
    step++;

// update start-of-step information for the new time step

    for(int i=0;i<2*ng;i++){x0.at(i)=x1[i];}
    for(int i=0;i<2*ng;i++){e0.at(i)=e1[i];}
    for(int i=0;i<2*ng;i++){u0.at(i)=u1[i];}
    for(int i=0;i<ng;i++){V0.at(i)=V1[i];}
    for(int i=0;i<ng;i++){ec0.at(i)=ec1[i];}

  }

  cout<<"Normal termination."<<endl;

  return 0;
}

// return pressure given the energy

double P(double d,double e){return (GAMMA-1.0)*d*e;}

// return energy given the pressure

double E(double d,double p){return p/((GAMMA-1.0)*d);}

// empty a vector

void vempty(vector<double>&v){vector<double> e;v.swap(e);return;}
