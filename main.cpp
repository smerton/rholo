// main program for rholo
// Riemann Hydro One-dimensional Low Order hydrodynamics test code

// Author S. R. Merton

#define DTSTART 0.01  // insert a macro for the first time step
#define ENDTIME 0.25 // insert a macro for the end time
#define GAMMA 1.4 // ration of specific heats for ideal gases

#include <iostream>
#include <vector>
#include <iomanip>

// eos lookups

double P(double d,double e); // return pressure as a function of energy
double E(double d,double p); // return energy as a function of pressure

using namespace std;

int main(){

  cout<<"main(): Starting up main loop..."<<endl;

// global data

  int const n(10),ng(n+2);            // no. ncells and ghosts
  vector<double> d(ng),p(ng),V0(ng),V1(ng),m(ng); // pressure, density, volume, mass
  vector<double> e0(2*ng),e1(2*ng);   // fe energy
  vector<double> ec0(2*ng),ec1(2*ng); // cell energy
  vector<double> u0(2*ng),u1(2*ng);   // velocity
  vector<double> x0(2*ng),x1(2*ng);   // coordinates
  double time(0.0),dt(DTSTART);       // start time and time step
  int step(0);                        // step number

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

  cout<<fixed<<setprecision(17);      // set output precision

// time integration

  while(time<ENDTIME+dt){

    cout<<"  step "<<step<<" time= "<<time<<" dt= "<<dt<<endl;

// move nodes to full step position

    for(int i=0;i<2*ng;i++){x1[i]=x0[i]+u0[i]*dt;}

// update cell volumes at the full step

    for(int i=0;i<ng;i++){V1.at(i)=x1[2*i+1]-x1[2*i];} 

// update cell density at the full step

    for(int i=0;i<ng;i++){d.at(i)=V1[i]*m[i];} 

// update cell energy at the full step

    for(int i=0;i<ng;i++){ec1[i]=ec0[i]-p[i]*(V1[i]-V0[i]);}

// update cell pressure at the full step

    for(int i=0;i<ng;i++){p[i]=P(d[i],ec1[i]);}

// update node velocity at the full step




// advance the time step

    time+=dt;
    step++;

  }

  cout<<"Normal termination."<<endl;

  return 0;
}

// return pressure given the energy

double P(double d,double e){return (GAMMA-1.0)*d*e;}

// return energy given the pressure

double E(double d,double p){return p/((GAMMA-1.0)*d);}
