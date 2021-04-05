// Main program for RhoLo (Riemann-based Hydro in One-dimension at Low Order - RhoLo)
// RhoLo is an ultra simple 1-D discontinuous finite element (DG) hydrodynamics test code
// This solves the Euler equations in their non-conservative form in the fluid frame (the Lagrangian frame)
// using Riemann boundary coniditions on each element, initial implementation is only first order in time

// Author S. R. Merton

#define DTSTART 0.00005  // insert a macro for the first time step
#define ENDTIME 0.25    // insert a macro for the end time
#define GAMMA 1.4       // ratio of specific heats for ideal gases
#define ECUT 1.0-8      // cut-off on the energy field

#include <iostream>
#include <vector>
#include <iomanip>
#include "riemann.h"
#include <cmath>
#include "matrix.h"

// sigantures for eos lookups

double P(double d,double e); // eos returns pressure as a function of energy
double E(double d,double p); // invert the eos to get energy if we only have pressure

// signature for emptying a vector

void vempty(vector<double>&v);

using namespace std;

int main(){

  cout<<"main(): Starting up main loop..."<<endl;

// global data

  int const n(500),ng(n+2);                               // no. ncells and ghosts
  vector<double> d(ng),p(ng),V0(ng),V1(ng),m(ng);       // pressure, density, volume & mass
  vector<double> e0(2*ng),e1(2*ng);                     // fe DG energy field
  vector<double> ec0(ng),ec1(ng);                       // cell energy field
  vector<double> u0(2*ng),u1(2*ng);                     // velocity field
  vector<double> x0(2*ng),x1(2*ng);                     // spatial coordinates
  double time(0.0),dt(DTSTART);                         // start time and time step
  int step(0);                                          // step number
  double l[3]={1.0,0.0,1.0},r[3]={0.125,0.0,0.1};       // left/right flux states for Riemann solver
  int const nloc(2),ngi(2);                             // no. of local nodes and integration points
  double N[nloc][ngi],NX[nloc][ngi];                    // fe shapes and their derivatives
  double NN[nloc][nloc],NXN[nloc][nloc],NNX[nloc][nloc],SN[nloc][nloc]; // mass matrix, divergence term and surface block
  vector<double> normal(2*ng);                          // normal to each face

// initialise the problem (Sod's shock tube) - hack this to run something else

  double dx(1.0/n);x0.at(0)=0.0-dx;x0.at(1)=0.0;
  for(int i=1;i<ng;i++){x0.at(2*i)=x0[2*(i-1)+1];x0.at(2*i+1)=x0[2*i]+dx;}
  for(int i=0;i<ng;i++){p.at(i)=(0.5*(x0[2*i]+x0[2*i+1])<=0.5)?1.0:0.1;}
  for(int i=0;i<ng;i++){d.at(i)=(0.5*(x0[2*i]+x0[2*i+1])<=0.5)?1.0:0.125;}
  for(int i=0;i<ng;i++){e0.at(2*i)=E(d[i],p[i]);e0.at(2*i+1)=E(d[i],p[i]);}
  for(int i=0;i<ng;i++){ec0.at(i)=E(d[i],p[i]);}
  for(int i=0;i<ng;i++){V0.at(i)=x0[2*i+1]-x0[2*i];}
  for(int i=0;i<ng;i++){m.at(i)=d[i]*V0[i];}
  for(int i=0;i<2*ng;i++){u0.at(i)=0.0;u1.at(i)=0.0;}
  for(int i=0;i<2*ng;i++){normal.at(i)=(i%2)?1.0:-1.0;}

// fe spaces for p1 DG element type - hack this for high order

  N[0][0]=0.5*(1.0+1.0/sqrt(3.0));N[0][1]=0.5*(1.0-1.0/sqrt(3.0));
  N[1][0]=0.5*(1.0-1.0/sqrt(3.0));N[1][1]=0.5*(1.0+1.0/sqrt(3.0));
  NX[0][0]=-0.5;NX[0][1]=-0.5;NX[1][0]=0.5;NX[1][1]=0.5;
  SN[0][0]=1.0;SN[0][1]=0.0;SN[1][0]=0.0;SN[1][1]=1.0;

// start the Riemann solver from initial flux states

  Riemann R(Riemann::exact,l,r);
  cout<<" pstar= "<<R.pstar<<endl;
  cout<<" ustar= "<<R.ustar<<endl;

// set output precision

  cout<<fixed<<setprecision(17);

// time integration

  while(time<ENDTIME+dt){

    cout<<"  step "<<step<<" time= "<<time<<" dt= "<<dt<<endl;

// evolve the Riemann problem to the current time level

    vector<double> rx;vempty(rx); // sample point coordinates
    for(long i=0;i<ng;i++){rx.push_back(x0[2*i]);rx.push_back(0.5*(x0[2*i]+x0[2*i+1]));rx.push_back(x0[2*i+1]);}
    R.profile(&rx,time);

// move the nodes to their full-step position

    for(long i=1;i<=n;i++){

// fluxes on left and right sides of face 0 (left boundary of cell)

      l[0]=d[i-1];l[1]=u0[2*(i-1)+1];l[2]=p[i-1];
      r[0]=d[i];r[1]=u0[2*i];r[2]=p[i];
      Riemann f0(Riemann::exact,l,r);

// fluxes on left and right sides of face 1 (right boundary of cell)

      l[0]=d[i];l[1]=u0[2*i+1];l[2]=p[i];
      r[0]=d[i+1];r[1]=u0[2*(i+1)];r[2]=p[i+1];
      Riemann f1(Riemann::exact,l,r);

      if(i==1){x1.at(1)=x0[1]+f0.ustar*dt;x1.at(0)=x0[0];} // move ghost cell on left mesh boundary

      x1.at(2*i)=x0[2*i]+f0.ustar*dt;x1.at(2*i+1)=x0[2*i+1]+f1.ustar*dt;

      if(i==n){x1.at(2*(n+1))=x0[2*(n+1)]+f1.ustar*dt;x1.at(2*(n+1)+1)=x0[2*(n+1)+1];} // move ghost cell on right mesh boundary

    }

// update cell volumes at the full-step

    for(int i=0;i<ng;i++){V1.at(i)=x1[2*i+1]-x1[2*i];if(V1[i]<0.0){cout<<"ERROR:  -'ve volume in cell "<<i<<endl;exit(1);}} 

// update cell density at the full-step

    for(int i=0;i<ng;i++){d.at(i)=m[i]/V1[i];} 

// update cell energy at the full-step

    for(int i=0;i<ng;i++){ec1.at(i)=max(ECUT,ec0[i]-(p[i]*(V1[i]-V0[i]))/m[i]);}

// construct the full-step DG energy field

    for(int i=1;i<=n;i++){

      double dx(x1[2*i+1]-x1[2*i]); // cell width for Jacobian

// fluxes on face 0 of i (left boundary of cell)

      l[0]=d[i-1];l[1]=u0[2*(i-1)+1];l[2]=p[i-1];
      r[0]=d[i];r[1]=u0[2*i];r[2]=p[i];
      Riemann f0(Riemann::exact,l,r);

// fluxes on face 1 of i (right boundary of cell)

      l[0]=d[i];l[1]=u0[2*i+1];l[2]=p[i];
      r[0]=d[i+1];r[1]=u0[2*(i+1)];r[2]=p[i+1];
      Riemann f1(Riemann::exact,l,r);

// pressure and velocity on each face

      double pstar[2]={f0.pstar,f1.pstar};
      double ustar[2]={f0.ustar,f1.ustar};

// matrix problem for one element

      Matrix A(nloc);double b[nloc],soln[nloc];

// assemble DG energy field for one element

      for(int iloc=0;iloc<nloc;iloc++){
        b[iloc]=0.0;
        for(int jloc=0;jloc<nloc;jloc++){
          NN[iloc][jloc]=0.0;NXN[iloc][jloc]=0.0;NNX[iloc][jloc]=0.0;
          for(int gi=0;gi<ngi;gi++){
            NN[iloc][jloc]+=N[iloc][gi]*N[jloc][gi]*dx/2.0; // mass matrix
            NNX[iloc][jloc]-=N[iloc][gi]*NX[jloc][gi];      // divergence term (for a continuous finite element energy field)
            NXN[iloc][jloc]+=NX[iloc][gi]*N[jloc][gi];      // divergence term (if by parts, use this for DG)
          }
          A.write(iloc,jloc,NN[iloc][jloc]);                // commit to address space in the matrix class
//          b[iloc]+=NNX[iloc][jloc]*ustar[jloc]*p[i]/d[i]; // source - for continuous finite element
          b[iloc]+=(NXN[iloc][jloc]*ustar[jloc]-normal[2*i+jloc]*SN[iloc][jloc]*ustar[jloc])*p[i]/d[i]; // source - discontinuous, for DG
        }
      }

      A.solve(soln,b);e1[2*i]=max(ECUT,e0[2*i]+soln[0]*dt); e1[2*i+1]=max(ECUT,e0[2*i+1]+soln[1]*dt);

    }

// update cell pressure at the full-step using PdV / DG energy field

    for(int i=0;i<ng;i++){p.at(i)=P(d[i],ec1[i]);} // use PdV

// update nodal DG velocities at the full step

    for(int i=1;i<=n;i++){

      double dx(x1[2*i+1]-x1[2*i]); // cell width for Jacobian

// fluxes on face 0 of i (left boundary of cell)

      l[0]=d[i-1];l[1]=u0[2*(i-1)+1];l[2]=p[i-1];
      r[0]=d[i];r[1]=u0[2*i];r[2]=p[i];

      Riemann f0(Riemann::exact,l,r);

// fluxes on face 1 of i (right boundary of cell)

      l[0]=d[i];l[1]=u0[2*i+1];l[2]=p[i];
      r[0]=d[i+1];r[1]=u0[2*(i+1)];r[2]=p[i+1];

      Riemann f1(Riemann::exact,l,r);

// pressure and velocity on each face

      double pstar[2]={f0.pstar,f1.pstar};
      double ustar[2]={f0.ustar,f1.ustar};

// matrix problem for one element

      Matrix A(nloc);double b[nloc],soln[nloc];

// assemble acceleration field for one element

      for(int iloc=0;iloc<nloc;iloc++){
        b[iloc]=0.0;
        for(int jloc=0;jloc<nloc;jloc++){
          NN[iloc][jloc]=0.0;NXN[iloc][jloc]=0.0;NNX[iloc][jloc]=0.0;
          for(int gi=0;gi<ngi;gi++){
            NN[iloc][jloc]+=N[iloc][gi]*N[jloc][gi]*dx/2.0; // mass matrix
            NNX[iloc][jloc]-=N[iloc][gi]*NX[jloc][gi];      // grad term (for a continuous finite element energy field)
            NXN[iloc][jloc]+=NX[iloc][gi]*N[jloc][gi];      // grad term (if by parts, use this for DG)
          }
          A.write(iloc,jloc,NN[iloc][jloc]);                // commit to address space in the matrix class
//          b[iloc]+=NNX[iloc][jloc]*pstar[jloc]/d[i]; // source - for continuous finite element
          b[iloc]+=(NXN[iloc][jloc]*pstar[jloc]-normal[2*i+jloc]*SN[iloc][jloc]*pstar[jloc])/d[i]; // source - discontinuous, for DG
        }
      }

      A.solve(soln,b);u1[2*i]=u0[2*i]+soln[0]*dt; u1[2*i+1]=u0[2*i+1]+soln[1]*dt;

    }

// impose constraint on acceleration field at domain boundaries

    u1[0]=u0[0];
    u1[2*ng]=u0[2*ng];

// some output - toggle this to output either the exact solutions from the Riemann solver or the finite element solution generated by the code

    for(int i=0;i<3*ng;i++){cout<<rx[i]<<" "<<R.density(i)<<" "<<R.pressure(i)<<" "<<R.velocity(i)<<" "<<R.energy(i)<<endl;} // exact solution from Riemann solver
//    for(long i=1;i<=n;i++){cout<<x1[2*i]  <<" "<<d[i]<<" "<<p[i]<<" "<<u1[2*i]<<" "<<e1[2*i]<<endl;cout<<x1[2*i+1]<<" "<<d[i]<<" "<<p[i]<<" "<<u1[2*i+1]<<" "<<e1[2*i+1]<<endl;} // DG solution

// advance the time step

    time+=dt;
    step++;

// advance the solution for the new time step

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

// invert the eos to return energy given the pressure

double E(double d,double p){return p/((GAMMA-1.0)*d);}

// empty a vector

void vempty(vector<double>&v){vector<double> e;v.swap(e);return;}
