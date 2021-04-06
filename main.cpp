// Main program for RhoLo (Riemann-based Hydro in One-dimension at Low Order - RhoLo)
// RhoLo is an ultra simple 1-D discontinuous finite element (DG) hydrodynamics test code
// This solves the Euler equations in their non-conservative form in the fluid frame (the Lagrangian frame)
// using Riemann boundary coniditions on each element, initial implementation is only first order in time

// Author S. R. Merton

#define DTSTART 0.0005   // insert a macro for the first time step
#define ENDTIME 0.25     // insert a macro for the end time
#define GAMMA 1.4        // ratio of specific heats for ideal gases
#define ECUT 1.0-8       // cut-off on the energy field
#define NSAMPLES 100    // number of sample points for the exact solution

#include <iostream>
#include <vector>
#include <iomanip>
#include "riemann.h"
#include <cmath>
#include "matrix.h"
#include "shape.h"

// sigantures for eos lookups

double P(double d,double e); // eos returns pressure as a function of energy
double E(double d,double p); // invert the eos to get energy if we only have pressure

// signature for emptying a vector

void vempty(vector<double>&v);

using namespace std;

int main(){

  cout<<"main(): Starting up main loop..."<<endl;

// global data

  int const n(100),ng(n+2),order(1);                     // no. ncells, ghosts and element order
  Shape S1(order),S2(order);                            // load FE stencils for energy/momentum equation
  vector<double> d(ng),p(ng),V0(ng),V1(ng),m(ng);       // pressure, density, volume & mass
  vector<double> e0(S1.nloc()*ng),e1(S1.nloc()*ng);     // fe DG energy field
  vector<double> ec0(ng),ec1(ng);                       // cell energy field
  vector<double> u0(S2.nloc()*ng),u1(S2.nloc()*ng);     // velocity field
  vector<double> x0(S2.nloc()*ng),x1(S2.nloc()*ng);     // spatial coordinates
  double time(0.0),dt(DTSTART);                         // start time and time step
  int step(0);                                          // step number
  double l[3]={1.0,0.0,1.0},r[3]={0.125,0.0,0.1};       // left/right flux states for Riemann solver
  double S1S[S1.nloc()][S1.nloc()]={};                  // empty surface block for S1 shape function
  double S2S[S2.nloc()][S2.nloc()]={};                  // empty surface block for S2 shape function
  vector<double> normal(S2.nloc()*ng);                  // normal to each face

// initialise the problem (Sod's shock tube) - hack this to run something else

  double dx(1.0/n);for(int i=0;i<S2.nloc();i++){x0.at(i)=-dx+i*dx/order;}
  for(int i=1;i<ng;i++){for(int j=0;j<S2.nloc();j++){x0.at(S2.nloc()*i+j)=x0[S2.nloc()*i+j-1]+j*dx/order;}}
  for(int i=0;i<ng;i++){p.at(i)=(0.5*(x0[S2.nloc()*i]+x0[S2.nloc()*(i+1)-1])<=0.5)?1.0:0.1;}
  for(int i=0;i<ng;i++){d.at(i)=(0.5*(x0[S2.nloc()*i]+x0[S2.nloc()*(i+1)-1])<=0.5)?1.0:0.125;}
  for(int i=0;i<ng;i++){for(int j=0;j<S1.nloc();j++){e0.at(S1.nloc()*i+j)=E(d[i],p[i]);}}
  for(int i=0;i<ng;i++){ec0.at(i)=E(d[i],p[i]);}
  for(int i=0;i<ng;i++){V0.at(i)=x0[S2.nloc()*(i+1)-1]-x0[S2.nloc()*i];}
  for(int i=0;i<ng;i++){m.at(i)=d[i]*V0[i];}
  for(int i=0;i<S2.nloc()*ng;i++){u0.at(i)=0.0;u1.at(i)=0.0;}
  for(int i=0;i<S2.nfaces()*ng;i++){normal.at(i)=(i%2)?1.0:-1.0;}

// surface blocks on S1 and S2 finite element stencils, these couple to the upwind/downwind element on each face

  S1S[0][0]=1.0;S1S[S1.nloc()-1][S1.nloc()-1]=1.0;
  S2S[0][0]=1.0;S2S[S2.nloc()-1][S2.nloc()-1]=1.0;

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
    for(long i=0;i<NSAMPLES;i++){rx.push_back(double(i)/double(NSAMPLES));}
    R.profile(&rx,time);

// move the nodes to their full-step position

    for(long i=1;i<=n;i++){

// fluxes on left and right sides of face 0 (left boundary of cell)

      l[0]=d[i-1];l[1]=u0[S2.nloc()*i-1];l[2]=p[i-1];
      r[0]=d[i];r[1]=u0[S2.nloc()*i];r[2]=p[i];
      Riemann f0(Riemann::exact,l,r);

// fluxes on left and right sides of face 1 (right boundary of cell)

      l[0]=d[i];l[1]=u0[S2.nloc()*(i+1)-1];l[2]=p[i];
      r[0]=d[i+1];r[1]=u0[S2.nloc()*(i+1)];r[2]=p[i+1];
      Riemann f1(Riemann::exact,l,r);

      if(i==1){x1.at(1)=x0[1]+f0.ustar*dt;x1.at(0)=x0[0];} // move ghost cell on left mesh boundary

      x1.at(S2.nloc()*i)=x0[S2.nloc()*i]+f0.ustar*dt;x1.at(S2.nloc()*(i+1)-1)=x0[S2.nloc()*(i+1)-1]+f1.ustar*dt;

      if(i==n){x1.at(S2.nloc()*(n+1))=x0[S2.nloc()*(n+1)]+f1.ustar*dt;x1.at(S2.nloc()*(n+2)-1)=x0[S2.nloc()*(n+2)-1];} // move ghost cell on right mesh boundary

    }

// update cell volumes at the full-step

    for(int i=0;i<ng;i++){V1.at(i)=x1[S2.nloc()*(i+1)-1]-x1[S2.nloc()*i];if(V1[i]<0.0){cout<<"ERROR:  -'ve volume in cell "<<i<<endl;exit(1);}} 

// update cell density at the full-step

    for(int i=0;i<ng;i++){d.at(i)=m[i]/V1[i];} 

// update cell energy at the full-step

    for(int i=0;i<ng;i++){ec1.at(i)=max(ECUT,ec0[i]-(p[i]*(V1[i]-V0[i]))/m[i]);}

// construct the full-step DG energy field

    for(int i=1;i<=n;i++){

      double dx(x1[S2.nloc()*(i+1)-1]-x1[S2.nloc()*i]); // cell width for Jacobian

// fluxes on face 0 of i (left boundary of cell)

      l[0]=d[i-1];l[1]=u0[S1.nloc()*i-1];l[2]=p[i-1];
      r[0]=d[i];r[1]=u0[S1.nloc()*i];r[2]=p[i];
      Riemann f0(Riemann::exact,l,r);

// fluxes on face 1 of i (right boundary of cell)

      l[0]=d[i];l[1]=u0[S1.nloc()*(i+1)-1];l[2]=p[i];
      r[0]=d[i+1];r[1]=u0[S1.nloc()*(i+1)];r[2]=p[i+1];
      Riemann f1(Riemann::exact,l,r);

// pressure and velocity on each face

//      double pstar[S1.nloc()]={f0.pstar,{},f1.pstar};
      double pstar[S1.nloc()]={};pstar[0]=f0.pstar;pstar[S1.nloc()-1]=f1.pstar;
//      double ustar[S1.nloc()]={f0.ustar,{},f1.ustar};
      double ustar[S1.nloc()]={};ustar[0]=f0.ustar;ustar[S1.nloc()-1]=f1.ustar;

// matrix problem for one element

      Matrix A(S1.nloc());double b[S1.nloc()],soln[S1.nloc()];

// assemble DG energy field for one element

      for(int iloc=0;iloc<S1.nloc();iloc++){
        b[iloc]=0.0;
        for(int jloc=0;jloc<S1.nloc();jloc++){
          double nn(0.0),nxn(0.0),nnx(0.0);
          for(int gi=0;gi<S1.ngi();gi++){
            nn+=S1.value(iloc,gi)*S1.value(jloc,gi)*S1.wgt(gi)*dx/2.0; // mass matrix
            nnx-=S1.value(iloc,gi)*S1.dvalue(jloc,gi)*S1.wgt(gi);      // divergence term (for continuous finite elements)
            nxn+=S1.dvalue(iloc,gi)*S1.value(jloc,gi)*S1.wgt(gi);      // divergence term (if by parts, use this for DG)
          }
          A.write(iloc,jloc,nn);                // commit to address space in the matrix class
//          b[iloc]+=nnx*ustar[jloc]*p[i]/d[i]; // source - for continuous finite elements
          b[iloc]+=(nxn*ustar[jloc]-normal[S1.nfaces()*i+jloc]*S1S[iloc][jloc]*ustar[jloc])*p[i]/d[i]; // source - discontinuous, for DG
        }
      }

      A.solve(soln,b);

// advance the solution

      for(int iloc=0;iloc<S1.nloc();iloc++){
        e1[S1.nloc()*i+iloc]=max(ECUT,e0[S1.nloc()*i+iloc]+soln[iloc]*dt);
      }

    }

// update cell pressure at the full-step using PdV / DG energy field

    for(int i=0;i<ng;i++){p.at(i)=P(d[i],ec1[i]);} // use PdV

// update nodal DG velocities at the full step

    for(int i=1;i<=n;i++){

      double dx(x1[S2.nloc()*(i+1)-1]-x1[S2.nloc()*i]); // cell width for Jacobian

// fluxes on face 0 of i (left boundary of cell)

      l[0]=d[i-1];l[1]=u0[S2.nloc()*i-1];l[2]=p[i-1];
      r[0]=d[i];r[1]=u0[S2.nloc()*i];r[2]=p[i];

      Riemann f0(Riemann::exact,l,r);

// fluxes on face 1 of i (right boundary of cell)

      l[0]=d[i];l[1]=u0[S2.nloc()*(i+1)-1];l[2]=p[i];
      r[0]=d[i+1];r[1]=u0[S2.nloc()*(i+1)];r[2]=p[i+1];

      Riemann f1(Riemann::exact,l,r);

// pressure and velocity on each face

//      double pstar[S2.nloc()]={f0.pstar,{},f1.pstar};
      double pstar[S2.nloc()]={};pstar[0]=f0.pstar;pstar[S2.nloc()-1]=f1.pstar;
//      double ustar[S2.nloc()]={f0.ustar,{},f1.ustar};
      double ustar[S2.nloc()]={};ustar[0]=f0.ustar;ustar[S2.nloc()-1]=f1.ustar;

// matrix problem for one element

      Matrix A(S2.nloc());double b[S2.nloc()],soln[S2.nloc()];

// assemble acceleration field for one element

      for(int iloc=0;iloc<S2.nloc();iloc++){
        b[iloc]=0.0;
        for(int jloc=0;jloc<S2.nloc();jloc++){
          double nn(0.0),nxn(0.0),nnx(0.0);
          for(int gi=0;gi<S2.ngi();gi++){
            nn+=S2.value(iloc,gi)*S2.value(jloc,gi)*S2.wgt(gi)*dx/2.0; // mass matrix
            nnx-=S2.value(iloc,gi)*S2.dvalue(jloc,gi)*S2.wgt(gi);      // grad term (for continuous finite elements)
            nxn+=S2.dvalue(iloc,gi)*S2.value(jloc,gi)*S2.wgt(gi);      // grad term (if by parts, use this for DG)
          }
          A.write(iloc,jloc,nn);                // commit to address space in the matrix class
//          b[iloc]+=nnx*pstar[jloc]/d[i]; // source - for continuous finite elements
          b[iloc]+=(nxn*pstar[jloc]-normal[S2.nfaces()*i+jloc]*S2S[iloc][jloc]*pstar[jloc])/d[i]; // source - discontinuous, for DG
        }
      }

      A.solve(soln,b);

// advance the solution

      for(int iloc=0;iloc<S2.nloc();iloc++){
        u1[S2.nloc()*i+iloc]=u0[S2.nloc()*i+iloc]+soln[iloc]*dt;
      }

    }

// impose a constraint on the acceleration field at domain boundaries to stop the mesh taking off

    for(int i=0;i<S2.nloc();i++){u1.at(i)=u0[i];u1.at(S2.nloc()*(n+1)+i)=u0[S2.nloc()*(n+1)+i];}

// some output - toggle this to output either the exact solutions from the Riemann solver or the finite element solution generated by the code

//    for(int i=0;i<NSAMPLES;i++){cout<<rx[i]<<" "<<R.density(i)<<" "<<R.pressure(i)<<" "<<R.velocity(i)<<" "<<R.energy(i)<<endl;} // exact solution from Riemann solver
    for(long i=1;i<=n;i++){cout<<x1[S2.nloc()*i]  <<" "<<d[i]<<" "<<p[i]<<" "<<u1[S2.nloc()*i]<<" "<<e1[S1.nloc()*i]<<endl;cout<<x1[S2.nloc()*i+S2.nloc()-1]<<" "<<d[i]<<" "<<p[i]<<" "<<u1[S2.nloc()*i+S2.nloc()-1]<<" "<<e1[S1.nloc()*i+S1.nloc()-1]<<endl;} // DG solution

// advance the time step

    time+=dt;
    step++;

// advance the solution for the new time step

    for(int i=0;i<S2.nloc()*ng;i++){x0.at(i)=x1[i];}
    for(int i=0;i<S1.nloc()*ng;i++){e0.at(i)=e1[i];}
    for(int i=0;i<S2.nloc()*ng;i++){u0.at(i)=u1[i];}
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
