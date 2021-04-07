// Main program for RhoLo (Riemann-based Hydro in One-dimension at Low Order - RhoLo)
// RhoLo is an ultra simple 1-D discontinuous finite element (DG) hydrodynamics test code
// This solves the Euler equations in their non-conservative form in the fluid frame (the Lagrangian frame)
// using Riemann boundary coniditions on each element, initial implementation is only first order in time

// Author S. R. Merton

#define DTSTART 0.0005    // insert a macro for the first time step
#define ENDTIME 0.25      // insert a macro for the end time
#define GAMMA 1.4         // ratio of specific heats for ideal gases
#define ECUT 1.0-8        // cut-off on the energy field
#define NSAMPLES 100      // number of sample points for the exact solution
#define VISFREQ 20       // frequency of the graphics dumps
#define VD vector<double> // vector of doubles

#include <iostream>
#include <vector>
#include <iomanip>
#include "riemann.h"
#include <cmath>
#include "matrix.h"
#include "shape.h"
#include "mesh.h"

// sigantures for eos lookups

double P(double d,double e); // eos returns pressure as a function of energy
double E(double d,double p); // invert the eos to get energy if we only have pressure

void vempty(vector<double>&v); // signature for emptying a vector
void silo(Mesh*M,VD*x0,VD*d,VD*p,VD*m,VD*ec0,VD*V0,VD*u0,VD*e0,Shape*S[],int cycle,double time);         // signature for silo output

using namespace std;

int main(){

  cout<<"main(): Starting up main loop..."<<endl;

// initialise a new mesh from file

  Mesh M((char*) "input.mesh");

// global data

  int const n(100),ng(n+2),order(1);                    // no. ncells, ghosts and element order
  Shape S1(order),S2(order),S3(order);                // load FE stencils for energy/momentum equation
  Shape*S[2]={&S3,&S3};                                 // pack the shapes into an array
  vector<double> d(ng),p(ng),V0(ng),V1(ng),m(ng);       // pressure, density, volume & mass
  vector<double> e0(S3.nloc()*ng),e1(S3.nloc()*ng);     // fe DG energy field
  vector<double> ec0(ng),ec1(ng);                       // cell-centred energy field
  vector<double> u0(S3.nloc()*ng),u1(S3.nloc()*ng);     // velocity field
  vector<double> x0(S3.nloc()*ng),x1(S3.nloc()*ng);     // spatial coordinates
  double time(0.0),dt(DTSTART);                         // start time and time step
  int step(0);                                          // step number
  double l[3]={1.0,0.0,1.0},r[3]={0.125,0.0,0.1};       // left/right flux states for Riemann solver
  double S1S[S1.nloc()][S1.nloc()]={};                  // empty surface block for S1 shape function
  double S2S[S2.nloc()][S2.nloc()]={};                  // empty surface block for S2 shape function
  double S3S[S3.nloc()][S3.nloc()]={};                  // empty surface block for S3 shape function

// initialise the problem

  double dx(1.0/n);
  for(int i=0;i<S3.nloc();i++){x0.at(i)=-dx+i*dx/(S3.order());}
  for(int i=1;i<ng;i++){for(int j=0;j<S3.nloc();j++){x0.at(S3.nloc()*i+j)=x0[S3.nloc()*i-1]+j*dx/S3.order();}}
  for(int i=0;i<ng;i++){p.at(i)=(0.5*(x0[S3.nloc()*i]+x0[S3.nloc()*(i+1)-1])<=0.5)?l[2]:r[2];}
  for(int i=0;i<ng;i++){d.at(i)=(0.5*(x0[S3.nloc()*i]+x0[S3.nloc()*(i+1)-1])<=0.5)?l[0]:r[0];}
  for(int i=0;i<ng;i++){for(int j=0;j<S3.nloc();j++){e0.at(S3.nloc()*i+j)=E(d[i],p[i]);}}
  for(int i=0;i<ng;i++){ec0.at(i)=E(d[i],p[i]);}
  for(int i=0;i<ng;i++){V0.at(i)=x0[S3.nloc()*(i+1)-1]-x0[S3.nloc()*i];}
  for(int i=0;i<ng;i++){m.at(i)=d[i]*V0[i];}
  for(int i=0;i<S3.nloc()*ng;i++){u1.at(i)=0.0;}
  for(int i=0;i<ng;i++){for(int j=0;j<S3.nloc();j++){u0.at(i*S3.nloc()+j)=(0.5*(x0[S3.nloc()*i]+x0[S3.nloc()*(i+1)-1])<=0.5)?l[1]:r[1];}}

// surface blocks on S1,S2 and S3 finite element stencils, these couple to the upwind/downwind element on each face

  S1S[0][0]=-1.0;S1S[S1.nloc()-1][S1.nloc()-1]=1.0; // contains the outward pointing unit normal on each face
  S2S[0][0]=-1.0;S2S[S2.nloc()-1][S2.nloc()-1]=1.0; // contains the outward pointing unit normal on each face
  S3S[0][0]=-1.0;S3S[S3.nloc()-1][S3.nloc()-1]=1.0; // contains the outward pointing unit normal on each face

// start the Riemann solver from initial flux states

  Riemann R(Riemann::exact,l,r);
  cout<<" pstar= "<<R.pstar<<endl;
  cout<<" ustar= "<<R.ustar<<endl;

// set output precision

  cout<<fixed<<setprecision(17);

// time integration

  while(time<ENDTIME+dt){

    cout<<"  step "<<step<<" time= "<<time<<" dt= "<<dt<<endl;

// graphics output

    if(step%VISFREQ==0){silo(&M,&x0,&d,&p,&m,&ec0,&V0,&u0,&e0,S,step,time);}

// evolve the Riemann problem to the current time level

    vector<double> rx;vempty(rx); // sample point coordinates
    for(long i=0;i<NSAMPLES;i++){rx.push_back(double(i)/double(NSAMPLES));}
    R.profile(&rx,time);

// move the nodes to their full-step position

    for(long i=1;i<=n;i++){

// fluxes on left and right sides of face 0 (left boundary of cell)

      l[0]=d[i-1];l[1]=u0[S3.nloc()*i-1];l[2]=p[i-1];
      r[0]=d[i];r[1]=u0[S3.nloc()*i];r[2]=p[i];
      Riemann f0(Riemann::exact,l,r);

// fluxes on left and right sides of face 1 (right boundary of cell)

      l[0]=d[i];l[1]=u0[S3.nloc()*(i+1)-1];l[2]=p[i];
      r[0]=d[i+1];r[1]=u0[S3.nloc()*(i+1)];r[2]=p[i+1];
      Riemann f1(Riemann::exact,l,r);

//      double ustar[S3.nloc()]={};ustar[0]=f0.ustar;ustar[1]=0.5*(f0.ustar+f1.ustar);ustar[S3.nloc()-1]=f1.ustar;
//      double ustar[S3.nloc()]={};ustar[0]=f0.ustar;ustar[1]=u0[S3.nloc()*i+1];ustar[S3.nloc()-1]=f1.ustar;
      double ustar[S3.nloc()]={};ustar[0]=f0.ustar;ustar[S3.nloc()-1]=f1.ustar;

      if(i==1){for(int iloc=0;iloc<S3.nloc();iloc++){x1.at(S3.reflect(iloc,0))=x0[S3.reflect(iloc,0)]+ustar[iloc]*dt;}} // move ghost cell on left mesh boundary

      for(int iloc=0;iloc<S3.nloc();iloc++){x1.at(S3.nloc()*i+iloc)=x0[S3.nloc()*i+iloc]+ustar[iloc]*dt;}

      if(i==n){for(int iloc=0;iloc<S3.nloc();iloc++){x1.at(S3.nloc()*(n+1)+S3.reflect(iloc,0))=x0[S3.nloc()*(n+1)+S3.reflect(iloc,0)]+ustar[iloc]*dt;}}// move ghost cell on right mesh boundary

    }

// update cell volumes at the full-step

    for(int i=0;i<ng;i++){V1.at(i)=x1[S3.nloc()*(i+1)-1]-x1[S3.nloc()*i];if(V1[i]<0.0){cout<<"ERROR:  -'ve volume in cell "<<i<<endl;exit(1);}}

// update cell density at the full-step

    for(int i=0;i<ng;i++){d.at(i)=m[i]/V1[i];} 

// update cell energy at the full-step

    for(int i=0;i<ng;i++){ec1.at(i)=max(ECUT,ec0[i]-(p[i]*(V1[i]-V0[i]))/m[i]);}

// construct the full-step DG energy field

    for(int i=1;i<=n;i++){

      double dx(x1[S3.nloc()*(i+1)-1]-x1[S3.nloc()*i]); // cell width for Jacobian

// fluxes on face 0 of i (left boundary of cell)

      l[0]=d[i-1];l[1]=u0[S3.nloc()*i-1];l[2]=p[i-1];
      r[0]=d[i];r[1]=u0[S3.nloc()*i];r[2]=p[i];
      Riemann f0(Riemann::exact,l,r);

// fluxes on face 1 of i (right boundary of cell)

      l[0]=d[i];l[1]=u0[S3.nloc()*(i+1)-1];l[2]=p[i];
      r[0]=d[i+1];r[1]=u0[S3.nloc()*(i+1)];r[2]=p[i+1];
      Riemann f1(Riemann::exact,l,r);

// pressure and velocity on each face

//      double pstar[S3.nloc()]={};pstar[0]=f0.pstar;pstar[1]=0.5*(f0.pstar+f1.pstar);pstar[S3.nloc()-1]=f1.pstar;
//      double ustar[S3.nloc()]={};ustar[0]=f0.ustar;ustar[1]=0.5*(f0.ustar+f1.ustar);ustar[S3.nloc()-1]=f1.ustar;
//      double u0vol[S3.nloc()]={};u0vol[0]=f0.ustar;u0vol[1]=0.5*(f0.ustar+f1.ustar);u0vol[S3.nloc()-1]=f1.ustar;

      double pstar[S3.nloc()]={};pstar[0]=f0.pstar;pstar[S3.nloc()-1]=f1.pstar;
      double ustar[S3.nloc()]={};ustar[0]=f0.ustar;;ustar[S3.nloc()-1]=f1.ustar;
      double u0vol[S3.nloc()]={};u0vol[0]=f0.ustar;u0vol[S3.nloc()-1]=f1.ustar;

// matrix problem for one element

      Matrix A(S3.nloc());double b[S3.nloc()],soln[S3.nloc()];

// assemble DG energy field for one element

      for(int iloc=0;iloc<S3.nloc();iloc++){
        b[iloc]=0.0;
        for(int jloc=0;jloc<S3.nloc();jloc++){
          double nn(0.0),nxn(0.0),nnx(0.0);
          for(int gi=0;gi<S3.ngi();gi++){
            nn+=S3.value(iloc,gi)*S3.value(jloc,gi)*S3.wgt(gi)*dx/2.0; // mass matrix
            nnx-=S3.value(iloc,gi)*S3.dvalue(jloc,gi)*S3.wgt(gi);      // divergence term (for continuous finite elements)
            nxn+=S3.dvalue(iloc,gi)*S3.value(jloc,gi)*S3.wgt(gi);      // divergence term (if by parts, use this for DG)
          }
          A.write(iloc,jloc,nn);                // commit to address space in the matrix class
//          b[iloc]+=nnx*ustar[jloc]*p[i]/d[i]; // source - for continuous finite elements
          b[iloc]+=(nxn*u0vol[jloc]-S3S[iloc][jloc]*ustar[jloc])*p[i]/d[i]; // source - discontinuous, for DG
        }
      }

      A.solve(soln,b);

// advance the solution

      for(int iloc=0;iloc<S3.nloc();iloc++){
        e1[S3.nloc()*i+iloc]=max(ECUT,e0[S3.nloc()*i+iloc]+soln[iloc]*dt);
      }

    }

// update cell pressure at the full-step using PdV / DG energy field

    for(int i=0;i<ng;i++){p.at(i)=P(d[i],ec1[i]);} // use PdV

// update nodal DG velocities at the full step

    for(int i=1;i<=n;i++){

      double dx(x1[S3.nloc()*(i+1)-1]-x1[S3.nloc()*i]); // cell width for Jacobian

// fluxes on face 0 of i (left boundary of cell)

      l[0]=d[i-1];l[1]=u0[S3.nloc()*i-1];l[2]=p[i-1];
      r[0]=d[i];r[1]=u0[S3.nloc()*i];r[2]=p[i];

      Riemann f0(Riemann::exact,l,r);

// fluxes on face 1 of i (right boundary of cell)

      l[0]=d[i];l[1]=u0[S3.nloc()*(i+1)-1];l[2]=p[i];
      r[0]=d[i+1];r[1]=u0[S3.nloc()*(i+1)];r[2]=p[i+1];

      Riemann f1(Riemann::exact,l,r);

// pressure and velocity on each face

//      double pstar[S3.nloc()]={};pstar[0]=f0.pstar;pstar[1]=0.5*(f0.pstar+f1.pstar);pstar[S3.nloc()-1]=f1.pstar;
//      double pvol[S3.nloc()]={};pvol[0]=f0.pstar;pvol[1]=0.5*(f0.pstar+f1.pstar);pvol[S3.nloc()-1]=f1.pstar;
//      double ustar[S3.nloc()]={};ustar[0]=f0.ustar;ustar[1]=u0[S3.nloc()*i+1];ustar[S3.nloc()-1]=f1.ustar;

      double pstar[S3.nloc()]={};pstar[0]=f0.pstar;pstar[S3.nloc()-1]=f1.pstar;
      double pvol[S3.nloc()]={};pvol[0]=f0.pstar;pvol[S3.nloc()-1]=f1.pstar;
      double ustar[S3.nloc()]={};ustar[0]=f0.ustar;ustar[S3.nloc()-1]=f1.ustar;

// matrix problem for one element

      Matrix A(S3.nloc());double b[S3.nloc()],soln[S3.nloc()];

// assemble acceleration field for one element

      for(int iloc=0;iloc<S3.nloc();iloc++){
        b[iloc]=0.0;
        for(int jloc=0;jloc<S3.nloc();jloc++){
          double nn(0.0),nxn(0.0),nnx(0.0);
          for(int gi=0;gi<S3.ngi();gi++){
            nn+=S3.value(iloc,gi)*S3.value(jloc,gi)*S3.wgt(gi)*dx/2.0; // mass matrix
//            nnx-=S3.value(iloc,gi)*S3.dvalue(jloc,gi)*S3.wgt(gi);      // grad term (for continuous finite elements)
            nxn+=S3.dvalue(iloc,gi)*S3.value(jloc,gi)*S3.wgt(gi);      // grad term (if by parts, use this for DG)
          }
          A.write(iloc,jloc,nn);                // commit to address space in the matrix class
//          b[iloc]+=nnx*pstar[jloc]/d[i]; // source - for continuous finite elements
          b[iloc]+=(nxn*pvol[jloc]-S3S[iloc][jloc]*pstar[jloc])/d[i]; // source - discontinuous, for DG
        }
      }

      A.solve(soln,b);

// advance the solution

      for(int iloc=0;iloc<S3.nloc();iloc++){
        u1[S3.nloc()*i+iloc]=u0[S3.nloc()*i+iloc]+soln[iloc]*dt;
      }

// impose a constraint on the acceleration field at domain boundaries to stop the mesh taking off

      for(int i=0;i<S3.nloc();i++){u1.at(i)=u0[i];u1.at(S3.nloc()*(n+1)+i)=u0[S3.nloc()*(n+1)+i];}

    }

// some output - toggle this to output either the exact solutions from the Riemann solver or the finite element solution generated by the code

//    for(int i=0;i<NSAMPLES;i++){cout<<rx[i]<<" "<<R.density(i)<<" "<<R.pressure(i)<<" "<<R.velocity(i)<<" "<<R.energy(i)<<endl;} // exact solution from Riemann solver
    for(long i=1;i<=n;i++){for(int iloc=0;iloc<S3.nloc();iloc++){cout<<x1[S3.nloc()*i+iloc]<<" "<<d[i]<<" "<<p[i]<<" "<<u1[S3.nloc()*i+iloc]<<" "<<e1[S3.nloc()*i+iloc]<<endl;}} // high-order DG

// advance the time step

    time+=dt;
    step++;

// advance the solution for the new time step

    for(int i=0;i<S3.nloc()*ng;i++){u0.at(i)=u1[i];}
    for(int i=0;i<S3.nloc()*ng;i++){e0.at(i)=e1[i];}
    for(int i=0;i<S3.nloc()*ng;i++){x0.at(i)=x1[i];}
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
