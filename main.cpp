// Main program for rholo - an ultra simple 1-D discontinuous finite element (DG) hydrodynamics test code
// Solves the Euler equations in their non-conservative form in the fluid frame (the Lagrangian frame)
// using Riemann boundary coniditions on each element, initial implementation is only first order in time
// Riemann-based Hydro in One-dimensional at Low Order (rholo)

// Author S. R. Merton

#define DTSTART 0.0005  // insert a macro for the first time step
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

  int const n(100),ng(n+2);                               // no. ncells and ghosts
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
  double b[nloc],soln[nloc];                            // matrix equation for one element
  vector<double> normal(2*ng);                          // normal to each face
  Matrix A(nloc);                                       // matrix for one element

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
    for(long i=0;i<ng;i++){rx.push_back(0.5*(x0[2*i]+x0[2*i+1]));}
    R.profile(&rx,time);

// move the nodes to their full-step position

//    for(long i=0;i<n;i++){x1.at(2*i+2)=x0[2*i+2]+R.velocity(2*i)*dt;x1.at(2*i+3)=x0[2*i+3]+R.velocity(2*i+1)*dt;}
    for(long i=1;i<=n;i++){

// fluxes on left and right sides of face 0 (left boundary of cell)

      l[0]=d[i-1];l[1]=R.velocity(i-1);l[2]=p[i-1];
      r[0]=d[i];r[1]=R.velocity(i);r[2]=p[i];
      Riemann f0(Riemann::exact,l,r);

// fluxes on left and right sides of face 1 (right boundary of cell)

      l[0]=d[i];l[1]=R.velocity(i);l[2]=p[i];
      r[0]=d[i+1];r[1]=R.velocity(i+1);r[2]=p[i+1];
      Riemann f1(Riemann::exact,l,r);

//      cout<<"el "<<i<<" u left= "<<f0.ustar<<" u right "<<f1.ustar;
//      cout<<" l= "<<l[0]<<" "<<l[1]<<" "<<l[2];
//      cout<<" r= "<<r[0]<<" "<<r[1]<<" "<<r[2]<<endl;

      if(i==1){x1.at(1)=x0[1]+f0.ustar*dt;x1.at(0)=x0[0];}

      x1.at(2*i)=x0[2*i]+f0.ustar*dt;x1.at(2*i+1)=x0[2*i+1]+f1.ustar*dt;

      if(i==n){x1.at(2*(n+1))=x0[2*(n+1)]+f1.ustar*dt;x1.at(2*(n+1)+1)=x0[2*(n+1)+1];}

    }

//    l[0]=R.density(2*i);l[1]=R.velocity(2*i);l[2]=R.pressure(2*i);
//    r[0]=R.density(2*i+1);r[1]=R.velocity(2*i+1);r[2]=R.pressure(2*i+1);
//    Riemann R1(Riemann::exact,l,r);
//    for(long i=0;i<n;i++){x1.at(2*i+2)=x0[2*i+2]+R.velocity(2*i)*dt;x1.at(2*i+3)=x0[2*i+3]+R1.ustar*dt;}

// update cell volumes at the full-step

    for(int i=0;i<ng;i++){V1.at(i)=x1[2*i+1]-x1[2*i];if(V1[i]<0.0){cout<<"ERROR:  -'ve volume in cell "<<i<<endl;exit(1);}} 

// update cell density at the full-step

    for(int i=0;i<ng;i++){d.at(i)=m[i]/V1[i];} 

// update cell energy at the full-step

    for(int i=0;i<ng;i++){ec1.at(i)=max(ECUT,ec0[i]-(p[i]*(V1[i]-V0[i]))/m[i]);}

// construct the DG energy field

    for(int i=1;i<=n;i++){

      double dx(x1[2*i+1]-x1[2*i]); // cell width for Jacobian

// fluxes on face 0 (left boundary of cell)

      l[0]=d[i-1];l[1]=R.velocity(i-1);l[2]=p[i-1];
      r[0]=d[i];r[1]=R.velocity(i);r[2]=p[i];
      Riemann f0(Riemann::exact,l,r);

// fluxes on face 1 (right boundary of cell)

      l[0]=d[i];l[1]=R.velocity(i);l[2]=p[i];
      r[0]=d[i+1];r[1]=R.velocity(i+1);r[2]=p[i+1];
      Riemann f1(Riemann::exact,l,r);

// pressure and velocity on each face

      double pstar[2]={f0.pstar,f1.pstar};
      double ustar[2]={f0.ustar,f1.ustar};

// assemble DG matrix equation for the energy field in one element

      for(int iloc=0;iloc<nloc;iloc++){
        b[iloc]=0.0;soln[iloc]=0.0;
        for(int jloc=0;jloc<nloc;jloc++){
          NN[iloc][jloc]=0.0;NXN[iloc][jloc]=0.0;NNX[iloc][jloc]=0.0;A.write(iloc,jloc,0.0);
          for(int gi=0;gi<ngi;gi++){
            NN[iloc][jloc]+=N[iloc][gi]*N[jloc][gi]*dx/2.0; // mass matrix
//            NXN[iloc][jloc]+=NX[iloc][gi]*N[jloc][gi];  // divergence term (for use if by parts)
            NNX[jloc][iloc]-=N[iloc][gi]*NX[jloc][gi];  // divergence term
          }
          A.write(iloc,jloc,NN[iloc][jloc]);                // commit to the matrix class
//          b[iloc]+=(NXN[iloc][jloc]*ustar[jloc]-normal[2*i+iloc]*SN[iloc][jloc]*pstar[jloc])/d[i]; // source (by parts)
          b[iloc]+=NNX[iloc][jloc]*ustar[jloc]/d[i]; // source
        }
      }

// solution for one element

      A.solve(soln,b);e1[2*i]=max(ECUT,e0[2*i]+soln[0]*dt); e1[2*i+1]=max(ECUT,e0[2*i+1]+soln[1]*dt);

// debug
//      cout<<"el "<<i<<" n= "<<normal[2*i]<<" "<<normal[2*i+1]<<" (xl,xr)= "<<x1[2*i]<<" "<<x1[2*i+1]<<endl;
//      cout<<"         ustar[]="<<ustar[0]<<" "<<ustar[1]<<" pstar[]="<<pstar[0]<<" "<<pstar[1]<<endl;
//      for(int iloc=0;iloc<nloc;iloc++){
//        cout<<iloc<<" ";
//        for(int jloc=0;jloc<nloc;jloc++){cout<<NN[iloc][jloc]<<" ";}
//        cout<<endl;
//      }
//      cout<<endl;
//      for(int iloc=0;iloc<nloc;iloc++){
//        cout<<iloc<<" ";
//        for(int jloc=0;jloc<nloc;jloc++){cout<<NXN[iloc][jloc]<<" ";}
//        cout<<endl;
//     }
//      cout<<endl;
//      exit(1);
// debug

    }

// update cell pressure at the full-step using either energy field

    for(int i=0;i<ng;i++){p.at(i)=P(d[i],ec1[i]);} // use PdV

// update nodal velocities at the full step

// ...

// some output

//    for(int i=0;i<ng;i++){cout<<rx[i]<<" "<<R.density(i)<<" "<<R.pressure(i)<<" "<<R.velocity(i)<<" "<<R.energy(i)<<endl;}
    for(long i=1;i<=n;i++){cout<<0.5*(x1[2*i]+x1[2*i+1])<<" "<<d[i]<<" "<<p[i]<<" "<<" "<<ec1[i]<<endl;}
//    for(long i=1;i<=n;i++){cout<<x1[2*i]<<" "<<e1[i]<<endl;cout<<x1[2*i+1]<<" "<<e1[2*i+1]<<endl;}

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
