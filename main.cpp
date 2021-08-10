// Finite element variant of RhoLo (Riemann-based Hydro in One-dimension at Low Order - RhoLo)
// RhoLo is an ultra simple 1-D discontinuous finite element (DG) hydrodynamics test code
// This finite element variant solves the Euler equations in their non-conservative form in the fluid frame (the Lagrangian frame)
// using a mixed continuous finite element method (cell-centred thermodynamic variable d,rho,e with node 
// centred kinematic variables u,a) and bulk viscosity q to increase entropy across element boundaries, initial 
// implementation is only first order in time
// for graphics: convert -density 300 filename.png filename.pdf

// Author S. R. Merton

#define DTSTART 0.0005    // insert a macro for the first time step
#define ENDTIME 0.25      // insert a macro for the end time
#define GAMMA 1.4         // ratio of specific heats for ideal gases
#define ECUT 1.0e-8        // cut-off on the energy field
#define NSAMPLES 500      // number of sample points for the exact solution
#define VISFREQ 10000     // frequency of the graphics dumps
#define VD vector<double> // vector of doubles
#define VTOL 1.0e-10      // threshold for volume errors
#define COURANT 0.333     // Courant number for CFL condition
#define DTSFACTOR 0.5     // safety factor on time-step control
#define GNOD i*(S.nloc()-1)+j // global FE node address
#define KNOD i*(K.nloc()-1)+j // global kinematic node address
#define TNOD i*(T.nloc()-1)+j // global thermodynamic node address
#define DX0 x0[i*(S.nloc()-1)+S.nloc()-1]-x0[i*(S.nloc()-1)] // width of cell i at t0
#define DX1 x1[i*(S.nloc()-1)+S.nloc()-1]-x1[i*(S.nloc()-1)] // width of cell i at t1
#define CENTROID 0.5*(x1[i*(S.nloc()-1)+S.nloc()-1]+x1[i*(S.nloc()-1)]) // centroid of cell i

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

  ofstream f1,f2,f3,f4,f5;                                 // files for output
  int const n(100),ng(n+4);                             // no. ncells, no. ghosts
  Shape S(1,2),T(1,3),K(2,3);                           // load FE stencils
  int nnodes(n*(S.nloc()-1)+1),ngnodes(ng*(S.nloc()-1)+1); // no. FE nodes
  int nknodes(n*(K.nloc()-1)+1),ngknodes(ng*(K.nloc()-1)+1); // no. kinematic nodes
  int ntnodes(n*(T.nloc()-1)+1),ngtnodes(ng*(T.nloc()-1)+1); // no. thermodynamic nodes
  double const cl(0.3),cq(1.0);                         // linear & quadratic coefficients for bulk viscosity
  vector<double> dinit(ng);                             // initial density field
  vector<double> d0(ng),d1(ng),V0(ng),V1(ng),m(ng);     // density, volume & mass
  vector<double> e0(ng),e1(ng);                         // cell-centred energy field
  vector<double> e2(ng*T.nloc()),e3(ng*T.nloc());       // discontinuous energy field
  vector<double> c(ng),p(ng),q(ng);                     // element sound speed, pressure and bulk viscosity
  vector<double> u0(ngnodes),u1(ngnodes);               // node velocity
  vector<double> u2(ngknodes),u3(ngknodes);             // high-order node velocity
  vector<double> x0(ngnodes),x1(ngnodes);               // node coordinates
  vector<double> x2(ngtnodes),x3(ngtnodes);             // thermodynamic node coordinates
  vector<double> x4(ngknodes),x5(ngknodes);             // kinematic node coordinates
  vector<double> xc(ng);                                // cell centroids
  vector<double> dt_cfl(ng);                            // element time-step
  vector<double> qd(T.ngi()),qp(T.ngi()),qq(T.ngi());   // data at each integration point
  vector<double> qc(T.ngi()),qe(T.ngi());               // data at each integration point
  double Vn[S.nloc()],mn[ngnodes];                      // volume and mass of an FE node
  double detJ[ng*K.ngi()],detJ0[ng*K.ngi()];            // determinant of the Jacobian at each Gauss point
  double ke(0.0),ie(0.0);                               // kinetic and internal energy for conservation checks
  double time(0.0),dt(DTSTART);                         // start time and time step
  int step(0);                                          // step number
  double l[3]={1.0,0.0,1.0},r[3]={0.125,0.0,0.1};       // left/right flux states for the problem
//  double l[3]={1.0,-2.0,0.4},r[3]={1.0,2.0,0.4};      // left/right flux states for the problem
//  double l[3]={1.0,0.0,1000.0},r[3]={1.0,0.0,0.01};   // left/right flux states for the problem

// initialise the problem

  double dx(1.0/n);x0.at(0)=-2.0*dx;x1.at(0)=x0[0];xc[0]=-1.5*dx;
  x2.at(0)=x0[0];x3.at(0)=x0[0];x4.at(0)=x0[0];x5.at(0)=x0[0];
  for(int i=1;i<ng;i++){xc.at(i)=xc[i-1]+dx;}
  for(int i=1;i<ngnodes;i++){x0.at(i)=x0[i-1]+1.0/(n*(S.nloc()-1));x1.at(i)=x0[i];}
  for(int i=1;i<ngtnodes;i++){x2.at(i)=x2[i-1]+1.0/(n*(T.nloc()-1));x3.at(i)=x2[i];}
  for(int i=1;i<ngknodes;i++){x4.at(i)=x4[i-1]+1.0/(n*(K.nloc()-1));x5.at(i)=x4[i];}
  for(int i=0;i<ng;i++){p.at(i)=(xc[i]<=0.5)?l[2]:r[2];}
  for(int i=0;i<ng;i++){dinit.at(i)=(xc[i]<=0.5)?l[0]:r[0];}
  for(int i=0;i<ng;i++){d0.at(i)=(xc[i]<=0.5)?l[0]:r[0];}
  for(int i=0;i<ng;i++){d1.at(i)=(xc[i]<=0.5)?l[0]:r[0];}
  for(int i=0;i<ng;i++){e0.at(i)=E(d0[i],p[i]);}
  for(int i=0;i<ng;i++){e1.at(i)=E(d1[i],p[i]);}
  for(int i=0;i<ng;i++){for(int j=0;j<T.nloc();j++){e2.at(i*T.nloc()+j)=E(d0[i],p[i]);}}
  for(int i=0;i<ng;i++){for(int j=0;j<T.nloc();j++){e3.at(i*T.nloc()+j)=E(d1[i],p[i]);}}
  for(int i=0;i<ng;i++){for(int j=0;j<S.nloc();j++){mn[GNOD]=0.0;}}

// Initialise jacobian for strong mass conservation

  for(int i=0;i<ng;i++){
    for(int gi=0;gi<K.ngi();gi++){
      detJ0[i*K.ngi()+gi]=0.0;detJ[i*K.ngi()+gi]=0.0;
      for(int j=0;j<K.nloc();j++){detJ0[i*K.ngi()+gi]+=K.dvalue(j,gi)*x4[KNOD];}
      for(int j=0;j<K.nloc();j++){detJ[i*K.ngi()+gi]+=K.dvalue(j,gi)*x5[KNOD];}
    }
  }

  for(int i=0;i<ng;i++){
    V0.at(i)=0.0;V1.at(i)=0.0;double detJ[S.ngi()];
    for(int gi=0;gi<S.ngi();gi++){detJ[gi]=0.0;for(int j=0;j<S.nloc();j++){detJ[gi]+=S.dvalue(j,gi)*x0[GNOD];}}
    for(int j=0;j<S.nloc();j++){Vn[j]=0.0;for(int gi=0;gi<S.ngi();gi++){Vn[j]+=S.value(j,gi)*detJ[gi]*S.wgt(gi);}}
    for(int j=0;j<S.nloc();j++){V0.at(i)+=Vn[j];V1.at(i)=V0[i];mn[GNOD]+=d0[i]*Vn[j];}
  }

  for(int i=0;i<ng;i++){m.at(i)=d0[i]*V0[i];}
  for(int i=0;i<ng;i++){for(int j=0;j<S.nloc();j++){u0.at(GNOD)=(x0[GNOD]<=0.5)?l[1]:r[1];u1.at(GNOD)=u0[GNOD];}}
  for(int i=0;i<ng;i++){for(int j=0;j<K.nloc();j++){u2.at(KNOD)=(x4[KNOD]<=0.5)?l[1]:r[1];u3.at(KNOD)=u2[KNOD];}}

  for(int i=0;i<ng;i++){q.at(i)=0.0;}
  for(int i=0;i<ng;i++){c.at(i)=sqrt(GAMMA*p[i]/d0[i]);}
  for(int i=2;i<ng-1;i++){for(int j=0;j<S.nloc();j++){ke+=0.5*mn[GNOD]*u0[GNOD]*u0[GNOD];}}
  for(int i=0;i<ng;i++){ie+=e1[i]*m[i];}

// start the Riemann solvers from initial flux states

  Riemann R0(Riemann::exact,l,r),R1(Riemann::exact,l,r),R2(Riemann::exact,l,r);

// set output precision

  cout<<fixed<<setprecision(17);

// time integration

  while(time<ENDTIME+dt){

// calculate a new stable time-step

    for(int i=0;i<ng;i++){double l(DX0);dt_cfl.at(i)=(COURANT*l/sqrt((c[i]*c[i])+2.0*q[i]/d0[i]));} // impose the CFL limit on each element

    double dt=DTSFACTOR*(*min_element(dt_cfl.begin(), dt_cfl.end())); // reduce across element and apply a saftey factor

    cout<<fixed<<setprecision(5)<<"  step "<<step<<" time= "<<time<<" dt= "<<dt<<fixed<<setprecision(5)<<" energy (i/k/tot)= "<<ie<<" "<<ke<<" "<<ie+ke<<endl;

// move the nodes to their full-step position

    for(int i=0;i<ng;i++){for(int j=0;j<S.nloc();j++){x1.at(GNOD)=x0[GNOD]+u0[GNOD]*dt;}xc[i]=CENTROID;}

// debug
//    if(step>0){for(int i=0;i<ng;i++){for(int j=0;j<S.nloc();j++){x1.at(GNOD)=x0[GNOD]+R2.velocity(GNOD)*dt;}xc[i]=CENTROID;}}
// debug

// update kinetic energy for conservation checks

    ke=0.0;for(int i=2;i<ng-1;i++){for(int j=0;j<S.nloc();j++){ke+=0.5*mn[GNOD]*u0[GNOD]*u0[GNOD];}}

// evolve the Riemann problems to the end of the time-step on the end of time-step meshes

    vector<double> r0x,rx;vempty(r0x); // sample point coordinates
//    for(long i=0;i<NSAMPLES;i++){r0x.push_back(x1[0]+(i*(x1[ng]-x1[0])/double(NSAMPLES)));} // sample points
    for(long i=0;i<NSAMPLES;i++){r0x.push_back(0.0+(i*(1.0-0.0)/double(NSAMPLES)));} // sample points
    R0.profile(&r0x,time+dt); // Riemann solution at the sample points along the mesh
    rx.clear();for(int i=0;i<ng;i++){rx.push_back(xc[i]);}
    R1.profile(&rx,time+dt);
    rx.clear();for(int i=0;i<ng;i++){for(int j=0;j<S.nloc()-1;j++){rx.push_back(x1[GNOD]);}}
    R2.profile(&rx,time+dt);

// update cell volumes at the full-step

    for(int i=0;i<ng;i++){
      V1.at(i)=0.0;double detJ[S.ngi()];
      for(int gi=0;gi<S.ngi();gi++){detJ[gi]=0.0;for(int j=0;j<S.nloc();j++){detJ[gi]+=S.dvalue(j,gi)*x1[GNOD];}}
      for(int j=0;j<S.nloc();j++){for(int gi=0;gi<S.ngi();gi++){V1.at(i)+=S.value(j,gi)*detJ[gi]*S.wgt(gi);}}
    }

// update cell density at the full-step

    for(int i=0;i<ng;i++){d1.at(i)=m[i]/V1[i];}

// debug
//    for(int i=0;i<ng;i++){d1.at(i)=R1.density(i);} // 3*i+1 is cell-centre address
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
//    for(int i=0;i<ng;i++){e1.at(i)=R1.energy(i);} // 3*i+1 is cell-centre address
// debug

// update cell pressure at the full-step

    for(int i=0;i<ng;i++){p.at(i)=P(d1[i],e1[i]);if(p[i]<0.0){cout<<"-'ve pressure detected in cell "<<i<<" e1= "<<e1[i]<<endl;exit(1);}}

// bulk q

    for(int i=0;i<ng;i++){
      c.at(i)=sqrt(GAMMA*p[i]/d1[i]);
      double l(DX1),divu((d0[i]-d1[i])/(d1[i]*dt));
      if(divu<0.0){
        q.at(i)=d0[i]*l*divu*((cq*l*divu)-cl*c[i]);
      }else{
        q.at(i)=0.0; // turn off q as cell divergence indicates expansion
      }
    }

// assemble acceleration field

  int matsiz((n+2)*(S.nloc()-1)+1);
  Matrix A(matsiz);double b[matsiz],x[matsiz];for(int i=0;i<matsiz;i++){b[i]=0.0;x[i]=0.0;}

//  double m[S.nloc()][S.nloc()]={};m[0][0]=0.5;m[S.nloc()-1][S.nloc()-1]=0.5;// for mass lumping

// next block codes for matrix assembly with a continuous Galerkin finite element type

  for(int iel=0;iel<ng;iel++){
    double detJ[S.ngi()];
    for(int gi=0;gi<S.ngi();gi++){detJ[gi]=0.0;for(int j=0;j<S.nloc();j++){detJ[gi]+=S.dvalue(j,gi)*x1[iel*(S.nloc()-1)+j];}}
    for(int iloc=0;iloc<S.nloc();iloc++){
      int i(iel*(S.nloc()-1)+iloc); // column address in the global matrix
      if((i>=S.nloc()-1&&i<=(ng-1)*(S.nloc()-1))){for(int gi=0;gi<S.ngi();gi++){b[i-(S.nloc()-1)]+=(p[iel]+q[iel])*S.dvalue(iloc,gi)*S.wgt(gi);}} // integrate the shape derivative for rhs
      for(int jloc=0;jloc<S.nloc();jloc++){
        double nn(0.0); // mass matrix
        int j(iel*(S.nloc()-1)+jloc); // row address in the global matrix
        for(int gi=0;gi<S.ngi();gi++){
          nn+=S.value(iloc,gi)*S.value(jloc,gi)*S.wgt(gi)*detJ[gi]; // DG & notes use this - double check ??
        }
        if((i>=S.nloc()-1&&i<=(ng-1)*(S.nloc()-1))&&(j>=S.nloc()-1&&j<=(ng-1)*(S.nloc()-1))){
          A.add(i-(S.nloc()-1),j-(S.nloc()-1),d1[iel]*nn);
//          A.add(i-(S.nloc()-1),i-(S.nloc()-1),d1[iel]*nn); // use mass lumping
        }
      }
    }
  }

// solve global system

  A.solve(x,b);

// debug
//  for(int i=0;i<matsiz;i++){
//    cout<<i<<" vel= "<<x[i]<<" "<<b[i]<<endl;
//  }
//  for(int iloc=0;iloc<S.nloc();iloc++){
//    cout<<"iloc "<<iloc<<" shape values: ";
//    for(int gi=0;gi<S.ngi();gi++){
//      cout<<S.dvalue(iloc,gi)<<" ";
//    }
//    cout<<endl;
//  }
//  exit(1);
// debug

// advance the solution

  for(int i=0;i<matsiz;i++){u1.at(i+S.nloc()-1)=u0[i+S.nloc()-1]+x[i]*dt;}

// impose boundary constraints on the acceleration field

//  u1.at(1)=u1[2];u1.at(ng-1)=u1[ng-2];
//  u1.at(0)=u1[1];u1.at(ng)=u1[ng-1];

// some output

    f1.open("exact.dat");f2.open("dpe.dat");f3.open("q.dat");f4.open("u.dat");f5.open("e.dat");
    f1<<fixed<<setprecision(17);f2<<fixed<<setprecision(17);f3<<fixed<<setprecision(17);f4<<fixed<<setprecision(17);f5<<fixed<<setprecision(17);

    for(int i=0;i<NSAMPLES;i++){f1<<r0x[i]<<" "<<R0.density(i)<<" "<<R0.pressure(i)<<" "<<R0.velocity(i)<<" "<<R0.energy(i)<<endl;}
    for(int i=2;i<n+2;i++){f2<<xc[i]<<" "<<d1[i]<<" "<<p[i]<<" "<<e1[i]<<" "<<endl;}
    for(int i=2;i<n+2;i++){f3<<xc[i]<<" "<<q[i]<<endl;}
    for(int i=2;i<n+2;i++){for(int j=0;j<S.nloc();j++){f4<<x1[GNOD]<<" "<<u1[GNOD]<<endl;}}
    for(int i=2;i<n+2;i++){for(int j=0;j<T.nloc();j++){f5<<x3[TNOD]<<" "<<e3[i*T.nloc()+j]<<endl;}}

    f1.close();f2.close();f3.close();f4.close();f5.close();

// advance the time step

    time+=dt;
    step++;

// debug
//      if(step==20){cout<<" debug stop 1..."<<endl;exit(1);}
//    for(int i=0;i<ng;i++){for(int j=0;j<S.nloc();j++){u1.at(GNOD)=R2.velocity(GNOD);}} // 3*i,3*i+2 are nodal address
// debug

// advance the solution for the new time step

    for(int i=0;i<ng;i++){for(int j=0;j<S.nloc();j++){u0.at(GNOD)=u1[GNOD];}}
    for(int i=0;i<ng;i++){for(int j=0;j<K.nloc();j++){u2.at(KNOD)=u3[KNOD];}}
    for(int i=0;i<ng;i++){for(int j=0;j<S.nloc();j++){x0.at(GNOD)=x1[GNOD];}}
    for(int i=0;i<ng*(T.nloc()-1)+1;i++){x2.at(i)=x3[i];}
    for(int i=0;i<ng*(K.nloc()-1)+1;i++){x4.at(i)=x5[i];}
    for(int i=0;i<ng*T.nloc();i++){e2.at(i)=e3[i];}
    for(int i=0;i<ng;i++){e0.at(i)=e1[i];}
    for(int i=0;i<ng;i++){V0.at(i)=V1[i];}
    for(int i=0;i<ng;i++){d0.at(i)=d1[i];}
//    for(int i=0;i<ng;i++){c.at(i)=sqrt(GAMMA*(p[i]+q[i])/d1[i]);}

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
