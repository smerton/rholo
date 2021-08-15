// Finite element variant of RhoLo (Riemann-based Hydro in One-dimension at Low Order - RhoLo)
// RhoLo is an ultra simple 1-D discontinuous finite element (DG) hydrodynamics test code
// This finite element variant solves the Euler equations in their non-conservative form in the fluid frame (the Lagrangian frame)
// using a mixed continuous finite element method (cell-centred thermodynamic variable d,rho,e with node 
// centred kinematic variables u,a) and bulk viscosity q to increase entropy across element boundaries, initial 
// implementation is only first order in time
// for graphics: convert -density 300 filename.png filename.pdf

// Author S. R. Merton

#define DTSTART 0.0005        // insert a macro for the first time step
#define ENDTIME 0.25          // insert a macro for the end time
#define GAMMA 1.4             // ratio of specific heats for ideal gases
#define ECUT 1.0e-8           // cut-off on the energy field
#define NSAMPLES 500          // number of sample points for the exact solution
#define VISFREQ 10000         // frequency of the graphics dumps
#define VD vector<double>     // vector of doubles
#define VTOL 1.0e-10          // threshold for volume errors
#define COURANT 0.333         // Courant number for CFL condition
#define DTSFACTOR 0.5         // safety factor on time-step control
#define KNOD i*(K.nloc()-1)+j // global node number on kinematic mesh
#define TNOD i*T.nloc()+j     // global node number on thermodynamic mesh
#define GPNT i*T.ngi()+gi     // global address of Gauss point gi in element i
#define DX0 x0[i*(K.nloc()-1)+K.nloc()-1]-x0[i*(K.nloc()-1)]            // cell width at start of step
#define DX1 x1[i*(K.nloc()-1)+K.nloc()-1]-x1[i*(K.nloc()-1)]            // cell width at end of step
#define CENTROID 0.5*(x1[i*(K.nloc()-1)+K.nloc()-1]+x1[i*(K.nloc()-1)]) // cell centroid

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

  ofstream f1,f2,f3,f4;                                 // files for output
  Shape K(2,3),T(1,3);                                  // p_n,p_n-1 shape functions
  int const n(10),ng(n+4);                              // no. ncells, no. ghosts
  int long nk(n*(K.nloc()-1)+1),nkg(ng*(K.nloc()-1)+1); // no. kinematic nodes, no. kinematic ghosts
  int long nt(n*T.nloc()),ntg(ng*T.nloc());             // no. thermodynamic nodes, no. thermodynamic ghosts
  double const cl(0.3),cq(1.0);                         // linear & quadratic coefficients for bulk viscosity
  vector<double> dinit(ng);                             // initial density field inside an element
  vector<double> d0(ng*T.ngi()),d1(ng*T.ngi());         // density at each Gauss point in each element
  vector<double> V0(ng),V1(ng),m(ng),xc(ng);            // volume, mass & centroid
  vector<double> e0(ntg),e1(ntg);                       // discontinuous FE energy field
  vector<double> c(ng*T.ngi()),p(ng*T.ngi());           // element sound speed & pressure at each Gauss point
  vector<double> q(ng*T.ngi());                         // bulk viscosity at each Gauss point
  vector<double> u0(nkg),u1(nkg);                       // node velocity
  vector<double> x0(nkg),x1(nkg),x2(ntg),x3(ntg);       // node coordinates
  vector<double> dt_cfl(ng*T.ngi());                    // element time-step at each Gauss point
  vector<double> detJ0(ng*T.ngi()),detJ(ng*T.ngi());    // determinant of the Jacobian
  double ke(0.0),ie(0.0);                               // kinetic and internal energy for conservation checks
  double time(0.0),dt(DTSTART);                         // start time and time step
  int step(0);                                          // step number
  double l[3]={1.0,0.0,1.0},r[3]={0.125,0.0,0.1};       // left/right flux states for the problem: Sod
//  double l[3]={1.0,-2.0,0.4},r[3]={1.0,2.0,0.4};      // left/right flux states for the problem: 123 (R2R)
//  double l[3]={1.0,0.0,1000.0},r[3]={1.0,0.0,0.01};   // left/right flux states for the problem: blast wave

// initialise the problem

  double dx(1.0/n);x0.at(0)=-2.0*dx;x1.at(0)=x0[0];
  for(long i=0;i<nkg;i++){x0.at(i)=x0[0]+i*dx/(K.nloc()-1);x1.at(i)=x0[i];}
  for(int i=0;i<ng;i++){xc.at(i)=CENTROID;}
  for(int i=0;i<ng;i++){for(int gi=0;gi<T.ngi();gi++){p.at(GPNT)=(xc[i]<=0.5)?l[2]:r[2];}}
  for(int i=0;i<ng;i++){for(int gi=0;gi<T.ngi();gi++){d0.at(GPNT)=(xc[i]<=0.5)?l[0]:r[0];}}
  for(int i=0;i<ng;i++){for(int gi=0;gi<T.ngi();gi++){d1.at(GPNT)=(xc[i]<=0.5)?l[0]:r[0];}}
  for(int i=0;i<ng;i++){for(int j=0;j<T.nloc();j++){e0.at(TNOD)=(xc[i]<=0.5)?E(l[0],l[2]):E(r[0],r[2]);}}
  for(int i=0;i<ng;i++){for(int j=0;j<T.nloc();j++){e1.at(TNOD)=(xc[i]<=0.5)?E(l[0],l[2]):E(r[0],r[2]);}}
  for(int i=0;i<ng;i++){V0.at(i)=DX0;V1.at(i)=DX1;}
  for(int i=0;i<ng;i++){dinit.at(i)=(xc[i]<=0.5)?l[0]:r[0];}
  for(int i=0;i<ng;i++){m.at(i)=(xc[i]<=0.5)?l[0]*V0[i]:r[0]*V0[i];}
  for(long i=0;i<nkg;i++){u0.at(i)=(x0[i]<=0.5)?l[1]:r[1];u1.at(i)=u0[i];}
  for(int i=0;i<ng;i++){for(int gi=0;gi<T.ngi();gi++){q.at(GPNT)=0.0;}}
  for(int i=0;i<ng;i++){for(int gi=0;gi<T.ngi();gi++){c.at(GPNT)=sqrt(GAMMA*p[GPNT]/d0[GPNT]);}}

// set thermodynamic node positions at time-0

  for(int i=0;i<ng;i++){
    for(int t=0;t<T.nloc();t++){
      double pos(-1.0+t*2.0/(T.nloc()-1));long tloc(i*T.nloc()+t);x2.at(tloc)=0.0;
      for(int j=0;j<K.nloc();j++){x2.at(tloc)+=K.value(j,pos)*x0[KNOD];}
    }
  }

// set time-0 Jacobian

  for(int i=0;i<ng;i++){
    for(int gi=0;gi<T.ngi();gi++){
      detJ0.at(GPNT)=0.0;
      for(int j=0;j<K.nloc();j++){
        detJ0.at(GPNT)+=K.dvalue(j,gi)*x0[KNOD];
      }
      if(detJ0.at(GPNT)<0.0){cout<<"-'ve determinant of J0 detected in cell "<<i<<endl;exit(1);}
    }
  }

// start the Riemann solvers from initial flux states

  Riemann R0(Riemann::exact,l,r),R1(Riemann::exact,l,r);

// set output precision

  cout<<fixed<<setprecision(17);

// time integration

  while(time<ENDTIME+dt){

// calculate a new stable time-step that will impose the CFL limit on each quadrature point

    for(int i=0;i<ng;i++){
      for(int gi=0;gi<K.ngi();gi++){
        dt_cfl.at(GPNT)=COURANT*(DX0/sqrt((c[GPNT]*c[GPNT])+2.0*q[GPNT]/d0[GPNT]));
      }
    }

// reduce across element and apply a saftey factor

    double dt=DTSFACTOR*(*min_element(dt_cfl.begin(), dt_cfl.end()));

    cout<<fixed<<setprecision(5)<<"  step "<<step<<" time= "<<time<<" dt= "<<dt;
    cout<<fixed<<setprecision(5)<<" energy (i/k/tot)= "<<ie<<" "<<ke<<" "<<ie+ke<<endl;

// move the nodes to their full-step position

    for(long i=0;i<nkg;i++){x1.at(i)=x0[i]+u0[i]*dt;}

    for(int i=0;i<ng;i++){
      for(int t=0;t<T.nloc();t++){
        double pos(-1.0+t*2.0/(T.nloc()-1));long tloc(i*T.nloc()+t);x3.at(tloc)=0.0;
        for(int j=0;j<K.nloc();j++){x3.at(tloc)+=K.value(j,pos)*x1[KNOD];}
      }
    }

// update mesh centroids

    for(int i=0;i<ng;i++){xc.at(i)=CENTROID;}

// update kinetic energy for conservation checks

    ke=0.0;for(int i=0;i<ng;i++){for(int j=0;j<K.nloc();j++){ke+=0.25*(m[i])*u0[KNOD]*u0[KNOD];}}

// evolve the Riemann problems to the end of the time-step on the end of time-step meshes

    vector<double> r0x,rx;vempty(r0x);vempty(rx); // sample point coordinates
    for(long i=0;i<NSAMPLES;i++){r0x.push_back(0.0+(i*(1.0-0.0)/double(NSAMPLES)));} // sample points
    R0.profile(&r0x,time+dt); // Riemann solution at the sample points along the mesh
    for(int i=0;i<nkg;i++){rx.push_back(x0[i]);}
    R1.profile(&rx,time+dt);

// update cell volumes at the full-step

    for(int i=0;i<ng;i++){V1.at(i)=DX1;if(V1[i]<VTOL){cout<<"-'ve volume detected in cell "<<i<<endl;exit(1);}}

// update Jacobian

    for(int i=0;i<ng;i++){
      for(int gi=0;gi<T.ngi();gi++){
        detJ.at(GPNT)=0.0;
        for(int j=0;j<K.nloc();j++){
          detJ.at(GPNT)+=K.dvalue(j,gi)*x1[KNOD];
        }
        if(detJ.at(GPNT)<0.0){cout<<"-'ve determinant of J detected in cell "<<i<<endl;exit(1);}
      }
    }

// update cell density at the full-step

    for(int i=0;i<ng;i++){for(int gi=0;gi<T.ngi();gi++){d1[GPNT]=dinit[i]*detJ0[GPNT]/detJ[GPNT];}}

// assemble finite element energy field on discontinuous thermodynamic grid

    {Matrix A(T.nloc());double b[T.nloc()],x[T.nloc()];
    for(int i=0;i<ng;i++){
      for(int iloc=0;iloc<T.nloc();iloc++){
        for(int jloc=0;jloc<T.nloc();jloc++){
          double nn(0.0); // DG mass matrix
           for(int gi=0;gi<T.ngi();gi++){
             nn+=d1[GPNT]*T.value(iloc,gi)*T.value(jloc,gi)*detJ[GPNT]*T.wgt(gi);
           }
           A.write(iloc,jloc,nn);
        }
      }

// solve local system

      A.solve(x,b);

// advance the solution

      for(int j=0;j<T.nloc();j++){e1.at(TNOD)=max(ECUT,e0[TNOD]-x[j]*dt);}

    }}







// debug
  cout<<"debug stop."<<endl;
  exit(1);
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
    for(int iloc=0;iloc<K.nloc();iloc++){
      int i(iel+iloc); // column address in the global matrix
      if((i>0&&i<ng)){for(int gi=0;gi<K.ngi();gi++){b[i-1]+=(p[iel]+q[iel])*K.dvalue(iloc,gi)*K.wgt(gi);}} // integrate the shape derivative for rhs
      for(int jloc=0;jloc<K.nloc();jloc++){
        double nn(0.0); // mass matrix
        int j(iel+jloc); // row address in the global matrix
        for(int gi=0;gi<K.ngi();gi++){
          nn+=K.value(iloc,gi)*K.value(jloc,gi)*K.wgt(gi)*0.5*(x1[iel+1]-x1[iel]); // DG & notes use this - double check ??
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

    f1.open("exact.dat");f2.open("dpe.dat");f3.open("q.dat");f4.open("u.dat");
    f1<<fixed<<setprecision(17);f2<<fixed<<setprecision(17);f3<<fixed<<setprecision(17);f4<<fixed<<setprecision(17);

    for(int i=0;i<NSAMPLES;i++){f1<<r0x[i]<<" "<<R0.density(i)<<" "<<R0.pressure(i)<<" "<<R0.velocity(i)<<" "<<R0.energy(i)<<endl;}
    for(int i=2;i<n+2;i++){f2<<0.5*(x1[i]+x1[i+1])<<" "<<d1[i]<<" "<<p[i]<<" "<<e1[i]<<" "<<endl;}
    for(int i=2;i<n+2;i++){f3<<0.5*(x1[i]+x1[i+1])<<" "<<q[i]<<endl;}
    for(int i=2;i<n+3;i++){f4<<x1[i]<<" "<<u1[i]<<endl;}

    f1.close();f2.close();f3.close();f4.close();

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

  cout<<"Normal termination."<<endl;

  return 0;
}

// return pressure given the energy

double P(double d,double e){return (GAMMA-1.0)*d*e;}

// invert the eos to return energy given the pressure

double E(double d,double p){return p/((GAMMA-1.0)*d);}

// empty a vector

void vempty(vector<double>&v){vector<double> e;v.swap(e);return;}
