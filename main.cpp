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
// git push -vu NAS high-order
//
// then check url's have been added with something like
//
// git remote -v
//
// and also check .git/config
//
// then git pull --all will pull from 1st url in location 1 (NAS)
// and 1st url in location 2 (github)
//
// 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
//
// for graphics: convert -density 300 filename.png filename.pdf

// Author S. R. Merton

#define DTSTART 0.0005          // insert a macro for the first time step
#define ENDTIME 0.20            // insert a macro for the end time
#define GAMMA 1.4               // ratio of specific heats for ideal gases
#define ECUT 1.0e-8             // cut-off on the energy field
#define NSAMPLES 500            // number of sample points for the exact solution
#define VISFREQ 10000           // frequency of the graphics dumps
#define VD vector<double>       // vector of doubles
#define VTOL 1.0e-10            // threshold for volume errors
#define COURANT 0.333           // Courant number for CFL condition
#define DTSFACTOR 0.75          // safety factor on time-step control
#define KNOD i*(K.nloc()-1)+k   // global node number on kinematic mesh
#define TNOD i*T.nloc()+j       // global node number on thermodynamic mesh
#define GPNT i*T.ngi()+gi       // global address of Gauss point gi in element i
#define DX0 (x0[i*(K.nloc()-1)+K.nloc()-1]-x0[i*(K.nloc()-1)])          // cell width at start of step
#define DX1 (x1[i*(K.nloc()-1)+K.nloc()-1]-x1[i*(K.nloc()-1)])          // cell width at end of step
#define CENTROID 0.5*(x1[i*(K.nloc()-1)+K.nloc()-1]+x1[i*(K.nloc()-1)]) // cell centroid
#define ROW (i-1)*(K.nloc()-1)+iloc                                     // row address in global matrix
#define COL (i-1)*(K.nloc()-1)+jloc                                     // column address in global matrix
#define XGI for(int j=0;j<T.nloc();j++){xgi+=T.value(j,gi)*x3[TNOD];}   // coordinates of integration point
#define NROWS (ng-2)*(K.nloc()-1)+1                                     // number of rows in the global matrix
#define NCOLS (ng-2)*(K.nloc()-1)+1                                     // number of columns in the global matrix (= no. rows)

#include <iostream>
#include <vector>
#include <iomanip>
#include <cmath>
#include <fstream>
#include "riemann.h"
#include "matrix.h"
#include "shape.h"
#include <bits/stdc++.h>
#include <chrono>

// sigantures for eos lookups

double P(double d,double e); // eos returns pressure as a function of energy
double E(double d,double p); // invert the eos to get energy if we only have pressure
void vempty(vector<double>&v); // signature for emptying a vector

using namespace std;
using namespace std::chrono;

int main(){

// get a timepoint

  auto start_main = high_resolution_clock::now();

  cout<<"main(): Starting up main loop..."<<endl;

// global data

  ofstream f1,f2,f3,f4,f5;                              // files for output
  Shape K(2,10),T(1,10);                                // p_n,q_n-1 shape functions
  int const n(50),ng(n+4);                              // no. ncells, no. ghosts
  int long nk(n*(K.nloc()-1)+1),nkg(ng*(K.nloc()-1)+1); // no. kinematic nodes, no. kinematic ghosts
  int long nt(n*T.nloc()),ntg(ng*T.nloc());             // no. thermodynamic nodes, no. thermodynamic ghosts
  double const cl(1.0),cq(1.0);                         // linear & quadratic coefficients for bulk viscosity
  vector<double> dinit(ng);                             // initial density field inside an element
  vector<double> d0_t(ng*T.ngi()),d1_t(ng*T.ngi());     // density at each Gauss point in each element
  vector<double> d0_k(ng*K.ngi()),d1_k(ng*K.ngi());     // density at each Gauss point in each element
  vector<double> V0(ng),V1(ng),m(ng),xc(ng);            // volume, mass & centroid
  vector<double> nodmass_t(ntg),nodmass_k(nkg);         // nodal mass
  vector<double> nodvol_t(ntg),nodvol_k(nkg);           // nodal volume
  vector<double> e0(ntg),e1(ntg);                       // discontinuous FE energy field
  vector<double> c(ng*T.ngi()),p(ng*T.ngi());           // element sound speed & pressure at each Gauss point
  vector<double> q(ng*T.ngi());                         // bulk viscosity at each Gauss point
  vector<double> u0(nkg),u1(nkg);                       // node velocity
  vector<double> x0(nkg),x1(nkg),x2(ntg),x3(ntg);       // node coordinates
  vector<double> dt_cfl(ng*T.ngi());                    // element time-step at each Gauss point
  vector<double> detJ0_t(ng*T.ngi()),detJ_t(ng*T.ngi());    // determinant of the Jacobian
  vector<double> detJ0_k(ng*K.ngi()),detJ_k(ng*K.ngi());    // determinant of the Jacobian
  vector<double> l0(ng);                                // initial length scale
//  double F[nkg][ntg];                                 // force matrix (located on stack)
  double** F=new double*[nkg];                          // force matrix (located on  heap)
  for(long i=0;i<nkg;i++){
    F[i]=new double[ntg];
  }

  double ke(0.0),ie(0.0);                               // kinetic and internal energy for conservation checks
  double time(0.0),dt(DTSTART);                         // start time and time step
  int step(0);                                          // step number
  high_resolution_clock::time_point start,stop;         // high resolution timers
  duration<double> time_span[10];                       // time accumulated in different parts of the code
  double l[3]={1.0,0.0,1.0},r[3]={0.125,0.0,0.1};       // left/right flux states for the problem: Sod
//  double l[3]={1.0,-2.0,0.4},r[3]={1.0,2.0,0.4};      // left/right flux states for the problem: 123 (R2R)
//  double l[3]={1.0,0.0,1000.0},r[3]={1.0,0.0,0.01};   // left/right flux states for the problem: blast wave

// initialise the problem

  double dx(1.0/n);x0.at(0)=-2.0*dx;x1.at(0)=x0[0];x2.at(0)=-2.0*dx;x3.at(0)=x2[0];
  for(long i=0;i<nkg;i++){x0.at(i)=x0[0]+i*dx/(K.nloc()-1);x1.at(i)=x0[i];}
  for(long i=0;i<ntg;i++){x2.at(i)=x2[0]+i*dx/(T.nloc()-1);x3.at(i)=x2[i];}
  for(int i=0;i<ng;i++){xc.at(i)=CENTROID;}
  for(int i=0;i<ng;i++){for(int gi=0;gi<K.ngi();gi++){p.at(GPNT)=(xc[i]<=0.5)?l[2]:r[2];}}
  for(int i=0;i<ng;i++){for(int gi=0;gi<T.ngi();gi++){d0_t.at(GPNT)=(xc[i]<=0.5)?l[0]:r[0];}}
  for(int i=0;i<ng;i++){for(int gi=0;gi<K.ngi();gi++){d0_k.at(GPNT)=(xc[i]<=0.5)?l[0]:r[0];}}
  for(int i=0;i<ng;i++){for(int gi=0;gi<T.ngi();gi++){d1_t.at(GPNT)=(xc[i]<=0.5)?l[0]:r[0];}}
  for(int i=0;i<ng;i++){for(int gi=0;gi<K.ngi();gi++){d1_k.at(GPNT)=(xc[i]<=0.5)?l[0]:r[0];}}
  for(int i=0;i<ng;i++){for(int j=0;j<T.nloc();j++){e0.at(TNOD)=(xc[i]<=0.5)?E(l[0],l[2]):E(r[0],r[2]);}}
  for(int i=0;i<ng;i++){for(int j=0;j<T.nloc();j++){e1.at(TNOD)=(xc[i]<=0.5)?E(l[0],l[2]):E(r[0],r[2]);}}
  for(int i=0;i<ng;i++){V0.at(i)=DX0;V1.at(i)=DX1;}
  for(int i=0;i<ng;i++){dinit.at(i)=(xc[i]<=0.5)?l[0]:r[0];}
  for(int i=0;i<ng;i++){m.at(i)=(xc[i]<=0.5)?l[0]*V0[i]:r[0]*V0[i];}
  for(long i=0;i<nkg;i++){u0.at(i)=(x0[i]<=0.5)?l[1]:r[1];u1.at(i)=u0[i];}
  for(int i=0;i<ng;i++){for(int gi=0;gi<T.ngi();gi++){q.at(GPNT)=0.0;}}
  for(int i=0;i<ng;i++){for(int gi=0;gi<T.ngi();gi++){c.at(GPNT)=sqrt(GAMMA*p[GPNT]/d0_k[GPNT]);}}
  for(long i=0;i<nkg;i++){for(long j=0;j<ntg;j++){F[i][j]=0.0;}}
  for(int i=0;i<ng;i++){l0.at(i)=DX0/(K.nloc()-1);} // initial nodal displacement

// initialise the timers

  for(int i=0;i<10;i++){
    start=high_resolution_clock::now();
    time_span[i]=duration_cast<duration<double>>(start-start);
  }

// check integration rules on thermodynamic and kinematic stencils look consistent, they need to match

  if(T.ngi()!=K.ngi()){
    cout<<"No. of integration points for thermodynamics ("<<T.ngi()<<") doesn't match the no. for kinematics ("<<K.ngi()<<")"<<endl;
    exit(1);
  }

// set thermodynamic node positions at time-0

  for(int i=0;i<ng;i++){
    for(int j=0;j<T.nloc();j++){
      double pos(-1.0+j*2.0/(T.nloc()-1));x2.at(TNOD)=0.0;
      for(int k=0;k<K.nloc();k++){x2.at(TNOD)+=K.value(k,pos)*x0[KNOD];}
    }
  }

// set time-0 Jacobians

  for(int i=0;i<ng;i++){

    for(int gi=0;gi<K.ngi();gi++){
      detJ0_k.at(GPNT)=0.0;
      for(int k=0;k<K.nloc();k++){
        detJ0_k.at(GPNT)+=K.dvalue(k,gi)*x0[KNOD];
      }
      if(detJ0_k.at(GPNT)<0.0){cout<<"-'ve determinant of J0_k detected in cell "<<i<<endl;exit(1);}
    }

    for(int gi=0;gi<T.ngi();gi++){
      detJ0_t.at(GPNT)=0.0;
      for(int j=0;j<T.nloc();j++){
        detJ0_t.at(GPNT)+=T.dvalue(j,gi)*x2[TNOD];
      }
      if(detJ0_t.at(GPNT)<0.0){cout<<"-'ve determinant of J0_t detected in cell "<<i<<endl;exit(1);}
    }

  }

// set nodal masses - these should not change with time

  for(long i=0;i<ntg;i++){nodmass_t[i]=0.0;}
  for(long i=0;i<nkg;i++){nodmass_k[i]=0.0;}

  for(int i=0;i<ng;i++){

    for(int j=0;j<T.nloc();j++){
      for(int gi=0;gi<T.ngi();gi++){
        nodmass_t[TNOD]+=dinit[i]*T.value(j,gi)*detJ0_t[GPNT]*T.wgt(gi);
      }
    }

    for(int k=0;k<K.nloc();k++){
      for(int gi=0;gi<K.ngi();gi++){
        nodmass_k[KNOD]+=dinit[i]*K.value(k,gi)*detJ0_k[GPNT]*K.wgt(gi);
      }
    }

  }

// start the Riemann solvers from initial flux states

  Riemann R0(Riemann::exact,l,r),R1(Riemann::exact,l,r),R2(Riemann::exact,l,r),R3(Riemann::exact,l,r);

// set output precision

  cout<<fixed<<setprecision(17);

// time integration

  while(time<ENDTIME+dt){

// calculate a new stable time-step that will impose the CFL limit on each quadrature point

    for(int i=0;i<ng;i++){
      for(int gi=0;gi<T.ngi();gi++){
        double l(l0[i]*detJ_k[GPNT]/detJ0_k[GPNT]);
        dt_cfl.at(GPNT)=COURANT*T.wgt(gi)*(l/sqrt((c[GPNT]*c[GPNT])+2.0*q[GPNT]/d0_k[GPNT]));
      }
    }

// reduce across element and apply a saftey factor

    double dt=DTSFACTOR*(*min_element(dt_cfl.begin(), dt_cfl.end()));
//    dt=DTSTART;cout<<"DT HARDWIRED !! "<<endl;

    cout<<fixed<<setprecision(5)<<"  step "<<step<<" time= "<<time<<" dt= "<<dt;
    cout<<fixed<<setprecision(5)<<" energy (i/k/tot)= "<<ie<<" "<<ke<<" "<<ie+ke<<endl;

// move nodes to their full-step position

    for(long i=0;i<nkg;i++){x1.at(i)=x0[i]+u0[i]*dt;}

// update thermodynamic node positions using a finite element method, these are a subset of kinematic node positions

    for(int i=0;i<ng;i++){
      for(int j=0;j<T.nloc();j++){
        double pos(-1.0+j*2.0/(T.nloc()-1));x3.at(TNOD)=0.0;
        for(int k=0;k<K.nloc();k++){x3.at(TNOD)+=K.value(k,pos)*x1[KNOD];}
      }
    }

// update mesh centroids

    for(int i=0;i<ng;i++){xc.at(i)=CENTROID;}

// update kinetic energy for conservation checks

    ke=0.0;for(int i=0;i<ng;i++){for(int k=0;k<K.nloc();k++){ke+=0.25*(m[i])*u0[KNOD]*u0[KNOD];}}

// evolve the Riemann problems to the end of the time-step on the end of time-step meshes

    start = high_resolution_clock::now();

    vector<double> r0x,rx,rx2,rx3;vempty(r0x);vempty(rx);vempty(rx2);vempty(rx3); // sample point coordinates
    for(long i=0;i<NSAMPLES;i++){r0x.push_back(0.0+(i*(1.0-0.0)/double(NSAMPLES)));} // sample points
//    for(long i=0;i<NSAMPLES;i++){r0x.push_back(x1[0]+(i*(x1[nkg-1]-x1[0])/double(NSAMPLES)));} // sample points
    R0.profile(&r0x,time+dt); // Riemann solution at the sample points along the mesh
    for(int i=0;i<nkg;i++){rx.push_back(x0[i]);}
    R1.profile(&rx,time+dt);
    for(int i=0;i<ng;i++){for(int gi=0;gi<T.ngi();gi++){double xgi(0.0);XGI;rx2.push_back(xgi);}}
    R2.profile(&rx2,time+dt);
    for(int i=0;i<ntg;i++){rx3.push_back(x3[i]);}
    R3.profile(&rx3,time+dt);

    stop = high_resolution_clock::now();

    time_span[5]+=duration_cast<duration<double>>(stop-start);

// update cell volumes at the full-step

    for(int i=0;i<ng;i++){V1.at(i)=DX1;if(V1[i]<VTOL){cout<<"-'ve volume detected in cell "<<i<<endl;exit(1);}}

// update Jacobians

    for(int i=0;i<ng;i++){

      for(int gi=0;gi<K.ngi();gi++){
        detJ_k.at(GPNT)=0.0;
        for(int k=0;k<K.nloc();k++){
          detJ_k.at(GPNT)+=K.dvalue(k,gi)*x1[KNOD];
        }
//        if(detJ_k.at(GPNT)<0.0){cout<<"-'ve determinant of Jk detected in cell "<<i<<endl;exit(1);}
      }

      for(int gi=0;gi<T.ngi();gi++){
        detJ_t.at(GPNT)=0.0;
        for(int j=0;j<T.nloc();j++){
          detJ_t.at(GPNT)+=T.dvalue(j,gi)*x3[TNOD];
        }
//        if(detJ_t.at(GPNT)<0.0){cout<<"-'ve determinant of Jt detected in cell "<<i<<endl;exit(1);}
      }
    }

// update cell density at the full-step

    for(int i=0;i<ng;i++){for(int gi=0;gi<T.ngi();gi++){d1_t[GPNT]=dinit[i]*detJ0_t[GPNT]/detJ_t[GPNT];}}
    for(int i=0;i<ng;i++){for(int gi=0;gi<K.ngi();gi++){d1_k[GPNT]=dinit[i]*detJ0_k[GPNT]/detJ_k[GPNT];}}

// set nodal volumes

    for(long i=0;i<ntg;i++){nodvol_t[i]=0.0;}
    for(long i=0;i<nkg;i++){nodvol_k[i]=0.0;}

    for(int i=0;i<ng;i++){
      for(int j=0;j<T.nloc();j++){
        for(int gi=0;gi<T.ngi();gi++){
          nodvol_t[TNOD]+=T.value(j,gi)*detJ_t[GPNT]*T.wgt(gi);
        }
      }
      for(int k=0;k<K.nloc();k++){
        for(int gi=0;gi<K.ngi();gi++){
          nodvol_k[KNOD]+=K.value(k,gi)*detJ_k[GPNT]*K.wgt(gi);
        }
      }
    }

// debug
//    for(int i=0;i<ng;i++){for(int gi=0;gi<T.ngi();gi++){d1_t.at(GPNT)=R2.density(GPNT);d1_k.at(GPNT)=R2.density(GPNT);}}
// debug

// assemble finite element energy field on discontinuous thermodynamic grid
// for( struct {int i; double j;} v = {0, 3.0}; v.i < 10; v.i++, v.j+=0.1)

    start = high_resolution_clock::now();

    {Matrix A(T.nloc());double b[T.nloc()],x[T.nloc()];
    for(int i=0;i<ng;i++){
      for(int j=0;j<T.nloc();j++){
        b[j]=0.0;for(int k=0;k<K.nloc();k++){b[j]+=F[KNOD][TNOD]*u1[KNOD];}
        for(int k=0;k<T.nloc();k++){
          double nn(0.0); // DG mass matrix
          for(int gi=0;gi<K.ngi();gi++){
            nn+=d1_k[GPNT]*T.value(j,gi)*T.value(k,gi)*detJ_k[GPNT]*K.wgt(gi);
          }
          A.write(j,k,nn);
        }
      }

    stop = high_resolution_clock::now();

    time_span[1]+=duration_cast<duration<double>>(stop-start);

// solve local system

      start = high_resolution_clock::now();

      A.solve(x,b);

      stop = high_resolution_clock::now();

      time_span[3]+=duration_cast<duration<double>>(stop-start);

// advance the solution

      for(int j=0;j<T.nloc();j++){e1.at(TNOD)=max(ECUT,e0[TNOD]-x[j]*dt);}
    }}

// debug
//    for(int i=0;i<ntg;i++){e1.at(i)=R3.energy(i);}
// debug

// update internal energy for conservation checks

    ie=0.0;for(int i=0;i<ng;i++){for(int j=0;j<T.nloc();j++){ie+=e1[TNOD]*0.5*m[i];}}

// update pressure at the full-step at the integration points

    for(int i=0;i<ng;i++){
      for(int gi=0;gi<K.ngi();gi++){
        double egi(0.0);
        for(int j=0;j<T.nloc();j++){
          egi+=T.value(j,gi)*e1[TNOD];
        }
        p.at(GPNT)=P(d1_k[GPNT],egi); // this is the EOS call
      }
    }

// debug
//    for(int i=0;i<ng;i++){for(int gi=0;gi<T.ngi();gi++){p.at(GPNT)=R2.pressure(GPNT);}}
// debug

// update sound speed

    for(int i=0;i<ng;i++){for(int gi=0;gi<K.ngi();gi++){c.at(GPNT)=sqrt(GAMMA*p[GPNT]/d1_k[GPNT]);}}

// bulk q

    for(int i=0;i<ng;i++){
//      double l(DX1);
      for(int gi=0;gi<K.ngi();gi++){
        double l(l0[i]*(detJ_k[GPNT]/detJ0_k[GPNT]));
        double divu((d0_k[GPNT]-d1_k[GPNT])/(d1_k[GPNT]*dt));
        if(divu<0.0){
          q.at(GPNT)=d0_k[GPNT]*l*divu*((cq*l*divu)-cl*c[GPNT]);
//          q.at(GPNT)=d0_t[GPNT]*l*divu*((cq*l*divu)-cl*c[GPNT]); // sometimes smoother (e.g. p1q2,40) ??
        }else{
          q.at(GPNT)=0.0; // turn off q as cell divergence indicates expansion
        }
      }
    }

// assemble force matrix to connect thermodynamic/kinematic spaces, this can be used as rhs of both e/u eqns

    start = high_resolution_clock::now();

    for(long i=0;i<nkg;i++){for(long j=0;j<ntg;j++){F[i][j]=0.0;}}
    for(int i=0;i<ng;i++){
      for(int k=0;k<K.nloc();k++){
        for(int j=0;j<T.nloc();j++){
          for(int gi=0;gi<K.ngi();gi++){
            F[KNOD][TNOD]+=(p[GPNT]+q[GPNT])*K.dvalue(k,gi)*T.value(j,gi)*K.wgt(gi);
          }
        }
      }
    }

    stop = high_resolution_clock::now();

    time_span[0]+=duration_cast<duration<double>>(stop-start);

// assemble acceleration field

    start = high_resolution_clock::now();

    {Matrix A(NROWS);double b[NROWS],x[NROWS];for(long i=0;i<NROWS;i++){b[i]=0.0;x[i]=0.0;}

// insert the boundary terms into start/end addresses of the source

//    {int i(0),k(K.nloc()-1);for(int j=0;j<T.nloc();j++){b[0]+=F[KNOD][TNOD]*1.0;}}
//    {int i(ng-1),k(0);for(int j=0;j<T.nloc();j++){b[NROWS-1]+=F[KNOD][TNOD]*1.0;}}
// debug
    {int i(0);
      for(int iloc=0;iloc<K.nloc();iloc++){
        for(int gi=0;gi<K.ngi();gi++){
          b[ROW]+=(p[GPNT]+q[GPNT])*K.dvalue(iloc,gi)*K.wgt(gi);
        }
      }
    }

    {int i(ng-1);
      for(int iloc=0;iloc<K.nloc();iloc++){
        for(int gi=0;gi<K.ngi();gi++){
          b[ROW]+=(p[GPNT]+q[GPNT])*K.dvalue(iloc,gi)*K.wgt(gi);
        }
      }
    }


// debug

    for(int i=1;i<ng-1;i++){int k(0);
      for(int iloc=0;iloc<K.nloc();iloc++,k++){
        for(int j=0;j<T.nloc();j++){b[ROW]+=F[KNOD][TNOD]*1.0;} // load vector
        for(int jloc=0;jloc<K.nloc();jloc++){
          double nn(0.0); // mass matrix
          for(int gi=0;gi<K.ngi();gi++){
            nn+=d1_k[GPNT]*K.value(iloc,gi)*K.value(jloc,gi)*detJ_k[GPNT]*K.wgt(gi);
          }
          A.add(ROW,COL,nn);
        }
      }

    }

    A.add(0,0,A.read(0,0));
    A.add((ng-2)*(K.nloc()-1),(ng-2)*(K.nloc()-1),A.read((ng-2)*(K.nloc()-1),(ng-2)*(K.nloc()-1)));

    stop = high_resolution_clock::now();

    time_span[2]+=duration_cast<duration<double>>(stop-start);

// solve global system

    start = high_resolution_clock::now();

    A.solve(x,b);

    stop = high_resolution_clock::now();

    time_span[4]+=duration_cast<duration<double>>(stop-start);

// update acceleration field

    for(int i=1;i<ng-1;i++){for(int k=0;k<K.nloc();k++){int iloc(k);u1.at(KNOD)=u0[KNOD]+x[ROW]*dt;}}

    }

// impose boundary constraints on the acceleration field

    for(int j=0;j<K.nloc()-1;j++){u1.at(j)=u1[K.nloc()];}
    for(int j=0;j<K.nloc()-1;j++){int i(ng-1);u1.at(i*(K.nloc()-1)+j+1)=u1[i*(K.nloc()-1)];}

// debug
//    for(long i=0;i<nkg;i++){u1.at(i)=R1.velocity(i);}
// debug

// some output

    start = high_resolution_clock::now();

    f1.open("exact.dat");f2.open("e.dat");f3.open("u.dat");f4.open("dp.dat");f5.open("mesh.dat");
    f1<<fixed<<setprecision(17);f2<<fixed<<setprecision(17);f3<<fixed<<setprecision(17);f4<<fixed<<setprecision(17);
    for(int i=0;i<NSAMPLES;i++){
      f1<<r0x[i]<<" "<<R0.density(i)<<" "<<R0.pressure(i)<<" "<<R0.velocity(i)<<" "<<R0.energy(i)<<endl;
    }
    for(long i=0;i<ntg;i++){f2<<x3[i]<<" "<<e1[i]<<endl;}
    for(long i=0;i<nkg;i++){f3<<x1[i]<<" "<<u1[i]<<endl;}
    for(int i=0;i<ng;i++){for(int gi=0;gi<T.ngi();gi++){double xgi(0.0);XGI;f4<<xgi<<" "<<d1_k[GPNT]<<" "<<p[GPNT]<<endl;}}
//    cout<<"NODE POS: "<<time<<" "<<x1[(ng/2)*(K.nloc()-1)]<<" "<<x1[(ng/2)*(K.nloc()-1)+1]<<endl;
    for(int i=0;i<ng;i++){
      int k;
      k=0;f5<<x1[KNOD]<<" 0.0"<<endl;f5<<x1[KNOD]<<" 9.0"<<endl;f5<<x1[KNOD]<<" 0.0"<<endl;
      k=K.nloc()-1;f5<<x1[KNOD]<<" 0.0"<<endl;f5<<x1[KNOD]<<" 9.0"<<endl;f5<<x1[KNOD]<<" 0.0"<<endl;
    }
    f1.close();f2.close();f3.close();f4.close();f5.close();

    stop = high_resolution_clock::now();

    time_span[6]+=duration_cast<duration<double>>(stop-start);

// advance the time step

    time+=dt;
    step++;

// debug
//    for(long i=0;i<nkg;i++){u1.at(i)=R1.velocity(i);}
// debug

// advance the solution for the new time step

    for(int i=0;i<nkg;i++){u0.at(i)=u1[i];}
    for(int i=0;i<nkg;i++){x0.at(i)=x1[i];}
    for(int i=0;i<ntg;i++){e0.at(i)=e1[i];}
    for(int i=0;i<ntg;i++){x2.at(i)=x3[i];}
    for(int i=0;i<ng;i++){V0.at(i)=V1[i];}
    for(int i=0;i<ng*T.ngi();i++){d0_t.at(i)=d1_t[i];}
    for(int i=0;i<ng*K.ngi();i++){d0_k.at(i)=d1_k[i];}

// debug
//  cout<<"debug stop."<<endl;
//  exit(1);
// debug

  }

// get a timepoint

  auto stop_main = high_resolution_clock::now();

// subtract timepoints to get a measure of runtime

  time_span[10]=duration_cast<duration<double>>(stop_main-start_main);

  cout<<endl<<"  Breakdown of time accumulated in each part of the calculation:"<<endl<<endl;

  cout<<"    Force Calculation                   "<<time_span[0].count()<<"s "<<time_span[0]*100.0/time_span[10]<<"%"<<endl;
  cout<<"    Energy Field Assembly               "<<time_span[1].count()<<"s "<<time_span[1]*100.0/time_span[10]<<"%"<<endl;
  cout<<"    Acceleration Field Assembly         "<<time_span[2].count()<<"s "<<time_span[2]*100.0/time_span[10]<<"%"<<endl;
  cout<<"    Matrix Inversion for Thermodynamics "<<time_span[3].count()<<"s "<<time_span[3]*100.0/time_span[10]<<"%"<<endl;
  cout<<"    Matrix Inversion for Kinematics     "<<time_span[4].count()<<"s "<<time_span[4]*100.0/time_span[10]<<"%"<<endl;
  cout<<"    Riemann Solvers                     "<<time_span[5].count()<<"s "<<time_span[5]*100.0/time_span[10]<<"%"<<endl;
  cout<<"    Output                              "<<time_span[6].count()<<"s "<<time_span[6]*100.0/time_span[10]<<"%"<<endl;
  cout<<"    Main                                "<<time_span[10].count()<<"s "<<time_span[10]*100.0/time_span[10]<<"%"<<endl;

// total that has been timed

  duration<double> total_time;
  start=high_resolution_clock::now();
  total_time=duration_cast<duration<double>>(start-start); // initialise for the summation

  for(int i=0;i<10;i++){total_time+=time_span[i];}

// subtracting from main gives us an estimate of what has not been timed, this is unknown

  duration<double> unknown;
  unknown=time_span[10]-total_time;

  cout<<endl;

  cout<<"    total timed                         "<<total_time.count()<<"s "<<total_time*100.0/time_span[10]<<"%"<<endl;
  cout<<"    unknown                             "<<unknown.count()<<"s "<<unknown*100.0/time_span[10]<<"%"<<endl;

  cout<<endl<<"  Run took "<<time_span[10].count()<<" seconds."<<endl;
  cout<<"  Normal termination."<<endl;

// release heap storage

  for(long i=0;i<nkg;i++){
    delete[] F[i];
  }
  delete[] F;

  return 0;
}

// return pressure given the energy

double P(double d,double e){return (GAMMA-1.0)*d*e;}

// invert the eos to return energy given the pressure

double E(double d,double p){return p/((GAMMA-1.0)*d);}

// empty a vector

void vempty(vector<double>&v){vector<double> e;v.swap(e);return;}

// debug - check nodal quantities, insert this soon after nodmass_k, nodvol_k etc
//
//  {
//  cout<<endl;
//  cout<<"nodal mass check:"<<endl;
//  double nmassk[K.nloc()],nmasst[T.nloc()];
//  double meshmassk(0.0),meshmasst(0.0);
//
//  for(int i=0;i<ng;i++){
//
//    double sum_k(0.0),sum_t(0.0);
//
//    for(int k=0;k<K.nloc();k++){
//      nmassk[k]=0.0;
//      for(int gi=0;gi<K.ngi();gi++){
//        nmassk[k]+=d1_k[GPNT]*K.value(k,gi)*detJ_k[GPNT]*K.wgt(gi);
//      }
//    }
//    for(int j=0;j<T.nloc();j++){
//      nmasst[j]=0.0;
//      for(int gi=0;gi<T.ngi();gi++){
//        nmasst[j]+=d1_t[GPNT]*T.value(j,gi)*detJ_t[GPNT]*T.wgt(gi);
//      }
//    }
//
//    if(i<10){
//      cout<<i<<"  m[i]= "<<m[i]<<" K: ";
//    }else{
//      cout<<i<<" m[i]= "<<m[i]<<" K: ";
//    }
//
//    for(int k=0;k<K.nloc();k++){sum_k+=nodmass_k[KNOD];cout<<nodmass_k[KNOD]<<" ";}
////    for(int k=0;k<K.nloc();k++){sum_k+=nmassk[k];cout<<nmassk[k]<<" ";}
//    cout<<" T: ";
//
//    for(int j=0;j<T.nloc();j++){sum_t+=nodmass_t[TNOD];cout<<nodmass_t[TNOD]<<" ";}
////    for(int j=0;j<T.nloc();j++){sum_t+=nmasst[j];cout<<nmasst[j]<<" ";}
//    cout<<" sum_k= "<<sum_k<<" sum_t= "<<sum_t<<endl;
//
//  }
//
//  for(long i=0;i<ntg;i++){meshmasst+=nodmass_t[i];}
//  for(long i=0;i<nkg;i++){meshmassk+=nodmass_k[i];}
//
//  cout<<"meshmassk= "<<meshmassk<<" meshmasst= "<<meshmasst<<endl;
//
//  cout<<endl;
//  cout<<"nodal volume check:"<<endl;
//  double nvolk[K.nloc()],nvolt[T.nloc()];
//  double meshvolk(0.0),meshvolt(0.0);
//
//  for(int i=0;i<ng;i++){
//
//    double sum_k(0.0),sum_t(0.0);
//
//    for(int k=0;k<K.nloc();k++){
//      nvolk[k]=0.0;
//      for(int gi=0;gi<K.ngi();gi++){
//        nvolk[k]+=K.value(k,gi)*detJ_k[GPNT]*K.wgt(gi);
//      }
//    }
//    for(int j=0;j<T.nloc();j++){
//      nvolt[j]=0.0;
//      for(int gi=0;gi<T.ngi();gi++){
//        nvolt[j]+=T.value(j,gi)*detJ_t[GPNT]*T.wgt(gi);
//      }
//    }
//
//    if(i<10){
//      cout<<i<<"  V1[i]= "<<V1[i]<<" K: ";
//    }else{
//      cout<<i<<" V1[i]= "<<V1[i]<<" K: ";
//    }
//
//    for(int k=0;k<K.nloc();k++){sum_k+=nodvol_k[KNOD];cout<<nodvol_k[KNOD]<<" ";}
////    for(int k=0;k<K.nloc();k++){sum_k+=nvolk[k];cout<<nvolk[k]<<" ";}
//    cout<<" T: ";
//
//    for(int j=0;j<T.nloc();j++){sum_t+=nodvol_t[TNOD];cout<<nodvol_t[TNOD]<<" ";}
////    for(int j=0;j<T.nloc();j++){sum_t+=nvolt[j];cout<<nvolt[j]<<" ";}
//    cout<<" sum_k= "<<sum_k<<" sum_t= "<<sum_t<<endl;
//
//  }
//
//  for(long i=0;i<ntg;i++){meshvolt+=nodvol_t[i];}
//  for(long i=0;i<nkg;i++){meshvolk+=nodvol_k[i];}
//
//  cout<<"meshvolk= "<<meshvolk<<" meshvolt= "<<meshvolt<<endl;
//
//  cout<<endl;
//
//  cout<<"x1[] nodvol_k"<<endl;
//  for(long i=0;i<nkg;i++){cout<<x1[i]<<" "<<nodvol_k[i]<<endl;}
//  cout<<endl;
//
//  cout<<"x3[] nodvol_t"<<endl;
//  for(long i=0;i<ntg;i++){cout<<x3[i]<<" "<<nodvol_t[i]<<endl;}
//  cout<<endl;
//
//  cout<<"x1[] nodmass_k"<<endl;
//  for(long i=0;i<nkg;i++){cout<<x1[i]<<" "<<nodmass_k[i]<<endl;}
//  cout<<endl;
//
//  cout<<"x3[] nodmass_t"<<endl;
//  for(long i=0;i<ntg;i++){cout<<x3[i]<<" "<<nodmass_t[i]<<endl;}
//  cout<<endl;
//
//  }

// debug - check nodal quantities
