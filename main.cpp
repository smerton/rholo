// High order finite element variant of RhoLo (Really High Order Lagrangian Operator - RhoLo)
// RhoLo is an ultra simple 1-D high order finite element hydrodynamics test code.
// This variant solves the Euler equations in their non-conservative form in the
// fluid frame (the Lagrangian frame) using a high order finite element method (discontinuous
// thermodynamic variables d,rho,e with node centred kinematic variables u,a) and bulk
// viscosity q to increase entropy across element boundaries, initial implementation is only
// first order in time
//
// Author S. R. Merton
//
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
// for graphics: convert -density 300 filename.png filename.pdf
//
//

#define DTSTART 0.0005          // insert a macro for the first time step
#define ENDTIME 0.20            // insert a macro for the end time
#define ECUT 1.0e-8             // cut-off on the energy field
#define NSAMPLES 1000           // number of sample points for the exact solution
#define VISFREQ 10000           // frequency of the graphics dumps
#define VD vector<double>       // vector of doubles
#define VTOL 1.0e-10            // threshold for volume errors
#define COURANT 0.333           // Courant number for CFL condition
#define DTSFACTOR 0.75          // safety factor on time step control
#define KNOD i*(K.nloc()-1)+k   // global node number on kinematic mesh
#define TNOD i*T.nloc()+j       // global node number on thermodynamic mesh
#define GPNT i*T.ngi()+gi       // global address of Gauss point gi in element i
#define DX0 (x0[i*(K.nloc()-1)+K.nloc()-1]-x0[i*(K.nloc()-1)])          // cell width at start of step
#define DX1 (x1[i*(K.nloc()-1)+K.nloc()-1]-x1[i*(K.nloc()-1)])          // cell width at end of step
#define CENTROID 0.5*(x1[i*(K.nloc()-1)+K.nloc()-1]+x1[i*(K.nloc()-1)]) // cell centroid
#define ROW i*(K.nloc()-1)+iloc                                         // row address in global matrix
#define COL i*(K.nloc()-1)+jloc                                         // column address in global matrix
#define XGI for(int j=0;j<T.nloc();j++){xgi+=T.value(j,gi)*x3[TNOD];}   // coordinates of integration point
#define NROWS n*(K.nloc()-1)+1                                          // number of rows in the global matrix
#define NCOLS n*(K.nloc()-1)+1                                          // number of columns in the global matrix (= no. rows)

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
#include "timer.h"
#include "eos.h"

// function signatures

void vempty(vector<double>&v);        // signature for emptying a vector
int iaddr(int iel,int iel1,int iel2); // signature for element address function
void get_exact(double t);             // exact solutions at time t
vector<double> r0x,rx,rx2,rx3;        // sample point coordinates for Riemann solver

using namespace std;
using namespace chrono;

int main(){

  cout<<"main(): Starting up main loop..."<<endl;

// global data

  ofstream f1,f2,f3,f4,f5,f6;                           // files for output
  ifstream f7;                                          // files for input
  Shape K(2,9),T(1,9);                                  // p_n,q_n-1 shape functions
  int const n(100);                                     // no. ncells
  int long nk(n*(K.nloc()-1)+1);                        // no. kinematic nodes
  int long nt(n*T.nloc());                              // no. thermodynamic nodes
  double const cl(1.0),cq(1.0);                         // linear & quadratic coefficients for bulk viscosity
  vector<double> dinit(n);                              // initial density field inside an element
  vector<double> d0_t(n*T.ngi()),d1_t(n*T.ngi());       // density at each Gauss point in each element
  vector<double> d0_k(n*K.ngi()),d1_k(n*K.ngi());       // density at each Gauss point in each element
  vector<double> V0(n),V1(n),m(n),xc(n);                // volume, mass & centroid
  vector<double> e0(nt),e1(nt);                         // discontinuous FE energy field
  vector<double> c(n*T.ngi()),p(n*T.ngi());             // element sound speed & pressure at each Gauss point
  vector<double> q(n*T.ngi());                          // bulk viscosity at each Gauss point
  vector<double> u0(nk),u1(nk);                         // node velocity
  vector<double> x0(nk),x1(nk),x2(nt),x3(nt);           // node coordinates
  vector<double> dt_cfl(n*T.ngi());                     // element time-step at each Gauss point
  vector<double> detJ0_k(n*K.ngi()),detJ_k(n*K.ngi());  // determinant of the Jacobian
  vector<double> l0(n);                                 // initial length scale
  vector<vector<long> > nzloc;                          // locations of non-zeroes in the inverse mass matrix
  double ke(0.0),ie(0.0);                               // kinetic and internal energy for conservation checks
  double time(0.0),dt(DTSTART);                         // start time and time step
  int step(0);                                          // step number
  Timer timers(20);                                     // time acccumulated in different parts of code
  double l[3]={1.0,0.0,1.0},r[3]={0.125,0.0,0.1};       // left/right flux states for the problem: Sod
//  double l[3]={1.0,-2.0,0.4},r[3]={1.0,2.0,0.4};      // left/right flux states for the problem: 123 (R2R)
//  double l[3]={1.0,0.0,1000.0},r[3]={1.0,0.0,0.01};   // left/right flux states for the problem: blast wave

//  heap storage

  double** F=new double*[nk];                           // force matrix (located on  heap)
  for(long i=0;i<nk;i++){
    F[i]=new double[nt];
  }

  Matrix KMASS(NROWS);                                  // mass matrix for kinematic/thermodynamic fields
  Matrix KMASSI(NROWS);                                 // inverse mass matrix for kinematics

// initialise the problem

  double dx(1.0/n);x0.at(0)=0.0;x1.at(0)=x0[0];x2.at(0)=0.0;x3.at(0)=x2[0];
  for(long i=0;i<nk;i++){x0.at(i)=x0[0]+i*dx/(K.nloc()-1);x1.at(i)=x0[i];}
  for(long i=0;i<nt;i++){x2.at(i)=x2[0]+i*dx/(T.nloc()-1);x3.at(i)=x2[i];}
  for(int i=0;i<n;i++){xc.at(i)=CENTROID;}
  for(int i=0;i<n;i++){for(int gi=0;gi<K.ngi();gi++){p.at(GPNT)=(xc[i]<=0.5)?l[2]:r[2];}}
  for(int i=0;i<n;i++){for(int gi=0;gi<T.ngi();gi++){d0_t.at(GPNT)=(xc[i]<=0.5)?l[0]:r[0];}}
  for(int i=0;i<n;i++){for(int gi=0;gi<K.ngi();gi++){d0_k.at(GPNT)=(xc[i]<=0.5)?l[0]:r[0];}}
  for(int i=0;i<n;i++){for(int gi=0;gi<T.ngi();gi++){d1_t.at(GPNT)=(xc[i]<=0.5)?l[0]:r[0];}}
  for(int i=0;i<n;i++){for(int gi=0;gi<K.ngi();gi++){d1_k.at(GPNT)=(xc[i]<=0.5)?l[0]:r[0];}}
  for(int i=0;i<n;i++){for(int j=0;j<T.nloc();j++){e0.at(TNOD)=(xc[i]<=0.5)?E(l[0],l[2]):E(r[0],r[2]);}}
  for(int i=0;i<n;i++){for(int j=0;j<T.nloc();j++){e1.at(TNOD)=(xc[i]<=0.5)?E(l[0],l[2]):E(r[0],r[2]);}}
  for(int i=0;i<n;i++){V0.at(i)=DX0;V1.at(i)=DX1;}
  for(int i=0;i<n;i++){dinit.at(i)=(xc[i]<=0.5)?l[0]:r[0];}
  for(int i=0;i<n;i++){m.at(i)=(xc[i]<=0.5)?l[0]*V0[i]:r[0]*V0[i];}
  for(long i=0;i<nk;i++){u0.at(i)=(x0[i]<0.5)?l[1]:(x0[i]==0.5)?0.0:r[1];u1.at(i)=u0[i];}
  for(int i=0;i<n;i++){for(int gi=0;gi<T.ngi();gi++){q.at(GPNT)=0.0;}}
  for(int i=0;i<n;i++){for(int gi=0;gi<T.ngi();gi++){c.at(GPNT)=C(p[GPNT],d0_k[GPNT]);}}
  for(long i=0;i<nk;i++){for(long j=0;j<nt;j++){F[i][j]=0.0;}}
  for(int i=0;i<n;i++){l0.at(i)=DX0/(K.nloc()-1);} // initial nodal displacement

// initialise the high res timers

  timers.Init();

// start a timer for main

  timers.Start(0);

// check integration rules on thermodynamic and kinematic stencils look consistent, they need to match

  if(T.ngi()!=K.ngi()){
    cout<<"No. of integration points for thermodynamics ("<<T.ngi()<<") doesn't match the no. for kinematics ("<<K.ngi()<<")"<<endl;
    exit(1);
  }

// set thermodynamic node positions at time-0

  for(int i=0;i<n;i++){
    for(int j=0;j<T.nloc();j++){
      double pos(-1.0+j*2.0/(T.nloc()-1));x2.at(TNOD)=0.0;
      for(int k=0;k<K.nloc();k++){x2.at(TNOD)+=K.value(k,pos)*x0[KNOD];}
    }
  }

// set time-0 Jacobians

  for(int i=0;i<n;i++){
    for(int gi=0;gi<K.ngi();gi++){
      detJ0_k.at(GPNT)=0.0;
      for(int k=0;k<K.nloc();k++){
        detJ0_k.at(GPNT)+=K.dvalue(k,gi)*x0[KNOD];
      }
      if(detJ0_k.at(GPNT)<0.0){cout<<"-'ve determinant of J0_k detected in cell "<<i<<endl;exit(1);}
    }
  }

// assemble mass matrix for acceleration field

  for(int iel=-1;iel<=n;iel++){
    int i(iaddr(iel,0,n-1)); // get element address
    for(int iloc=0;iloc<K.nloc();iloc++){
      for(int jloc=0;jloc<K.nloc();jloc++){
        double nn(0.0); // mass matrix
        for(int gi=0;gi<K.ngi();gi++){
          nn+=d0_k[GPNT]*K.value(iloc,gi)*K.value(jloc,gi)*detJ0_k[GPNT]*K.wgt(gi);
        }
        KMASS.add(ROW,COL,nn);
      }
    }
  }

// invert the mass matrices

  cout<<"Inverting mass matrix for the acceleration field..."<<endl;

  timers.Start(1);

  KMASSI.inverse2(&KMASS); // lapack

  timers.Stop(1);

  cout<<"Done."<<endl;

// non-zeroes in the inverse mass matrix

  for(long i=0;i<NROWS;i++){
    vector<long> nzrow;
    for(long j=0;j<NCOLS;j++){
      if(abs(KMASSI.read(i,j))>1.0e-12){
        nzrow.push_back(j); // append location of the non-zero
      }
    }
    nzloc.push_back(nzrow); // append all the non-zeroes on this row
  }

// start the Riemann solvers from initial flux states

  Riemann R0(Riemann::exact,l,r),R1(Riemann::exact,l,r),R2(Riemann::exact,l,r),R3(Riemann::exact,l,r);

// set output precision

  cout<<fixed<<setprecision(17);

// time integration

  while(time<=ENDTIME){

// calculate a new stable time step that will impose the CFL limit on each quadrature point

    for(int i=0;i<n;i++){
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

    for(long i=0;i<nk;i++){x1.at(i)=x0[i]+u0[i]*dt;}

// update mesh centroids

    for(int i=0;i<n;i++){xc.at(i)=CENTROID;}

// update kinetic energy for conservation checks

    ke=0.0;for(int i=0;i<n;i++){for(int k=0;k<K.nloc();k++){ke+=0.25*(m[i])*u0[KNOD]*u0[KNOD];}}

// exact solutions at end of current time step

    timers.Start(2);

//    get_exact(time+dt);

    timers.Stop(2);

// update cell volumes at the full-step

    for(int i=0;i<n;i++){V1.at(i)=DX1;if(V1[i]<VTOL){cout<<"-'ve volume detected in cell "<<i<<endl;exit(1);}}

// update Jacobians

    for(int i=0;i<n;i++){

      for(int gi=0;gi<K.ngi();gi++){
        detJ_k.at(GPNT)=0.0;
        for(int k=0;k<K.nloc();k++){
          detJ_k.at(GPNT)+=K.dvalue(k,gi)*x1[KNOD];
        }
//        if(detJ_k.at(GPNT)<0.0){cout<<"-'ve determinant of Jk detected in cell "<<i<<endl;exit(1);}
      }

    }

// update cell density at the full-step

    for(int i=0;i<n;i++){for(int gi=0;gi<K.ngi();gi++){d1_k[GPNT]=dinit[i]*detJ0_k[GPNT]/detJ_k[GPNT];}}

// debug
//    for(int i=0;i<n;i++){for(int gi=0;gi<K.ngi();gi++){d1_k.at(GPNT)=R2.density(GPNT);}}
// debug

// assemble finite element energy field on the discontinuous thermodynamic grid
// for( struct {int i; double j;} v = {0, 3.0}; v.i < 10; v.i++, v.j+=0.1)

    timers.Start(3);

    {Matrix A(T.nloc());double b[T.nloc()],x[T.nloc()];
    for(int i=0;i<n;i++){
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

// solve local system

      A.solve(x,b);

// advance the solution

      for(int j=0;j<T.nloc();j++){e1.at(TNOD)=max(ECUT,e0[TNOD]-x[j]*dt);}
    }}

    timers.Stop(3);

// debug
//    for(int i=0;i<nt;i++){e1.at(i)=R3.energy(i);}
// debug

// update internal energy for conservation checks

    ie=0.0;for(int i=0;i<n;i++){for(int j=0;j<T.nloc();j++){ie+=e1[TNOD]*0.5*m[i];}}

// update pressure at the full-step at the integration points

    for(int i=0;i<n;i++){
      for(int gi=0;gi<K.ngi();gi++){
        double egi(0.0);
        for(int j=0;j<T.nloc();j++){
          egi+=T.value(j,gi)*e1[TNOD];
        }
        p.at(GPNT)=P(d1_k[GPNT],egi); // this is the EOS call
      }
    }

// debug
//    for(int i=0;i<n;i++){for(int gi=0;gi<T.ngi();gi++){p.at(GPNT)=R2.pressure(GPNT);}}
// debug

// update sound speed

    for(int i=0;i<n;i++){for(int gi=0;gi<K.ngi();gi++){c.at(GPNT)=C(p[GPNT],d1_k[GPNT]);}}

// bulk q

    for(int i=0;i<n;i++){
      for(int gi=0;gi<K.ngi();gi++){
        double l(l0[i]*(detJ_k[GPNT]/detJ0_k[GPNT]));
        double divu((d0_k[GPNT]-d1_k[GPNT])/(d1_k[GPNT]*dt));
        if(divu<0.0){
          q.at(GPNT)=d0_k[GPNT]*l*divu*((cq*l*divu)-cl*c[GPNT]);
        }else{
          q.at(GPNT)=0.0; // compression switch to turn off q as cell divergence indicates expansion
        }
      }
    }

// assemble force matrix to connect thermodynamic/kinematic spaces, this can be used as rhs of both e/u eqns

    timers.Start(4);

    for(long i=0;i<nk;i++){for(long j=0;j<nt;j++){F[i][j]=0.0;}}
    for(int i=0;i<n;i++){
      for(int k=0;k<K.nloc();k++){
        for(int j=0;j<T.nloc();j++){
          for(int gi=0;gi<K.ngi();gi++){
            F[KNOD][TNOD]+=(p[GPNT]+q[GPNT])*K.dvalue(k,gi)*T.value(j,gi)*K.wgt(gi);
          }
        }
      }
    }

    timers.Stop(4);

// assemble acceleration field

    timers.Start(5);

    {double b[NROWS],x[NROWS];for(long i=0;i<NROWS;i++){b[i]=0.0;x[i]=0.0;}

    for(int i=0;i<n;i++){
      for(int k=0;k<K.nloc();k++){
        for(int j=0;j<T.nloc();j++){
          b[KNOD]+=F[KNOD][TNOD];
        }
      }
    }

// sourcing on boundary of acceleration field

    for(int j=0;j<T.nloc();j++){int i(0);b[0]+=F[K.nloc()-1][TNOD];}
    for(int j=0;j<T.nloc();j++){int i(n-1);b[nk-1]+=F[nk-K.nloc()][TNOD];}

// solve global system

    for(long i=0;i<NROWS;i++){
      x[i]=0.0;
      for(long j=0;j<NCOLS;j++){
        x[i]+=KMASSI.read(i,j)*b[j];
      }
    }

// new CSR version - avoids zeroes in the inverse mass matrix for a faster matvec

//    for(long irow=0;irow<NROWS;irow++){
//      x[irow]=0.0;
//      for(int j=0;j<nzloc[irow].size();j++){
//        x[irow]+=KMASSI.read(irow,nzloc[irow][j])*b[nzloc[irow][j]];
//      }
//    }

// advance solution

    for(int i=0;i<n;i++){for(int k=0;k<K.nloc();k++){u1.at(KNOD)=u0[KNOD]+x[KNOD]*dt;}}

    }

    timers.Stop(5);

// debug
//    for(long i=0;i<nk;i++){u1.at(i)=R1.velocity(i);}
// debug

// advance the time step

    time+=dt;
    step++;

// debug
//    for(long i=0;i<nk;i++){u1.at(i)=R1.velocity(i);}
// debug

// advance the solution for the new time step

    for(int i=0;i<nk;i++){u0.at(i)=u1[i];}
    for(int i=0;i<nk;i++){x0.at(i)=x1[i];}
    for(int i=0;i<nt;i++){e0.at(i)=e1[i];}
    for(int i=0;i<nt;i++){x2.at(i)=x3[i];}
    for(int i=0;i<n;i++){V0.at(i)=V1[i];}
    for(int i=0;i<n*T.ngi();i++){d0_t.at(i)=d1_t[i];}
    for(int i=0;i<n*K.ngi();i++){d0_k.at(i)=d1_k[i];}

// debug
//  cout<<"debug stop."<<endl;
//  output(); // might want this ??
//  exit(1);
// debug

  }

// some output

    timers.Start(6);

    f1.open("exact.dat");f2.open("e.dat");f3.open("u.dat");f4.open("dp.dat");f5.open("mesh.dat");f6.open("r.dat");
    f1<<fixed<<setprecision(17);f2<<fixed<<setprecision(17);f3<<fixed<<setprecision(17);f4<<fixed<<setprecision(17);
    vempty(r0x);
    for(long i=0;i<NSAMPLES;i++){r0x.push_back(0.0+(i*(1.0-0.0)/double(NSAMPLES)));} // sample points
//    for(long i=0;i<NSAMPLES;i++){r0x.push_back(x1[0]+(i*(x1[nkg-1]-x1[0])/double(NSAMPLES)));} // sample points
    R0.profile(&r0x,time); // Riemann solution at the sample points along the mesh
    for(int i=0;i<NSAMPLES;i++){
      f1<<r0x[i]<<" "<<R0.density(i)<<" "<<R0.pressure(i)<<" "<<R0.velocity(i)<<" "<<R0.energy(i)<<endl;
    }
    for(int i=1;i<NSAMPLES;i++){if(R0.region(i)!=R0.region(i-1)){f6<<r0x[i]<<" -20000.0"<<endl<<r0x[i]<<" 20000.0"<<endl<<r0x[i]<<" -20000.0"<<endl;}} // sample region
    for(int i=0;i<n;i++){
      for(int j=0;j<T.nloc();j++){
        double pos(-1.0+j*2.0/(T.nloc()-1));x3.at(TNOD)=0.0;
        for(int k=0;k<K.nloc();k++){x3.at(TNOD)+=K.value(k,pos)*x1[KNOD];}
      }
    }

    for(long i=0;i<nt;i++){f2<<x3[i]<<" "<<e1[i]<<endl;}
    for(long i=0;i<nk;i++){f3<<x1[i]<<" "<<u1[i]<<endl;}
    for(int i=0;i<n;i++){for(int gi=0;gi<T.ngi();gi++){double xgi(0.0);XGI;f4<<xgi<<" "<<d1_k[GPNT]<<" "<<p[GPNT]<<endl;}}
//    cout<<"NODE POS: "<<time<<" "<<x1[(ng/2)*(K.nloc()-1)]<<" "<<x1[(ng/2)*(K.nloc()-1)+1]<<endl;
    for(int i=0;i<n;i++){
      int k;
      k=0;f5<<x1[KNOD]<<" -10.0"<<endl;f5<<x1[KNOD]<<" 10.0"<<endl;f5<<x1[KNOD]<<" -10.0"<<endl;
      k=K.nloc()-1;f5<<x1[KNOD]<<" -10.0"<<endl;f5<<x1[KNOD]<<" 10.0"<<endl;f5<<x1[KNOD]<<" -10.0"<<endl;
    }
    f1.close();f2.close();f3.close();f4.close();f5.close();f6.close();

    timers.Stop(6);

// stop timer for main

  timers.Stop(0);

// estimate convergence error in the L1/L2 norms using a Riemann solution as the exact solution

  vector<double> rx;vempty(rx);  // sample points
  double xstart(0.0),xstop(1.0); // sample range - should be whole mesh though bc.s may artificially reduce convergence
  for(long i=0;i<n;i++){for(int k=0;k<K.nloc();k++){if(x1[KNOD]>=xstart&&x1[KNOD]<=xstop){rx.push_back(x1[KNOD]);}}}
  R0.profile(&rx,ENDTIME); // evolve a Riemann problem on the sample range to final time level

// this is just to get the exact solution from file (a high res calculation):
//  vempty(rx);for(long i=0;i<nkg;i++){rx.push_back(x1[i]);};R0.profile(&rx,ENDTIME);
//  f6.open("runs/p1q1h2000/u.dat");
//  vector<double> xf,uf;
//  while(!f6.eof()){double a,b;f6>>a>>b;xf.push_back(a);uf.push_back(b);}
//  f7.close();

  double l1(0.0),l2(0.0),l1r(0.0),l2r(0.0),l1d(0.0),l2d(0.0),l1n(0.0),l2n(0.0);
  int ii(0);

// numerator and denominator for each norm

  for(int i=0;i<n;i++){
    for(int k=0;k<K.nloc();k++){
      if(x1[KNOD]>=xstart&&x1[KNOD]<=xstop){
        double err(abs(R0.velocity(ii)-u1[KNOD])); // absolute error

        double kvol(0.0); // volume of node k
        for(int gi=0;gi<K.ngi();gi++){
          kvol+=K.value(k,gi)*detJ_k[GPNT]*K.wgt(gi);
        }

        l1d+=abs(R0.velocity(ii)); // l1 denominator
        l2d+=R0.velocity(ii)*R0.velocity(ii); // l2 denominator
        l1n+=err*kvol; // l1 numerator
        l2n+=err*err*kvol; // l2 numerator
        ii++;
      }
    }
  }

// construct norms and relative errors

  l1=l1n/rx.size(); // L1 error norm
  l1r=l1n/l1d; // L1 relative error norm
  l2=sqrt(l2n/rx.size()); // L2 error norm
  l2r=sqrt(l2n/l2d); // L2 relative error norm

  cout<<endl<<fixed<<setprecision(10)<<"  Error Estimates (grid spacing h= "<<1.0/n<<"):"<<endl;

  cout<<"  L1 norm= "<<l1<<" (relative error= "<<l1r<<")"<<endl;
  cout<<"  L2 norm= "<<l2<<" (relative error= "<<l2r<<")"<<endl;
  cout<<"  No. points sampled= "<<rx.size()<<endl;
  cout<<"  Range sampled= "<<xstart<<","<<xstop<<endl;

// output high resolution timings

  cout<<endl<<"  Breakdown of time accumulated in each part of the calculation:"<<endl<<endl;

  cout<<"    Matrix Inverter                     "<<timers.Span(1)<<"s "<<timers.Span(1)*100.0/timers.Span(0)<<"%"<<endl;
  cout<<"    Riemann Solvers                     "<<timers.Span(2)<<"s "<<timers.Span(2)*100.0/timers.Span(0)<<"%"<<endl;
  cout<<"    Energy Field Assembly               "<<timers.Span(3)<<"s "<<timers.Span(3)*100.0/timers.Span(0)<<"%"<<endl;
  cout<<"    Force Calculation                   "<<timers.Span(4)<<"s "<<timers.Span(4)*100.0/timers.Span(0)<<"%"<<endl;
  cout<<"    Acceleration Field Assembly         "<<timers.Span(5)<<"s "<<timers.Span(5)*100.0/timers.Span(0)<<"%"<<endl;
  cout<<"    Output                              "<<timers.Span(6)<<"s "<<timers.Span(6)*100.0/timers.Span(0)<<"%"<<endl;
  cout<<"    Main                                "<<timers.Span(0)<<"s "<<timers.Span(0)*100.0/timers.Span(0)<<"%"<<endl;

// total timed excluding main

  double total(timers.Total()-timers.Span(0));

// subtracting total timed from main gives us an estimate of what has not been timed, this is unknown

  double unknown(timers.Span(0)-total);

  cout<<endl;

  cout<<"    total timed                         "<<total<<"s "<<total*100.0/timers.Span(0)<<"%"<<endl;
  cout<<"    unknown                             "<<unknown<<"s "<<unknown*100.0/timers.Span(0)<<"%"<<endl;

  cout<<endl<<"  Run took "<<timers.Span(0)<<" seconds."<<endl;
  cout<<"  Normal termination."<<endl;

// release heap storage

  for(long i=0;i<nk;i++){
    delete[] F[i];
    F[i]=NULL;
  }
  delete[] F;
  F=NULL;

  return 0;
}

// empty a vector

void vempty(vector<double>&v){vector<double> e;v.swap(e);return;}

// element address

int iaddr(int iel,int iel1,int iel2){

  int i(iel);

  if(iel<iel1){

    i=iel1;

  }else if(iel>iel2){

    i=iel2;

  }

  return i;

}

// evolve the Riemann problems to the given time to find the exact solutions

  void get_exact(double t){

// evolve the Riemann problems to the given time to find the exact solutions at that time

// this was set up as a debugging tool:
// R0 is the exact solution at NSAMPLES sample points for comparing against solutions from the code
// R1.velocity used as exact solution to overide finite element velocity field
// R2.pressure used as exact solution to overide Gauss point pressures
// R2.density used as exact solution to overide Gauss point densities
// R3.energy used as exact solution to overide finite element energy field

// reset sample point coordinates for various meshes

//    vempty(r0x);vempty(rx);vempty(rx2);vempty(rx3);

//    for(long i=0;i<NSAMPLES;i++){r0x.push_back(0.0+(i*(1.0-0.0)/double(NSAMPLES)));}
////    for(long i=0;i<NSAMPLES;i++){r0x.push_back(x1[0]+(i*(x1[nkg-1]-x1[0])/double(NSAMPLES)));}
//    R0.profile(&r0x,t); // Riemann solution at the sample points along the mesh
//    for(int i=0;i<nkg;i++){rx.push_back(x0[i]);}
//    R1.profile(&rx,t);
//    for(int i=0;i<ng;i++){for(int gi=0;gi<T.ngi();gi++){double xgi(0.0);XGI;rx2.push_back(xgi);}}
//    R2.profile(&rx2,t);

// thermodynamic node coordinates using a finite element interpolation of kinematic node coordinates
//
//    for(int i=0;i<ng;i++){
//      for(int j=0;j<T.nloc();j++){
//        double pos(-1.0+j*2.0/(T.nloc()-1));x3.at(TNOD)=0.0;
//        for(int k=0;k<K.nloc();k++){x3.at(TNOD)+=K.value(k,pos)*x1[KNOD];}
//      }
//    }
//
//    for(int i=0;i<ntg;i++){rx3.push_back(x3[i]);}
//    R3.profile(&rx3,t);
//
    return;

  }

// debug - nodal masses, insert this above time step loop

// set nodal masses - these should not change with time

//  vector<double> nodmass_t(ntg),nodmass_k(nkg);         // nodal mass
//
//  for(long i=0;i<ntg;i++){nodmass_t[i]=0.0;}
//  for(long i=0;i<nkg;i++){nodmass_k[i]=0.0;}
//
//
//    for(int j=0;j<T.nloc();j++){
//      for(int gi=0;gi<T.ngi();gi++){
//        nodmass_t[TNOD]+=dinit[i]*T.value(j,gi)*detJ0_t[GPNT]*T.wgt(gi);
//      }
//    }
//
//    for(int k=0;k<K.nloc();k++){
//      for(int gi=0;gi<K.ngi();gi++){
//        nodmass_k[KNOD]+=dinit[i]*K.value(k,gi)*detJ0_k[GPNT]*K.wgt(gi);
//      }
//    }
//
//  }
//



// debug - nodal volumes, insert this INSIDE time step loop

// set nodal volumes

//  vector<double> nodvol_t(ntg),nodvol_k(nkg);           // nodal volume

//    for(long i=0;i<ntg;i++){nodvol_t[i]=0.0;}
//    for(long i=0;i<nkg;i++){nodvol_k[i]=0.0;}
//
//    for(int i=0;i<ng;i++){
//      for(int j=0;j<T.nloc();j++){
//        for(int gi=0;gi<T.ngi();gi++){
//          nodvol_t[TNOD]+=T.value(j,gi)*detJ_t[GPNT]*T.wgt(gi);
//        }
//      }
//      for(int k=0;k<K.nloc();k++){
//        for(int gi=0;gi<K.ngi();gi++){
//          nodvol_k[KNOD]+=K.value(k,gi)*detJ_k[GPNT]*K.wgt(gi);
//        }
//      }
//    }

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
