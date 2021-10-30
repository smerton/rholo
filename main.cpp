// Finite element variant of RhoLo (Really High Order Lagrangian Operator - RhoLo)
// RhoLo is an ultra simple finite element (DG) hydrodynamics test code
// This finite element variant solves the Euler equations in their non-conservative form in the fluid frame (the Lagrangian frame)
// using a mixed continuous finite element method (cell-centred thermodynamic variable d,rho,e with node 
// centred kinematic variables u,a) and bulk viscosity q to increase entropy across element boundaries, initial 
// implementation is only first order in time
//
// Author S. R. Merton
//
// how to set up push and pull from multiple remotes (e.g. home NAS and github):
//
// git remote add NAS ssh://smerton@192.168.1.79/shares/git/rholo.git
// git remote set-url --add NAS git@github.com:smerton/rholo.git
// git remote add github git@github.com:smerton/rholo.git
// git push -vu NAS finite-element
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

#define DTSTART 0.0005    // insert a macro for the first time step
#define ENDTIME 0.20      // insert a macro for the end time
#define ECUT 1.0e-8       // cut-off on the energy field
#define NSAMPLES 1000     // number of sample points for the exact solution
//#define VISFREQ 200     // frequency of the graphics dump steps
#define VISFREQ 0.05      // frequency of the graphics dump times
#define VD vector<double> // vector of doubles
#define VVD vector<VD>    // vector of VD
#define VI vector<int>    // vector of ints
#define VTOL 1.0e-10      // threshold for volume errors
#define COURANT 0.333     // Courant number for CFL condition
#define DTSFACTOR 0.5     // safety factor on time-step control

#include <iostream>
#include <vector>
#include <iomanip>
#include <cmath>
#include <fstream>
#include "riemann.h"
#include "matrix.h"
#include "shape.h"
#include <bits/stdc++.h>
#include "eos.h"
#include "mesh.h"
#include <ctime>

// function signatures

string date();                                                                         // function to return the date and time
void vempty(vector<double>&v);                                                         // signature for emptying a vector
void silo(VD*x0,VD*d0,VD*p,VD*e0,VD*q,VD*c,VD*u0,VI*mat,int step,double time,Shape*S); // silo graphics output
void initial_data(int const n,int const ndims, int const nmats, VI const &mat,         // output material states
                  VD const &d, VD const &V, VD const &m, VD const &e, VD const &p, 
                  VVD const &x, VVD const &u, int const &s, double const &t);

using namespace std;

int main(){

  cout<<"main(): Starting up main loop..."<<endl;

// global data

  Mesh M("mesh/input.mesh");                              // load a new mesh from file
  Shape S(1);                                             // load a p1 shape function
  ofstream f1,f2,f3,f4,f5;                                // files for output
  int const n(M.NCells()),ndims(M.NDims());               // no. ncells and no. dimensions
  int const nmats(M.NMaterials());                        // number of materials
  double const cl(0.3),cq(1.0);                           // linear & quadratic coefficients for bulk viscosity
  vector<double> d0(n),d1(n),V0(n),V1(n),m(n);            // density, volume & mass
  vector<double> e0(n),e1(n);                             // cell-centred energy field
  vector<double> c(n),p(n),q(n);                          // element sound speed, pressure and bulk viscosity
  vector<vector<double> > u0(ndims),u1(ndims);            // node velocity
  vector<vector<double> > x0(ndims),x1(ndims);            // node coordinates
  vector<double> dt_cfl(n);                               // element time-step
  vector<double> dts(2);                                  // time-step for each condition (0=CFL, 1=graphics hits)
  vector<int> mat(n);                                     // element material numbers
  double ke(0.0),ie(0.0);                                 // kinetic and internal energy for conservation checks
  double time(0.0),dt(DTSTART);                           // start time and time step
  int step(0);                                            // step number

  vector<vector<double> > state={{1.000,0.000,1.000},     // initial flux state in each material for Sod's shock tube 
                                 {0.125,0.000,0.100}};

//  vector<vector<double> > state={{1.000,-2.000,0.400},  // initial flux state in each material for the 123 problem 
//                                 {1.000, 2.000,0.400}};

//  vector<vector<double> > state={{1.000,0.000,1000.0},  // initial flux state in each material for the blast wave
//                                 {1.000,0.000,0.0100}};

  double l[3]={state[0][0],state[0][1],state[0][2]};      // left flux state for input to the Riemann solvers
  double r[3]={state[1][0],state[1][1],state[1][2]};      // right flux state for input to the Riemann solvers

// initialise the problem

  M.InitCoords(x0); // set initial coordinates
  M.InitCoords(x1); // set initial coordinates
  for(int i=0;i<n;i++){mat.at(i)=M.Material(i);}
  for(int i=0;i<n;i++){V0.at(i)=M.Volume(i);}
  for(int i=0;i<n;i++){V1.at(i)=M.Volume(i);}
  for(int i=0;i<n;i++){int imat(mat[i]);p.at(i)=state[imat-1][2];}  // load pressure field from intial flux state
  for(int i=0;i<n;i++){int imat(mat[i]);d0.at(i)=state[imat-1][0];} // load density field from initial flux state
  for(int i=0;i<n;i++){int imat(mat[i]);d1.at(i)=state[imat-1][0];}
  for(int i=0;i<n;i++){m.at(i)=d0[i]*V0[i];}                        // calculate the initial mass field
  for(int i=0;i<n;i++){e0.at(i)=E(d0[i],p[i]);}                     // invert the eos to start the energy field
  for(int i=0;i<n;i++){e1.at(i)=E(d1[i],p[i]);}
  for(int i=0;i<n;i++){q.at(i)=0.0;}
  for(int i=0;i<n;i++){c.at(i)=sqrt(GAMMA*p[i]/d0[i]);}

  for(int idim=0;idim<ndims;idim++){
    for(int i=0;i<n;i++){
for(int iloc=0;iloc<S.nloc();iloc++){
  long j(M.Vertex(i,iloc));
u0.at(idim).
}
    }

    cout<<" v size= "<<u0.at(idim).size()<<endl;
  }





//  for(int i=0;i<ng+1;i++){u0.at(i)=(x0[i]<=0.5)?l[1]:r[1];}
//  for(int i=0;i<ng+1;i++){u1.at(i)=u0[i];}

// debug
// output initial data
  initial_data(n,ndims,nmats,mat,d0,V0,m,e0,p,x0,u0,step,time);
  cout<<"main(): 2D not yet operational, stopping."<<endl;
  exit(1);
// debug

// start the Riemann solvers from initial flux states

  Riemann R0(Riemann::exact,l,r),R1(Riemann::exact,l,r);

// set output precision

  cout<<fixed<<setprecision(17);

// time integration

  while(time<=ENDTIME){

// calculate a new stable time-step

//    for(int i=0;i<ng;i++){double l(x0[i+1]-x0[i]);dt_cfl.at(i)=(COURANT*l/sqrt((c[i]*c[i])+2.0*q[i]/d0[i]));} // impose the CFL limit on each element

// reduce across element and apply a saftey factor

//    dts.at(0)=DTSFACTOR*(*min_element(dt_cfl.begin(), dt_cfl.end()));

// calculate a time step so we hit I/O event times smoothly and accurately

//    double dt_g(-1.0*remainder(time,VISFREQ)); //  time until next graphics event
//    dts.at(1)=1.0e6; // next I/O event too far away to care

//    if( dt_g/(4.0*dt)<1.0 && dt_g>1.0e-12){
//      if(dt_g/dt<4.000000001) {dts.at(1)=dt_g/4.0;}
//      if(dt_g/dt<3.000000001) {dts.at(1)=dt_g/3.0;}
//      if(dt_g/dt<2.000000001) {dts.at(1)=dt_g/2.0;}
//      if(dt_g/dt<1.000000001) {dts.at(1)=dt_g;}
//      if(abs(1.0-(dts[1]/dt))>1.0e-5){
//        cout<<"       time-step cut on I/O event approach, dt lowered by "<<(1.0-dts[1]/dt)*100.00<<"% to "<<dts[1]<<endl;
//      }
//    }

// choose the smallest time step

//    dt=(*min_element(dts.begin(), dts.end()));
////    dt=DTSTART;cout<<"DT HARDWIRED !! "<<endl;

    cout<<fixed<<setprecision(5)<<"  step "<<step<<" time= "<<time<<" dt= "<<dt<<fixed<<setprecision(5)<<" energy (i/k/tot)= "<<ie<<" "<<ke<<" "<<ie+ke<<endl;

// graphics output

//      if(abs(remainder(time,VISFREQ))<1.0e-12){
//        silo(&x0,&d0,&p,&e0,&q,&c,&u0,&mat,step,time,&S);
//      }else{
//        if(abs(remainder(step,VISFREQ))==0){silo(&x0,&d0,&p,&e0,&q,&c,&u0,&mat,step,time,&S);}
//      }

// move the nodes to their full-step position

//    for(int i=0;i<ng+1;i++){x1.at(i)=x0[i]+u0[i]*dt;}

// update kinetic energy for conservation checks

//    ke=0.0;for(int i=1;i<ng;i++){ke+=0.25*(m[i-1]+m[i])*u0[i]*u0[i];} // include level 1 halo

// evolve the Riemann problems to the end of the time-step on the end of time-step meshes

//    vector<double> r0x,rx;vempty(r0x); // sample point coordinates
////    for(long i=0;i<NSAMPLES;i++){r0x.push_back(x1[0]+(i*(x1[ng]-x1[0])/double(NSAMPLES)));} // sample points
//    for(long i=0;i<NSAMPLES;i++){r0x.push_back(0.0+(i*(1.0-0.0)/double(NSAMPLES)));} // sample points
//    R0.profile(&r0x,time+dt); // Riemann solution at the sample points along the mesh
//    rx.clear();for(int i=0;i<ng;i++){rx.push_back(x1[i]);rx.push_back(0.5*(x1[i]+x1[i+1]));rx.push_back(x1[i+1]);}
//    R1.profile(&rx,time+dt);

// update cell volumes at the full-step

//    for(int i=0;i<ng;i++){V1.at(i)=x1[i+1]-x1[i];if(V1[i]<VTOL){cout<<"-'ve volume detected in cell "<<i<<endl;exit(1);}}

// update cell density at the full-step

//    for(int i=0;i<ng;i++){d1.at(i)=m[i]/V1[i];}

// debug
//    for(int i=0;i<ng;i++){d1.at(i)=R1.density(3*i+1);} // 3*i+1 is cell-centre address
// debug

// update cell energy at the full-step

//    for(int i=0;i<ng;i++){e1.at(i)=max(ECUT,e0[i]-((p[i]+q[i])*(V1[i]-V0[i]))/m[i]);}

// update internal energy for conservation checks

//    ie=0.0;for(int i=0;i<ng;i++){ie+=e1[i]*m[i];}

// debug
//    for(int i=0;i<ng;i++){e1.at(i)=R1.energy(3*i+1);} // 3*i+1 is cell-centre address
// debug

// update cell pressure at the full-step

//    for(int i=0;i<ng;i++){p.at(i)=P(d1[i],e1[i]);if(p[i]<0.0){cout<<"-'ve pressure detected in cell "<<i<<" e1= "<<e1[i]<<endl;exit(1);}}

// bulk q

//    for(int i=0;i<ng;i++){
//      c.at(i)=sqrt(GAMMA*p[i]/d1[i]);
//      double l(x1[i+1]-x1[i]),divu((d0[i]-d1[i])/(d1[i]*dt));
//      if(divu<0.0){
//        q.at(i)=d0[i]*l*divu*((cq*l*divu)-cl*c[i]);
//      }else{
//        q.at(i)=0.0; // turn off q as cell divergence indicates expansion
//      }
//    }

// assemble acceleration field

//  Matrix A(n+3);double b[n+3],x[n+3];for(int i=0;i<n+3;i++){b[i]=0.0;x[i]=0.0;}
//  double m[2][2]={};m[0][0]=0.5;m[1][1]=0.5;// for mass lumping

// next block codes for matrix assembly with a continuous Galerkin finite element type

//  for(int iel=0;iel<ng;iel++){
//    for(int iloc=0;iloc<S.nloc();iloc++){
//      int i(iel+iloc); // column address in the global matrix
//      if((i>0&&i<ng)){for(int gi=0;gi<S.ngi();gi++){b[i-1]+=(p[iel]+q[iel])*S.dvalue(0,iloc,gi)*S.wgt(gi);}} // integrate the shape derivative for rhs
//      for(int jloc=0;jloc<S.nloc();jloc++){
//        double nn(0.0); // mass matrix
//        int j(iel+jloc); // row address in the global matrix
//        for(int gi=0;gi<S.ngi();gi++){
//          nn+=S.value(iloc,gi)*S.value(jloc,gi)*S.wgt(gi)*0.5*(x1[iel+1]-x1[iel]); // DG & notes use this - double check ??
//        }
//        if((i>0&&i<ng)&&(j>0&&j<ng)){
//          A.add(i-1,j-1,d1[iel]*nn);
////          A.add(i-1,j-1,d1[iel]*(x1[iel+1]-x1[iel])*m[iloc][jloc]); // use mass lumping
//        }
//      }
//    }
//
//  }

// solve global system

//  A.solve(x,b);

// advance the solution

//  for(int i=0;i<n+3;i++){u1.at(i+1)=u0[i+1]+x[i]*dt;}

// impose boundary constraints on the acceleration field

//  u1.at(1)=u1[2];u1.at(ng-1)=u1[ng-2];
//  u1.at(0)=u1[1];u1.at(ng)=u1[ng-1];

// some output

//    f1.open("exact.dat");f2.open("dpe.dat");f3.open("q.dat");f4.open("u.dat");f5.open("r.dat");
//    f1<<fixed<<setprecision(17);f2<<fixed<<setprecision(17);f3<<fixed<<setprecision(17);f4<<fixed<<setprecision(17);
//    for(int i=0;i<NSAMPLES;i++){f1<<r0x[i]<<" "<<R0.density(i)<<" "<<R0.pressure(i)<<" "<<R0.velocity(i)<<" "<<R0.energy(i)<<endl;}
//    for(int i=1;i<NSAMPLES;i++){if(R0.region(i)!=R0.region(i-1)){f5<<r0x[i]<<" -20000.0"<<endl<<r0x[i]<<" 20000.0"<<endl<<r0x[i]<<" -20000.0"<<endl;}} // sample region
//    for(int i=2;i<n+2;i++){f2<<0.5*(x1[i]+x1[i+1])<<" "<<d1[i]<<" "<<p[i]<<" "<<e1[i]<<" "<<endl;}
//    for(int i=2;i<n+2;i++){f3<<0.5*(x1[i]+x1[i+1])<<" "<<q[i]<<endl;}
//    for(int i=2;i<n+3;i++){f4<<x1[i]<<" "<<u1[i]<<endl;}
//    f1.close();f2.close();f3.close();f4.close();f5.close();

// advance the time step

    time+=dt;
    step++;

// debug
//    for(int i=0;i<ng;i++){u1.at(i)=R1.velocity(3*i);u1.at(i+1)=R1.velocity(3*i+2);} // 3*i,3*i+2 are nodal address
// debug

// advance the solution for the new time step

//    for(int i=0;i<ng+1;i++){u0.at(i)=u1[i];}
//    for(int i=0;i<ng+1;i++){x0.at(i)=x1[i];}
//    for(int i=0;i<ng;i++){e0.at(i)=e1[i];}
//    for(int i=0;i<ng;i++){V0.at(i)=V1[i];}
//    for(int i=0;i<ng;i++){d0.at(i)=d1[i];}
//    for(int i=0;i<ng;i++){c.at(i)=sqrt(GAMMA*(p[i]+q[i])/d1[i]);}

  }

// estimate convergence rate in the L1/L2 norms using a Riemann solution as the exact solution

//  vector<double> rx;vempty(rx);  // sample points
//  double xstart(0.0),xstop(1.0); // sample range - should be whole mesh though bc.s may artificially reduce convergence
//  for(long i=0;i<ng+1;i++){if(x1[i]>=xstart&&x1[i]<=xstop){rx.push_back(x1[i]);}}
//  R0.profile(&rx,ENDTIME); // evolve a Riemann problem on the sample range to final time level
//
//  double l1(0.0),l2(0.0),l1r(0.0),l2r(0.0),l1d(0.0),l2d(0.0),l1n(0.0),l2n(0.0);
//  int ii(0);

// numerator and denominator for each norm

//  for(int i=0;i<ng+1;i++){
//    if(x1[i]>=xstart&&x1[i]<=xstop){
//      double err(abs(R0.velocity(ii)-u1[i])); // absolute error
//      l1d+=abs(R0.velocity(ii)); // l1 denominator
//      l2d+=R0.velocity(ii)*R0.velocity(ii); // l2 denominator
//      l1n+=err*V1[i]; // l1 numerator
//      l2n+=err*err*V1[i]; // l2 numerator
//      ii++;
//    }
//  }

// construct norms and relative errors

//  l1=l1n/ii; // L1 error norm
//  l1r=l1n/l1d; // L1 relative error norm
//  l2=sqrt(l2n/ii); // L2 error norm
//  l2r=sqrt(l2n/l2d); // L2 relative error norm

//  cout<<endl<<fixed<<setprecision(10)<<"  Error Estimates (grid spacing h= "<<(xstop-xstart)/ii<<"):"<<endl;

//  cout<<"  L1 norm= "<<l1<<" (relative error= "<<l1r<<")"<<endl;
//  cout<<"  L2 norm= "<<l2<<" (relative error= "<<l2r<<")"<<endl;
//  cout<<"  No. points sampled= "<<ii<<endl;
//  cout<<"  Range sampled= "<<xstart<<","<<xstop<<endl;

//  cout<<"Normal termination."<<endl;

  return 0;
}

// empty a vector

void vempty(vector<double>&v){vector<double> e;v.swap(e);return;}

// return today's date

string date(){time_t now = time(0);return ctime(&now);}

// output initial data - args are passed by ref to prevent a copy, and are const to prevent changes

void initial_data(int const n,int const ndims,int const nmats,VI const &mat,VD const &d,VD const &V,
                  VD const &m,VD const &e,VD const &p,VVD const &x, VVD const &u, int const &s, double const &t){

  cout<<endl;
  cout<<"             RhoLo"<<endl;
  cout<<endl;
  cout<<"Date: "<<date()<<endl;

  cout<<"Initial data for the simulation"<<endl;

  cout<<"Number of dimensions:             "<<ndims<<endl;
  cout<<"Number of cells:                  "<<n<<endl;
  cout<<"Number of materials:              "<<nmats<<endl;

// material average data

  vector<double> volm(nmats);            // material volume
  vector<double> pm(nmats);              // material pressure
  vector<double> mm(nmats);              // material mass
  vector<double> em(nmats);              // total energy
  vector<double> um(nmats);              // material velocity
  vector<double> iem(nmats);             // internal energy
  vector<double> kem(nmats);             // kinetic energy

// totals

  double dtot(0.0),mtot(0.0),vtot(0.0),ptot(0.0),etot(0.0),ietot(0.0),ketot(0.0);

// state of each material

  for(int i=0;i<n;i++){
    int imat(mat[i]);
    volm.at(imat-1)+=V[i];               // subtract 1 to convert to address in the vector
    pm.at(imat-1)+=p[i]*V[i];
    mm.at(imat-1)+=m[i];
    iem.at(imat-1)+=e[i]*m[i];
    kem.at(imat-1)+=0.0;
    mtot+=m[i];
    vtot+=V[i];
    ietot+=e[i]*m[i];
    ketot+=0.0;
  }


// print out the material averages in a table

  cout<<endl;
  cout<<fixed<<setprecision(6)<<"Material states at time "<<t<<" on step "<<s<<endl;
  cout<<"~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"<<endl;
  cout<<"Mat.  "<<"Density  "<<"Mass     "<<"Volume   "<<"Pressure "<<"Energy (Spec. I/I/K/Total)"<<endl;
  cout<<"~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"<<endl;
  for(int i=0;i<nmats;i++){
    int matnum(i+1);
    em[i]=iem[i]+kem[i];
    cout<<fixed<<setprecision(6);
    cout<<matnum<<"     "<<mm[i]/volm[i]<<" "<<mm[i]<<" "<<volm[i]<<" "<<pm[i]/volm[i]<<" "<<iem[i]/mm[i]<<" "<<kem[i]<<" "<<em[i]<<endl;
  }

  dtot=mtot/vtot;
  ietot=ietot/mtot;

  cout<<"total "<<dtot<<" "<<mtot<<" "<<vtot<<" "<<P(dtot,ietot)<<" "<<ietot<<" "<<ketot<<" "<<etot<<endl;

  cout<<"~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"<<endl;

  cout<<endl;

  return;

}
