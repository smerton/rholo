// Finite element variant of RhoLo (Riemann-based Hydro in One-dimension at Low Order - RhoLo)
// RhoLo is an ultra simple 1-D discontinuous finite element (DG) hydrodynamics test code
// This finite element variant solves the Euler equations in their non-conservative form in the fluid frame (the Lagrangian frame)
// using a mixed continuous finite element method (cell-centred thermodynamic variable d,rho,e with node 
// centred kinematic variables u,a) and bulk viscosity q to increase entropy across element boundaries, initial 
// implementation is only first order in time
// for graphics: convert -density 300 filename.png filename.pdf

// Author S. R. Merton

#define DTSTART 0.0005          // insert a macro for the first time step
#define ENDTIME 0.25            // insert a macro for the end time
#define GAMMA 1.4               // ratio of specific heats for ideal gases
#define ECUT 1.0e-8             // cut-off on the energy field
#define NSAMPLES 500            // number of sample points for the exact solution
#define VISFREQ 10000           // frequency of the graphics dumps
#define VD vector<double>       // vector of doubles
#define VTOL 1.0e-10            // threshold for volume errors
#define COURANT 0.333           // Courant number for CFL condition
#define DTSFACTOR 0.5           // safety factor on time-step control
#define KNOD i*(K.nloc()-1)+k   // global node number on kinematic mesh
#define TNOD i*T.nloc()+j       // global node number on thermodynamic mesh
#define GPNT i*T.ngi()+gi       // global address of Gauss point gi in element i
#define DX0 x0[i*(K.nloc()-1)+K.nloc()-1]-x0[i*(K.nloc()-1)]            // cell width at start of step
#define DX1 x1[i*(K.nloc()-1)+K.nloc()-1]-x1[i*(K.nloc()-1)]            // cell width at end of step
#define CENTROID 0.5*(x1[i*(K.nloc()-1)+K.nloc()-1]+x1[i*(K.nloc()-1)]) // cell centroid
#define ROW (i-1)*(K.nloc()-1)+iloc                                     // row address in global matrix
#define COL (i-1)*(K.nloc()-1)+jloc                                     // column address in global matrix
#define XGI for(int j=0;j<T.nloc();j++){xgi+=T.value(j,gi)*x3[TNOD];}   // coordinates of integration point

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
  int const n(100),ng(n+4);                             // no. ncells, no. ghosts
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
  double F[nkg][ntg],FT[ntg][nkg];                      // force matrix
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
    for(int j=0;j<T.nloc();j++){
      double pos(-1.0+j*2.0/(T.nloc()-1));x2.at(TNOD)=0.0;
      for(int k=0;k<K.nloc();k++){x2.at(TNOD)+=K.value(k,pos)*x0[KNOD];}
    }
  }

// set time-0 Jacobian

  for(int i=0;i<ng;i++){
    for(int gi=0;gi<T.ngi();gi++){
      detJ0.at(GPNT)=0.0;
      for(int k=0;k<K.nloc();k++){
        detJ0.at(GPNT)+=K.dvalue(k,gi)*x0[KNOD];
      }
      if(detJ0.at(GPNT)<0.0){cout<<"-'ve determinant of J0 detected in cell "<<i<<endl;exit(1);}
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
      for(int gi=0;gi<K.ngi();gi++){
        dt_cfl.at(GPNT)=COURANT*(DX0/sqrt((c[GPNT]*c[GPNT])+2.0*q[GPNT]/d0[GPNT]));
      }
    }

// reduce across element and apply a saftey factor

    double dt=DTSFACTOR*(*min_element(dt_cfl.begin(), dt_cfl.end()));

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

    vector<double> r0x,rx,rx2,rx3;vempty(r0x);vempty(rx);vempty(rx2);vempty(rx3); // sample point coordinates
//    for(long i=0;i<NSAMPLES;i++){r0x.push_back(0.0+(i*(1.0-0.0)/double(NSAMPLES)));} // sample points
    for(long i=0;i<NSAMPLES;i++){r0x.push_back(x1[0]+(i*(x1[nkg-1]-x1[0])/double(NSAMPLES)));} // sample points
    R0.profile(&r0x,time+dt); // Riemann solution at the sample points along the mesh
    for(int i=0;i<nkg;i++){rx.push_back(x0[i]);}
    R1.profile(&rx,time+dt);
    for(int i=0;i<ng;i++){for(int gi=0;gi<T.ngi();gi++){double xgi(0.0);XGI;rx2.push_back(xgi);}}
    R2.profile(&rx2,time+dt);
    for(int i=0;i<ntg;i++){rx3.push_back(x3[i]);}
    R3.profile(&rx3,time+dt);

// update cell volumes at the full-step

    for(int i=0;i<ng;i++){V1.at(i)=DX1;if(V1[i]<VTOL){cout<<"-'ve volume detected in cell "<<i<<endl;exit(1);}}

// update Jacobian

    for(int i=0;i<ng;i++){
      for(int gi=0;gi<T.ngi();gi++){
        detJ.at(GPNT)=0.0;
        for(int k=0;k<K.nloc();k++){
          detJ.at(GPNT)+=K.dvalue(k,gi)*x1[KNOD];
        }
        if(detJ.at(GPNT)<0.0){cout<<"-'ve determinant of J detected in cell "<<i<<endl;exit(1);}
      }
    }

// assemble force matrix to connect thermodynamic/kinematic spaces, this can be used as rhs of both e/u eqns

    for(long i=0;i<nkg;i++){for(long j=0;j<ntg;j++){F[i][j]=0.0;}}
    for(int i=0;i<ng;i++){
      for(int k=0;k<K.nloc();k++){
        for(int j=0;j<T.nloc();j++){
          double f(0.0);
          for(int gi=0;gi<T.ngi();gi++){
//            f+=(p[GPNT]+q[GPNT])*K.dvalue(k,gi)*T.value(j,gi)*detJ[gi]*T.wgt(gi)*JI;
            f+=(p[GPNT]+q[GPNT])*K.dvalue(k,gi)*T.value(j,gi)*T.wgt(gi);
          }
          F[KNOD][TNOD]=+f;
        }
      }
    }
    for(long i=0;i<nkg;i++){for(long j=0;j<ntg;j++){FT[j][i]=F[i][j];}} // write transpose

// update cell density at the full-step

    for(int i=0;i<ng;i++){for(int gi=0;gi<T.ngi();gi++){d1[GPNT]=dinit[i]*detJ0[GPNT]/detJ[GPNT];}}

// debug
    for(int i=0;i<ng;i++){for(int gi=0;gi<T.ngi();gi++){d1.at(GPNT)=R2.density(GPNT);}}
// debug

// update Jacobian for energy solve

    for(int i=0;i<ng;i++){
      for(int gi=0;gi<T.ngi();gi++){
        detJ.at(GPNT)=0.0;
        for(int j=0;j<T.nloc();j++){
          detJ.at(GPNT)+=T.dvalue(j,gi)*x3[TNOD];
        }
        if(detJ.at(GPNT)<0.0){cout<<"-'ve determinant of J detected in cell "<<i<<endl;exit(1);}
      }
    }

// assemble finite element energy field on discontinuous thermodynamic grid
// for( struct {int i; double j;} v = {0, 3.0}; v.i < 10; v.i++, v.j+=0.1)

    {Matrix A(T.nloc());double b[T.nloc()],x[T.nloc()];
    for(int i=0;i<ng;i++){int k(0);
      for(int iloc=0;iloc<T.nloc();iloc++,k++){
        b[iloc]=0.0;int j(0);
        for(int jloc=0;jloc<T.nloc();jloc++,j++){
          b[iloc]+=FT[TNOD][KNOD]*u1[KNOD];
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
    for(int i=0;i<ntg;i++){e1.at(i)=R3.energy(i);}
// debug

// update internal energy for conservation checks

    ie=0.0;for(int i=0;i<ng;i++){for(int j=0;j<T.nloc();j++){ie+=e1[TNOD]*0.5*m[i];}}

// update pressure at the full-step at the integration points

    for(int i=0;i<ng;i++){
      for(int gi=0;gi<T.ngi();gi++){
        double egi(0.0);
        for(int j=0;j<T.nloc();j++){
          egi+=T.value(j,gi)*e1[TNOD];
        }
        p.at(GPNT)=P(d1[GPNT],egi); // this is the EOS call
      }
    }

// debug
    for(int i=0;i<ng;i++){for(int gi=0;gi<T.ngi();gi++){p.at(GPNT)=R2.pressure(GPNT);}}
// debug

// update sound speed

    for(int i=0;i<ng;i++){for(int gi=0;gi<T.ngi();gi++){c.at(GPNT)=sqrt(GAMMA*p[GPNT]/d1[GPNT]);}}

// bulk q

    for(int i=0;i<ng;i++){
      double l(DX1);
      for(int gi=0;gi<T.ngi();gi++){
        double divu((d0[GPNT]-d1[GPNT])/(d1[GPNT]*dt));
        if(divu<0.0){
          q.at(GPNT)=d0[GPNT]*l*divu*((cq*l*divu)-cl*c[GPNT]);
        }else{
          q.at(GPNT)=0.0; // turn off q as cell divergence indicates expansion
        }
      }
    }

// update Jacobian for momentum solve

    for(int i=0;i<ng;i++){
      for(int gi=0;gi<T.ngi();gi++){
        detJ.at(GPNT)=0.0;
        for(int k=0;k<K.nloc();k++){
          detJ.at(GPNT)+=K.dvalue(k,gi)*x1[KNOD];
        }
        if(detJ.at(GPNT)<0.0){cout<<"-'ve determinant of J detected in cell "<<i<<endl;exit(1);}
      }
    }

// assemble acceleration field

    {Matrix A((ng-2)*(K.nloc()-1)+1);double b[(ng-2)*(K.nloc()-1)+1],x[(ng-2)*(K.nloc()-1)+1];
    for(int i=1;i<ng-1;i++){int k(0);
      for(int iloc=0;iloc<K.nloc();iloc++,k++){
        b[ROW]=0.0;x[ROW]=0.0;int j(0);for(int jloc=0;jloc<T.nloc();jloc++,j++){b[ROW]+=F[KNOD][TNOD]*1.0;}
        for(int jloc=0;jloc<K.nloc();jloc++){
          double nn(0.0); // mass matrix
          for(int gi=0;gi<K.ngi();gi++){
            nn+=d1[GPNT]*K.value(iloc,gi)*K.value(jloc,gi)*detJ[GPNT]*K.wgt(gi);
          }
          A.add(ROW,COL,nn);
        }
      }

    }

// solve global system

    A.solve(x,b);

// update acceleration field

    for(int i=1;i<ng-1;i++){for(int k=0;k<K.nloc();k++){int iloc(k);u1.at(KNOD)=u1[KNOD]+x[ROW]*dt;}}

    }

// impose boundary constraints on the acceleration field

    for(int j=0;j<K.nloc()-1;j++){u1.at(j)=u1[K.nloc()];}
    for(int j=0;j<K.nloc()-1;j++){int i(ng-1);u1.at(i*(K.nloc()-1)+j+1)=u1[i*(K.nloc()-1)];}

// debug
    for(long i=0;i<nkg;i++){u1.at(i)=R1.velocity(i);}
// debug

// some output

    f1.open("exact.dat");f2.open("e.dat");f3.open("u.dat");f4.open("dp.dat");
    f1<<fixed<<setprecision(17);f2<<fixed<<setprecision(17);f3<<fixed<<setprecision(17);f4<<fixed<<setprecision(17);
    for(int i=0;i<NSAMPLES;i++){
      f1<<r0x[i]<<" "<<R0.density(i)<<" "<<R0.pressure(i)<<" "<<R0.velocity(i)<<" "<<R0.energy(i)<<endl;
    }
    for(long i=0;i<ntg;i++){f2<<x3[i]<<" "<<e1[i]<<endl;}
    for(long i=0;i<nkg;i++){f3<<x1[i]<<" "<<u1[i]<<endl;}
    for(int i=0;i<ng;i++){for(int gi=0;gi<T.ngi();gi++){double xgi(0.0);XGI;f4<<xgi<<" "<<d1[GPNT]<<" "<<p[GPNT]<<endl;}}
//    cout<<"NODE POS: "<<time<<" "<<x1[(ng/2)*(K.nloc()-1)]<<" "<<x1[(ng/2)*(K.nloc()-1)+1]<<endl;
    f1.close();f2.close();f3.close();f4.close();

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
    for(int i=0;i<ng*T.ngi();i++){d0.at(i)=d1[i];}

// debug
//  cout<<"debug stop."<<endl;
//  exit(1);
// debug

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
