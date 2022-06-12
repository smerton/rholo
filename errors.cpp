// estimate convergence rate in the L1/L2 norms using an exact solution

// Author S. R. Merton

#include <iostream>
#include <vector>
#include <cmath>       // sin/cos
#include "globals.h"   // defines
#include "mesh.h"      // mesh class
#include "shape.h"     // shape class
#include "riemann.h"   // riemann solver
#include "jacobian.h"  // jacobian
#include "tests.h"     // test problems

using namespace std;

void errors(Mesh const &M,Shape const &S,Shape const &T,VVD const &xk,VVD const &xt,VVD const &u,VD const &V,int const &test_problem){

  double l1(0.0),l2(0.0),h(0.0);
  double dpi(4.0*atan(1.0));

  vector<vector<double> > u_exact(M.NDims());

// exact solution for the test problem

  switch(test_problem){

    case(TAYLOR):

// Taylor Green vortex

      for(long i=0;i<xk.at(0).size();i++){
        u_exact.at(0).push_back(sin(dpi*xk.at(0).at(i))*cos(dpi*xk.at(1).at(i)));
        u_exact.at(1).push_back(-cos(dpi*xk.at(0).at(i))*sin(dpi*xk.at(1).at(i)));
      }

      break;

    case(RAYLEIGH):

// Rayleigh-Taylor instability

      return;

      break;

    case(NOH):

// Noh stagnation shock

      break;

    case(SEDOV):

// Sedov expanding shock

      break;

    case(SOD):

// Sod's shock tube

      break;

    case(R2R):

// 123 problem

      break;

    case(SALTZMANN):

// Saltzmann piston

      return;

      break;

  }

// compute error between computational and exact solution

  for(int i=0;i<M.NCells();i++){

    double l1_i(0.0),l2_i(0.0);
    VD detJ(S.ngi(),0.0);
    VVVD detDJ(M.NDims(),VVD(S.nloc(),VD(S.ngi(),0.0)));

    jacobian(i,xk,M,S,detJ,detDJ);

    h+=V.at(i);

// solution and exact solution at each quadrature point

    vector<double> ugi_exact(M.NDims(),0.0),ugi(M.NDims(),0.0);
    for(int gi=0;gi<S.ngi();gi++){
      for(int iloc=0;iloc<S.nloc();iloc++){
        ugi_exact.at(0)+=u_exact.at(0).at(M.GlobalNode_CFEM(i,iloc))*S.value(iloc,gi);
        ugi_exact.at(1)+=u_exact.at(1).at(M.GlobalNode_CFEM(i,iloc))*S.value(iloc,gi);
        ugi.at(0)+=u.at(0).at(M.GlobalNode_CFEM(i,iloc))*S.value(iloc,gi);
        ugi.at(1)+=u.at(1).at(M.GlobalNode_CFEM(i,iloc))*S.value(iloc,gi);

      }

// integrate error across the element

      l1_i+=abs((ugi_exact.at(0)-ugi.at(0)))*S.wgt(gi);
      l1_i+=abs((ugi_exact.at(1)-ugi.at(1)))*S.wgt(gi);

      l2_i+=(ugi_exact.at(0)-ugi.at(0))*(ugi_exact.at(0)-ugi.at(0))*S.wgt(gi);
      l2_i+=(ugi_exact.at(1)-ugi.at(1))*(ugi_exact.at(1)-ugi.at(1))*S.wgt(gi);

    }

// contribution to norm

    l1+=l1_i*V.at(i);
    l2+=l2_i*V.at(i);

  }

  h=sqrt(h/M.NCells());
  l1/=M.NCells();
  l2=sqrt(l2/M.NCells());


  cout<<endl;
  cout<<"  errors(): Error Estimators (average grid spacing h= "<<h<<")"<<endl;
  cout<<"  errors(): L1 norm= "<<l1<<" (relative error= "<<l1<<")"<<endl;
  cout<<"  errors(): L2 norm= "<<l2<<" (relative error= "<<l2<<")"<<endl;

  return;

}
