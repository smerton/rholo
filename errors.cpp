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

// compute error between computational and exact solution in the L1/L2 norms

  for(int i=0;i<M.NCells();i++){

    double l1_i(0.0),l2_i(0.0);

    VD detJ(S.ngi(),0.0);
    VVVD detDJ(M.NDims(),VVD(S.nloc(),VD(S.ngi(),0.0)));

    jacobian(i,xk,M,S,detJ,detDJ);

    h+=V.at(i);

// solution and exact solution at each quadrature point

//    vector<double> ugi_exact(M.NDims(),0.0),ugi(M.NDims(),0.0);

    for(int gi=0;gi<S.ngi();gi++){
      double xgi(0.0),ygi(0.0);
      vector<double> ugi_exact(M.NDims(),0.0),ugi(M.NDims(),0.0);
      for(int iloc=0;iloc<S.nloc();iloc++){
        xgi+=xk.at(0).at(M.GlobalNode_CFEM(i,iloc))*S.value(iloc,gi);
        ygi+=xk.at(1).at(M.GlobalNode_CFEM(i,iloc))*S.value(iloc,gi);
//        ugi_exact.at(0)+=u_exact.at(0).at(M.GlobalNode_CFEM(i,iloc))*S.value(iloc,gi); // better to calculate exact
//        ugi_exact.at(1)+=u_exact.at(1).at(M.GlobalNode_CFEM(i,iloc))*S.value(iloc,gi); // at gi surely ??
        ugi.at(0)+=u.at(0).at(M.GlobalNode_CFEM(i,iloc))*S.value(iloc,gi);
        ugi.at(1)+=u.at(1).at(M.GlobalNode_CFEM(i,iloc))*S.value(iloc,gi);
      }

// exact solution for the test problem calculated at the integration point gi

      switch(test_problem){

        case(TAYLOR):

// Taylor Green vortex

          ugi_exact.at(0)=sin(dpi*xgi)*cos(dpi*ygi);
          ugi_exact.at(1)=-cos(dpi*xgi)*sin(dpi*ygi);

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

// integrate errors across the element

      l1_i+=abs((ugi_exact.at(0)-ugi.at(0)))*detJ.at(gi)*S.wgt(gi);
      l1_i+=abs((ugi_exact.at(1)-ugi.at(1)))*detJ.at(gi)*S.wgt(gi);

      l2_i+=(ugi_exact.at(0)-ugi.at(0))*(ugi_exact.at(0)-ugi.at(0))*detJ.at(gi)*S.wgt(gi);
      l2_i+=(ugi_exact.at(1)-ugi.at(1))*(ugi_exact.at(1)-ugi.at(1))*detJ.at(gi)*S.wgt(gi);

    }

// contribution in the L1/L2 norms

    l1+=l1_i;
    l2+=l2_i;

  }

// normalise to the number of cells on the mesh

  h=sqrt(h/M.NCells());

  l1/=M.NCells();
  l2=sqrt(l2/M.NCells());

// output errors in the L1/L2 norms

  cout<<endl;
  cout<<"  Error Estimators on average grid spacing h= "<<h<<":"<<endl;
  cout<<"    L1 norm= "<<l1<<" (relative error= "<<l1<<")"<<endl;
  cout<<"    L2 norm= "<<l2<<" (relative error= "<<l2<<")"<<endl;

  return;

}
