// This function codes for the jacobian of the lagrangian motion

// Author S. R. Merton

#include <iostream>
#include <vector>
#include <iomanip>
#include "globals.h"   // defines
#include "shape.h"     // signature of the shape class
#include "mesh.h"      // signature of the mesh class

using namespace std;

// calculate a jacobian and the determinant

void jacobian(int const &i,VVD const &x,Mesh const &M,Shape const &S,VD &detJ,VVVD &detDJ){

// loop over quadrature points and calculate the jacobian

  for(int gi=0;gi<S.ngi();gi++){

    double dxdu(0.0),dydu(0.0),dxdv(0.0),dydv(0.0);

// derivatives of the physical coordinates at the quadrature points

    for(int iloc=0;iloc<S.nloc();iloc++){
      long gloc=(S.type()==CONTINUOUS)?M.GlobalNode_CFEM(i,iloc):M.GlobalNode_DFEM(i,iloc);
      dxdu+=x.at(0).at(gloc)*S.dvalue(0,iloc,gi); // dx/du
      dxdv+=x.at(0).at(gloc)*S.dvalue(1,iloc,gi); // dx/dv
      dydu+=x.at(1).at(gloc)*S.dvalue(0,iloc,gi); // dy/du
      dydv+=x.at(1).at(gloc)*S.dvalue(1,iloc,gi); // dy/dv
    }

// calculate the determinant at the quadrature point and commit to the vector

    detJ.at(gi)=dxdu*dydv-dxdv*dydu;

// determinants for the deriavtives at the quadrature points

    for(int iloc=0;iloc<S.nloc();iloc++){
      detDJ.at(0).at(iloc).at(gi)=(dydv*S.dvalue(0,iloc,gi)-dydu*S.dvalue(1,iloc,gi))/detJ[gi];
      detDJ.at(1).at(iloc).at(gi)=(-dxdv*S.dvalue(0,iloc,gi)+dxdu*S.dvalue(1,iloc,gi))/detJ[gi];
    }

  }

  return;

}

// calculate a jacobian for the Lagrangian motion and return the determinant

void jacobian(int const &i,VVD const &x0,VVD const &x,Mesh const &M,Shape const &S,VD &detJs,VVVD &Js){

// loop over quadrature points and calculate jacobians J0 and J

  for(int gi=0;gi<S.ngi();gi++){

    double dx0du(0.0),dy0du(0.0),dx0dv(0.0),dy0dv(0.0);
    double dxdu(0.0),dydu(0.0),dxdv(0.0),dydv(0.0);

// derivatives of the physical coordinates at the quadrature points

    for(int iloc=0;iloc<S.nloc();iloc++){

      long gloc=(S.type()==CONTINUOUS)?M.GlobalNode_CFEM(i,iloc):M.GlobalNode_DFEM(i,iloc);

// at time 0

//      dx0du+=x0.at(0).at(gloc)*S.dvalue(0,iloc,gi); // dx0/du
//      dx0dv+=x0.at(0).at(gloc)*S.dvalue(1,iloc,gi); // dx0/dv
//      dy0du+=x0.at(1).at(gloc)*S.dvalue(0,iloc,gi); // dy0/du
//      dy0dv+=x0.at(1).at(gloc)*S.dvalue(1,iloc,gi); // dy0/dv

//      dx0du=1.0; // mod for direction s
//      dx0dv=1.0; // mod for direction s
//      dy0du+=x0.at(1).at(gloc)*S.dvalue(0,iloc,gi); // dy0/du
//      dy0dv+=x0.at(1).at(gloc)*S.dvalue(1,iloc,gi); // dy0/dv

//      dx0du+=x0.at(0).at(gloc)*S.dvalue(0,iloc,gi); // dx0/du
//      dx0dv=1.0; // mod for direction s
//      dy0du=1.0; // mod for direction s
//      dy0dv+=x0.at(1).at(gloc)*S.dvalue(1,iloc,gi); // dy0/dv

//      dx0du=1.0; // mod for direction s
//      dx0dv=1.0; // mod for direction s
//      dy0du=1.0;
//      dy0dv+=x0.at(1).at(gloc)*S.dvalue(1,iloc,gi); // dy0/dv

      dx0du=1.0; // mod for direction s
      dx0dv+=x0.at(0).at(gloc)*S.dvalue(1,iloc,gi); // dx0/dv
      dy0du+=x0.at(1).at(gloc)*S.dvalue(0,iloc,gi); // dy0/du
      dy0dv=1.0; // mod for direction s

// at time t

//      dxdu+=x.at(0).at(gloc)*S.dvalue(0,iloc,gi); // dx/du
//      dxdv+=x.at(0).at(gloc)*S.dvalue(1,iloc,gi); // dx/dv
//      dydu+=x.at(1).at(gloc)*S.dvalue(0,iloc,gi); // dy/du
//      dydv+=x.at(1).at(gloc)*S.dvalue(1,iloc,gi); // dy/dv

//      dxdu=1.0; // mod for direction s
//      dxdv=1.0; // mod for direction s
//      dydu+=x.at(1).at(gloc)*S.dvalue(0,iloc,gi); // dy/du
//      dydv+=x.at(1).at(gloc)*S.dvalue(1,iloc,gi); // dy/dv

//      dxdu+=x.at(0).at(gloc)*S.dvalue(0,iloc,gi); // dx/du
//      dxdv=1.0; // mod for direction s
//      dydu=1.0; // mod for direction s
//      dydv+=x.at(1).at(gloc)*S.dvalue(1,iloc,gi); // dy/dv

//      dxdu=1.0; // mod for direction s
//      dxdv=1.0; // mod for direction s
//      dydu=1.0;
//      dydv+=x.at(1).at(gloc)*S.dvalue(1,iloc,gi); // dy/dv

      dxdu=1.0; // mod for direction s
      dxdv+=x.at(0).at(gloc)*S.dvalue(1,iloc,gi); // dx/dv
      dydu+=x.at(1).at(gloc)*S.dvalue(0,iloc,gi); // dy/du
      dydv=1.0; // mod for direction s

    }

// define the jacobian for time-0, this maps to the isoparametric element from the time-0 element

    vector<vector<double> > J0{{dx0du,dx0dv},
                               {dy0du,dy0dv}};

// determinant of J0

    double detJ0(dx0du*dy0dv-dx0dv*dy0du);

// the inverse of J0 is simply the adjugate divided by the determinant

    vector<vector<double> > invJ0{{ dy0dv/detJ0,-dx0dv/detJ0},
                                  {-dy0du/detJ0, dx0du/detJ0}};

// define the jacobian for time-t, this maps to the isoparametric element from the time-t element

    vector<vector<double> > J{{dxdu,dxdv},
                              {dydu,dydv}};

// determinant of J

    double detJ(dxdu*dydv-dxdv*dydu);

// define a jacobian of the lagrangian motion, this maps to the time-0 element from the time-t element

//    vector<vector<double> > Js{{invJ0[0][0]*J[0][0]+invJ0[0][1]*J[1][0],invJ0[0][0]*J[0][1]+invJ0[0][1]*J[1][1]},
//                               {invJ0[1][0]*J[0][0]+invJ0[1][1]*J[1][0],invJ0[1][0]*J[0][1]+invJ0[1][1]*J[1][1]}};

    Js.at(gi)[0][0]=invJ0[0][0]*J[0][0]+invJ0[0][1]*J[1][0];
    Js.at(gi)[0][1]=invJ0[0][0]*J[0][1]+invJ0[0][1]*J[1][1];
    Js.at(gi)[1][0]=invJ0[1][0]*J[0][0]+invJ0[1][1]*J[1][0];
    Js.at(gi)[1][1]=invJ0[1][0]*J[0][1]+invJ0[1][1]*J[1][1];

// determinant of Js

    detJs.at(gi)=(Js.at(gi)[0][0]*Js.at(gi)[1][1]-Js.at(gi)[0][1]*Js.at(gi)[1][0]);

  }

  return;

}

// calculate a jacobian at the local nodes

void jacobian(int const &i,VVD const &x,Mesh const &M,Shape const &S,VD &detJ){

// node positions and displacement in local coordinates

  double xpos[S.nloc()],ypos[S.nloc()],disp(2.0/S.order());

  for(int jsloc=0,kloc=0;jsloc<S.sloc();jsloc++){
    for(int isloc=0;isloc<S.sloc();isloc++,kloc++){
      xpos[kloc]=-1.0+isloc*disp;
      ypos[kloc]=-1.0+jsloc*disp;
    }
  }

// loop over local nodes and calculate the jacobian

  for(int iloc=0;iloc<S.nloc();iloc++){

    double dxdu(0.0),dydu(0.0),dxdv(0.0),dydv(0.0);

// derivatives of the physical coordinates at the nodes

    for(int jloc=0;jloc<S.nloc();jloc++){
      long gloc=(S.type()==CONTINUOUS)?M.GlobalNode_CFEM(i,jloc):M.GlobalNode_DFEM(i,jloc);
      dxdu+=x.at(0).at(gloc)*S.dvalue(0,jloc,xpos[iloc],ypos[iloc]); // dx/du
      dxdv+=x.at(0).at(gloc)*S.dvalue(1,jloc,xpos[iloc],ypos[iloc]); // dx/dv
      dydu+=x.at(1).at(gloc)*S.dvalue(0,jloc,xpos[iloc],ypos[iloc]); // dy/du
      dydv+=x.at(1).at(gloc)*S.dvalue(1,jloc,xpos[iloc],ypos[iloc]); // dy/dv
    }

// calculate the determinant at the quadrature point and commit to the vector

    detJ.at(iloc)=dxdu*dydv-dxdv*dydu;

  }

  return;

}
