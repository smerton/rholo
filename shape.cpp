// Function definitions for the shape class

// Author S. R. Merton

#include <iostream>
#include <vector>
#include "shape.h"
#include "quadrature.h"
#include "polynomial.h"

using namespace std;

// constructor to instantiate a new shape object and all its attributes

Shape::Shape(int n){

// set the polyhedral order

  morder=n;

// set number of local nodes

  mnloc=this->order()+1;

// set number of surface nodes on a face

  msloc=1;

// set number of integration points adequate for the shape

  mngi=this->order()+1;

// set number of faces

  mnfaces=2;

// set up reflection across each face

  vector<int> reflections[this->nfaces()];

  for(int i=0;i<this->nloc();i++){
    reflections[0].push_back(this->nloc()-i-1);
    reflections[1].push_back(this->nloc()-i-1);
  }

  for(int iface=0;iface<this->nfaces();iface++){
    mreflect.push_back(reflections[iface]);
  }

// load a quadrature rule for the numerical integration of the polynomials Px

  QuadratureRule Q(this->order()+1);

// shape value at each integration point

  for(int i=0;i<this->nloc();i++){
    vector<double> ivalue,idvalue;
    for(int gi=0;gi<this->ngi();gi++){
      ivalue.push_back(this->value(i,Q.x[gi]));
      idvalue.push_back(this->dvalue(i,Q.x[gi]));
    }
    mvalue.push_back(ivalue);
    mdvalue.push_back(idvalue);
  }

// quadrature weights

  for(int gi=0;gi<this->ngi();gi++){
    mwgt.push_back(Q.w[gi]);
  }

}

// destructor function to release storage associated with a Shape class object

Shape::~Shape(){}

// accessor functions to member data

int Shape::order(){return morder;} // returns the polyhedral order
int Shape::nloc(){return mnloc;} // returns number of local nodes
int Shape::sloc(){return msloc;} // returns number of nodes on the surface
int Shape::ngi(){return mngi;} // number of Gauss integration points
int Shape::nfaces(){return mnfaces;} // number of element faces
int Shape::reflect(int iloc,int iface){return mreflect[iface][iloc];} // refelct iloc across face iface

double Shape::wgt(int gi){return mwgt[gi];} // quadrature weight of integration point gi
double Shape::value(int i,int gi){return mvalue[i][gi];} // shape i value at Gauss point gi
double Shape::value(int i,double x){Polynomial P(this->order(),i,-1.0,1.0);return P.value(x);} // shape i value at coordinate x
double Shape::dvalue(int i,int gi){return mdvalue[i][gi];} // derivative i value at Gauss point gi
double Shape::dvalue(int i,double x){Polynomial P(this->order(),i,-1.0,1.0);return P.dvalue(x);} // shape i derivative value at coordinate x
