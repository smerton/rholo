// Function definitions for the shape class

// Author S. R. Merton

#include <iostream>
#include <vector>
#include "shape.h"
#include "quadrature.h"
#include "polynomial.h"
#include "matrix.h"  // matrix operations

double delta(int i,int j); // kronecker delta function

using namespace std;

// constructor to instantiate a new shape object and all its attributes in local coordinates

Shape::Shape(int n){

// set the polyhedral order

  morder=n;

// set the number of dimensions

  mndims=2; // shouldn't we be getting this from the mesh class ?

// set number of local nodes

  mnloc=(this->order()+1)*(this->order()+1);

// set number of surface nodes on a face

  msloc=this->order()+1;

// set number of integration points adequate for the shape

  mngi=(this->order()+1)*(this->order()+1);

// set number of faces

  mnfaces=4;

// set up reflection across each face

  vector<int> reflections[this->nfaces()];

  reflections[0].push_back(2);
  reflections[1].push_back(1);
  reflections[2].push_back(2);
  reflections[3].push_back(3);

  for(int iface=0;iface<this->nfaces();iface++){
    mreflect.push_back(reflections[iface]);
  }

// node positions

  for(int j=0,k=0;j<sloc();j++){
    for(int i=0;i<sloc();i++,k++){
      mpos[0].push_back(i);
      mpos[1].push_back(j);
    }
  }

// load a quadrature rule for the numerical integration of the polynomials Px

  QuadratureRule Q(this->order()+1);

// shape and derivative value at each integration point

  for(int i=0;i<nloc();i++){
    vector<double> values;
    vector<vector<double> > dvalues(mndims);
    for(int jpt=0;jpt<Q.npoints;jpt++){
      for(int ipt=0;ipt<Q.npoints;ipt++){
        values.push_back(value(i,Q.x[ipt],Q.x[jpt]));
        for(int idim=0;idim<mndims;idim++){
          dvalues.at(idim).push_back(dvalue(idim,i,Q.x[ipt],Q.x[jpt]));
        }
      }
    }
    mvalue.push_back(values);

    for(int idim=0;idim<ndims();idim++){
      mdvalue[idim].push_back(dvalues.at(idim));
    }

  }

// quadrature weights

  for(int ipt=0,gi=0;ipt<Q.npoints;ipt++){
    for(int jpt=0;jpt<Q.npoints;jpt++,gi++){
      mwgt.push_back(Q.w[ipt]*Q.w[jpt]);
    }
  }

}

// constructor to instantiate a new shape object and all its attributes in global coordinates

Shape::Shape(int n,vector<vector<double> > x){

// this is only coded for bilinear case

  if(n!=1){
    cout<<"Shape::Shape(): Only bilinear shapes can be constructed in global coordinates, n= "<<n<<endl;
    cout<<"Shape::Shape(): Stopping"<<endl;
    exit(1);
  }

// set the polyhedral order

  morder=n;

// set the number of dimensions

  mndims=2; // shouldn't we be getting this from the mesh class ?

// set number of local nodes

  mnloc=(this->order()+1)*(this->order()+1);

// set number of surface nodes on a face

  msloc=this->order()+1;

// set number of integration points adequate for the shape

  mngi=(this->order()+1)*(this->order()+1);

// set number of faces

  mnfaces=4;

// set up reflection across each face

  vector<int> reflections[this->nfaces()];

  reflections[0].push_back(2);
  reflections[1].push_back(1);
  reflections[2].push_back(2);
  reflections[3].push_back(3);

  for(int iface=0;iface<this->nfaces();iface++){
    mreflect.push_back(reflections[iface]);
  }

// node positions

  for(int j=0,k=0;j<sloc();j++){
    for(int i=0;i<sloc();i++,k++){
      mpos[0].push_back(i);
      mpos[1].push_back(j);
    }
  }

  cout<<"Setting up a shape in global coordinates:"<<endl;
  for(int iloc=0;iloc<mnloc;iloc++){
    cout<<" node "<<iloc<<" x,y= "<<x.at(0).at(iloc)<<","<<x.at(1).at(iloc)<<endl;
  }

// assemble a matrix equation for the coefficients

  Matrix A(mnloc),AI(mnloc);

  for(int iloc=0;iloc<mnloc;iloc++){
    A.write(iloc,0,1.0);
    A.write(iloc,1,x.at(0).at(iloc));
    A.write(iloc,2,x.at(1).at(iloc));
    A.write(iloc,3,x.at(0).at(iloc)*x.at(1).at(iloc));
  }

// check for singularities which are expected if x1=x3 and y1=y4

//  if(A.singular()){
//    cout<<"Shape::Shape(): Element has a singularity, unable to construct."<<endl;
//    cout<<"Shape::Shape(): Stopping."<<endl;
//    exit(1);
//  }

// solve to get the coefficients

  for(int iloc=0;iloc<mnloc;iloc++){

// rhs

    double b[mnloc];

    for(int jloc=0;jloc<mnloc;jloc++){
      b[jloc]=delta(iloc,jloc);
    }

// inverse

    AI.inverse2(&A); // lapack drivers dgetrf_ and dgetri_

// solve for the coefficient for shape function iloc

    vector<double> coeffs_iloc;
    for(int jloc=0;jloc<mnloc;jloc++){
      double xsoln(0.0);
      for(int kloc=0;kloc<mnloc;kloc++){
        xsoln+=AI.read(jloc,kloc)*b[kloc];
      }
      coeffs_iloc.push_back(xsoln);
    }

// store coefficients for shape function iloc

    mcoeff.push_back(coeffs_iloc);

  }

}

// destructor function to release storage associated with a Shape class object

Shape::~Shape(){}

/// member function to prolongate the vector u[] to a new vector v[] in the destination element

void Shape::prolongate(double*u,double*v,int p){

// create the destination element

  Shape M(p);double P[this->nloc()][M.nloc()]={};

// form a matrix to map u to M

  for(int i=0;i<M.nloc();i++){
    v[i]=0.0;
    for(int j=0;j<this->nloc();j++){
      for(int gi=0;gi<M.ngi();gi++){
        P[i][j]+=M.value(i,gi)*this->value(j,gi)*M.wgt(gi);
      }
      v[i]+=P[i][j]*u[j];
    }
  }

  return;

}

// accessor functions to member data

int Shape::order() const {return morder;} // returns the polyhedral order
int Shape::ndims() const {return mndims;} // returns the number of dimensions
int Shape::nloc() const {return mnloc;} // returns number of local nodes
int Shape::sloc() const {return msloc;} // returns number of nodes on the surface
int Shape::ngi() const {return mngi;} // number of Gauss integration points
int Shape::nfaces() const {return mnfaces;} // number of element faces
int Shape::reflect(int iloc,int iface) const {return mreflect[iface][iloc];} // refelct iloc across face iface
int Shape::pos(int idim,int iloc) const {return mpos[idim][iloc];} // node position in dimension idim

double Shape::wgt(int gi) const {return mwgt[gi];} // quadrature weight of integration point gi
double Shape::value(int i,int gi) const {return mvalue[i][gi];} // shape i value at Gauss point gi
double Shape::value(int i,double u,double v) const {Polynomial P1(order(),pos(0,i),-1.0,1.0),P2(order(),pos(1,i),-1.0,1.0);return P1.value(u)*P2.value(v);} // shape i value at coordinate x,y
double Shape::value(int i,vector<double> x) const {return(mcoeff.at(i).at(0)+mcoeff.at(i).at(1)*x.at(0)+mcoeff.at(i).at(2)*x.at(1)+mcoeff.at(i).at(3)*x.at(0)*x.at(1));}
double Shape::dvalue(int idim,int i,int gi) const {return mdvalue[idim][i][gi];} // derivative i value at Gauss point gi
double Shape::dvalue(int idim,int i,double u,double v) const {Polynomial P1(order(),pos(0,i),-1.0,1.0),P2(order(),pos(1,i),-1.0,1.0);double dval[2];dval[0]=P1.dvalue(u)*P2.value(v);dval[1]=P1.value(u)*P2.dvalue(v);return dval[idim];}// shape i derivative value at coordinate x
double Shape::dvalue(int idim,int i,vector<double> x) const {double dx(mcoeff.at(i).at(1)+mcoeff.at(i).at(3)*x.at(1)),dy(mcoeff.at(i).at(2)+mcoeff.at(i).at(3)*x.at(0));return((idim==0)?dx:dy);} // shape i derivative value at global coordinate x

// kronecker delta function

double delta(int i,int j){return (i==j)?1.0:0.0;}
