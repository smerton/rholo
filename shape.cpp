// Function definitions for the shape class

// Author S. R. Merton

#include <iostream>
#include <vector>
#include <algorithm> // min_element, max_element
#include "shape.h"
#include "quadrature.h"
#include "polynomial.h"
#include "matrix.h"  // matrix operations

double delta(int i,int j); // kronecker delta function
template <typename T> int sgn(T val); // return type safe sign of the argument

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

// node coordinates

  mr=x;

//  cout<<"Setting up a shape in global coordinates:"<<endl;
//  for(int iloc=0;iloc<mnloc;iloc++){
//    cout<<" node "<<iloc<<" x,y= "<<x.at(0).at(iloc)<<","<<x.at(1).at(iloc)<<endl;
//  }

// assemble a matrix equation for the coefficients

  Matrix A(mnloc),AI(mnloc);

  for(int iloc=0;iloc<mnloc;iloc++){
    A.write(iloc,0,1.0);
    A.write(iloc,1,mr.at(0).at(iloc));
    A.write(iloc,2,mr.at(1).at(iloc));
    A.write(iloc,3,mr.at(0).at(iloc)*x.at(1).at(iloc));
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
//double Shape::dvalue(int idim,int i,vector<double> x) const {double dx(mcoeff.at(1).at(i)+mcoeff.at(3).at(i)*x.at(1)),dy(mcoeff.at(2).at(i)+mcoeff.at(3).at(i)*x.at(0));return((idim==0)?dx:dy);} // shape i derivative value at global coordinate x
double Shape::coeff(int i,int j) const {return mcoeff.at(i).at(j);} // coefficient j in the polynomial expansion of shape i

// integrate shape i on the range (x1,x2),(y1,y2)

//double Shape::integrate(int i,double x1,double x2,double y1,double y2) const {double x(x2-x1),y(y2-y1);return(coeff(i,0)*x*y+0.5*coeff(i,1)*x*x*y+0.5*coeff(i,2)*y*y*x+0.5*coeff(i,3)*x*y*y);}

// integrate shape i on the range (x1,x2) (y1,y2)

double Shape::integrate(int i,double x1,double x2,double y1,double y2) const {

  double t1(coeff(i,0)*(x2-x1)*(y2-y1));
  double t2(0.5*coeff(i,1)*(x2*x2-x1*x1)*(y2-y1));
  double t3(0.5*coeff(i,2)*(y2*y2-y1*y1)*(x2-x1));
  double t4(0.25*coeff(i,3)*(x2*x2-x1*x1)*(y2*y2-y1*y1));

  return(t1+t2+t3+t4);

}

// integrate derivative i of shape j on the range (x1,x2) (y1,y2)

double Shape::integrate(int i,int j,double x1,double x2,double y1,double y2) const {

  double t1(coeff(j,1)*(x2-x1)*(y2-y1));
  double t2(0.5*coeff(j,3)*(x2-x1)*(y2*y2-y1*y1));
  double t3(coeff(j,2)*(x2-x1)*(y2-y1));
  double t4(0.5*coeff(j,3)*(x2*x2-x1*x1)*(y2-y1));
  double Idx(t1+t2);
  double Idy(t3+t4);

  return((i==0)?Idx:Idy);

}

// integrate shape i on the range r in global coordinates

double Shape::integrate(int i) const {

// this function will integrate the shape function i global coordinates
// without having to use isoparametrics

// rtmp.at(0) is a vector containing the x-coordiates of the 4 nodes in physical space
// rtmp.at(1) is a vector containing the y-coordiates of the 4 nodes in physical space

// In global coordinates the integral is calcualted by dividing the element x1,y1,x2,y2,x3,y3,x4,y4
// into three domains comprising a left triangle (I1), middle trapezium (I2) and right triangle (I3).
// The contribution from each component is determined using functions of x for
// the limits of the inner integration (on y) and summed to get total
// for example:

//                       (x2,y2)
//
//                         o   .
//                      .  .     .
//                   .     .      o.
//                .        .      .  .
//             .           .      .    .
//          .              .      .      .
//(x1,y1).                 .      .        . (x4,y4)
//    o                    .      .         o
//                         .      .        .
//       .                 .      .       .
//                I1       .  I2  .  I3  .
//                         .      .     .
//            .            .      .    .
//                         .      .   .
//                         .      .  .
//                  .      .      . .
//                         .      ..
//                         o      .
//                             .  o
//
//                             (x3,y3)

  vector<vector<double> > rtmp=mr;

  int imin,imax;
  double x1,x2,x3,x4;
  double y1,y2,y3,y4;

// acquire x1,y1

  imin=(distance(begin(rtmp.at(0)),min_element(begin(rtmp.at(0)),end(rtmp.at(0)))));

  x1=rtmp.at(0).at(imin);
  y1=rtmp.at(1).at(imin);

// remove from coordinate vector so we don't find it again when looking for x4,y4

  rtmp.at(0).erase(rtmp.at(0).begin()+imin);
  rtmp.at(1).erase(rtmp.at(1).begin()+imin);

// acquire x4,y4

  imax=(distance(begin(rtmp.at(0)),max_element(begin(rtmp.at(0)),end(rtmp.at(0)))));

  x4=rtmp.at(0).at(imax);
  y4=rtmp.at(1).at(imax);

// remove from coordinate vector so we don't find it again when looking for x2,y2

  rtmp.at(0).erase(rtmp.at(0).begin()+imax);
  rtmp.at(1).erase(rtmp.at(1).begin()+imax);

// acquire x2,y2

  imin=(distance(begin(rtmp.at(0)),min_element(begin(rtmp.at(0)),end(rtmp.at(0)))));

  x2=rtmp.at(0).at(imin);
  y2=rtmp.at(1).at(imin);

// acquire x3,y3

  imax=(distance(begin(rtmp.at(0)),max_element(begin(rtmp.at(0)),end(rtmp.at(0)))));

  x3=rtmp.at(0).at(imax);
  y3=rtmp.at(1).at(imax);

// split the element into 3 integration domains

  double m1,m2,m3,m4,m5,m6; // gradients of cell sides on each division
  double c1,c2,c3,c4,c5,c6; // y intercept for each cell side on each division
  double o1,o2,o3,o4,o5,o6; // x off-sets for each division

// set parameters for the functions that form the integration limits of the inner integral on each domain

  if(y2>y3){
     cout<<"y2>y3"<<endl;
//    m1=(y3-y1)/(x3-x1); // gradient in function f1(x) for lower integration limit
    m1=(y3-y1)/(sgn(x3-x1)*max(1.0e-10,abs(x3-x1))); // gradient in function f1(x) for lower integration limit
//    m2=(y2-y1)/(x2-x1); // gradient in function f2(x) for upper integration limit
    m2=(y2-y1)/(sgn(x2-x1)*max(1.0e-10,abs(x2-x1))); // gradient in function f2(x) for upper integration limit
//    m3=(y4-y2)/(x4-x2); // gradient in function g1(x) for lower integration limit
    m3=(y4-y2)/(sgn(x4-x2)*max(1.0e-10,abs(x4-x2))); // gradient in function g1(x) for lower integration limit
//    m4=(y4-y2)/(x4-x2); // gradient in function g2(x) for upper integration limit
    m4=(y4-y2)/(sgn(x4-x2)*max(1.0e-10,abs(x4-x2))); // gradient in function g2(x) for upper integration limit
//    m5=(y4-y3)/(x4-x3); // gradient in function h1(x) for lower integration limit
    m5=(y4-y3)/(sgn(x4-x3)*max(1.0e-10,abs(x4-x3))); // gradient in function h1(x) for lower integration limit
//    m6=(y4-y2)/(x4-x2); // gradient in function h2(x) for upper integration limit
    m6=(y4-y2)/(sgn(x4-x2)*max(1.0e-10,abs(x4-x2))); // gradient in function h2(x) for upper integration limit
    c1=y1; // intercept on y-axis in function f1(x) for lower integration limit
    c2=y1; // intercept on y-axis in function f2(x) for upper integration limit
    c3=y2; // intercept on y-axis in function g1(x) for lower integration limit
    c4=y2; // intercept on y-axis in function g2(x) for upper integration limit
    c5=y3; // intercept on y-axis in function h1(x) for lower integration limit
    c6=y2; // intercept on y-axis in function h2(x) for upper integration limit
    o1=x1; // off-set along x-axis in function f1(x) for lower integration limit
    o2=x1; // off-set along x-axis in function f2(x) for upper integration limit
    o3=x2; // off-set along x-axis in function g1(x) for lower integration limit
    o4=x2; // off-set along x-axis in function g2(x) for upper integration limit
    o5=x3; // off-set along x-axis in function h1(x) for lower integration limit
    o6=x2; // off-set along x-axis in function h2(x) for upper integration limit
  }else{
     cout<<"y2<=y3"<<endl;
//    m1=(y2-y1)/(x2-x1); // gradient in function f1(x) for lower integration limit
    m1=(y2-y1)/(sgn(x2-x1)*max(1.0e-10,abs(x2-x1))); // gradient in function f1(x) for lower integration limit
//    m2=(y3-y1)/(x3-x1); // gradient in function f2(x) for upper integration limit
    m2=(y3-y1)/(sgn(x3-x1)*max(1.0e-10,abs(x3-x1))); // gradient in function f2(x) for upper integration limit
//    m3=(y3-y1)/(x3-x1); // gradient in function g1(x) for lower integration limit
    m3=(y3-y1)/(sgn(x3-x1)*max(1.0e-10,abs(x3-x1))); // gradient in function g1(x) for lower integration limit
//    m4=(y3-y1)/(x3-x1); // gradient in function g2(x) for upper integration limit
    m4=(y3-y1)/(sgn(x3-x1)*max(1.0e-10,abs(x3-x1))); // gradient in function g2(x) for upper integration limit
//    m5=(y4-y2)/(x4-x2); // gradient in function h1(x) for lower integration limit
    m5=(y4-y2)/(sgn(x4-x2)*max(1.0e-10,abs(x4-x2))); // gradient in function h1(x) for lower integration limit
//    m6=(y4-y3)/(x4-x3); // gradient in function h2(x) for upper integration limit
    m6=(y4-y3)/(sgn(x4-x3)*max(1.0e-10,abs(x4-x3))); // gradient in function h2(x) for upper integration limit
    c1=y1; // intercept on y-axis in function f1(x) for lower integration limit
    c2=y1; // intercept on y-axis in function f2(x) for upper integration limit
    c3=y1; // intercept on y-axis in function g1(x) for lower integration limit
    c4=y1; // intercept on y-axis in function g2(x) for upper integration limit
    c5=y2; // intercept on y-axis in function h1(x) for lower integration limit
    c6=y3; // intercept on y-axis in function h2(x) for upper integration limit
    o1=x1; // off-set along x-axis in function f1(x) for lower integration limit
    o2=x1; // off-set along x-axis in function f2(x) for upper integration limit
    o3=x1; // off-set along x-axis in function g1(x) for lower integration limit
    o4=x1; // off-set along x-axis in function g2(x) for upper integration limit
    o5=x2; // off-set along x-axis in function h1(x) for lower integration limit
    o6=x3; // off-set along x-axis in function g2(x) for upper integration limit
  }
// debug
  cout<<endl;
  cout<<" x1,y1= "<<x1<<" "<<y1<<endl;
  cout<<" x2,y2= "<<x2<<" "<<y2<<endl;
  cout<<" x3,y3= "<<x3<<" "<<y3<<endl;
  cout<<" x4,y4= "<<x4<<" "<<y4<<endl;
  cout<<" m1,m2= "<<m1<<" "<<m2<<endl;
  cout<<" m3,m4= "<<m3<<" "<<m4<<endl;
  cout<<" m5,m6= "<<m5<<" "<<m6<<endl;
  cout<<" c1,c2= "<<c1<<" "<<c2<<endl;
  cout<<" c3,c4= "<<c3<<" "<<c4<<endl;
  cout<<" c5,c6= "<<c5<<" "<<c6<<endl;
  cout<<" o1,o2= "<<o1<<" "<<o2<<endl;
  cout<<" o3,o4= "<<o3<<" "<<o4<<endl;
  cout<<" o5,o6= "<<o5<<" "<<o6<<endl;
  cout<<" a0= "<<coeff(i,0)<<endl;
  cout<<" a1= "<<coeff(i,1)<<endl;
  cout<<" a2= "<<coeff(i,2)<<endl;
  cout<<" a3= "<<coeff(i,3)<<endl;



// node 0, cartesian element
x1=0.0000000000;
x2=0.0833333333;
x3=0.0000000000;
x4=0.0833333333;

y1=0.0000000000;
y2=0.0000000000;
y3=0.0833333333;
y4=0.0833333333;

m1=(y2-y1)/(sgn(x2-x1)*max(1.0e-10,abs(x2-x1)));
m2=(y4-y3)/(sgn(x4-x3)*max(1.0e-10,abs(x4-x3)));
c1=y1;
c2=y3;
o1=0.0;
o2=0.0;

  cout<<endl;
  cout<<" x1,y1 should be = "<<x1<<" "<<y1<<endl;
  cout<<" x2,y2 should be = "<<x2<<" "<<y2<<endl;
  cout<<" x3,y3 should be = "<<x3<<" "<<y3<<endl;
  cout<<" x4,y4 should be = "<<x4<<" "<<y4<<endl;
  cout<<" m1,m2 should be = "<<m1<<" "<<m2<<endl;
  cout<<" c1,c2 should be = "<<c1<<" "<<c2<<endl;
  cout<<" o1,o2 should be = "<<o1<<" "<<o2<<endl;
  cout<<" a0 should be = "<<coeff(i,0)<<endl;
  cout<<" a1 should be = "<<coeff(i,1)<<endl;
  cout<<" a2 should be = "<<coeff(i,2)<<endl;
  cout<<" a3 should be = "<<coeff(i,3)<<endl;
// debug

// evaluate integral on first range triangle (x1,y1) -> (x2,y2)

  double a(coeff(i,0)),b(coeff(i,1)),c(coeff(i,2)),d(coeff(i,3));

  double term1((x2*x2-x1*x1)*(a*m2-a*m1)+(x2-x1)*(-a*m2*o2+a*c2+a*m1*o1+a*c1));
  double term2((x2*x2*x2-x1*x1*x1)*(0.5*b*m2-0.5*b*m1)+(x2*x2-x1*x1)*(-0.5*b*m2*o2+0.5*b*c2+0.5*b*m1*o1-0.5*b*c1));
  double term31((x2*x2*x2-x1*x1*x1)*(m2*m2-m1*m1)/3.0);
  double term32((x2*x2-x1*x1)*(-m2*m2*o2+m2*c2-m1*c1-m1*m1*o1));
  double term33((x2-x1)*(m2*m2*o2*o2-2.0*m2*o2*c2+c2*c2+m1*m1*o1*o1+2.0*m1*c1*o1+c1*c1));
  double term3(0.5*c*(term31+term32+term33));
  double term41(0.25*(x2*x2*x2*x2-x1*x1*x1*x1)*m2*m2);
  double term42((x2*x2*x2-x1*x1*x1)*(2.0*m2*m2*o2+2.0*m2*c2)/3.0);
  double term43(0.5*(x2*x2-x1*x1)*(m2*m2*o2*o2+c2*c2-2.0*m2*o2*c2));
  double term44(0.25*(x2*x2*x2*x2-x1*x1*x1*x1)*m1*m1);
  double term45((x2*x2*x2-x1*x1*x1)*(2.0*m1*m1*o1+2.0*m1*c1)/3.0);
  double term46(0.5*(x2*x2-x1*x1)*(m1*m1*o1*o1+c1*c1-2.0*m1*o1*c1));
  double term4(0.5*d*(term41+term42+term43-term44-term45-term46));

// contribution to the integral from first range triangle
cout<<"term1 "<<term1<<endl;
cout<<"term2 "<<term2<<endl;
cout<<"term3 "<<term3<<endl;
cout<<"term4 "<<term4<<endl;

//term1=a*(x2-x1)*(y3-y1);
//term2=0.5*b*(x2*x2-x1*x1)*(y3-y1);
//term3=0.5*c*(x2-x1)*(y3*y3-y1*y1);
//term4=0.25*d*(x2*x2-x1*x1)*(y3*y3-y1*y1);

cout<<"term1 should be "<<term1<<endl;
cout<<"term2 should be "<<term2<<endl;
cout<<"term3 should be "<<term3<<endl;
cout<<"term4 should be "<<term4<<endl;

  double I1(term1+term2+term3+term4);

//  I1=a*(x2-x1)*(y2-y1)+b*(x2-x1)*(y2-y1)+0.5*c*(x2-x1)*(y2*y2-y1*y1)+0.25*d*(x2*x2-x1*x1)*(y2*y2-y1*y1);
  cout<<" I1= "<<I1<<endl;

// evaluate integral on middle range trapezium (x2,y2) -> (x3,x3)

  double term5((x3*x3-x2*x2)*(a*m3-a*m2)+(x3-x2)*(-a*m3*o3+a*c3+a*m2*o2+a*c2));
  double term6((x3*x3*x3-x2*x2*x2)*(0.5*b*m3-0.5*b*m2)+(x3*x3-x2*x2)*(-0.5*b*m3*o3+0.5*b*c3+0.5*b*m2*o2-0.5*b*c2));
  double term71((x3*x3*x3-x2*x2*x2)*(m3*m3-m2*m2)/3.0);
  double term72((x3*x3-x2*x2)*(-m3*m3*o3+m3*c3-m2*c2-m2*m2*o2));
  double term73((x3-x2)*(m3*m3*o3*o3-2.0*m3*o3*c3+c3*c3+m2*m2*o2*o2+2.0*m2*c2*o2+c2*c2));
  double term7(0.5*c*(term71+term72+term73));
  double term81(0.25*(x3*x3*x3*x3-x2*x2*x2*x2)*m3*m3);
  double term82((x3*x3*x3-x2*x2*x2)*(2.0*m3*m3*o3+2.0*m3*c3)/3.0);
  double term83(0.5*(x3*x3-x2*x2)*(m3*m3*o3*o3+c3*c3-2.0*m3*o3*c3));
  double term84(0.25*(x3*x3*x3*x3-x2*x2*x2*x2)*m2*m2);
  double term85((x3*x3*x3-x2*x2*x2)*(2.0*m2*m2*o2+2.0*m2*c2)/3.0);
  double term86(0.5*(x3*x3-x2*x2)*(m2*m2*o2*o2+c2*c2-2.0*m2*o2*c2));
  double term8(0.5*d*(term81+term82+term83-term84-term85-term86));

// contribution to the integral from middle range trapezium

  double I2(term5+term6+term7+term8);

// evaluate integral on third range triangle (x3,y3) -> (x4,x4)

  double term9((x4*x4-x3*x3)*(a*m4-a*m3)+(x4-x3)*(-a*m4*o4+a*c4+a*m3*o3+a*c3));
  double term10((x4*x4*x4-x3*x3*x3)*(0.5*b*m4-0.5*b*m3)+(x4*x4-x3*x3)*(-0.5*b*m4*o4+0.5*b*c4+0.5*b*m3*o3-0.5*b*c3));
  double term111((x4*x4*x4-x3*x3*x3)*(m4*m4-m3*m3)/3.0);
  double term112((x4*x4-x3*x3)*(-m4*m4*o4+m4*c4-m3*c3-m3*m3*o3));
  double term113((x4-x3)*(m4*m4*o4*o4-2.0*m4*o4*c4+c4*c4+m3*m3*o3*o3+2.0*m3*c3*o3+c3*c3));
  double term11(0.5*c*(term111+term112+term113));
  double term121(0.25*(x4*x4*x4*x4-x3*x3*x3*x3)*m4*m4);
  double term122((x4*x4*x4-x3*x3*x3)*(2.0*m4*m4*o4+2.0*m4*c4)/3.0);
  double term123(0.5*(x4*x4-x3*x3)*(m4*m4*o4*o4+c4*c4-2.0*m4*o4*c4));
  double term124(0.25*(x4*x4*x4*x4-x3*x3*x3*x3)*m3*m3);
  double term125((x4*x4*x4-x3*x3*x3)*(2.0*m3*m3*o3+2.0*m3*c3)/3.0);
  double term126(0.5*(x4*x4-x3*x3)*(m3*m3*o3*o3+c3*c3-2.0*m3*o3*c3));
  double term12(0.5*d*(term121+term122+term123-term124-term125-term126));

// contribution to the integral from third range tiangle

  double I3(term9+term10+term11+term12);

// sum contribution from all 3 ranges to assemble the entire element

  return(I1+I2+I3);

}

// kronecker delta function

double delta(int i,int j){return (i==j)?1.0:0.0;}

// type safe function to return the sign of the argument

template <typename T> int sgn(T val) {return( (val>=T(0))?T(1):T(-1));} // -1 or 1
