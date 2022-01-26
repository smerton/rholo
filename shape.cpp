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

// this function will integrate the shape function i in global coordinates
// without having to use isoparametrics or quadrature rules.

// The integral is calculated by dividing the element ABDC into
// three domains comprising a first-range triangle (I1), mid-range 
// trapezium (I2) and third-range triangle (I3).
// Each range is integrated separately and the resulting integrals
// I1,I2 and I3 and summed to obtain the total integral.
// The contribution from each component is determined using functions of x for
// the limits of the inner integration on y, which will lie along two sides of
// the cell.
// For example:

//                       (x2,y2)
//                         C
//                         o   .
//                      .  .     .
//                   .     .      o.
//                .        .      .  .
//             .           .      .    .
//          .              .      .      .
//(x0,y0).                 .      .        . (x3,y3)
//    o                    .      .         oD
//   A                     .      .        .
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
//                                B
//                             (x1,y1)

// test case
// node 0 x,y= -0.0000010400,-0.0000010400
// node 1 x,y= 0.0703333333,0.0000000000
// node 2 x,y= 0.0000000000,0.0703333333
// node 3 x,y= 0.0741409452,0.0741409452
// node 0 I1x= -0.0354838810 I1y= -0.0354838810
//       I2x= -0.0351666667 I2y= -0.0351666667
//       Gx= -0.0352599367
//       IS= 0.0012813525
//       IG= -0.1079455423

// integation limits

  double xl,xu;

// coefficients in the shape function polynomial N_i(x,y)=a+bx+cy+dxy

//  double a(coeff(i,0)),b(coeff(i,1)),c(coeff(i,2)),d(coeff(i,3));
  double a(1.0),b(0.0),c(0.0),d(0.0); // should recover the element volume

// terms in the integral of the polynomial a+bx+cy+dxy

  double t1; // axy
  double t2; // 0.5*bx^2y
  double t3; // 0.5*cy^2x
  double t4; // 0.25*dx^2y^2
  double xuxl2,xuxl3,xuxl4;     // extra terms used to make up t1,t2,t3,t4

// parameters of the two lines that make the lower and upper limits of the y integral

  double ml,mu; // gradients
  double cl,cu; // y intercepts

// coordinates of the nodes

  double x0(mr.at(0).at(0));
  double x1(mr.at(0).at(1));
  double x2(mr.at(0).at(2));
  double x3(mr.at(0).at(3));

  double y0(mr.at(1).at(0));
  double y1(mr.at(1).at(1));
  double y2(mr.at(1).at(2));
  double y3(mr.at(1).at(3));

// split the element into 3 integration domains

// acquire integration limits for I1 triangle

  if(abs(x1-x0)<=abs(x2-x0)){

// x integral is along x0<=x<=x1, set curve parameters for lower limit of y integral yl(x) = ml*(x-ol)+cl

    xl=x0;
    ml=(y1-y0)/(sgn(x1-x0)*max(1.0e-10,abs(x1-x0))); // line AB is the lower limit of the integration
    cl=y0-ml*x0;

 // set curve parameters for upper limit of y integral yu(x)=mu*(x-ou)+cu

    xu=x1;
    mu=(y2-y0)/(sgn(x2-x0)*max(1.0e-10,abs(x2-x0))); // line AC is the upper limit of the integration
    cu=y0-mu*x0;

  }else{

// x integral is along x0<=x<=x2, set curve parameters for lower limit of y integral yl(x) = ml*(x-ol)+cl

    xl=x0;
    ml=(y1-y0)/(sgn(x1-x0)*max(1.0e-10,abs(x1-x0))); // line AB is the lower limit of the integration
    cl=y0-ml*x0;

 // set curve parameters for upper limit of y integral yu(x)=mu*(x-ou)+cu

    xu=x2;
    mu=(y2-y0)/(sgn(x2-x0)*max(1.0e-10,abs(x2-x0))); // line AC is the upper limit of the integration
    cu=y0-mu*x0;

  }

// debug
  cout<<endl;
  cout<<"Shape::integrate(): coefficients:"<<endl;
  cout<<"Shape::integrate():   a= "<<a<<endl;
  cout<<"Shape::integrate():   b= "<<b<<endl;
  cout<<"Shape::integrate():   c= "<<c<<endl;
  cout<<"Shape::integrate():   d= "<<d<<endl;
  cout<<"Shape::integrate(): coords:"<<endl;
  cout<<"Shape::integrate():   x0 y0= "<<x0<<" "<<y0<<endl;
  cout<<"Shape::integrate():   x1 y2= "<<x1<<" "<<y1<<endl;
  cout<<"Shape::integrate():   x2 y3= "<<x2<<" "<<y2<<endl;
  cout<<"Shape::integrate():   x3 y4= "<<x3<<" "<<y3<<endl;
// debug

// compute I1 integral components on the range xl<x<xu, yl(x)<y<yu(x)

  xuxl2=xu*xu-xl*xl;
  xuxl3=xu*xu*xu-xl*xl*xl;
  xuxl4=xu*xu*xu*xu-xl*xl*xl*xl;

  t1=a*(xuxl2*0.5*(mu-ml)+(xu-xl)*(cu-cl));
  t2=0.5*b*(xuxl3*(mu-ml)+xuxl2*(cu-cl));
  t3=0.5*c*((xuxl3*(mu*mu-ml*ml)/3.0)+xuxl2*(mu*cu-ml*cl)+(xu-xl)*(cu*cu+cl*cl));
  t4=0.5*d*(0.25*xuxl4*mu*mu+(xuxl3*2.0*mu*cu/3.0)+0.5*xuxl2*cu*cu-0.25*xuxl4*ml*ml-(xuxl3*2.0*ml*cl/3.0)-0.5*xuxl2*cl*cl);

  cout<<"Shape::integrate(): terms in I1= "<<t1<<" "<<t2<<" "<<t3<<" "<<t4<<endl;

// contribution to element from first range triangle

  double I1(t1+t2+t3+t4);

// debug
  cout<<"Shape::integrate(): yl(x),yu(x) parameters for I1:"<<endl;
  cout<<"Shape::integrate():   ml= "<<ml<<endl;
  cout<<"Shape::integrate():   mu= "<<mu<<endl;
  cout<<"Shape::integrate():   xl= "<<xl<<endl;
  cout<<"Shape::integrate():   xu= "<<xu<<endl;
  cout<<"Shape::integrate():   cl= "<<cl<<endl;
  cout<<"Shape::integrate():   cu= "<<cu<<endl;
  cout<<"Shape::integrate(): I1= "<<I1<<endl;
// debug

// acquire integration limits for I2 trapezium

  if(abs(x3-x1)<=abs(x3-x2)){

// x integral is along x2<=x<=x1, set curve parameters for lower limit of y integral yl(x) = ml*(x-ol)+cl

    xl=x2;
    ml=(y1-y0)/(sgn(x1-x0)*max(1.0e-10,abs(x1-x0))); // line AB is the lower limit of the integration
    cl=y0-ml*x0;

 // set curve parameters for upper limit of y integral yu(x)=mu*(x-ou)+cu

    xu=x1;
    mu=(y3-y2)/(sgn(x3-x2)*max(1.0e-10,abs(x3-x2))); // line CD is the upper limit of the integration
    cu=y2-mu*x2;

  }else{

// x integral is along x1<=x<=x2, set curve parameters for lower limit of y integral yl(x) = ml*(x-ol)+cl

    xl=x1;
    ml=(y3-y1)/(sgn(x3-x1)*max(1.0e-10,abs(x3-x1))); // line BD is the lower limit of the integration
    cl=y1-ml*x1;

 // set curve parameters for upper limit of y integral yu(x)=mu*(x-ou)+cu

    xu=x2;
    mu=(y2-y0)/(sgn(x2-x0)*max(1.0e-10,abs(x2-x0))); // line AC is the upper limit of the integration
    cu=y0-mu*x0;

  }

// compute I2 integral components on the range xl<x<xu, yl(x)<y<yu(x)

  xuxl2=xu*xu-xl*xl;
  xuxl3=xu*xu*xu-xl*xl*xl;
  xuxl4=xu*xu*xu*xu-xl*xl*xl*xl;

  t1=a*(xuxl2*0.5*(mu-ml)+(xu-xl)*(cu-cl));
  t2=0.5*b*(xuxl3*(mu-ml)+xuxl2*(cu-cl));
  t3=0.5*c*((xuxl3*(mu*mu-ml*ml)/3.0)+xuxl2*(mu*cu-ml*cl)+(xu-xl)*(cu*cu+cl*cl));
  t4=0.5*d*(0.25*xuxl4*mu*mu+(xuxl3*2.0*mu*cu/3.0)+0.5*xuxl2*cu*cu-0.25*xuxl4*ml*ml-(xuxl3*2.0*ml*cl/3.0)-0.5*xuxl2*cl*cl);

  cout<<"Shape::integrate(): terms in I2= "<<t1<<" "<<t2<<" "<<t3<<" "<<t4<<endl;

// contribution to element from mid-range trapezium

  double I2(t1+t2+t3+t4);

// debug
  cout<<"Shape::integrate(): yl(x),yu(x) parameters for I2:"<<endl;
  cout<<"Shape::integrate():   ml= "<<ml<<endl;
  cout<<"Shape::integrate():   mu= "<<mu<<endl;
  cout<<"Shape::integrate():   xl= "<<xl<<endl;
  cout<<"Shape::integrate():   xu= "<<xu<<endl;
  cout<<"Shape::integrate():   cl= "<<cl<<endl;
  cout<<"Shape::integrate():   cu= "<<cu<<endl;
  cout<<"Shape::integrate(): I2= "<<I2<<endl;
// debug

// acquire integration limits for I3 triangle

  if(abs(x3-x1)<=abs(x3-x2)){

// x integral is along x1<=x<=x3, set curve parameters for lower limit of y integral yl(x) = ml*(x-ol)+cl

    xl=x1;
    ml=(y3-y1)/(sgn(x3-x1)*max(1.0e-10,abs(x3-x1))); // line BD is the lower limit of the integration
    cl=y1-ml*x1;

 // set curve parameters for upper limit of y integral yu(x)=mu*(x-ou)+cu

    xu=x3;
    mu=(y3-y2)/(sgn(x3-x2)*max(1.0e-10,abs(x3-x2))); // line CD is the upper limit of the integration
    cu=y2-mu*x2;

  }else{

// x integral is along x2<=x<=x3, set curve parameters for lower limit of y integral yl(x) = ml*(x-ol)+cl

    xl=x2;
    ml=(y3-y1)/(sgn(x3-x1)*max(1.0e-10,abs(x3-x1))); // line BD is the lower limit of the integration
    cl=y1-ml*x1;

 // set curve parameters for upper limit of y integral yu(x)=mu*(x-ou)+cu

    xu=x3;
    mu=(y3-y2)/(sgn(x3-x2)*max(1.0e-10,abs(x3-x2))); // line CD is the upper limit of the integration
    cu=y2-mu*x2;

  }

// compute I3 integral components on the range xl<x<xu, yl(x)<y<yu(x)

  xuxl2=xu*xu-xl*xl;
  xuxl3=xu*xu*xu-xl*xl*xl;
  xuxl4=xu*xu*xu*xu-xl*xl*xl*xl;

  t1=a*(xuxl2*0.5*(mu-ml)+(xu-xl)*(cu-cl));
  t2=0.5*b*(xuxl3*(mu-ml)+xuxl2*(cu-cl));
  t3=0.5*c*((xuxl3*(mu*mu-ml*ml)/3.0)+xuxl2*(mu*cu-ml*cl)+(xu-xl)*(cu*cu+cl*cl));
  t4=0.5*d*(0.25*xuxl4*mu*mu+(xuxl3*2.0*mu*cu/3.0)+0.5*xuxl2*cu*cu-0.25*xuxl4*ml*ml-(xuxl3*2.0*ml*cl/3.0)-0.5*xuxl2*cl*cl);

  cout<<"Shape::integrate(): terms in I3= "<<t1<<" "<<t2<<" "<<t3<<" "<<t4<<endl;

// contribution to element from third range triangle

  double I3(t1+t2+t3+t4);


// debug
  cout<<"Shape::integrate(): yl(x),yu(x) parameters for I3:"<<endl;
  cout<<"Shape::integrate():   ml= "<<ml<<endl;
  cout<<"Shape::integrate():   mu= "<<mu<<endl;
  cout<<"Shape::integrate():   xl= "<<xl<<endl;
  cout<<"Shape::integrate():   xu= "<<xu<<endl;
  cout<<"Shape::integrate():   cl= "<<cl<<endl;
  cout<<"Shape::integrate():   cu= "<<cu<<endl;
  cout<<"Shape::integrate(): I3= "<<I3<<endl;
// debug

// sum contribution from all 3 ranges to assemble the entire element

// debug
  cout<<"Shape::integrate():    I1+I2+I3= "<<I1+I2+I3<<endl;
//  exit(1);
// debug

  return(I1+I2+I3);

}

// kronecker delta function

double delta(int i,int j){return (i==j)?1.0:0.0;}

// type safe function to return the sign of the argument

template <typename T> int sgn(T val) {return( (val>=T(0))?T(1):T(-1));} // -1 or 1
