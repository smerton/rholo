// Function definitions for the shape class

// Author S. R. Merton

#include <iostream>
#include <vector>
#include <algorithm> // min_element, max_element
#include "shape.h"
#include "quadrature.h"
#include "polynomial.h"
#include "matrix.h"  // matrix operations

int factorial(int n); // factorial function
double delta(int i,int j); // kronecker delta function
template <typename T> int sgn(T val); // return type safe sign of the argument

using namespace std;

// constructor to instantiate a new shape object and all its attributes in local coordinates

Shape::Shape(int n){

// set the polyhedral order

  morder=n;

// set the type of shape, assume CONITNUOUS

  mtype=CONTINUOUS;

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

  mr.resize(ndims());

  double dx(2.0/morder),dy(2.0/morder);
  for(int isuby=0,k=0;isuby<morder+1;isuby++){
    for(int isubx=0;isubx<morder+1;isubx++,k++){
      mr.at(0).push_back(-1.0+isubx*dx);
      mr.at(1).push_back(-1.0+isuby*dy);
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

// constructor to instantiate a new shape object of type t using ngi integration point and to set all its attributes in local coordinates

Shape::Shape(int n,int ngi,int t){

// set the polyhedral order

  morder=n;

// set the type of shape, assume CONITNUOUS

  mtype=t;

// set the number of dimensions

  mndims=2; // shouldn't we be getting this from the mesh class ?

// set number of local nodes

  mnloc=(this->order()+1)*(this->order()+1);

// set number of surface nodes on a face

  msloc=this->order()+1;

// set number of integration points adequate for the shape

  int ngi_min((this->order()+1)); // min number needed

  if(ngi<ngi_min){
    cout<<"Shape::Shape(): "<<ngi_min<<" quadrature points needed, "<<ngi<<" requested."<<endl;
    exit(1);
  }

// ngi is the number of points on-axis so we square it to obtain correct number for a 2D element

  mngi=ngi*ngi;

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

  QuadratureRule Q(ngi);

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

// loop over powers of x and y for each column of the matrix

  for(int iloc=0;iloc<mnloc;iloc++){
    for(int i=0,jloc=0;i<=morder;i++){
      for(int j=0;j<=morder;j++,jloc++){
        A.write(iloc,jloc,pow(mr.at(0).at(iloc),i)*pow(mr.at(1).at(iloc),j));
      }
    }

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

// member function to prolongate the vector u[] of nodal values to a new vector v[] in a destination element of order p
// this will map nodal data in vector u from the "this->" element to a new element of order p

void Shape::prolongate(double*u,double*v,int p) const {

// create the destination element of polyhedral order p

  Shape M(p,sqrt(this->ngi()),CONTINUOUS);

// form a prolongation operator to map u[] onto destination stencil M

  Matrix P(M.nloc(),this->nloc());

  for(int iloc=0;iloc<M.nloc();iloc++){
    for(int jloc=0;jloc<this->nloc();jloc++){
      double pij(0.0);
      for(int gi=0;gi<M.ngi();gi++){
        pij+=M.value(iloc,gi)*this->value(jloc,gi)*M.wgt(gi);
      }
      P.write(iloc,jloc,pij);
    }
  }

// form a mass matrix and its inverse on destination stencil M

  Matrix MMASS(M.nloc()),MMASSI(M.nloc());

  for(int iloc=0;iloc<M.nloc();iloc++){
    for(int jloc=0;jloc<M.nloc();jloc++){
      double nmass(0.0);
      for(int gi=0;gi<M.ngi();gi++){
        nmass+=M.value(iloc,gi)*M.value(jloc,gi)*M.wgt(gi);
      }
      MMASS.write(iloc,jloc,nmass);
    }
  }

// invert the stencil M mass matrix

  MMASSI.inverse2(&MMASS);

// normalise prolongation operator to respect Galerkin orthogonality

  Matrix PNORM(M.nloc(),this->nloc());

  PNORM.product(&MMASSI,&P);

// prolongate the source vector u[] of length this>nloc() to destination vector v[] which will be of length M.nloc()

  for(int iloc=0;iloc<M.nloc();iloc++){
    v[iloc]=0.0;
    for(int jloc=0;jloc<this->nloc();jloc++){
      v[iloc]+=PNORM.read(iloc,jloc)*u[jloc];
    }
  }

  return;

}

// accessor functions to member data

int Shape::order() const {return morder;} // returns the polyhedral order
int Shape::type() const {return mtype;} // returns the shape type as either CONTINUOUS or DISCONTINUOUS
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

// scatter values at the integration points to the nodes
// example: vector<double> arr_iloc(S.values(arr_gi));

vector<double> Shape::values(vector<double> const &v) const{

// no. integrations must match the no. nodes as we need a square matrix to invert

  if(nloc()!=ngi()){
    cout<<"Shape::node_values(): Number of nodes ("<<nloc()<<") does not match the number of integration points("<<ngi()<<"), stopping."<<endl;
    exit(1);
  }

  vector<double> u(v);

// store shape value at integration points in a square matrix

  Matrix NMAT(nloc()),NMATI(nloc());
  for(int gi=0;gi<ngi();gi++){
    for(int iloc=0;iloc<nloc();iloc++){
      NMAT.write(gi,iloc,value(iloc,gi));
    }
  }

// invert the matrix

  NMATI.inverse2(&NMAT);

// perform a mat-vec to recover the nodal values

  for(int iloc=0;iloc<nloc();iloc++){
    u.at(iloc)=0.0;
    for(int gi=0;gi<ngi();gi++){
      u.at(iloc)+=v.at(gi)*NMATI.read(iloc,gi);
    }
  }

  return u;

}

//double Shape::value(int i,vector<double> x) const {return(mcoeff.at(i).at(0)+mcoeff.at(i).at(1)*x.at(0)+mcoeff.at(i).at(2)*x.at(1)+mcoeff.at(i).at(3)*x.at(0)*x.at(1));}

double Shape::value(int iloc,vector<double> x) const {

  double pval(0.0);

  for(int i=0,k=0;i<=morder;i++){
    for(int j=0;j<=morder;j++,k++){
      pval+=mcoeff.at(iloc).at(k)*pow(x.at(0),i)*pow(x.at(1),j);
    }
  }

  return pval;

}

double Shape::dvalue(int idim,int i,int gi) const {return mdvalue[idim][i][gi];} // derivative i value at Gauss point gi
double Shape::dvalue(int idim,int i,double u,double v) const {Polynomial P1(order(),pos(0,i),-1.0,1.0),P2(order(),pos(1,i),-1.0,1.0);double dval[2];dval[0]=P1.dvalue(u)*P2.value(v);dval[1]=P1.value(u)*P2.dvalue(v);return dval[idim];}// shape i derivative value at coordinate x
double Shape::dvalue(int idim,int i,vector<double> x) const {double dx(mcoeff.at(i).at(1)+mcoeff.at(i).at(3)*x.at(1)),dy(mcoeff.at(i).at(2)+mcoeff.at(i).at(3)*x.at(0));return((idim==0)?dx:dy);} // shape i derivative value at global coordinate x
//double Shape::dvalue(int idim,int i,vector<double> x) const {double dx(mcoeff.at(1).at(i)+mcoeff.at(3).at(i)*x.at(1)),dy(mcoeff.at(2).at(i)+mcoeff.at(3).at(i)*x.at(0));return((idim==0)?dx:dy);} // shape i derivative value at global coordinate x
double Shape::coeff(int i,int j) const {return mcoeff.at(i).at(j);} // coefficient j in the polynomial expansion of shape i

// integrate derivative idim of shape j in global coordinates

double Shape::integrate(int idim,int i) const {

// integation limits

  double xl,xu;

// coefficients in the shape function polynomial N_i(x,y)=a+bx+cy+dxy

  double a(coeff(i,0)),b(coeff(i,1)),c(coeff(i,2)),d(coeff(i,3));

// integrals of the derivatives

  double Ix(0.0),Iy(0.0),I1(0.0),I2(0.0),I3(0.0);

// terms in the integral of the derivatives of the polynomial a+bx+cy+dxy

  double t1; // axy
  double t2; // 0.5*bx^2y
  double t3; // 0.5*cy^2x
  double t4; // 0.25*dx^2y^2
  double xxxu,xxu;                        // extra terms used when inserting the upper integration limits
  double xxxl,xxl;                        // extra terms used when inserting the lower integration limits
  double xuxl2,xuxl3;                     // extra terms used to make up t1,t2,t3,t4
  double tt1,tt2,tt3,tt4,tt5,tt6,tt7,tt8; // extra terms used to make up t1,t2,t3,t4

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
//  if(idim==1){
//    cout<<endl;
//    cout<<"Shape::integrate(): coefficients:"<<endl;
//    cout<<"Shape::integrate():   a= "<<a<<endl;
//    cout<<"Shape::integrate():   b= "<<b<<endl;
//    cout<<"Shape::integrate():   c= "<<c<<endl;
//    cout<<"Shape::integrate():   d= "<<d<<endl;
//    cout<<"Shape::integrate(): coords:"<<endl;
//    cout<<"Shape::integrate():   x0 y0= "<<x0<<" "<<y0<<endl;
//    cout<<"Shape::integrate():   x1 y2= "<<x1<<" "<<y1<<endl;
//    cout<<"Shape::integrate():   x2 y3= "<<x2<<" "<<y2<<endl;
//    cout<<"Shape::integrate():   x3 y4= "<<x3<<" "<<y3<<endl;
//  }
// debug

// compute I1 integral components on the range xl<x<xu, yl(x)<y<yu(x)

  xxu=xu*xu;
  xxxu=xxu*xu;
  xxl=xl*xl;
  xxxl=xxl*xl;

  xuxl2=xxu-xxl;
  xuxl3=xxxu-xxxl;

// I1 contributions to the derivatives

  tt1=0.5*xuxl2*(mu-ml);
  tt2=(xu-xl)*(cu-cl);
  tt3=xuxl3*mu*mu/3.0;
  tt4=xuxl2*mu*cu;
  tt5=(xu-xl)*cu*cu;
  tt6=xuxl3*ml*ml/3.0;
  tt7=xuxl2*ml*cl;
  tt8=(xu-xl)*cl*cl;

// sum in I1 contribution to the x derivative

  Ix+=b*(tt1+tt2);
  Ix+=0.5*d*(tt3+tt4+tt5-tt6-tt7-tt8);

  tt1=0.5*xuxl2*(mu-ml);
  tt2=(xu-xl)*(cu-cl);
  tt3=xuxl3*(mu-ml)/3.0;
  tt4=0.5*xuxl2*(cu-cl);

// sum in I1 contribution to the y derivative

  I1=c*(tt1+tt2)+d*(tt3+tt4);
  Iy+=c*(tt1+tt2);
  Iy+=d*(tt3+tt4);

//  if(idim==1){
//    cout<<"Shape::integrate(): terms in I1= "<<tt1<<" "<<tt2<<" "<<tt3<<" "<<tt4<<" "<<tt5<<" "<<tt6<<" "<<tt7<<" "<<tt8<<endl;
//  }

// debug
//  if(idim==1){
//    cout<<"Shape::integrate(): yl(x),yu(x) parameters for I1:"<<endl;
//    cout<<"Shape::integrate():   ml= "<<ml<<endl;
//    cout<<"Shape::integrate():   mu= "<<mu<<endl;
//    cout<<"Shape::integrate():   xl= "<<xl<<endl;
//    cout<<"Shape::integrate():   xu= "<<xu<<endl;
//    cout<<"Shape::integrate():   cl= "<<cl<<endl;
//    cout<<"Shape::integrate():   cu= "<<cu<<endl;
//    cout<<"Shape::integrate(): Iy(I1)= "<<I1<<endl;
//  }
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

  xxu=xu*xu;
  xxxu=xxu*xu;
  xxl=xl*xl;
  xxxl=xxl*xl;

  xuxl2=xxu-xxl;
  xuxl3=xxxu-xxxl;

// I2 contributions to the derivatives

  tt1=0.5*xuxl2*(mu-ml);
  tt2=(xu-xl)*(cu-cl);
  tt3=xuxl3*mu*mu/3.0;
  tt4=xuxl2*mu*cu;
  tt5=(xu-xl)*cu*cu;
  tt6=xuxl3*ml*ml/3.0;
  tt7=xuxl2*ml*cl;
  tt8=(xu-xl)*cl*cl;

// sum in I2 contribution to the x derivative

  Ix+=b*(tt1+tt2);
  Ix+=0.5*d*(tt3+tt4+tt5-tt6-tt7-tt8);

  tt1=0.5*xuxl2*(mu-ml);
  tt2=(xu-xl)*(cu-cl);
  tt3=xuxl3*(mu-ml)/3.0;
  tt4=0.5*xuxl2*(cu-cl);

// sum in I2 contribution to the y derivative

  I2=c*(tt1+tt2)+d*(tt3+tt4);
  Iy+=c*(tt1+tt2);
  Iy+=d*(tt3+tt4);

//  if(idim==1){
//    cout<<"Shape::integrate(): terms in I2= "<<tt1<<" "<<tt2<<" "<<tt3<<" "<<tt4<<" "<<tt5<<" "<<tt6<<" "<<tt7<<" "<<tt8<<endl;
//  }

// debug
//  if(idim==1){
//    cout<<"Shape::integrate(): yl(x),yu(x) parameters for I2:"<<endl;
//    cout<<"Shape::integrate():   ml= "<<ml<<endl;
//    cout<<"Shape::integrate():   mu= "<<mu<<endl;
//    cout<<"Shape::integrate():   xl= "<<xl<<endl;
//    cout<<"Shape::integrate():   xu= "<<xu<<endl;
//    cout<<"Shape::integrate():   cl= "<<cl<<endl;
//    cout<<"Shape::integrate():   cu= "<<cu<<endl;
//    cout<<"Shape::integrate(): Iy(I2)= "<<Iy<<endl;
//  }
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

  xxu=xu*xu;
  xxxu=xxu*xu;
  xxl=xl*xl;
  xxxl=xxl*xl;

  xuxl2=xxu-xxl;
  xuxl3=xxxu-xxxl;

// I3 contributions to the derivatives

  tt1=0.5*xuxl2*(mu-ml);
  tt2=(xu-xl)*(cu-cl);
  tt3=xuxl3*mu*mu/3.0;
  tt4=xuxl2*mu*cu;
  tt5=(xu-xl)*cu*cu;
  tt6=xuxl3*ml*ml/3.0;
  tt7=xuxl2*ml*cl;
  tt8=(xu-xl)*cl*cl;

// sum in I3 contribution to the x derivative

  Ix+=b*(tt1+tt2);
  Ix+=0.5*d*(tt3+tt4+tt5-tt6-tt7-tt8);

  tt1=0.5*xuxl2*(mu-ml);
  tt2=(xu-xl)*(cu-cl);
  tt3=xuxl3*(mu-ml)/3.0;
  tt4=0.5*xuxl2*(cu-cl);

// sum in I3 contribution to the y derivative

  I3=c*(tt1+tt2)+d*(tt3+tt4);
  Iy+=c*(tt1+tt2);
  Iy+=d*(tt3+tt4);

//  if(idim==1){
//    cout<<"Shape::integrate(): terms in I3= "<<tt1<<" "<<tt2<<" "<<tt3<<" "<<tt4<<" "<<tt5<<" "<<tt6<<" "<<tt7<<" "<<tt8<<endl;
//  }

// debug
//  if(idim==1){
//    cout<<"Shape::integrate(): yl(x),yu(x) parameters for I3:"<<endl;
//    cout<<"Shape::integrate():   ml= "<<ml<<endl;
//    cout<<"Shape::integrate():   mu= "<<mu<<endl;
//    cout<<"Shape::integrate():   xl= "<<xl<<endl;
//    cout<<"Shape::integrate():   xu= "<<xu<<endl;
//    cout<<"Shape::integrate():   cl= "<<cl<<endl;
//    cout<<"Shape::integrate():   cu= "<<cu<<endl;
//    cout<<"Shape::integrate(): Iy(I3)= "<<I3<<endl;
//  }
// debug

  return((idim==0)?Ix:Iy);

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

  double a(coeff(i,0)),b(coeff(i,1)),c(coeff(i,2)),d(coeff(i,3));
//  double a(1.0),b(1.0),c(1.0),d(1.0); // should recover the element volume

// terms in the integral of the polynomial a+bx+cy+dxy

  double t1; // axy
  double t2; // 0.5*bx^2y
  double t3; // 0.5*cy^2x
  double t4; // 0.25*dx^2y^2
  double xxxxu,xxxu,xxu;            // extra terms used when inserting the upper integration limits
  double xxxxl,xxxl,xxl;            // extra terms used when inserting the lower integration limits
  double xuxl2,xuxl3,xuxl4;             // extra terms used to make up t1,t2,t3,t4
  double tt1,tt2,tt3,tt4,tt5,tt6,tt7,tt8,tt9,tt10,tt11,tt12,tt13; // extra terms used to make up t1,t2,t3,t4

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
//  cout<<endl;
//  cout<<"Shape::integrate(): coefficients:"<<endl;
//  cout<<"Shape::integrate():   a= "<<a<<endl;
//  cout<<"Shape::integrate():   b= "<<b<<endl;
//  cout<<"Shape::integrate():   c= "<<c<<endl;
//  cout<<"Shape::integrate():   d= "<<d<<endl;
//  cout<<"Shape::integrate(): coords:"<<endl;
//  cout<<"Shape::integrate():   x0 y0= "<<x0<<" "<<y0<<endl;
//  cout<<"Shape::integrate():   x1 y2= "<<x1<<" "<<y1<<endl;
//  cout<<"Shape::integrate():   x2 y3= "<<x2<<" "<<y2<<endl;
//  cout<<"Shape::integrate():   x3 y4= "<<x3<<" "<<y3<<endl;
// debug

// compute I1 integral components on the range xl<x<xu, yl(x)<y<yu(x)

  xxu=xu*xu;
  xxxu=xxu*xu;
  xxxxu=xxxu*xu;
  xxl=xl*xl;
  xxxl=xxl*xl;
  xxxxl=xxxl*xl;

  xuxl2=xxu-xxl;
  xuxl3=xxxu-xxxl;
  xuxl4=xxxxu-xxxxl;

  tt1=xuxl2*0.5*(mu-ml);
  tt2=(xu-xl)*(cu-cl);
  tt3=xuxl3*((mu-ml)/3.0);
  tt4=xuxl2*0.5*(cu-cl);
  tt5=(xuxl3*(mu*mu-ml*ml)/3.0);
  tt6=xuxl2*(mu*cu-ml*cl);
  tt7=(xu-xl)*(cu*cu-cl*cl);
  tt8=0.25*xuxl4*mu*mu;
  tt9=(xuxl3*2.0*mu*cu/3.0);
  tt10=0.5*xuxl2*cu*cu;
  tt11=0.25*xuxl4*ml*ml;
  tt12=(xuxl3*2.0*ml*cl/3.0);
  tt13=0.5*xuxl2*cl*cl;

  t1=a*(tt1+tt2);
  t2=b*(tt3+tt4);
  t3=0.5*c*(tt5+tt6+tt7);
  t4=0.5*d*(tt8+tt9+tt10-tt11-tt12-tt13);

//  cout<<"Shape::integrate(): terms in I1= "<<t1<<" "<<t2<<" "<<t3<<" "<<t4<<endl;

// contribution to element from first range triangle

  double I1(t1+t2+t3+t4);

// debug
//  cout<<"Shape::integrate(): yl(x),yu(x) parameters for I1:"<<endl;
//  cout<<"Shape::integrate():   ml= "<<ml<<endl;
//  cout<<"Shape::integrate():   mu= "<<mu<<endl;
//  cout<<"Shape::integrate():   xl= "<<xl<<endl;
//  cout<<"Shape::integrate():   xu= "<<xu<<endl;
//  cout<<"Shape::integrate():   cl= "<<cl<<endl;
//  cout<<"Shape::integrate():   cu= "<<cu<<endl;
//  cout<<"Shape::integrate(): I1= "<<I1<<endl;
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

  xxu=xu*xu;
  xxxu=xxu*xu;
  xxxxu=xxxu*xu;
  xxl=xl*xl;
  xxxl=xxl*xl;
  xxxxl=xxxl*xl;

  xuxl2=xxu-xxl;
  xuxl3=xxxu-xxxl;
  xuxl4=xxxxu-xxxxl;

  tt1=xuxl2*0.5*(mu-ml);
  tt2=(xu-xl)*(cu-cl);
  tt3=xuxl3*((mu-ml)/3.0);
  tt4=xuxl2*0.5*(cu-cl);
  tt5=(xuxl3*(mu*mu-ml*ml)/3.0);
  tt6=xuxl2*(mu*cu-ml*cl);
  tt7=(xu-xl)*(cu*cu-cl*cl);
  tt8=0.25*xuxl4*mu*mu;
  tt9=(xuxl3*2.0*mu*cu/3.0);
  tt10=0.5*xuxl2*cu*cu;
  tt11=0.25*xuxl4*ml*ml;
  tt12=(xuxl3*2.0*ml*cl/3.0);
  tt13=0.5*xuxl2*cl*cl;

  t1=a*(tt1+tt2);
  t2=b*(tt3+tt4);
  t3=0.5*c*(tt5+tt6+tt7);
  t4=0.5*d*(tt8+tt9+tt10-tt11-tt12-tt13);

//  cout<<"Shape::integrate(): terms in I2= "<<t1<<" "<<t2<<" "<<t3<<" "<<t4<<endl;

// contribution to element from mid-range trapezium

  double I2(t1+t2+t3+t4);

// debug
//  cout<<"Shape::integrate(): yl(x),yu(x) parameters for I2:"<<endl;
//  cout<<"Shape::integrate():   ml= "<<ml<<endl;
//  cout<<"Shape::integrate():   mu= "<<mu<<endl;
//  cout<<"Shape::integrate():   xl= "<<xl<<endl;
//  cout<<"Shape::integrate():   xu= "<<xu<<endl;
//  cout<<"Shape::integrate():   cl= "<<cl<<endl;
//  cout<<"Shape::integrate():   cu= "<<cu<<endl;
//  cout<<"Shape::integrate(): I2= "<<I2<<endl;
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

  xxu=xu*xu;
  xxxu=xxu*xu;
  xxxxu=xxxu*xu;
  xxl=xl*xl;
  xxxl=xxl*xl;
  xxxxl=xxxl*xl;

  xuxl2=xxu-xxl;
  xuxl3=xxxu-xxxl;
  xuxl4=xxxxu-xxxxl;

  tt1=xuxl2*0.5*(mu-ml);
  tt2=(xu-xl)*(cu-cl);
  tt3=xuxl3*((mu-ml)/3.0);
  tt4=xuxl2*0.5*(cu-cl);
  tt5=(xuxl3*(mu*mu-ml*ml)/3.0);
  tt6=xuxl2*(mu*cu-ml*cl);
  tt7=(xu-xl)*(cu*cu-cl*cl);
  tt8=0.25*xuxl4*mu*mu;
  tt9=(xuxl3*2.0*mu*cu/3.0);
  tt10=0.5*xuxl2*cu*cu;
  tt11=0.25*xuxl4*ml*ml;
  tt12=(xuxl3*2.0*ml*cl/3.0);
  tt13=0.5*xuxl2*cl*cl;

  t1=a*(tt1+tt2);
  t2=b*(tt3+tt4);
  t3=0.5*c*(tt5+tt6+tt7);
  t4=0.5*d*(tt8+tt9+tt10-tt11-tt12-tt13);

//  cout<<"Shape::integrate(): terms in I3= "<<t1<<" "<<t2<<" "<<t3<<" "<<t4<<endl;

// contribution to element from third range triangle

  double I3(t1+t2+t3+t4);


// debug
//  cout<<"Shape::integrate(): yl(x),yu(x) parameters for I3:"<<endl;
//  cout<<"Shape::integrate():   ml= "<<ml<<endl;
//  cout<<"Shape::integrate():   mu= "<<mu<<endl;
//  cout<<"Shape::integrate():   xl= "<<xl<<endl;
//  cout<<"Shape::integrate():   xu= "<<xu<<endl;
//  cout<<"Shape::integrate():   cl= "<<cl<<endl;
//  cout<<"Shape::integrate():   cu= "<<cu<<endl;
//  cout<<"Shape::integrate(): I3= "<<I3<<endl;
// debug

// sum contribution from all 3 ranges to assemble the entire element

// debug
//  cout<<"Shape::integrate():    I1+I2+I3= "<<I1+I2+I3<<endl;
//  exit(1);
// debug

  return(I1+I2+I3);

}

// given three collinear points (p,q,r) this function checks if point q lies on line pr

bool Shape::onSegment(Point p,Point q,Point r) const {return (q.x <= max(p.x, r.x) && q.x >= min(p.x, r.x) &&q.y <= max(p.y, r.y) && q.y >= min(p.y, r.y));}

// orientation of ordered triplet (p,q,r) to return 0 (p,q,r colinear), 1 (clockwise), 2 (counter-clockwise)

int Shape::orientation(Point p,Point q,Point r) const{

  double val((q.y-p.y)*(r.x-q.x)-(q.x-p.x)*(r.y-q.y));

  return (val==0.0)?0:((val>0.0)?1:2);

}

// returns true if line segment 'p1q1' and 'p2q2' intersect

bool Shape::doIntersect(Point p1,Point q1,Point p2,Point q2) const {

// four orientations needed for general and special cases

  int o1 = orientation(p1, q1, p2);
  int o2 = orientation(p1, q1, q2);
  int o3 = orientation(p2, q2, p1);
  int o4 = orientation(p2, q2, q1);
 
// general case

  if (o1 != o2 && o3 != o4) return true;
 
// special case: p1, q1 and p2 are collinear and p2 lies on segment p1q1

  if (o1 == 0 && onSegment(p1, p2, q1)) return true;
 
// special case: p1, q1 and p2 are collinear and q2 lies on segment p1q1

  if (o2 == 0 && onSegment(p1, q2, q1)) return true;
 
// special case: p2, q2 and p1 are collinear and p1 lies on segment p2q2

  if (o3 == 0 && onSegment(p2, p1, q2)) return true;
 
// special case: p2, q2 and q1 are collinear and q1 lies on segment p2q2

  if (o4 == 0 && onSegment(p2, q1, q2)) return true;

// none of the above

  return false;

}


// returns true if point p lies inside the element
// ref: geeksforgeeks.org/how-to-check-if-a-given-point-lies-inside-a-polygon

bool Shape::isInside(Point p) const {

// form a polygon from the node positions

  Point polygon[nfaces()*(order()+1)];
  int k=0;

  for(int iloc=0;iloc<=order();iloc++,k++){polygon[k]={mr.at(0).at(iloc),mr.at(1).at(iloc)};}
  for(int iloc=order();iloc<nloc();iloc+=(order()+1),k++){polygon[k]={mr.at(0).at(iloc),mr.at(1).at(iloc)};}
  for(int iloc=nloc()-1;iloc>=nloc()-order()-1;iloc--,k++){polygon[k]={mr.at(0).at(iloc),mr.at(1).at(iloc)};}
  for(int iloc=nloc()-order()-1;iloc>=0;iloc-=(order()+1),k++){polygon[k]={mr.at(0).at(iloc),mr.at(1).at(iloc)};}

// create a point for line segment from p to infinity

  Point extreme={10000.0, p.y};
 
// sum intersections of the above line with sides of the element

  int count=0,i=0;
  do{

    int next=(i+1)%mnloc;
 
// check if the line segment from 'p' to 'extreme' intersects the line segment from 'polygon[i]' to 'polygon[next]'

    if (doIntersect(polygon[i], polygon[next], p, extreme)){

// if the point 'p' is collinear with line segment 'i-next' check if it lies on segment

      if (orientation(polygon[i], p, polygon[next]) == 0)
      return onSegment(polygon[i], p, polygon[next]);
 
      count++;

    }

    i=next;

  } while(i!=0);
 
// return true if count is odd, false otherwise

  return count&1; // Same as (count%2 == 1)

}

// tests if point v lies inside the element

bool Shape::contains(vector<double> const &v) const{

  Point pp{v.at(0),v.at(1)};

  return isInside(pp);

}

// factorial n

int factorial(int n){return((n>1)?n*factorial(n-1):1);}

// kronecker delta function

double delta(int i,int j){return (i==j)?1.0:0.0;}

// type safe function to return the sign of the argument

template <typename T> int sgn(T val) {return( (val>=T(0))?T(1):T(-1));} // -1 or 1
