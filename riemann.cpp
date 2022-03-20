// Function definitions for members of the Riemann class

// Author S. R. Merton

#include <iostream>
#include"riemann.h"
#include <cmath>
#include <iomanip>
#include <algorithm>

#define NITER 20 // number iterations of the pressure function
#define TOL 1.0e-6 // convergence criteria for the pressure iteration
#define QMAX 2.0 // maximum permitted pressure ratio

using namespace std;

void empty(vector<double> &v);

Riemann::Riemann(solver s,double*l,double*r){

// Constructor to start a new Riemann problem and load all of its attributes
// within the scope of the calling function
//
// Arguments:
// s of type solver contains the name of the Riemann approximation to use, choices are:
// pvrs for the primitive variable Riemann solver
// exact for the exact Riemann solver
// l of type double* containg the left state of the contact
// r of type double* containg the right state of the contact

//  cout<<"Riemann::Riemann(): Starting up a Riemann solver of type "<<solvername[s]<<"..."<<endl;

// set up the left state

  mDl=l[0];
  mul=l[1];
  mPl=l[2];
  mgl=l[3];

// set up the right state

  mDr=r[0];
  mur=r[1];
  mPr=r[2];
  mgr=l[3];

// gamma values

  if(mgl!=mgr){
    cout<<"Riemann::Riemann(): left and right states have a different gamma, stopping."<<endl;
    exit(1);
  }

  g0=mgl;
  g1=(mgl-1.0)/(2.0*mgl);
  g2=(mgl+1.0)/(2.0*mgl);
  g3=2.0*mgl/(mgl-1.0);
  g4=2.0/(mgl-1.0);
  g5=2.0/(mgl+1.0);
  g6=(mgl-1.0)/(mgl+1.0);
  g7=0.5*(mgl-1.0);
  g8=1.0/mgl;
  g9=mgl-1.0;

//  cout<<"Riemann::Riemann(): left state (d,u,p) = "<<mDl<<","<<mul<<","<<mPl<<endl;
//  cout<<"Riemann::Riemann(): right state (d,u,p) = "<<mDr<<","<<mur<<","<<mPr<<endl;

// load the selevcted solver

  switch(s){
    case pvrs:{
      pvrs_solver();
      break;}
    case exact:{;
      exact_solver();
      break;}
    default:
      cout<<"Riemann::Riemann(): ERROR: solver "<<s<<"is not implemented."<<endl;
  }

}

// primitive variable Riemann solver
// this one is based on a simple acoustic wave structure
// reference: pp 280, Ch. 9 Toro

void Riemann::pvrs_solver(){

// sound speed

  cl=sqrt(g0*mPl/mDl);
  cr=sqrt(g0*mPr/mDr);

// acoustic impedance

  double zl(mPl*cl),zr(mPr*cr);

// approximate pressure along the discontinuity between states

  pstar=(zl*mPr+zr*mPl)/(zl+zr)+(zl*zr)*(mul-mur)/(zl+zr);

// velocity at the discontinuity

  ustar=(mPl-mPr)/(zl+zr)+(mur*zr+mul*zl)/(zl+zr);

  return;

}

// exact Riemann solver for ideal gases
// loads the two-shock approximation if a shock wave is presenting
// loads the two-rarefaction approximation where rarefactions are presenting

void Riemann::exact_solver(){

// sounds speeds left and right

  cl=sqrt(g0*mPl/mDl);
  cr=sqrt(g0*mPr/mDr);

// compute critical velocity

  double ducrit=g4*(cl+cr)-(mur-mul);

// stop if this has generated a vacuum

  if(ducrit<=0.0){
    cout<<"riemann(): Error - vacuum generated."<<endl;
    exit(1);
  }

// compute a guess value

  pstar=starte();

  double p0(pstar),fl(0.0),fld(0.0),fr(0.0),frd(0.0),du(mur-mul),cha(1.0e10);

// iterate on the pressure function

  while(cha>TOL){
    prefun(&fl,&fld,pstar,mDl,mPl,cl);
    prefun(&fr,&frd,pstar,mDr,mPr,cr);
    pstar-=(fl+fr+du)/(fld+frd);
    cha=2.0*abs((pstar-p0)/(pstar+p0));
    if(pstar<0.0)pstar=TOL;
    p0=pstar;
  }

// compute u

  ustar=0.5*(mul+mur+fr-fl);

  return;

}

// function to code for the initial guess in the pressure iteration

double Riemann::starte(){

// guess values from pvrs riemann solver

  double pv(0.5*(mPl+mPr)-0.125*(mur-mul)*(mDl+mDr)*(cl+cr));
  double pmin(min(mPl,mPr)),pmax(max(mPl,mPr));
  double qrat(pmax/pmin),p(0.0),pnu(0.0),pde(0.0);

  if(qrat<=QMAX&&(pmin<=pv&&pv<=pmax)){

// use pvrs solution as the guess

    p=max(TOL,pv);

  }else{

    if(pv<pmin){

// use trrs solution as the guess

      pnu=(cl+cr-g7*(mur-mul)); // numerator
      pde=(cl/pow(mPl,g1)+cr/pow(mPr,g1)); // denominator
      p=pow((pnu/pde),g3);

    }else{

// use tsrs solution as the guess

      double gel(sqrt((g5/mDl)/(g6*mPl+max(TOL,pv))));
      double ger(sqrt((g5/mDr)/(g6*mPr+max(TOL,pv))));
      p=(gel*mPl+ger*mPr-(mur-mul))/(gel+ger);
      p=max(TOL,p);

    }

  }

  return p;

}

// this codes for the pressure function

void Riemann::prefun(double*f,double*fd,double p,double dk,double pk,double ck){

  if(p<=pk){

// rarefection wave

    double prat(p/pk);
    (*f)=g4*ck*(pow(prat,g1)-1.0);
    (*fd)=(1.0/(dk*ck))*pow(prat,-g2);

  }else{

// shock wave

    double ak(g5/dk),bk(g6*pk),qrt(sqrt(ak/(bk+p)));
    (*f)=(p-pk)*qrt;
    (*fd)=(1.0-0.5*(p-pk)/(bk+p))*qrt;

  }

  return;

}

// function to implement the action to empty a vector

void empty(vector<double> &v){

  vector<double> e;

  v.assign(e.begin(),e.end());

  return;

}

// codes for the profile at each sample point (x,t)

void Riemann::profile(vector<double>*x,double time){

// clear member address space

  empty(mdensity);
  empty(mpressure);
  empty(menergy);
  empty(mvelocity);
  empty(mregion);

// load coordinates

  mSCoords.assign((*x).begin(),(*x).end());
  cout<<fixed<<setprecision(5);

// set domain boundaries

  double xr=*max_element(mSCoords.begin(), mSCoords.end());
  double xl=*min_element(mSCoords.begin(), mSCoords.end());

// obtain the profiles

  for(long i=0;i<mSCoords.size();i++){
    double s((mSCoords.at(i)-0.5*(xr-xl))/time),ds(0.0),us(0.0),ps(0.0),es(0.0),rs(0.0);
    sample(s,ps,us,ds,es,rs); // solution at point (s,t)
    mdensity.push_back(ds);
    mpressure.push_back(ps);
    menergy.push_back(es);
    mvelocity.push_back(us);
    mregion.push_back(rs);
  }

  return;

}

// this function codes for the density field at position i

double Riemann::density(long i){return mdensity[i];}

// codes for the pressure field at position i

double Riemann::pressure(long i){return mpressure[i];}

// codes energy field at position i

double Riemann::energy(long i){return menergy[i];}

// codes for velocity field at position i

double Riemann::velocity(long i){return mvelocity[i];}

// codes for region field at position i

double Riemann::region(long i){return mregion[i];}

// sample the solution at position s according to wave patterns
// p,u,d,e are pressure, velocity, density and energy solutions

void Riemann::sample(double s,double &p,double &u,double &d,double &e,double &r){

  if(s<=ustar){

// sample point is to the left of the contact

    if(pstar<=mPl){

// left fan

      double shl(mul-cl);

      if(s<=shl){

// left data state

        d=mDl;
        u=mul;
        p=mPl;
        r=1.0;

      }else{

        double cml(cl*pow((pstar/mPl),g1));
        double stl(ustar-cml);

        if(s>stl){

// middle left state

          d=mDl*pow((pstar/mPl),g8);
          u=ustar;
          p=pstar;
          r=3.0;

        }else{

// a left state (inside fan)

          u=g5*(cl+g7*mul+s);
          double c(g5*(cl+g7*(mul-s)));
          d=mDl*pow((c/cl),g4);
          p=mPl*pow((c/cl),g3);
          r=2.0;

        }

      }

    }else{

// left shock

      double pml(pstar/mPl);
      double sl(mul-cl*sqrt(g2*pml+g1));

      if(s<=sl){

// left data state

        d=mDl;
        u=mul;
        p=mPl;
        r=1.0;

      }else{

// middle left state (behind shock)

        d=mDl*(pml+g6)/(pml*g6+1.0);
        u=ustar;
        p=pstar;
        r=3.0;

      }

    }

  }else{

// right of contact

    if(pstar>mPr){

// right shock

      double pmr(pstar/mPr);
      double sr(mur+cr*sqrt(g2*pmr+g1));

      if(s>=sr){

// right data state

        d=mDr;
        u=mur;
        p=mPr;
        r=5.0;

      }else{

// middle right state (behind shock)

        d=mDr*(pmr+g6)/(pmr*g6+1.0);
        u=ustar;
        p=pstar;
        r=4.0;

      }

    }else{

// right fan

      double shr(mur+cr);

      if(s>=shr){

// right data state

        d=mDr;
        u=mur;
        p=mPr;
        r=5.0;

      }else{

        double cmr(cr*pow((pstar/mPr),g1));
        double str(ustar+cmr);

        if(s<=str){

//middle right state

          d=mDr*pow((pstar/mPr),g8);
          u=ustar;
          p=pstar;
          r=4.0;

        }else{

// fan right state (inside fan)

          u=g5*(-cr+g7*mur+s);
          double c(g5*(cr-g7*(mur-s)));
          d=mDr*pow((c/cr),g4);
          p=mPr*pow((c/cr),g3);
          r=3.0;

        }

      }

    }

  }

// set energy value at sample point

  e=p/d/(g0-1.0);

  return;

}

// Destructor function to release storage associated with a Shape class object

Riemann::~Riemann(){}
