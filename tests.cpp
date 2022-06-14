// Function definitions to code for various test problems

// S. R. Merton

#include <iostream>
#include <vector>
#include <iomanip>     // floating point precision
#include <cmath>       // for sqrt
#include <fstream>     // for file io
#include <algorithm>   // min_element, max_element
#include "globals.h"   // defines
#include "eos.h"       // eos lookups
#include "mesh.h"      // mesh class
#include "matrix.h"    // matrix class (needed to include bcs.h)
#include "shape.h"     // shape class
#include "riemann.h"   // riemann solver
#include "tests.h"     // test problems
#include "line.h"      // line class
#include "utilities.h" // vempty
#include "jacobian.h"  // jacobian
#include "bcs.h"       // VELOCITY bc definition

using namespace std;

// input overides for the Taylor-Green vortex

void init_TAYLOR(Mesh const &M,Shape const &S,Shape const &T,double const &dpi,VD &dinit,VVD &u0,VVD &u1,VD &e0,VD &e1,VVD const &xk,VVD const &xt,VD const &gamma,vector<int> const &mat,VD &detJ0,VVVD &detDJ0,VD &detJ,VVVD &detDJ,VD const &m){

  cout<<"init_TAYLOR(): Input overides for Taylor-Green..."<<endl;

// start x component of velocity field

  for(long i=0;i<u0.at(0).size();i++){
    u0.at(0).at(i)=sin(dpi*xk.at(0).at(i))*cos(dpi*xk.at(1).at(i));
    u1.at(0).at(i)=u0.at(0).at(i);
  }

// start y component of velocity field

  for(long i=0;i<u1.at(1).size();i++){
    u0.at(1).at(i)=-cos(dpi*xk.at(0).at(i))*sin(dpi*xk.at(1).at(i));
    u1.at(1).at(i)=u0.at(1).at(i);
  }

// start the energy field

  for(long i=0;i<M.NCells();i++){
    for(int iloc=0;iloc<T.nloc();iloc++){
      double xval(xt.at(0).at(M.GlobalNode_DFEM(i,iloc)));
      double yval(xt.at(1).at(M.GlobalNode_DFEM(i,iloc)));
      double rho(1.0);
      double p(0.25*rho*(cos(2.0*dpi*xval)+cos(2.0*dpi*yval))+1.0);

// invert the eos to obtain internal energy

      e0.at(M.GlobalNode_DFEM(i,iloc))=E(rho,p,gamma.at(mat.at(i)-1));
      e1.at(M.GlobalNode_DFEM(i,iloc))=e0.at(M.GlobalNode_DFEM(i,iloc));

    }
  }

  return;

}

// input overides for the Rayleigh-Taylor instability

void init_RAYLEIGH(Mesh const &M,Shape const &S,Shape const &T,double const &dpi,VD &dinit,VVD &u0,VVD &u1,VD &e0,VD &e1,VVD const &xk,VVD const &xt,VD const &gamma,vector<int> const &mat,VD &detJ0,VVVD &detDJ0,VD &detJ,VVVD &detDJ,VD const &m){

  cout<<"init_RAYLEIGH(): Input overides for the Rayleigh-Taylor instability test not coded yet."<<endl;

  exit(1);

  return;

}

// input overides for the Noh stagnation shock

void init_NOH(Mesh const &M,Shape const &S,Shape const &T,double const &dpi,VD &dinit,VVD &u0,VVD &u1,VD &e0,VD &e1,VVD const &xk,VVD const &xt,VD const &gamma,vector<int> const &mat,VD &detJ0,VVVD &detDJ0,VD &detJ,VVVD &detDJ,VD const &m){

  cout<<"init_NOH(): Input overides for Noh..."<<endl;

// set origin

  double xorig(0.5*(M.Min(0)+M.Max(0))),yorig(0.5*(M.Min(1)+M.Max(1)));

// correct origin for reflection

  if((M.bc_edge(0)==VELOCITY)){yorig=xk.at(1).at(0);} // ymin forced reflective
  if((M.bc_edge(3)==VELOCITY)){xorig=xk.at(0).at(0);} // xmin forced reflective

  double origin[2]={xorig,yorig}; // origin coordinates

// start velocity field

  for(long i=0;i<u0.at(0).size();i++){

    double rx(xk.at(0).at(i)-origin[0]);     // radial vector component from domain origin to node
    double ry(xk.at(1).at(i)-origin[1]);     // radial vector component from domain origin to node

// length of radial vector from domain origin to node

    double rnorm(sqrt(rx*rx+ry*ry));

// velocity is a radial vector from degree of freedom towards the domain origin

    u0.at(0).at(i)=-rx/max(rnorm,1.0e-12);
    u1.at(0).at(i)=u0.at(0).at(i);

    u0.at(1).at(i)=-ry/max(rnorm,1.0e-12);
    u1.at(1).at(i)=u0.at(1).at(i);

  }

  return;

}

// input overides for the Sedov explosion

void init_SEDOV(Mesh const &M,Shape const &S,Shape const &T,double const &dpi,VD &dinit,VVD &u0,VVD &u1,VD &e0,VD &e1,VVD const &xk,VVD const &xt,VD const &gamma,vector<int> const &mat,VD &detJ0,VVVD &detDJ0,VD &detJ,VVVD &detDJ,VD const &m){

  cout<<"init_SEDOV(): Input overides for Sedov..."<<endl;

// load energy field

  for(long i=0;i<M.NCells();i++){

// find cells to be sourced by checking for nodes at the origin

    bool zcentre(false);

    for(int iloc=0;iloc<S.nloc();iloc++){
      if( abs(xk.at(0).at(M.GlobalNode_CFEM(i,iloc)))<1.0e-7 && abs(xk.at(1).at(M.GlobalNode_CFEM(i,iloc)))<1.0e-7  ){
        zcentre=true;
      }
    }


    for(int iloc=0;iloc<T.nloc();iloc++){

    e0.at(M.GlobalNode_DFEM(i,iloc))=0.0;
    e1.at(M.GlobalNode_DFEM(i,iloc))=0.0;

// delta function at domain origin

      if(zcentre){
//        e0.at(i)=0.3014676/0.025; // drive ~ 12.058704, from another code see ref. paper for numbers
        e0.at(M.GlobalNode_DFEM(i,iloc))=0.25/m[i]; // place 1/4 of the drive in each of the 4 cells at the origin per unit mass
        e1.at(M.GlobalNode_DFEM(i,iloc))=e0.at(M.GlobalNode_DFEM(i,iloc));
      }

    }

  }

  return;

}

// input overides for the Saltzmann piston

void init_SALTZMANN(Mesh const &M,Shape const &S,Shape const &T,double const &dpi,VD &dinit,VVD &u0,VVD &u1,VD &e0,VD &e1,VVD &xk,VVD &xt,VD const &gamma,vector<int> const &mat,VD &detJ0,VVVD &detDJ0,VD &detJ,VVVD &detDJ,VD const &m){

  cout<<"init_SALTZMANN(): Input overides for Saltzmann..."<<endl;

// set some tolerances to be used here

  double tol(1.0e-6);
  double dx1(0.01);
  double dy1(0.01);

// perturb the mesh

  for(long i=0;i<xk.at(0).size();i++){

// original node position

    double xorig(xk.at(0).at(i));
    double yorig(xk.at(1).at(i));

// perturbation

    double w1(100.0*xorig+tol);
    double w2(100.0*yorig+tol);

// set new node position

// mesh distortion handled by generator so we don't need to do this
// also detJ0 has already been used so will need to move a few things
// around if we were to distort the mesh here...

//    xk.at(0).at(i)=w1*dx1+(10.0-w2)*dy1*sin(dpi*w1*0.01);

  }

// set initial energy field

  for(int i=0;i<M.NCells();i++){
    for(int iloc=0;iloc<T.nloc();iloc++){
      e0.at(M.GlobalNode_DFEM(i,iloc))=0.0001;
      e1.at(M.GlobalNode_DFEM(i,iloc))=0.0001;
    }
  }

  return;

}

// functions to produce lineouts in 1D and 2D

void lineouts(Mesh const &M,Shape const &S,Shape const &T,VD const &dinit,VD const &e,VVD const &xinit,VVD const &x,VVD const &xt,VVD const &u, int const &test_problem,vector<int> const &mat,VD const &g){

// select 2D lineouts for 2D problems and 1D lineouts for 1D problems

  switch(test_problem){

    case(TAYLOR):

// Taylor Green vortex

      lineouts_2d(M,S,T,dinit,e,xinit,x,xt,u,test_problem,mat,g);

      return;;

      break;

    case(RAYLEIGH):

// Rayleigh-Taylor instability

      lineouts_2d(M,S,T,dinit,e,xinit,x,xt,u,test_problem,mat,g);

      return;

      break;

    case(NOH):

// Noh stagnation shock

      lineouts_2d(M,S,T,dinit,e,xinit,x,xt,u,test_problem,mat,g);

      break;

    case(SEDOV):

// Sedov expanding shock

      lineouts_2d(M,S,T,dinit,e,xinit,x,xt,u,test_problem,mat,g);

      break;

    case(SOD):

// Sod's shock tube

      lineouts_1d(M,S,T,dinit,e,xinit,x,xt,u,test_problem,mat,g);

      break;

    case(R2R):

// 123 problem

      lineouts_1d(M,S,T,dinit,e,xinit,x,xt,u,test_problem,mat,g);

      break;

    case(SALTZMANN):

// Saltzmann piston

      return;

      break;

  }

  return;

}

// this function codes for some 1D lineouts

void lineouts_1d(Mesh const &M,Shape const &S,Shape const &T,VD const &dinit,VD const &e,VVD const &xinit,VVD const &x,VVD const &xt,VVD const &u, int const &test_problem,vector<int> const &mat,VD const &g){

// file handle for output

  ofstream f1,f2;

// establish the mesh limits

  double xmin(*min_element(x.at(0).begin(),x.at(0).end()));
  double xmax(*max_element(x.at(0).begin(),x.at(0).end()));
  double ymin(*min_element(x.at(1).begin(),x.at(1).end()));
  double ymax(*max_element(x.at(1).begin(),x.at(1).end()));
  double delta(0.001);

// decalre the lineout structure

  struct lineout_type {
    double x1,y1; // start point of each line
    double x2,y2; // end point of each line
    string filename; //filename to output
    string filehead; // file header
    int nsamples; // number of sample points per cell
    int lcentre; // cell left of mesh centre
    int rcentre; // cell right of mesh centre
    double xoffset; // allow an offset so we can coincide with the exact solution
    bool zskip; // skip sample points near the centre of the mesh
  } lineout;

  vector<lineout_type> Lineout;

// set up different lineouts for different problems

  switch(test_problem){

    case(SOD):

// Sod's shock tube

      lineout.x1=xmin;
      lineout.x2=xmax;
      lineout.y1=0.5*(ymin+ymax);
      lineout.y2=lineout.y1;
      lineout.filename="lineout_1.dat";
      lineout.filehead="# Sod lineout from (0.0,0.5) to (1.0,0.5) : Columns are x d e p ux uy";
      lineout.nsamples=5;
      lineout.rcentre=-10000;
      lineout.lcentre=-10000;
      lineout.xoffset=0.0;
      lineout.zskip=false;
      Lineout.push_back(lineout);

      break;

    case(R2R):

// 123 problem

      lineout.x1=xmin;
      lineout.x2=xmax;
      lineout.y1=0.5*(ymin+ymax);
      lineout.y2=lineout.y1;
      lineout.filename="lineout_1.dat";
      lineout.filehead="# 123 problem lineout from (xmin,0.5) to (xmax,0.5) : Columns are x d e p ux uy";
      lineout.nsamples=5;
      lineout.rcentre=M.NCells()/2;
      lineout.lcentre=lineout.rcentre-1;
      lineout.xoffset=0.8;
      lineout.zskip=true;
      Lineout.push_back(lineout);

      break;

  }

// loop over the lineouts and produce the output

  int iline(0);

  cout<<"lineouts_1d(): Lineout "<<iline+1<<" of "<<Lineout.size()<<" writing to file "<<Lineout.at(iline).filename<<" ..."<<endl;

// open the output file for the lineout and write the header part

  f1.open(Lineout.at(iline).filename);
  f1<<Lineout.at(iline).filehead<<endl;

// open a file to store the cell boundaries as these might be useful to inspect for 1D problems

  f2.open("cell_boundaries.dat");
  f2<<"# This file contains cell boundaries."<<endl;
  f2<<"# To unlink the cells and remove the unwanted diagonal lines after plotting:"<<endl;
  f2<<"# Use Data->Data set operations->Operation type:Split"<<endl;
  f2<<"# Set Length to 2"<<endl;

// split each cell into nsamples divisions for sampling

  for(int i=0;i<M.NCells(0);i++){

    double dxn(2.0/Lineout.at(iline).nsamples);

// set local coordinates of the sample points

    vector<double> xsample(Lineout.at(iline).nsamples+1),ysample(Lineout.at(iline).nsamples+1);

    xsample.at(0)=-1.0;
    ysample.at(0)=0.0;
    for(int isample=0;isample<Lineout.at(iline).nsamples;isample++){
      xsample.at(isample+1)=xsample.at(0)+(isample+1)*dxn;
      ysample.at(isample+1)=ysample.at(0);
    }

// jacobian at each sample point

    for(int isample=0;isample<=Lineout.at(iline).nsamples;isample++){
      double dxdu(0.0),dxdv(0.0),dydu(0.0),dydv(0.0),dxdu0(0.0),dxdv0(0.0),dydu0(0.0),dydv0(0.0);
      for(int iloc=0;iloc<S.nloc();iloc++){
        dxdu0+=S.dvalue(0,iloc,xsample.at(isample),ysample.at(isample))*xinit.at(0).at(M.GlobalNode_CFEM(i,iloc));
        dxdv0+=S.dvalue(1,iloc,xsample.at(isample),ysample.at(isample))*xinit.at(0).at(M.GlobalNode_CFEM(i,iloc));
        dydu0+=S.dvalue(0,iloc,xsample.at(isample),ysample.at(isample))*xinit.at(1).at(M.GlobalNode_CFEM(i,iloc));
        dydv0+=S.dvalue(1,iloc,xsample.at(isample),ysample.at(isample))*xinit.at(1).at(M.GlobalNode_CFEM(i,iloc));
        dxdu+=S.dvalue(0,iloc,xsample.at(isample),ysample.at(isample))*x.at(0).at(M.GlobalNode_CFEM(i,iloc));
        dxdv+=S.dvalue(1,iloc,xsample.at(isample),ysample.at(isample))*x.at(0).at(M.GlobalNode_CFEM(i,iloc));
        dydu+=S.dvalue(0,iloc,xsample.at(isample),ysample.at(isample))*x.at(1).at(M.GlobalNode_CFEM(i,iloc));
        dydv+=S.dvalue(1,iloc,xsample.at(isample),ysample.at(isample))*x.at(1).at(M.GlobalNode_CFEM(i,iloc));
      }

      double detJ0(dxdu0*dydv0-dxdv0*dydu0),detJ(dxdu*dydv-dxdv*dydu);

// reject sample points close to the cell edges as this helps to avoid sharp density spikes

      if(Lineout.at(iline).zskip){
        if((((isample==0)&&(i!=Lineout.at(iline).rcentre)) || ((isample==Lineout.at(iline).nsamples) &&(i!=Lineout.at(iline).lcentre)) )){continue;} // 123
//        bool zskip((i!=M.NCells()/2)||(i+1!=M.NCells()/2));
//        if(!zskip&&((abs(xsample.at(isample)+1.0)<1.0e-2)||(abs(xsample.at(isample)-1.0)<1.0e-2))){continue;}
      }else{
        if(((abs(xsample.at(isample)+1.0)<1.0e-2)||(abs(xsample.at(isample)-1.0)<1.0e-2))){continue;} // sod
      }

// initialise values for interpolation

      double interpolated_value[6]={Lineout.at(iline).xoffset,0.0,0.0,0.0,0.0,0.0}; // ordering is x d e p ux uy

// coordinate of sample point along the lineout

      for(int iloc=0;iloc<S.nloc();iloc++){interpolated_value[0]+=S.value(iloc,xsample.at(isample),ysample.at(isample))*x.at(0).at(M.GlobalNode_CFEM(i,iloc));}

// sample density field at the interpolation point

      interpolated_value[1]=dinit.at(i)*detJ0/detJ;

// sample energy field at the interpolation point

      for(int iloc=0;iloc<T.nloc();iloc++){interpolated_value[2]+=T.value(iloc,xsample.at(isample),ysample.at(isample))*e.at(M.GlobalNode_DFEM(i,iloc));}

// sample pressure field at the interpolation point

      interpolated_value[3]=P(interpolated_value[1],interpolated_value[2],g.at(mat.at(i)-1));

// sample velocity field x-component at interpolation point

      for(int iloc=0;iloc<S.nloc();iloc++){interpolated_value[4]+=S.value(iloc,xsample.at(isample),ysample.at(isample))*u.at(0).at(M.GlobalNode_CFEM(i,iloc));}

// sample velocity field y-component at interpolation point

      for(int iloc=0;iloc<S.nloc();iloc++){interpolated_value[5]+=S.value(iloc,xsample.at(isample),ysample.at(isample))*u.at(1).at(M.GlobalNode_CFEM(i,iloc));}

// output interpolated data along the lineout

      for(int j=0;j<6;j++){f1<<fixed<<setprecision(10)<<interpolated_value[j]<<" ";}
      f1<<endl;

    }

  }

// output cell boundaries

  for(int i=0;i<M.NCells(0);i++){
    f2<<setprecision(10)<<Lineout.at(iline).xoffset+x.at(0).at(M.GlobalNode_CFEM(i,0))<<" -1000000.0"<<endl;
    f2<<setprecision(10)<<Lineout.at(iline).xoffset+x.at(0).at(M.GlobalNode_CFEM(i,0))<<" 1000000.0"<<endl;
  }
  f2<<setprecision(10)<<Lineout.at(iline).xoffset+x.at(0).at(M.GlobalNode_CFEM(M.NCells()-1,S.order()))<<" -1000000.0"<<endl;
  f2<<setprecision(10)<<Lineout.at(iline).xoffset+x.at(0).at(M.GlobalNode_CFEM(M.NCells()-1,S.order()))<<" 1000000.0"<<endl;

// close the output files

  f1.close();
  f2.close();

  return;

}

// this function codes for some 2D lineouts

void lineouts_2d(Mesh const &M,Shape const &S,Shape const &T,VD const &dinit,VD const &e,VVD const &xinit,VVD const &x,VVD const &xt,VVD const &u, int const &test_problem,vector<int> const &mat,VD const &g){

// file handle for output

  ofstream f1;

// establish the mesh limits

  double xmin(*min_element(x.at(0).begin(),x.at(0).end()));
  double xmax(*max_element(x.at(0).begin(),x.at(0).end()));
  double ymin(*min_element(x.at(1).begin(),x.at(1).end()));
  double ymax(*max_element(x.at(1).begin(),x.at(1).end()));
  double delta(0.001);

// decalre the lineout structure

  struct lineout_type {
    double x1,y1; // start point of each line
    double x2,y2; // end point of each line
    string filename; //filename to output
    string filehead; // file header
    int nsamples; // number of sample points on each line
  } lineout;

  vector<lineout_type> Lineout;

// set up different lineouts for different problems

  switch(test_problem){

    case(TAYLOR):

// Taylor Green vortex

      return;;

      break;

    case(RAYLEIGH):

// Rayleigh-Taylor instability

      return;

      break;

    case(NOH):

// Noh stagnation shock

      for(int iline=0;iline<4;iline++){

        switch(iline){

          case(0):

            lineout.x1=0.5*(xmin+xmax);
            lineout.x2=0.5;
            lineout.y1=0.5*(ymin+ymax);
            lineout.y2=0.5;
            lineout.filename="lineout_1.dat";
            lineout.filehead="# Noh lineout from (0.0,0.0) to (0.5,0.05) : Columns are x d e p ux uy";
            lineout.nsamples=100;

            break;

          case(1):

            lineout.x1=0.5*(xmin+xmax);
            lineout.x2=-0.5;
            lineout.y1=0.5*(ymin+ymax);
            lineout.y2=0.5;
            lineout.filename="lineout_2.dat";
            lineout.filehead="# Noh lineout from (0.0,0.0) to (-0.5,0.5) : Columns are x d e p ux uy";
            lineout.nsamples=100;

            break;

          case(2):

            lineout.x1=0.5*(xmin+xmax);
            lineout.x2=-0.5;
            lineout.y1=0.5*(ymin+ymax);
            lineout.y2=-0.5;
            lineout.filename="lineout_3.dat";
            lineout.filehead="# Noh lineout from (0.0,0.0) to (-0.5,-0.5) : Columns are x d e p ux uy";
            lineout.nsamples=100;

            break;

          case(3):

            lineout.x1=0.5*(xmin+xmax);
            lineout.x2=0.5;
            lineout.y1=0.5*(ymin+ymax);
            lineout.y2=-0.5;
            lineout.filename="lineout_4.dat";
            lineout.filehead="# Noh lineout from (0.0,0.0) to (0.5,-0.5) : Columns are x d e p ux uy";
            lineout.nsamples=100;

            break;

        }

        Lineout.push_back(lineout);

      }

      break;

    case(SEDOV):

// Sedov expanding shock

      for(int iline=0;iline<4;iline++){

        switch(iline){

          case(0):

            lineout.x1=0.0+delta;
            lineout.x2=xmax;
            lineout.y1=0.0+delta;
            lineout.y2=ymax-delta;
            lineout.filename="lineout_1.dat";
            lineout.filehead="# Sedov lineout from (0.0,0.0) to (1.2,1.2) : Columns are x d e p ux uy";
            lineout.nsamples=100;

            break;

          case(1):

            lineout.x1=0.0-delta;
            lineout.x2=xmin;
            lineout.y1=0.0+delta;
            lineout.y2=ymax-delta;
            lineout.filename="lineout_2.dat";
            lineout.filehead="# Sedov lineout from (0.0,0.0) to (-1.2,1.2) : Columns are x d e p ux uy";
            lineout.nsamples=100;

            break;

          case(2):

            lineout.x1=0.0-delta;
            lineout.x2=xmin;
            lineout.y1=0.0-delta;
            lineout.y2=ymin+delta;
            lineout.filename="lineout_3.dat";
            lineout.filehead="# Sedov lineout from (0.0,0.0) to (-1.2,-1.2) : Columns are x d e p ux uy";
            lineout.nsamples=100;

            break;

          case(3):

            lineout.x1=0.0+delta;
            lineout.x2=xmax;
            lineout.y1=0.0-delta;
            lineout.y2=ymin+delta;
            lineout.filename="lineout_4.dat";
            lineout.filehead="# Sedov lineout from (0.0,0.0) to (1.2,-1.2) : Columns are x d e p ux uy";
            lineout.nsamples=100;

            break;

        }

        Lineout.push_back(lineout);

      }

      break;

    case(SOD):

// Sod's shock tube

      lineout.x1=-0.001;
      lineout.x2=1.0001;
      lineout.y1=0.5*(ymin+ymax);
      lineout.y2=lineout.y1;
      lineout.filename="lineout_1.dat";
      lineout.filehead="# Sod lineout from (0.0,0.5) to (1.0,0.5) : Columns are x d e p ux uy";
      lineout.nsamples=100;

      Lineout.push_back(lineout);

      break;

    case(R2R):

// 123 problem

      lineout.x1=xmin;
      lineout.x2=xmax;
      lineout.y1=0.5*(ymin+ymax);
      lineout.y2=lineout.y1;
      lineout.filename="lineout_1.dat";
      lineout.filehead="# 123 problem lineout from (xmin,0.5) to (xmax,0.5) : Columns are x d e p ux uy";
      lineout.nsamples=101; // odd avoids a sample at the stationary centre node which would give zero determinant

      Lineout.push_back(lineout);

      break;

    case(SALTZMANN):

// Saltzmann piston

      return;

      break;

  }

// loop over the lineouts and produce the output

  for(int iline=0;iline<Lineout.size();iline++){

    cout<<"lineouts_2d(): Lineout "<<iline+1<<" of "<<Lineout.size()<<" writing to file "<<Lineout.at(iline).filename<<" ..."<<endl;

// open the output file for the lineout and write the header part

    f1.open(Lineout.at(iline).filename);
    f1<<Lineout.at(iline).filehead<<endl;

// local node positions

    double xpos[S.nloc()],ypos[S.nloc()],dx(2.0/S.order()),dy(2.0/S.order());

    for(int isuby=0,k=0;isuby<S.order()+1;isuby++){
      for(int isubx=0;isubx<S.order()+1;isubx++,k++){
        xpos[k]=-1.0+isubx*dx;
        ypos[k]=-1.0+isuby*dy;
      }
    }

// set up line AB to sample along

    vector<double> A(2),B(2); // A and B are the two end points of the line

    A.at(0)=Lineout.at(iline).x1;A.at(1)=Lineout.at(iline).y1;
    B.at(0)=Lineout.at(iline).x2;B.at(1)=Lineout.at(iline).y2;

    Line AB(A,B);
    AB.divide(Lineout.at(iline).nsamples);

// place segment end points in a vector

    vector<vector<double> > end_point(AB.nsegments()+1);

    end_point.at(0).push_back(A.at(0));end_point.at(0).push_back(A.at(1));
    for(int iseg=0;iseg<AB.nsegments();iseg++){
      end_point.at(iseg+1).push_back(AB.coord(0,iseg));
      end_point.at(iseg+1).push_back(AB.coord(1,iseg));
    }

// search mesh for each point along the line

    for(int ipt=0;ipt<end_point.size();ipt++){

// set sample point coordinates and distance along the lineout from lineout origin

      double xpt(end_point.at(ipt).at(0)),ypt(end_point.at(ipt).at(1));
      double xdist(xpt-end_point.at(0).at(0)),ydist(ypt-end_point.at(0).at(1));

// sweep mesh and search for a cell that contains this point

      for(int i=0;i<M.NCells();i++){

// cell vertices

        vector<vector<double> > rk(2);
        for(int j=0;j<S.nloc();j++){
          rk.at(0).push_back(x.at(0).at(M.GlobalNode_CFEM(i,j)));
          rk.at(1).push_back(x.at(1).at(M.GlobalNode_CFEM(i,j)));
        }

// declare a shape function in global coordinates

        Shape Gk(S.order(),rk);

// if point ipt is not inside this cell cycle the loop

        if(!Gk.contains(end_point.at(ipt))){continue;}

// use global shape to map global coordinates to local coordinates

        double xloc(0.0),yloc(0.0);

        for(int iloc=0;iloc<S.nloc();iloc++){
          xloc+=Gk.value(iloc,end_point.at(ipt))*xpos[iloc];
          yloc+=Gk.value(iloc,end_point.at(ipt))*ypos[iloc];
        }

// jacobian and determinant at local coordinates

        double dxdu(0.0),dxdv(0.0),dydu(0.0),dydv(0.0),dxdu0(0.0),dxdv0(0.0),dydu0(0.0),dydv0(0.0);
        for(int iloc=0;iloc<S.nloc();iloc++){
          dxdu0+=S.dvalue(0,iloc,xloc,yloc)*xinit.at(0).at(M.GlobalNode_CFEM(i,iloc));
          dxdv0+=S.dvalue(1,iloc,xloc,yloc)*xinit.at(0).at(M.GlobalNode_CFEM(i,iloc));
          dydu0+=S.dvalue(0,iloc,xloc,yloc)*xinit.at(1).at(M.GlobalNode_CFEM(i,iloc));
          dydv0+=S.dvalue(1,iloc,xloc,yloc)*xinit.at(1).at(M.GlobalNode_CFEM(i,iloc));
          dxdu+=S.dvalue(0,iloc,xloc,yloc)*x.at(0).at(M.GlobalNode_CFEM(i,iloc));
          dxdv+=S.dvalue(1,iloc,xloc,yloc)*x.at(0).at(M.GlobalNode_CFEM(i,iloc));
          dydu+=S.dvalue(0,iloc,xloc,yloc)*x.at(1).at(M.GlobalNode_CFEM(i,iloc));
          dydv+=S.dvalue(1,iloc,xloc,yloc)*x.at(1).at(M.GlobalNode_CFEM(i,iloc));
        }

        double detJ0(dxdu0*dydv0-dxdv0*dydu0),detJ(dxdu*dydv-dxdv*dydu);

// initialise values for interpolation

        double interpolated_value[6]={0.0,0.0,0.0,0.0,0.0,0.0}; // ordering is x d e p ux uy

// coordinate of sample point along the lineout

        interpolated_value[0]=sqrt(xdist*xdist+ydist*ydist);

// sample density field at the interpolation point

        interpolated_value[1]=dinit.at(i)*detJ0/detJ;

// sample energy field at the interpolation point

        for(int iloc=0;iloc<T.nloc();iloc++){interpolated_value[2]+=T.value(iloc,xloc,yloc)*e.at(M.GlobalNode_DFEM(i,iloc));}

// sample pressure field at the interpolation point

        interpolated_value[3]=P(interpolated_value[1],interpolated_value[2],g.at(mat.at(i)-1));

// sample velocity field x-component at global coordinate ri(x,y)

        for(int iloc=0;iloc<S.nloc();iloc++){interpolated_value[4]+=S.value(iloc,xloc,yloc)*u.at(0).at(M.GlobalNode_CFEM(i,iloc));}

// sample velocity field y-component at global coordinate ri(x,y)

        for(int iloc=0;iloc<S.nloc();iloc++){interpolated_value[5]+=S.value(iloc,xloc,yloc)*u.at(1).at(M.GlobalNode_CFEM(i,iloc));}

// output interpolated data along the lineout

        for(int j=0;j<6;j++){f1<<fixed<<setprecision(10)<<interpolated_value[j]<<" ";}
        f1<<endl;

      }

    }

// close the output file

    f1.close();

  }

  if(Lineout.size()==0){cout<<"lineouts(): No lineouts defined for this problem."<<endl;}

  return;

}

// function to code for the exact solution where applicable

void exact(VVD const &s,VVD const &x,int const &test_problem,double const &time){

// file handle for output

  ofstream f1;

// establish the mesh limits

  double xmin(*min_element(x.at(0).begin(),x.at(0).end()));
  double xmax(*max_element(x.at(0).begin(),x.at(0).end()));
  double ymin(*min_element(x.at(1).begin(),x.at(1).end()));
  double ymax(*max_element(x.at(1).begin(),x.at(1).end()));

// decalre the lineout structure

  struct lineout_type {
    double x1,y1; // start point of each line
    double x2,y2; // end point of each line
    string filename; //filename to output
    string filehead; // file header
    int nsamples; // number of sample points on each line
  } lineout;

  vector<lineout_type> Lineout;

// set up an exact solution for each problems that has one

  switch(test_problem){

    case(TAYLOR):

// Taylor Green vortex

      return;;

      break;

    case(RAYLEIGH):

// Rayleigh-Taylor instability

      return;

      break;

    case(NOH):

// Noh stagnation shock

      return;;

      break;

    case(SEDOV):

// Sedov expanding shock

      return;;

      break;

    case(SOD):

// Sod's shock tube

      lineout.x1=xmin;
      lineout.x2=xmax;
      lineout.y1=0.5*(ymin+ymax);
      lineout.y2=lineout.y1;
      lineout.filename="exact.dat";
      lineout.filehead="# Sod exact solution from (0.0,0.5) to (1.0,0.5) : Columns are x d p e u";
      lineout.nsamples=NSAMPLES;

      {

        double l[4]={s.at(0)[0],s.at(0)[1],s.at(0)[3],s.at(0)[4]}; // left flux state
        double r[4]={s.at(1)[0],s.at(1)[1],s.at(1)[3],s.at(1)[4]}; // right flux state
        Riemann R(Riemann::exact,l,r); // start Riemann solver from initial flux states
        vector<double> rx;vempty(rx); // sample point coordinates for the Riemann solver
//        for(int i=0;i<NSAMPLES;i++){rx.push_back(0.0+(i*(1.0-0.0)/double(NSAMPLES)));} // sample points for the Riemann solver
        for(int i=0;i<NSAMPLES;i++){rx.push_back(i*(xmax-xmin)/NSAMPLES);} // sample points for the Riemann solver
        R.profile(&rx,time); // Riemann solution at NSAMPLE sample points

// append to the data structure

        Lineout.push_back(lineout);

        cout<<"exact(): exact solution writing to file "<<Lineout.at(0).filename<<" ..."<<endl;

// open the output file for the exact solution and write the header part

        f1.open(Lineout.at(0).filename);
        f1<<Lineout.at(0).filehead<<endl;

// output the exact solution

        f1<<fixed<<setprecision(10);
        for(int j=0;j<rx.size();j++){
          f1<<rx.at(j)<<" "<<R.density(j)<<" "<<R.pressure(j)<<" "<<R.energy(j)<<" "<<R.velocity(j)<<endl;
        }

// close the output file

        f1.close();

      }

      break;


    case(R2R):

// 123 problem

      lineout.x1=xmin;
      lineout.x2=xmax;
      lineout.y1=0.5*(ymin+ymax);
      lineout.y2=lineout.y1;
      lineout.filename="exact.dat";
      lineout.filehead="# 123 problem solution from (0.0,0.5) to (1.0,0.5) : Columns are x d p e u";
      lineout.nsamples=NSAMPLES;

      {

        double l[4]={s.at(0)[0],s.at(0)[1],s.at(0)[3],s.at(0)[4]}; // left flux state
        double r[4]={s.at(1)[0],s.at(1)[1],s.at(1)[3],s.at(1)[4]}; // right flux state
        Riemann R(Riemann::exact,l,r); // start Riemann solver from initial flux states
        vector<double> rx;vempty(rx); // sample point coordinates for the Riemann solver
        for(int i=0;i<NSAMPLES;i++){rx.push_back(i*(xmax-xmin)/NSAMPLES);} // sample points for the Riemann solver
        R.profile(&rx,time); // Riemann solution at NSAMPLE sample points

// append to the data structure

        Lineout.push_back(lineout);

        cout<<"exact(): exact solution writing to file "<<Lineout.at(0).filename<<" ..."<<endl;

// open the output file for the exact solution and write the header part

        f1.open(Lineout.at(0).filename);
        f1<<Lineout.at(0).filehead<<endl;

// output the exact solution

        f1<<fixed<<setprecision(10);
        for(int j=0;j<rx.size();j++){
          f1<<rx.at(j)<<" "<<R.density(j)<<" "<<R.pressure(j)<<" "<<R.energy(j)<<" "<<R.velocity(j)<<endl;
        }

// close the output file

        f1.close();

      }

      break;

    case(SALTZMANN):

// Saltzmann piston

      return;

      break;

  }

  return;

}

// function to code for a manufactured solution on the internal energy field

void manufactured_soln_ie(Mesh const &M,Shape const &S,Shape const &T,VVD const &x1,int const &test_problem,VD &detJ,VVVD &detDJ,VD &b){

// define pi

  double dpi(4.0*atan(1.0));

  switch(test_problem){

    case(TAYLOR):

// Taylor Green vortex

      for(int i=0;i<M.NCells();i++){

        jacobian(i,x1,M,S,detJ,detDJ);

// coordinates of the integration points

        vector<double> egi(S.ngi(),0.0);
        for(int gi=0;gi<S.ngi();gi++){
          double xgi(0.0),ygi(0.0);
          for(int iloc=0;iloc<S.nloc();iloc++){
            xgi+=x1.at(0).at(M.GlobalNode_CFEM(i,iloc))*S.value(iloc,gi);
            ygi+=x1.at(1).at(M.GlobalNode_CFEM(i,iloc))*S.value(iloc,gi);
          }

// compute manufactured solution at the integration points

         egi.at(gi)=(3.0*dpi/8.0)*((cos(3.0*dpi*xgi)*cos(dpi*ygi))-(cos(dpi*xgi)*cos(3.0*dpi*ygi)));

        }

// place manufactured solution on rhs of energy equation

        for(int iloc=0;iloc<T.nloc();iloc++){
          double bsum(0.0);
          for(int gi=0;gi<S.ngi();gi++){
            bsum-=T.value(iloc,gi)*egi.at(gi)*detJ.at(gi)*S.wgt(gi);
          }
          b.at(M.GlobalNode_DFEM(i,iloc))=bsum;
        }
      }

      break;

    case(RAYLEIGH):

// Rayleigh-Taylor instability

      return;

      break;

    case(NOH):

// Noh stagnation shock

      return;;

      break;

    case(SEDOV):

// Sedov expanding shock

      return;

      break;

    case(SOD):

// Sod's shock tube

      return;

      break;


    case(R2R):

// 123 problem

      return;

      break;

    case(SALTZMANN):

// Saltzmann piston

      return;

      break;

  }

  return;

}
