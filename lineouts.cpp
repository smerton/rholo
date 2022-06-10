// Functions to produce lineouts in 1D and 2D

// Author S. R. Merton

#define VD vector<double> // vector of doubles
#define VVD vector<VD>    // vector of VD
#define TAYLOR 1              // Taylor-Green vortex problem
#define RAYLEIGH 2            // Rayleigh-Taylor instability problem
#define NOH 3                 // Noh stagnation shock problem
#define SEDOV 4               // Sedov expanding shock problem
#define TRIPLE 5              // Triple point problem
#define SOD 6                 // Sod's shock tube problem
#define R2R 7                 // 123 problem
#define SALTZMANN 8           // Saltzmann piston problem

#include <iostream>
#include <vector>
#include <iomanip>    // floating point precision
#include <cmath>      // for sqrt
#include <fstream>    // for file io
#include <algorithm>  // min_element, max_element
#include "eos.h"      // eos lookups
#include "line.h"
#include "mesh.h"
#include "shape.h"
#include "lineouts.h" // function signatures

using namespace std;

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
      lineout.xoffset=0.9;
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

// reject sample points closest to the cell edges as this helps to avoid sharp density spikes

      if(Lineout.at(iline).zskip){
        if((((isample==0)&&(i!=Lineout.at(iline).rcentre)) || ((isample==Lineout.at(iline).nsamples) &&(i!=Lineout.at(iline).lcentre)) )){continue;}
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
  f2<<setprecision(10)<<Lineout.at(iline).xoffset+x.at(0).at(M.GlobalNode_CFEM(M.NCells()-1,S.order()))<<" -100000.0"<<endl;
  f2<<setprecision(10)<<Lineout.at(iline).xoffset+x.at(0).at(M.GlobalNode_CFEM(M.NCells()-1,S.order()))<<" 100000.0"<<endl;

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
