// This function writes members of the input class to the supplied filename

// Author S. R. Merton

#include <iostream>
#include "input.h"
#include <fstream>
#include <iomanip>
#include <ctime>     // date and time
#include <stdlib.h>  // getenv
#include <unistd.h>  // gethostname
#include <cmath>     // atan

using namespace std;

string gethost();
string date(){time_t now = time(0);return ctime(&now);} // function to return the date and time

void Output(char* outputfile,Input*I){

  cout<<"Output(): Writing input pipeline to meshfile "<<outputfile<<endl;

// write the input pipeline

  ofstream file(outputfile);

  file<<"MFEM mesh v1.0"<<endl;
  file<<"#"<<endl;
  file<<"#"<<I->Title()<<endl;

  for(int i=0;i<I->NDescriptions();i++){
    file<<"#"<<I->Description(i)<<endl;
  }

// write header

  file<<"#"<<endl;
  file<<"# Input file for RhoLo based on MFEM mesh file format"<<endl;
  file<<"# This file was created by gen on "<<date();
  file<<"# and can be visualised using GLVis"<<endl;
  file<<"# See https://mfem.org/mesh-format-v1.x for description of this format"<<endl;
  file<<"#"<<endl;
  file<<"# Stamp revealed the following header information:"<<endl;
  file<<"#     User: "<<getenv("USER")<<endl;
  file<<"#     Host: "<<gethost()<<endl;
  file<<"#     Home: "<<getenv("HOME")<<endl;
  file<<"#     pwd:  "<<getenv("PWD")<<endl;
  file<<"#"<<endl;
  file<<"# RhoLo uses the following which may be inconsistent with MFEM:"<<endl;
  file<<"#   element attribute for material number (column 1 in elements block)"<<endl;
  file<<"#   boundary attribute for mesh side, =0 for bottom, 1 for right, =2 for top, =3 for left (column 1 in boundary block)"<<endl;
  file<<"#   boundary type not used - it is not known what MFEM uses this for (column 2 in boundary block)"<<endl;
  file<<"#   FiniteElementCollection not used - but MFEM may use this to set polyhedral order (in nodes block)"<<endl;
  file<<"#"<<endl;
  file<<"# The following format is interpreted by RhoLo"<<endl;
  file<<"# dimension"<<endl;
  file<<"# <number of spatial dimensions of the mesh>"<<endl;
  file<<"#"<<endl;
  file<<"# elements"<<endl;
  file<<"# <number of elements>"<<endl;
  file<<"# <material number> <element type from MFEM geometry type list, see below> <node number of vertices>"<<endl;
  file<<"#"<<endl;
  file<<"# boundary"<<endl;
  file<<"# <number of element sides on the mesh boundary>"<<endl;
  file<<"# <element side coincident with mesh edge> <boundary condition, 1=free surface, 2= forced-reflective> <node number of vertices on this side>"<<endl;
  file<<"#"<<endl;
  file<<"# vertices"<<endl;
  file<<"# <number of mesh vertices>"<<endl;
  file<<"# <serialisation, =2 for x,y pairs>"<<endl;
  file<<"# <coordinate of vertex in each dimension>"<<endl;
  file<<"#"<<endl;
  file<<"#"<<endl;
  file<<"# MFEM Geometry Types (see mesh/geom.hpp):"<<endl;
  file<<"#"<<endl;
  file<<"# POINT       = 0"<<endl;
  file<<"# SEGMENT     = 1"<<endl;
  file<<"# TRIANGLE    = 2"<<endl;
  file<<"# SQUARE      = 3"<<endl;
  file<<"# TETRAHEDRON = 4"<<endl;
  file<<"# CUBE        = 5"<<endl;
  file<<"# PRISM       = 6"<<endl;
  file<<"#"<<endl;
  file<<endl;

  file<<"dimension"<<endl;
  file<<I->NDims()-1<<endl;
  file<<endl;

  int nx(I->Nodes(0)),ny(I->Nodes(1));
  int ncellsx(nx-1),ncellsy(ny-1);

// set up cell widths

  double dx[ncellsx],dy[ncellsy];

  for(int i=0;i<ncellsx;i++){dx[i]=(I->XMax()-I->XMin())/ncellsx;}
  for(int j=0;j<ncellsy;j++){dy[j]=(I->YMax()-I->YMin())/ncellsy;}

// some simluations use a distrorted grid

  if(I->zNoh()){
    cout<<"Output():: distorting mesh for Noh problem..."<<endl;

    double dx1(3*(0.0-I->XMin())/(2*ncellsx)); // high resolution quadrant
    double dx2(3*(0.0-I->XMin())/(1*ncellsx)); // low resolution quadrant
    double dy1(3*(0.0-I->YMin())/(2*ncellsy)); // high resolution quadrant
    double dy2(3*(0.0-I->YMin())/(1*ncellsy)); // low resolution quadrant

    for(int i=0;i<ncellsx;i++){
      if(i<2*ncellsx/3){
// high res quadrant
        dx[i]=dx1;
      }else{
        dx[i]=dx2;
      }
    }

    for(int j=0;j<ncellsy;j++){
      if(j>=ncellsy/3){
// high res quadrant
        dy[j]=dy1;
      }else{
        dy[j]=dy2;
      }
    }

    cout<<"Output():: done."<<endl;

  }

  if(I->zSedov()){
    cout<<"Output():: distorting mesh for Sedov problem..."<<endl;

    double dx1(3*(0.0-I->XMin())/(2*ncellsx)); // high resolution quadrant
    double dx2(3*(0.0-I->XMin())/(1*ncellsx)); // low resolution quadrant
    double dy1(3*(0.0-I->YMin())/(2*ncellsy)); // high resolution quadrant
    double dy2(3*(0.0-I->YMin())/(1*ncellsy)); // low resolution quadrant

    for(int i=0;i<ncellsx;i++){
      if(i<2*ncellsx/3){
// high res quadrant
        dx[i]=dx1;
      }else{
        dx[i]=dx2;
      }
    }

    for(int j=0;j<ncellsy;j++){
      if(j>=ncellsy/3){
// high res quadrant
        dy[j]=dy1;
      }else{
        dy[j]=dy2;
      }
    }

    cout<<"Output():: done."<<endl;

  }

// set up vertices

  double vx[nx*ny],vy[nx*ny];

  for(long i=0;i<nx*ny;i+=nx){vx[i]=I->XMin();}
  for(long i=0;i<nx;i++){vy[i]=I->YMin();}

  for(long j=0;j<ny;j++){
    for(long i=0;i<ncellsx;i++){
      long k(j*nx+i+1),k0(k-1);
      vx[k]=vx[k0]+dx[i];
    }
  }

  for(long j=0;j<ncellsy;j++){
    for(long i=0;i<nx;i++){
      long k((j+1)*nx+i),k0(k-nx);
      vy[k]=vy[k0]+dy[j];
    }
  }

// distort the mesh for Saltzmann, this requires vx[] and vy[] to be set

  if(I->zSaltzmann()){
    cout<<"Output():: distorting mesh for Saltzmann problem..."<<endl;

// set some tolerances to be used here

    double tol(1.0e-6);
    double dx1(0.01);
    double dy1(0.01);
    double dpi(4.0*atan(1.0));

// perturb the mesh

    for(long i=0;i<nx*ny;i++){

// original node position

        double xorig(vx[i]);
        double yorig(vy[i]);

// perturbation

        double w1(100.0*xorig+tol);
        double w2(100.0*yorig+tol);

// set new node position

        vx[i]=w1*dx1+(10.0-w2)*dy1*sin(dpi*w1*0.01);

    }

    cout<<"Output():: done."<<endl;

  }

// write elements block

  file<<"elements"<<endl;
  file<<I->NCells()<<endl;

  for(int j=0,iel=0;j<ncellsy;j++){
    for(int i=0;i<ncellsx;i++,iel++){

// corner nodes

      int i1(j*nx+i),i2(j*nx+i+1),i3((j+1)*nx+i),i4((j+1)*nx+i+1);

// calculate centroid

      double xc[2];
      xc[0]=0.25*(vx[i1]+vx[i2]+vx[i3]+vx[i4]);
      xc[1]=0.25*(vy[i1]+vy[i2]+vy[i3]+vy[i4]);

// assign material to element

      int mat;
      for(int imat=0;imat<I->NMaterials();imat++){
        if((xc[0]>I->Range(0,imat)&&xc[0]<I->Range(1,imat))&&(xc[1]>I->Range(2,imat)&&xc[1]<I->Range(3,imat))){
          mat=I->Material(imat);
          break;
        }
      }

      file<<mat<<" 3 "<<i1<<" "<<i2<<" "<<i4<<" "<<i3<<" "<<endl;

    }
  }
  file<<endl;

  file<<"boundary"<<endl;
  file<<2*ncellsx+2*ncellsy;
  file<<endl;

  for(int i=0;i<ncellsx;i++){file<<"0 1 "<<i<<" "<<i+1<<endl;}
  for(int j=0;j<ncellsy;j++){file<<"1 1 "<<(j+1)*nx-1<<" "<<(j+2)*nx-1<<endl;}
  for(int i=0;i<ncellsx;i++){file<<"2 1 "<<nx*ny-i-1<<" "<<nx*ny-i-2<<endl;}
  for(int j=0;j<ncellsy;j++){file<<"3 1 "<<nx*(ny-1)-(j*nx)<<" "<<nx*(ny-2)-(j*nx)<<endl;}

  file<<endl;
  file<<"vertices"<<endl;
  file<<nx*ny<<endl;
  file<<"2"<<endl;

  file<<fixed<<setprecision(17);

  for(int j=0,k=0;j<ny;j++){
    for(int i=0;i<nx;i++,k++){
      file<<vx[k]<<" "<<vy[k]<<endl;
    }
  }
  file<<endl;

  cout<<"Output(): Done."<<endl;

  return;

}

// old version of output function - RhoLo can't read this !

void Output_old(char* outputfile,Input*I){

  cout<<"Output_old(): Writing input pipeline to meshfile "<<outputfile<<endl;

// write the input pipeline

  ofstream file(outputfile);

  file<<"#"<<I->Title()<<endl;

  for(int i=0;i<I->NDescriptions();i++){
    file<<"#"<<I->Description(i)<<endl;
  }

// set up vertices

  double vx[I->Nodes(0)],vy[I->Nodes(1)],vz[I->Nodes(2)];
  double dx((I->XMax()-I->XMin())/(I->Nodes(0)-1));
  double dy((I->YMax()-I->YMin())/(I->Nodes(1)-1));
  double dz((I->ZMax()-I->ZMin())/(I->Nodes(2)-1));

  vx[0]=I->XMin();
  for(long i=1;i<I->Nodes(0);i++){vx[i]=vx[i-1]+dx;}

  vy[0]=I->YMin();
  for(long i=1;i<I->Nodes(1);i++){vy[i]=vy[i-1]+dy;}

  vz[0]=I->ZMin();
  for(long i=1;i<I->Nodes(2);i++){vz[i]=vz[i-1]+dz;}

// write vertices

  file<<"# z-coordinate of each mesh vertex"<<endl;

  file<<fixed<<setprecision(17);

  for(long k=0;k<I->Nodes(2);k++){
    for(long j=0;j<I->Nodes(1);j++){
      for(long i=0;i<I->Nodes(0);i++){
        file<<vz[k]<<endl;
      }
    }
  }

  file<<"# y-coordinate of each mesh vertex"<<endl;

  file<<fixed<<setprecision(17);

  for(long k=0;k<I->Nodes(2);k++){
    for(long j=0;j<I->Nodes(1);j++){
      for(long i=0;i<I->Nodes(0);i++){
        file<<vy[j]<<endl;
      }
    }
  }

  file<<"# x-coordinate of each mesh vertex"<<endl;

  file<<fixed<<setprecision(17);

  for(long k=0;k<I->Nodes(2);k++){
    for(long j=0;j<I->Nodes(1);j++){
      for(long i=0;i<I->Nodes(0);i++){
        file<<vx[i]<<endl;
      }
    }
  }

// write vertex numbers

  file<<"# mesh vertex number in each element corner"<<endl;

  long nx(max(1,I->Nodes(0)-1)),ny(max(1,I->Nodes(1)-1)),nz(max(1,I->Nodes(2)-1));
  for(long k=0;k<nz;k++){
    for(long j=0;j<ny;j++){
      for(long i=0;i<nx;i++){
        file<<(j*I->Nodes(0))+i<<endl;
        file<<(j*I->Nodes(0))+i+1<<endl;
        file<<(j*I->Nodes(0))+i+1+nx+1<<endl;
        file<<(j*I->Nodes(0))+i+1+nx<<endl;
      }
    }
  }

// number of vertices

  file<<"# number of vertices in each element (the sum of this is the length of the above list)"<<endl;

  for(long i=0;i<nx;i++){
    for(long j=0;j<ny;j++){
      for(long k=0;k<nz;k++){
        file<<"4"<<endl;
      }
    }
  }

// materials

  file<<"# material number in each element"<<endl;

// set up element-wise material array

  long ncellsx(max(1,I->Nodes(0)-1)),ncellsy(max(1,I->Nodes(1)-1)),ncellsz(max(1,I->Nodes(2)-1));
  long totele(ncellsx*ncellsy*ncellsz);
  int elmat[ncellsx*ncellsy*ncellsz];

  long iel(0);
  for(long k=0;k<ncellsz;k++){
    for(long j=0;j<ncellsy;j++){
      for(long i=0;i<ncellsx;i++){

// place ambient material in this cell first

        elmat[iel]=I->Ambient();

// compute element centre coordinates

        double xc=vx[i]+0.5*(vx[i+1]-vx[i]);
        double yc=vy[j]+0.5*(vy[j+1]-vy[j]);
        double zc=vz[k]+0.5*(vz[k+1]-vz[k]);

// find which material the point (xc,yc,zc) intersects and place this in the cell

        for(int imat=0;imat<I->NMaterials();imat++){
          if((xc>I->Range(0,imat)&&xc<I->Range(1,imat))&&(yc>I->Range(2,imat)&&yc<I->Range(3,imat))&&(zc>I->Range(4,imat)&&zc<I->Range(5,imat))){
            elmat[iel]=I->Material(imat);
            break;
          }
        }

        iel++;
      }
    }
  }

// write out material in each element

  for(long iel=0;iel<totele;iel++){
    file<<elmat[iel]<<endl;
  }

// output polyhedral types

  file<<"# polyhedral type to insert in each element (address on type list below)"<<endl;

  for(long iel=0;iel<totele;iel++){
    file<<"1"<<endl;
  }

// output polyhedral order for the kinematics

  file<<"# polyhedral order to insert in each element for kinematics (1 = linear, 2 = quadratic, 3= cubic, etc...)"<<endl;

  for(long iel=0;iel<totele;iel++){
    file<<I->Pk()<<endl;
  }

// output polyhedral order for the thermodynamics

  file<<"# polyhedral order to insert in each element for thermodynamics (1 = linear, 2 = quadratic, 3= cubic, etc...)"<<endl;

  for(long iel=0;iel<totele;iel++){
    file<<I->Pt()<<endl;
  }

  file<<"# polyhderal orders present for kinematics"<<endl;
  file<<I->Pk()<<endl;

  file<<"# polyhderal orders present for thermodynamics"<<endl;
  file<<I->Pt()<<endl;

  file<<"# number of polyhedral orders present"<<endl;
  file<<"1"<<endl;
  file<<"# name of each polyhedral type"<<endl;
  file<<"QUAD"<<endl;
  file<<"# polyhedral element types present (2001 = QUAD, 2002 = TRI, 3001 = HEX, 3002 = TET)"<<endl;
  file<<"2001"<<endl;
  file<<"# number of polyhedral element types present"<<endl;
  file<<"1"<<endl;

// output material data

  file<<"# initial energy in each material"<<endl;
  for(int imat=0;imat<I->NMaterials();imat++){
    file<<"0.0001"<<endl; // code starts with pressure and inverts EoS to get energy, so values here should be irrelevant
  }

  file<<"# initial density in each material"<<endl;
  for(int imat=0;imat<I->NMaterials();imat++){
    file<<I->Density(imat)<<endl;
  }

  file<<"# initial pressure in each material"<<endl;
  for(int imat=0;imat<I->NMaterials();imat++){
    file<<I->Pressure(imat)<<endl;
  }

  file<<"# initial velocity in x-direction in each material"<<endl;
  for(int imat=0;imat<I->NMaterials();imat++){
    file<<I->Velocity(0,imat)<<endl;
  }

  file<<"# initial velocity in y-direction in each material"<<endl;
  for(int imat=0;imat<I->NMaterials();imat++){
    file<<I->Velocity(1,imat)<<endl;
  }

  file<<"# initial velocity in z-direction in each material"<<endl;
  for(int imat=0;imat<I->NMaterials();imat++){
    file<<I->Velocity(2,imat)<<endl;
  }

  file<<"# material numbers present in this mesh"<<endl;
  for(int imat=0;imat<I->NMaterials();imat++){
    file<<I->Material(imat)<<endl;
  }

  file<<"# number of materials"<<endl;
  file<<I->NMaterials()<<endl;

  file<<"# number of cells"<<endl;
  file<<totele<<endl;

  long nvertices=I->Nodes(0)*I->Nodes(1)*I->Nodes(2);
  file<<"# number of vertices"<<endl;
  file<<nvertices<<endl;

  file<<"# boundary condition on each face of the mesh (anticlockwise from bottom face)"<<endl;

  for(int iface=0;iface<6;iface++){
    if(I->Boundary(iface).compare("R")==0){
      file<<"REFLECTIVE"<<endl;
    }else{
      file<<"TRANSMISSIVE"<<endl;
    }
  }

  file<<"# length of each mesh dimension"<<endl;

  for(int idim=0;idim<I->NDims();idim++){
    file<<I->Nodes(idim)<<endl;
  }

  file<<"# number of mesh dimensions"<<endl;
  file<<I->NDims()<<endl;

  cout<<"Output_old(): Done."<<endl;

  return;

}

// return user name

string gethost(){

  char*hostname(getenv("HOSTNAME"));
  string computername;

  if(hostname!=NULL) {
    computername=hostname;
    hostname=NULL;
  }else{
    hostname=new char[512];
    if(gethostname(hostname,512)==0){ // success = 0, failure = -1
      computername=hostname;
    }
    delete[]hostname;
    hostname=NULL;
  }

  return computername;
}
