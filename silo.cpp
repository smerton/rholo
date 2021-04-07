// silo graphics IO function

#include <iostream>
#include <vector>
#include "mesh.h"
#include "shape.h"
#include "silo.h"
#include <fstream>
#include <cmath>
#include <iomanip>

#define FILENAME "graphics.silo" // define the output filename
#define FILEINFO "This is a silo file created by silo.cpp and it contains polyhedral meshes." // define info for reader app
#define MESHNAME "MyFirstMesh" // name of a test mesh
#define TOL 1.0e-8 // tolerance for adding an element to the profile

void profiler(string*str,Mesh*M,Shape*S[]); // profiler

using namespace std;

void silo(Mesh*M,Shape*S[],int cycle,double time){

  DBfile*dbfile;
  int dberr;
  DBoptlist *optlist=NULL;

//  const char*filename(FILENAME);
  const char*fileinfo(FILEINFO);
  const char*meshname(MESHNAME);

// vector<char> filename(str.c_str(), str.c_str() + str.size() + 1);
  string str1("step_"+to_string(cycle)),str2("profile_"+to_string(cycle));
  const char*filename(str1.append(".silo").c_str());

// start profiler

  profiler(&str2,M,S);

// only one polyhedral order on the mesh is coded for here

  if(M->NOrders()>1){
    cout<<"silo(): only handling 1 polyhedral order, but there are "<<M->NOrders()<<" present."<<endl;
    return;
  }

//  only 2D meshes have been set up here

  if(M->Dim(M->NDims()-1)>1){
    cout<<"silo(): 3D meshes are not coded for, "<<M->Dim(M->NDims()-1)<<" nodes detected in 3rd dimension."<<endl;
    return;
  }

// load mesh vertices on kinematic grid

  double xcoords_k[M->KNodes()];double ycoords_k[M->KNodes()];double zcoords_k[M->KNodes()];
  for(long i=0;i<M->KNodes();i++){
    xcoords_k[i]=M->KCoord0(0,i);
    ycoords_k[i]=M->KCoord0(1,i);
    zcoords_k[i]=M->KCoord0(2,i);
  }

// load mesh vertices on thermodynamic grid

  double xcoords_t[M->TNodes()];double ycoords_t[M->TNodes()];double zcoords_t[M->TNodes()];
  for(long i=0;i<M->TNodes();i++){
    xcoords_t[i]=M->TCoord0(0,i);
    ycoords_t[i]=M->TCoord0(1,i);
    zcoords_t[i]=M->TCoord0(2,i);
  }

// place coordinates into silo output format

  double *coords_k[]={(double*)xcoords_k,(double*)ycoords_k,(double*)zcoords_k};
  double *coords_t[]={(double*)xcoords_t,(double*)ycoords_t,(double*)zcoords_t};

// set up the number of nodes on the meshes

  int nnodes_k=M->KNodes();
  int nnodes_t=M->TNodes();

  int nzones(int((M->Dim(0)-1)*(M->Dim(1)-1)));
  int nzones_k(((M->Dim(0)-1)*(S[1]->Order()+1)-1)*((M->Dim(1)-1)*(S[1]->Order()+1)-1));
  int nzones_t(((M->Dim(0)-1)*(S[0]->Order()+1)-1)*((M->Dim(1)-1)*(S[0]->Order()+1)-1));

  int ndims(3);
  int lnodelist_k(nzones_k*4); // includes dummies
  int lnodelist_t(nzones_t*4); // includes dummies
  int nshapes(1);

  int shapesize[nshapes];
  int shapecnt[nshapes];
  int shapecnt_k[nshapes];
  int shapecnt_t[nshapes];
  int shapetype[nshapes];

// setup data output arrays

  double var1[nzones_k]; // zone centred kinematic fields
  double var2[nnodes_k]; // node centred kinematic fields
  long var3[nnodes_k];   // node centred longs

// set up the material data structure

  int nmat(M->NMaterials());
  int matnos[nmat];
  int matlist[nzones_k];
  int nmatdims(1);
  int lmatlist[nmatdims];
  int mixlen(0);
  int mix_next[mixlen];
  int mix_mat[mixlen];
  int mix_vf[mixlen];

  lmatlist[0]=nzones;

// store material numbers present on the mesh

  for(int imat=0;imat<nmat;imat++){matnos[imat]=M->Materials(imat);}

// store material number in each zone

  for(long i=0;i<nzones_k;i++){matlist[i]=(M->mElement[i]>=0)?M->Material(M->mElement[i]):M->Material(0);}

// create a nodelist to display the kinematic data

  int nodelist_k[lnodelist_k];

  for(long ilist=0;ilist<M->lKNodeList();ilist++){
    nodelist_k[ilist]=M->KNodeList(ilist);
  }

// create a nodelist to display thethermodynamic data

  int nodelist_t[lnodelist_t];

  for(long ilist=0;ilist<M->lTNodeList();ilist++){
    nodelist_t[ilist]=M->TNodeList(ilist);
  }

// set the number of nodes per shape

  shapesize[0]=4;

// set the number of elements containing each shape

  shapecnt[0]=nzones;shapecnt_k[0]=nzones_k;shapecnt_t[0]=nzones_t;

// set the shape types

  shapetype[0]=DB_ZONETYPE_QUAD;

// create the silo database - this opens it also

  dbfile=DBCreate(filename,DB_CLOBBER,DB_LOCAL,fileinfo,DB_HDF5);

// create a zonelist for the unstructured grid ucd meshes

  dberr=DBPutZonelist2(dbfile,"zlist_k",nzones_k,M->NDims(),nodelist_k,lnodelist_k,0,0,0,shapetype,shapesize,shapecnt_k,nshapes,NULL);
  dberr=DBPutZonelist2(dbfile,"zlist_t",nzones_t,M->NDims(),nodelist_t,lnodelist_t,0,0,0,shapetype,shapesize,shapecnt_t,nshapes,NULL);

// write the meshes

  optlist = DBMakeOptlist(2);
  dberr=DBAddOption(optlist, DBOPT_DTIME, &time);
  dberr=DBAddOption(optlist, DBOPT_CYCLE, &cycle);
  dberr=DBPutPointmesh (dbfile,"Kinematics",M->NDims(),coords_k,nnodes_k,DB_DOUBLE,optlist);
  dberr=DBPutPointmesh (dbfile,"Thermodynamics",M->NDims(),coords_t,nnodes_t,DB_DOUBLE,optlist);
  dberr=DBPutUcdmesh   (dbfile,"Elements",M->NDims(),NULL,coords_k,nnodes_k,nzones_k,"zlist_k",NULL,DB_DOUBLE,optlist);
  dberr=DBFreeOptlist(optlist);

// write the materials

  dberr=DBPutMaterial(dbfile,"Materials","Elements",nmat,matnos,matlist,&nzones_k,1,NULL,NULL,NULL,NULL,mixlen,DB_INT,NULL);

// write the element numbers

  for(long i=0;i<nzones_k;i++){matlist[i]=M->mElement[i];}
  for(long i=0;i<nzones_k;i++){matlist[i]=(M->mElement[i]>=0)?M->mElement[i]:-1;}
  optlist = DBMakeOptlist(1);
  dberr=DBAddOption(optlist, DBOPT_UNITS, (void*)"");
  dberr=DBPutUcdvar1(dbfile,"elementID","Elements",matlist,nzones_k,NULL,0,DB_INT,DB_ZONECENT,optlist);
  dberr=DBFreeOptlist(optlist);

// write the node numbers

  for(long i=0;i<nnodes_k;i++){var3[i]=i;}
  optlist = DBMakeOptlist(1);
  dberr=DBAddOption(optlist, DBOPT_UNITS, (void*)"");
  dberr=DBPutUcdvar1(dbfile,"nodeID","Elements",var3,nnodes_k,NULL,0,DB_LONG,DB_NODECENT,optlist);
  dberr=DBFreeOptlist(optlist);

// write scalar data fields

  for(long i=0;i<nzones_k;i++){var1[i]=(M->mElement[i]>=0)?M->Density(M->mElement[i]):M->Density(0);}
  optlist = DBMakeOptlist(1);
  dberr=DBAddOption(optlist, DBOPT_UNITS, (void*)"g/cc");
  dberr=DBPutUcdvar1(dbfile,"density","Elements",var1,nzones_k,NULL,0,DB_DOUBLE,DB_ZONECENT,optlist);
  dberr=DBFreeOptlist(optlist);

  for(long i=0;i<nzones_k;i++){var1[i]=(M->mElement[i]>=0)?M->Pressure(M->mElement[i]):M->Pressure(0);}
  optlist = DBMakeOptlist(1);
  dberr=DBAddOption(optlist, DBOPT_UNITS, (void*)"Mb");
  dberr=DBPutUcdvar1(dbfile,"pressure","Elements",var1,nzones_k,NULL,0,DB_DOUBLE,DB_ZONECENT,optlist);
  dberr=DBFreeOptlist(optlist);

  for(long i=0;i<nzones_k;i++){var1[i]=(M->mElement[i]>=0)?M->Mass(M->mElement[i]):M->Mass(0);}
  optlist = DBMakeOptlist(1);
  dberr=DBAddOption(optlist, DBOPT_UNITS, (void*)"g");
  dberr=DBPutUcdvar1(dbfile,"mass","Elements",var1,nzones_k,NULL,0,DB_DOUBLE,DB_ZONECENT,optlist);
  dberr=DBFreeOptlist(optlist);

  for(long i=0;i<nzones_k;i++){var1[i]=(M->mElement[i]>=0)?M->CCEnergy0(M->mElement[i]):M->CCEnergy0(0);}
  optlist = DBMakeOptlist(1);
  dberr=DBAddOption(optlist, DBOPT_UNITS, (void*)"Mbcc");
  dberr=DBPutUcdvar1(dbfile,"pdv","Elements",var1,nzones_k,NULL,0,DB_DOUBLE,DB_ZONECENT,optlist);
  dberr=DBFreeOptlist(optlist);

  for(long i=0;i<nzones_k;i++){var1[i]=(M->mElement[i]>=0)?M->Volume0(M->mElement[i]):M->Volume0(0);}
  optlist = DBMakeOptlist(1);
  dberr=DBAddOption(optlist, DBOPT_UNITS, (void*)"cc");
  dberr=DBPutUcdvar1(dbfile,"volume","Elements",var1,nzones_k,NULL,0,DB_DOUBLE,DB_ZONECENT,optlist);
  dberr=DBFreeOptlist(optlist);

// map energy field onto kinematic mesh for the vis

  double data_t[S[0]->Nloc()];double data_k[S[1]->Nloc()];
  for(int i=0;i<nzones;i++){
    for(int iloc=0;iloc<S[0]->Nloc();iloc++){data_t[iloc]=M->Energy0(i*S[0]->Nloc()+iloc);}
    S[0]->Prolongate(data_t,data_k);
    for(int iloc=0;iloc<S[1]->Nloc();iloc++){var2[i*S[1]->Nloc()+iloc]=data_k[iloc];}
  }

  optlist = DBMakeOptlist(1);
  dberr=DBAddOption(optlist, DBOPT_UNITS, (void*)"Mbcc");
  dberr=DBPutUcdvar1(dbfile,"energy","Elements",var2,nnodes_k,NULL,0,DB_DOUBLE,DB_NODECENT,optlist);
  dberr=DBFreeOptlist(optlist);

// write vector data fields

  for(int i=0;i<nnodes_k;i++){var2[i]=M->Velocity0(0,i);}
  optlist = DBMakeOptlist(1);
  dberr=DBAddOption(optlist, DBOPT_UNITS, (void*)"cm/us");
  dberr=DBPutUcdvar1(dbfile,"velocity_x","Elements",var2,nnodes_k,NULL,0,DB_DOUBLE,DB_NODECENT,optlist); 
  dberr=DBFreeOptlist(optlist);

  for(int i=0;i<nnodes_k;i++){var2[i]=M->Velocity0(1,i);}
  optlist = DBMakeOptlist(1);
  dberr=DBAddOption(optlist, DBOPT_UNITS, (void*)"cm/us");
  dberr=DBPutUcdvar1(dbfile,"velocity_y","Elements",var2,nnodes_k,NULL,0,DB_DOUBLE,DB_NODECENT,optlist); 
  dberr=DBFreeOptlist(optlist);

// close the file

  dberr=DBClose(dbfile);

//  cout<<"silo(): done."<<endl;

  return;

}

// this peripheral function profiles the denisty, pressure and energy fields

void profiler(string*str,Mesh*M,Shape*S[]){

  const char*filename((*str).append(".dat").c_str());

  ofstream profile(filename);

// define the profile

  double yp=0.01875,zp=0.0;

// header

  profile<<"# This file contains a profile along the line (x,"<<yp<<")"<<endl;
  profile<<"# CELL CENTRE X     DENSITY             PRESSURE            ENERGY"<<endl;

// set required precision

  profile<<fixed<<setprecision(17);

// sweep the mesh

  for(long iel=0;iel<M->NCells();iel++){

// FE interpolation to get the coordinates of the element centre

    double xc(0.0),yc(0.0),zc(0.0);
    for(int iloc=0;iloc<M->KLoc(iel);iloc++){
      xc+=S[1]->Value(iloc,0.0,0.0)*M->KCoord(0,M->KNode(iel,iloc));
      yc+=S[1]->Value(iloc,0.0,0.0)*M->KCoord(1,M->KNode(iel,iloc));
      zc+=S[1]->Value(iloc,0.0,0.0)*M->KCoord(2,M->KNode(iel,iloc));
    }

// add element to the profile

    if(abs(yp-yc)<TOL&&abs(zp-zc)<TOL){
      profile<<xc<<" "<<M->Density(iel)<<" "<<M->Pressure(iel)<<" "<<M->CCEnergy(iel)<<endl;
    }

  }

  return;

}
