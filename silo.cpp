// silo graphics IO function

#include <iostream>
#include <vector>
#include "mesh.h"
#include "shape.h"
#include "silo.h"
#include <fstream>
#include <cmath>
#include <iomanip>
#include "riemann.h"

#define VD vector<double> // vector of doubles
#define FILENAME "graphics.silo" // define the output filename
#define FILEINFO "This is a silo file created by silo.cpp and it contains polyhedral meshes." // define info for reader app
#define MESHNAME "MyFirstMesh" // name of a test mesh
#define NSAMPLES 2000 // number of sample points

using namespace std;

void silo(Mesh*M,VD*x0,VD*d,VD*p,VD*m,VD*ec0,VD*V0,VD*u0,VD*e0,Shape*S[],int cycle,double time){

  DBfile*dbfile;
  int dberr;
  DBoptlist *optlist=NULL;

// number of dummy elements in each direction

  long ni( (M->Dim(0)-1)*(S[1]->order()+1)-1);
  long nj( (M->Dim(1)-1)*(S[1]->order()+1)-1);

//  const char*filename(FILENAME);
  const char*fileinfo(FILEINFO);
  const char*meshname(MESHNAME);

// vector<char> filename(str.c_str(), str.c_str() + str.size() + 1);
  string str1("step_"+to_string(cycle)),str2("profile_"+to_string(cycle));
  const char*filename(str1.append(".silo").c_str());

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

  double xcoords_k[M->KNodes()],ycoords_k[M->KNodes()],zcoords_k[M->KNodes()],xcoords_s[NSAMPLES],ycoords_s[NSAMPLES],zcoords_s[NSAMPLES];VD rx;

  for(long i=0;i<M->KNodes();i++){
//    xcoords_k[i]=M->KCoord0(0,i);
    ycoords_k[i]=M->KCoord0(1,i);
    zcoords_k[i]=M->KCoord0(2,i);
  }

// overide x-cordinate

  long ki(0);
  for(int j=0;j<(M->Dim(1)-1);j++){
    for(int i=1;i<M->Dim(0);i++){
      for(int jloc=0;jloc<S[1]->order()+1;jloc++){
        for(int iloc=0;iloc<S[1]->order()+1;iloc++){
          long ix(i*S[1]->nloc()+iloc);
          xcoords_k[ki]=(*x0)[ix];
          ki++;
        }
      }
    }
  }

// load mesh vertices on thermodynamic grid

  double xcoords_t[M->TNodes()];double ycoords_t[M->TNodes()];double zcoords_t[M->TNodes()];
  for(long i=0;i<M->TNodes();i++){
//    xcoords_t[i]=M->TCoord0(0,i);
    ycoords_t[i]=M->TCoord0(1,i);
    zcoords_t[i]=M->TCoord0(2,i);
  }

// overide x-cordinate

  long ti(0);
  for(int j=0;j<(M->Dim(1)-1);j++){
    for(int i=1;i<M->Dim(0);i++){
      for(int jloc=0;jloc<S[0]->order()+1;jloc++){
        for(int iloc=0;iloc<S[0]->order()+1;iloc++){
          long ix(i*S[0]->nloc()+iloc);
          xcoords_t[ti]=(*x0)[ix];
          ti++;
        }
      }
    }
  }

// sample point coordinates

  for(int i=0;i<NSAMPLES;i++){xcoords_s[i]=0.0+i/double(NSAMPLES);rx.push_back(xcoords_s[i]);ycoords_s[i]=0.0;zcoords_s[i]=0.0;}

// place coordinates into silo output format

  double *coords_k[]={(double*)xcoords_k,(double*)ycoords_k,(double*)zcoords_k};
  double *coords_t[]={(double*)xcoords_t,(double*)ycoords_t,(double*)zcoords_t};
  double *coords_s[]={(double*)xcoords_s,(double*)ycoords_s,(double*)zcoords_s};

// set up the number of nodes on the meshes

  int nnodes_k=M->KNodes();
  int nnodes_t=M->TNodes();

  int nzones(int((M->Dim(0)-1)*(M->Dim(1)-1)));
  int nzones_k(((M->Dim(0)-1)*(S[1]->order()+1)-1)*((M->Dim(1)-1)*(S[1]->order()+1)-1));
  int nzones_t(((M->Dim(0)-1)*(S[0]->order()+1)-1)*((M->Dim(1)-1)*(S[0]->order()+1)-1));

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
  double var4[NSAMPLES]; // sample point data

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
  dberr=DBPutPointmesh (dbfile,"Sample",1,coords_s,NSAMPLES,DB_DOUBLE,optlist);
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

//  for(long i=0;i<nzones_k;i++){var1[i]=(M->mElement[i]>=0)?M->Density(M->mElement[i]):M->Density(0);}
  for(long j=0;j<nj;j++){for(long i=0;i<ni;i++){var1[j*ni+i]=(M->mElement[i]>=0)?(*d)[M->mElement[i]+1]:(*d)[0];}}
  optlist = DBMakeOptlist(1);
  dberr=DBAddOption(optlist, DBOPT_UNITS, (void*)"g/cc");
  dberr=DBPutUcdvar1(dbfile,"density","Elements",var1,nzones_k,NULL,0,DB_DOUBLE,DB_ZONECENT,optlist);
  dberr=DBFreeOptlist(optlist);

//  for(long i=0;i<nzones_k;i++){var1[i]=(M->mElement[i]>=0)?M->Pressure(M->mElement[i]):M->Pressure(0);}
  for(long j=0;j<nj;j++){for(long i=0;i<ni;i++){var1[j*ni+i]=(M->mElement[i]>=0)?(*p)[M->mElement[i]+1]:(*p)[0];}}
  optlist = DBMakeOptlist(1);
  dberr=DBAddOption(optlist, DBOPT_UNITS, (void*)"Mb");
  dberr=DBPutUcdvar1(dbfile,"pressure","Elements",var1,nzones_k,NULL,0,DB_DOUBLE,DB_ZONECENT,optlist);
  dberr=DBFreeOptlist(optlist);

//  for(long i=0;i<nzones_k;i++){var1[i]=(M->mElement[i]>=0)?M->Mass(M->mElement[i]):M->Mass(0);}
  for(long j=0;j<nj;j++){for(long i=0;i<ni;i++){var1[j*ni+i]=(M->mElement[i]>=0)?(*m)[M->mElement[i]+1]:(*m)[0];}}
  optlist = DBMakeOptlist(1);
  dberr=DBAddOption(optlist, DBOPT_UNITS, (void*)"g");
  dberr=DBPutUcdvar1(dbfile,"mass","Elements",var1,nzones_k,NULL,0,DB_DOUBLE,DB_ZONECENT,optlist);
  dberr=DBFreeOptlist(optlist);

//  for(long i=0;i<nzones_k;i++){var1[i]=(M->mElement[i]>=0)?M->CCEnergy0(M->mElement[i]):M->CCEnergy0(0);}
  for(long j=0;j<nj;j++){for(long i=0;i<ni;i++){var1[j*ni+i]=(M->mElement[i]>=0)?(*ec0)[M->mElement[i]+1]:(*ec0)[0];}}
  optlist = DBMakeOptlist(1);
  dberr=DBAddOption(optlist, DBOPT_UNITS, (void*)"Mbcc");
  dberr=DBPutUcdvar1(dbfile,"pdv","Elements",var1,nzones_k,NULL,0,DB_DOUBLE,DB_ZONECENT,optlist);
  dberr=DBFreeOptlist(optlist);

//  for(long i=0;i<nzones_k;i++){var1[i]=(M->mElement[i]>=0)?M->Volume0(M->mElement[i]):M->Volume0(0);}
  for(long j=0;j<nj;j++){for(long i=0;i<ni;i++){var1[j*ni+i]=(M->mElement[i]>=0)?(*V0)[M->mElement[i]+1]:(*V0)[0];}}
  optlist = DBMakeOptlist(1);
  dberr=DBAddOption(optlist, DBOPT_UNITS, (void*)"cc");
  dberr=DBPutUcdvar1(dbfile,"volume","Elements",var1,nzones_k,NULL,0,DB_DOUBLE,DB_ZONECENT,optlist);
  dberr=DBFreeOptlist(optlist);

// map energy field onto kinematic mesh for the vis

//  double data_t[S[0]->nloc()];double data_k[S[1]->nloc()];
//  for(int i=0;i<nzones;i++){
//    for(int iloc=0;iloc<S[0]->nloc();iloc++){data_t[iloc]=M->Energy0(i*S[0]->nloc()+iloc);}
//    S[0]->prolongate(data_t,data_k,S[0]->order());
//    for(int iloc=0;iloc<S[1]->nloc();iloc++){var2[i*S[1]->nloc()+iloc]=data_k[iloc];}
//    for(int iloc=0;iloc<S[0]->nloc();iloc++){var2[i*S[0]->nloc()+iloc]=M->Energy0(i*S[0]->nloc()+iloc);} // ->prolongate not working correctly ??
//  }

  {long ti(0);for(int j=0;j<(M->Dim(1)-1);j++){for(int i=1;i<M->Dim(0);i++){for(int jloc=0;jloc<S[0]->order()+1;jloc++){for(int iloc=0;iloc<S[0]->order()+1;iloc++){long ix(i*S[0]->nloc()+iloc);var2[ti]=(*e0)[ix];ti++;}}}}}
  optlist = DBMakeOptlist(1);
  dberr=DBAddOption(optlist, DBOPT_UNITS, (void*)"Mbcc");
  dberr=DBPutUcdvar1(dbfile,"energy","Elements",var2,nnodes_k,NULL,0,DB_DOUBLE,DB_NODECENT,optlist);
  dberr=DBFreeOptlist(optlist);

// write vector data fields

//  for(int i=0;i<nnodes_k;i++){var2[i]=M->Velocitx0(0,i);}
  {long ki(0);for(int j=0;j<(M->Dim(1)-1);j++){for(int i=1;i<M->Dim(0);i++){for(int jloc=0;jloc<S[1]->order()+1;jloc++){for(int iloc=0;iloc<S[1]->order()+1;iloc++){long ix(i*S[1]->nloc()+iloc);var2[ki]=(*u0)[ix];ki++;}}}}}
  optlist = DBMakeOptlist(1);
  dberr=DBAddOption(optlist, DBOPT_UNITS, (void*)"cm/us");
  dberr=DBPutUcdvar1(dbfile,"velocity_x","Elements",var2,nnodes_k,NULL,0,DB_DOUBLE,DB_NODECENT,optlist); 
  dberr=DBFreeOptlist(optlist);

//  for(int i=0;i<nnodes_k;i++){var2[i]=M->Velocity0(1,i);}
  {long ki(0);for(int j=0;j<(M->Dim(1)-1);j++){for(int i=1;i<M->Dim(0);i++){for(int jloc=0;jloc<S[1]->order()+1;jloc++){for(int iloc=0;iloc<S[1]->order()+1;iloc++){long ix(i*S[1]->nloc()+iloc);var2[ki]=0.0;ki++;}}}}}
  optlist = DBMakeOptlist(1);
  dberr=DBAddOption(optlist, DBOPT_UNITS, (void*)"cm/us");
  dberr=DBPutUcdvar1(dbfile,"velocity_y","Elements",var2,nnodes_k,NULL,0,DB_DOUBLE,DB_NODECENT,optlist); 
  dberr=DBFreeOptlist(optlist);

// sample the exact solution

  double l[3]={1.0,0.0,1.0},r[3]={0.125,0.0,0.1};       // left/right flux states for Riemann solver
  Riemann R(Riemann::exact,l,r);R.profile(&rx,time);

  for(int i=0;i<NSAMPLES;i++){var4[i]=R.energy(i);}
  optlist = DBMakeOptlist(1);
  dberr=DBAddOption(optlist, DBOPT_UNITS, (void*)"Mbcc");
  dberr=DBPutPointvar1(dbfile,"exact_energy","Sample",var4,NSAMPLES,DB_DOUBLE,optlist);
  dberr=DBFreeOptlist(optlist);

  for(int i=0;i<NSAMPLES;i++){var4[i]=R.density(i);}
  optlist = DBMakeOptlist(1);
  dberr=DBAddOption(optlist, DBOPT_UNITS, (void*)"g/cc");
  dberr=DBPutPointvar1(dbfile,"exact_density","Sample",var4,NSAMPLES,DB_DOUBLE,optlist);
  dberr=DBFreeOptlist(optlist);

  for(int i=0;i<NSAMPLES;i++){var4[i]=R.pressure(i);}
  optlist = DBMakeOptlist(1);
  dberr=DBAddOption(optlist, DBOPT_UNITS, (void*)"Mb");
  dberr=DBPutPointvar1(dbfile,"exact_pressure","Sample",var4,NSAMPLES,DB_DOUBLE,optlist);
  dberr=DBFreeOptlist(optlist);

  for(int i=0;i<NSAMPLES;i++){var4[i]=R.velocity(i);}
  optlist = DBMakeOptlist(1);
  dberr=DBAddOption(optlist, DBOPT_UNITS, (void*)"cm/us");
  dberr=DBPutPointvar1(dbfile,"exact_velocity","Sample",var4,NSAMPLES,DB_DOUBLE,optlist);
  dberr=DBFreeOptlist(optlist);

// close the file

  dberr=DBClose(dbfile);

  return;

}
