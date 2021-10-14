// silo graphics IO function

// Author S. R. Merton

#define FILENAME "silo/graphics.silo"
#define CODENAME "rholo"
#define VERSION "1.0"
#define TITLE "MyFirstTitle"
#define FILEINFO "This is a silo file created by rholo and it contains polyhedral meshes." // define info for reader app
#define MESHNAME "MyFirstMesh"   // name of a test mesh
#define KNOD i*(K->nloc()-1)+k   // global node number on kinematic mesh
#define TNOD i*T->nloc()+j       // global node number on thermodynamic mesh
#define VD vector<double>        // vector of doubles
#define VI vector<int>           // vector of ints
#define GPNT i*K->ngi()+gi       // global address of Gauss point gi in element i

#include <iostream>
#include <vector>
#include <ctime>
#include <filesystem>
#include <algorithm>
#include "silo.h"
#include "shape.h"

// function signatures

std::string date();

using namespace std;

void silo(VD*x,VD*d,VD*p,VD*e,VD*q,VD*c,VD*u,VI*m,int step,double time,Shape*K,Shape*T){

  DBfile*dbfile;
  DBoptlist *optlist=NULL;
  int dberr;
  string str1("silo/step_"+to_string(step)+".silo");

  const char*filename(str1.c_str());
  const char*codename(CODENAME);
  const char*version(VERSION);
  const char*title(TITLE);
  const char*fileinfo(FILEINFO);
  const char*meshname(MESHNAME);

  int nx((x->size()-1)/K->order());    // number of cells in x direction
  int ny(1);                           // number of cells in y direction
  int nzones = nx*ny;                  // number of zones
  int nzonesg = nzones*K->ngi();       // number of gauss point zones
  int ndims = 2;                       // number of dimensions
  int nnodesk = 2*x->size();           // number of nodes on kinematic grid
  int nnodest = 2*nzones*T->nloc();    // number of nodes on thermodynamic grid
  int nnodesg = 2*(nzonesg+1);         // number of gauss points
  int origin = 0;                      // first address in nodelist arrays
  int lnodelist=(2*K->nloc()*nzones);  // length of the node list
  int lnodelistg=(4*nzonesg);          // length of the gauss point node list
  int nodelist[lnodelist];             // the node list
  int nodelistg[lnodelistg];           // the gauss point node list
  int nshapetypes(1);                  // number of different shape types on the mesh
  int shapesize[nshapetypes];          // number of nodes defining each shape
  int shapesizeg[nshapetypes];         // number of nodes defining each gauss shape
  int shapecounts[nshapetypes];        // number of zones of each shape type
  int shapecountsg[nshapetypes];       // number of gauss zones of each shape type

// set up material data structure

  int nmat(*max_element(m->begin(),m->end())); // number of materials
  int matdims[]={nx,ny};                       // material dimensions
  char*matname[nmat];                          // material names
  int matnos[nmat];                            // materials numbers present
  int matlist[nzones];                         // material number in each zone
  int mixlen(0);                               // number of mixed cells
  int mix_next[mixlen];                        // indices into mixed data arrays
  int mix_mat[mixlen];                         // material numbers for mixed zones
  int mix_vf[mixlen];                          // volunme fractions

// output arrays

  int elnos[nzones];                           // element numbers
  long nknos[nnodesk];                         // node numbers
  double var1[nnodesk];                        // zone centred scalars
  double varg[nzonesg];                        // gauss centred scalars

  cout<<"Writing a silo graphics dump to file "<<filename<<endl;

// store coordinates in correct format for silo and repeat for each mesh

  double xcoordsk[nnodesk];for(int i=0;i<x->size();i++){xcoordsk[i]=x->at(i);};for(int i=0;i<x->size();i++){xcoordsk[x->size()+i]=x->at(i);}
  double ycoordsk[nnodesk];for(int i=0;i<x->size();i++){ycoordsk[i]=0.0;};for(int i=0;i<x->size();i++){ycoordsk[x->size()+i]=0.25;}
  double *coordsk[]={xcoordsk,ycoordsk};

// thermodynamic coordinates

  double xcoordst[nnodest],ycoordst[nnodest];
  for(int i=0;i<nzones;i++){
      for(int j=0;j<T->nloc();j++){
        double pos(-1.0+j*2.0/(T->nloc()-1));xcoordst[i*T->nloc()+j]=0.0;
        for(int k=0;k<K->nloc();k++){
          xcoordst[i*T->nloc()+j]+=K->value(k,pos)*xcoordsk[i*(K->nloc()-1)+k];
        }
        xcoordst[(nzones*T->nloc())+i*T->nloc()+j]=xcoordst[i*T->nloc()+j];
      }
  }
  for(int i=0;i<nzones*T->nloc();i++){ycoordst[i]=0.0;};for(int i=0;i<nzones*T->nloc();i++){ycoordst[nzones*T->nloc()+i]=0.25;}
  double *coordst[]={xcoordst,ycoordst};

// Gauss point coordinates

  double xcoordsg[nnodesg]={},ycoordsg[nnodesg]={};
  for(int i=0,k=0;i<nzones;i++){

    k=0;xcoordsg[i*(K->ngi())]=x->at(KNOD);
    k=K->nloc()-1;xcoordsg[(i+1)*(K->ngi())]=x->at(KNOD);

// coordinates of the Gauss points

    double xgi[K->ngi()]={};
    for(int gi=0;gi<K->ngi();gi++){
      for(int k=0;k<K->nloc();k++){
        xgi[gi]+=K->value(k,gi)*x->at(KNOD);
      }
    }

// centre the node between adjacent Gauss points

    double xpos[K->ngi()-1]={};
    for(int gi=1;gi<K->ngi();gi++){
      xcoordsg[i*(K->ngi())+gi]=0.5*(xgi[gi]+xgi[gi-1]);
    }

  }

  for(int i=0;i<=nzonesg;i++){xcoordsg[nzonesg+1+i]=xcoordsg[i];}
  for(int i=0;i<=nzonesg;i++){ycoordsg[i]=0.0;};for(int i=0;i<=nzonesg;i++){ycoordsg[nzonesg+1+i]=0.25;}
  double *coordsg[]={xcoordsg,ycoordsg};

// connectivities

  for(int i=0;i<nx;i++){
    for(int j=0;j<2;j++){
      for(int k=0;k<K->nloc();k++){
        if(j==0){
          nodelist[2*i*K->nloc()+(j*K->nloc())+k]=i*(K->nloc()-1)+k;
        }else{
          nodelist[2*i*K->nloc()+(j*K->nloc())+k]=(i+nx+1)*(K->nloc()-1)-k+1;
        }
      }
    }
  }

  for(int i=0;i<nzonesg;i++){
    for(int j=0;j<2;j++){
      for(int k=0;k<2;k++){
        if(j==0){
          nodelistg[4*i+(2*j)+k]=i+k;
        }else{
          nodelistg[4*i+(2*j)+k]=(i+nzonesg+1)-k+1;
        }
      }
    }
  }

// zone shapes

  for(int i=0;i<nshapetypes;i++){shapesize[i]=2*K->nloc();shapesizeg[i]=4;}
  for(int i=0;i<nshapetypes;i++){shapecounts[i]=nzones;shapecountsg[i]=nzonesg;}

// material numbers

  for(int i=0;i<nmat;i++){matnos[i]=i+1;}
  for(int i=0;i<nzones;i++){matlist[i]=m->at(i);}
  for(int i=0;i<nmat;i++){matname[i]="Air";}

// disengage deprecation signalling

  dberr=DBSetDeprecateWarnings(0);

// create the silo database - this opens it also

  dbfile=DBCreate(filename,DB_CLOBBER,DB_LOCAL,fileinfo,DB_HDF5);

// write out connectivity information. //

  dberr=DBPutZonelist(dbfile, "zonelist", nzones, ndims, nodelist, lnodelist, origin ,shapesize, shapecounts, nshapetypes);
  dberr=DBPutZonelist(dbfile, "zonelistg", nzonesg, ndims, nodelistg, lnodelistg, origin ,shapesizeg, shapecountsg, nshapetypes);

// write out the meshes //

  optlist=DBMakeOptlist(2);
  dberr=DBAddOption(optlist,DBOPT_DTIME,&time);
  dberr=DBAddOption(optlist,DBOPT_CYCLE,&step);
  dberr=DBPutUcdmesh(dbfile,"Elements",ndims,NULL,coordsk,nnodesk,nzones,"zonelist",NULL,DB_DOUBLE,optlist);
  dberr=DBPutPointmesh(dbfile,"Kinematics",ndims,coordsk,nnodesk,DB_DOUBLE,optlist);
  dberr=DBPutPointmesh(dbfile,"Thermodynamics",ndims,coordst,nnodest,DB_DOUBLE,optlist);
  dberr=DBPutUcdmesh(dbfile,"Gauss",ndims,NULL,coordsg,nnodesg,nzonesg,"zonelistg",NULL,DB_DOUBLE,optlist);
  dberr=DBFreeOptlist(optlist);

// write out the material

  optlist=DBMakeOptlist(1);
  dberr=DBAddOption(optlist,DBOPT_MATNAMES,matname);
  dberr=DBPutMaterial(dbfile,"Materials","Elements",nmat,matnos,matlist,matdims,ndims,NULL,NULL,NULL,NULL,mixlen,DB_INT,optlist);
  dberr=DBFreeOptlist(optlist);

// write the element numbers

  for(int i=0;i<nzones;i++){elnos[i]=i;}
  optlist = DBMakeOptlist(1);
  dberr=DBAddOption(optlist, DBOPT_UNITS, (void*)"");
  dberr=DBPutUcdvar1(dbfile,"elementID","Elements",elnos,nzones,NULL,0,DB_INT,DB_ZONECENT,optlist);
  dberr=DBFreeOptlist(optlist);

// write the node numbers

  for(long i=0;i<nnodesk;i++){nknos[i]=i;}
  optlist = DBMakeOptlist(1);
  dberr=DBAddOption(optlist, DBOPT_UNITS, (void*)"");
  dberr=DBPutUcdvar1(dbfile,"nodeID","Elements",nknos,nnodesk,NULL,0,DB_LONG,DB_NODECENT,optlist);
  dberr=DBFreeOptlist(optlist);

// write scalar data fields

  for(long i=0;i<nzonesg;i++){varg[i]=d->at(i);}
  optlist = DBMakeOptlist(1);
  dberr=DBAddOption(optlist, DBOPT_UNITS, (void*)"g/cc");
  dberr=DBPutUcdvar1(dbfile,"density","Gauss",varg,nzonesg,NULL,0,DB_DOUBLE,DB_ZONECENT,optlist);
  dberr=DBFreeOptlist(optlist);

  for(long i=0;i<nzonesg;i++){varg[i]=p->at(i);}
  optlist = DBMakeOptlist(1);
  dberr=DBAddOption(optlist, DBOPT_UNITS, (void*)"Mb");
  dberr=DBPutUcdvar1(dbfile,"pressure","Gauss",varg,nzonesg,NULL,0,DB_DOUBLE,DB_ZONECENT,optlist);
  dberr=DBFreeOptlist(optlist);

  for(long i=0;i<nzonesg;i++){varg[i]=q->at(i);}
  optlist = DBMakeOptlist(1);
  dberr=DBAddOption(optlist, DBOPT_UNITS, (void*)"Mb");
  dberr=DBPutUcdvar1(dbfile,"bulk_q","Gauss",varg,nzonesg,NULL,0,DB_DOUBLE,DB_ZONECENT,optlist);
  dberr=DBFreeOptlist(optlist);

  for(long i=0;i<nzonesg;i++){varg[i]=c->at(i);}
  optlist = DBMakeOptlist(1);
  dberr=DBAddOption(optlist, DBOPT_UNITS, (void*)"cm/s");
  dberr=DBPutUcdvar1(dbfile,"sound_speed","Gauss",varg,nzonesg,NULL,0,DB_DOUBLE,DB_ZONECENT,optlist);
  dberr=DBFreeOptlist(optlist);

  for(int i=0;i<nzones;i++){for(int gi=0;gi<K->ngi();gi++){varg[GPNT]=0.0;for(int j=0;j<T->nloc();j++){varg[GPNT]+=T->value(j,gi)*e->at(TNOD);}}}
  optlist = DBMakeOptlist(1);
  dberr=DBAddOption(optlist, DBOPT_UNITS, (void*)"Mbcc");
  dberr=DBPutUcdvar1(dbfile,"energy","Gauss",varg,nzonesg,NULL,0,DB_DOUBLE,DB_ZONECENT,optlist);
  dberr=DBFreeOptlist(optlist);

// write vector data fields

  for(long i=0;i<x->size();i++){var1[i]=u->at(i);};for(int i=0;i<x->size();i++){var1[x->size()+i]=u->at(i);}
  optlist = DBMakeOptlist(1);
  dberr=DBAddOption(optlist, DBOPT_UNITS, (void*)"cm/s");
  dberr=DBPutUcdvar1(dbfile,"velocity_x","Elements",var1,nnodesk,NULL,0,DB_DOUBLE,DB_NODECENT,optlist); 
  dberr=DBFreeOptlist(optlist);

  for(long i=0;i<x->size();i++){var1[i]=0.0;};for(int i=0;i<x->size();i++){var1[x->size()+i]=0.0;}
  optlist = DBMakeOptlist(1);
  dberr=DBAddOption(optlist, DBOPT_UNITS, (void*)"cm/s");
  dberr=DBPutUcdvar1(dbfile,"velocity_y","Elements",var1,nnodesk,NULL,0,DB_DOUBLE,DB_NODECENT,optlist); 
  dberr=DBFreeOptlist(optlist);

// close the database

  dberr=DBClose(dbfile);

  return;

}

// return today's date

std::string date(){

  time_t now = time(0);

  return ctime(&now);

}
