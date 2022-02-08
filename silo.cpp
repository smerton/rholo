// silo graphics IO function

// Author S. R. Merton

#define FILENAME "silo/graphics.silo"
#define CODENAME "rholo"
#define VERSION "1.0"
#define TITLE "MyFirstTitle"
#define FILEINFO "This is a silo file created by rholo and it contains polyhedral meshes." // define info for reader app
#define MESHNAME "MyFirstMesh"      // name of a test mesh
#define KNOD i*(S->nloc()-1)+k      // global node number on kinematic mesh
#define TNOD i*S->nloc()+j          // global node number on thermodynamic mesh
#define VD vector<double>           // vector of doubles
#define VVD vector<vector<double> > // vector of vector of doubles
#define VI vector<int>              // vector of ints

#include <iostream>
#include <vector>
//#include <ctime>
#include <filesystem>
#include <algorithm>
#include "silo.h"
#include "shape.h"
#include "mesh.h"
#include <cmath>

// function signatures

std::string date();

using namespace std;

// local function to return the number of sample points and set up their global numbers

long NSampleNodes(Mesh const &M, int const &nsubs, vector<vector<long> > &SampleNode);

// local function to interpolate the vector vin on to the vector vout at sample points

void sample(Mesh const &M,int const &nsubs,Shape const &T, vector<double> const &vin,double vout[]);

void silo(VVD const &x,VVD const &xt,VD const &d,VD const &p,VD const &e,VD const &q,VD const &c,
          VVD const &u,VI const &m,int step,double time,Mesh const &M,VD const &g,Shape const &S,Shape const &T){

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

  cout<<"       silo(): Writing a silo graphics dump to file "<<filename<<endl;

// set up material data structure

  int nmat(M.NMaterials());             // number of materials
  int matdims[]={M.NCells(),1};         // material dimensions
  char*matname[nmat];                   // material names
  int matnos[nmat];                     // materials numbers present
  int *matlist;                         // material number in each zone
  int mixlen(0);                        // number of mixed cells
  int mix_next[mixlen];                 // indices into mixed data arrays
  int mix_mat[mixlen];                  // material numbers for mixed zones
  int mix_vf[mixlen];                   // volunme fractions

// reuse the following, several different meshes are required to support the different lengths of data
// some of these meshes subdivide the elements and so the following data can vary for the different meshes

  int nzones;                       // number of cells
  int nsubs;                        // number of subcells
  int ndims(M.NDims());             // number of dimensions
  long lnodelist;                   // length of the nodelist array
  int *nodelist;                    // nodelist array
  long nnodes;                      // number of nodes in the mesh
  int origin(0);                    // first address in nodelist arrays
  int nshapetypes(1);               // number of different shape types on the mesh
  int shapesize[nshapetypes];       // number of nodes defining each shape
  int shapecounts[nshapetypes];     // number of zones of each shape type
  double *xcoords;                  // stores coordinates in correct format for silo
  double *ycoords;                  // stores coordinates in correct format for silo
  double *coords[ndims];            // array of length ndims pointing to the coordinate arrays
  int *elnos;                       // element numnbers
  long *nodnos;                     // node numbers
  double *var1;                     // nodal variable to output against the mesh
  vector<vector<long> > SampleNode; // global numbers of the sample points
  vector<double> vin;               // vector to sample at the mesh sample points
  double *vout;              // data at the mesh sample points

// disengage deprecation signalling

  dberr=DBSetDeprecateWarnings(0);

// create the silo database - this opens it also

#ifdef USE_HDF
  dbfile=DBCreate(filename,DB_CLOBBER,DB_LOCAL,fileinfo,DB_HDF5);
#else
  dbfile=DBCreate(filename,DB_CLOBBER,DB_LOCAL,fileinfo,DB_PDB);
#endif

// draw a mesh called "Elements" to support the material numbers and element numbers

  Shape N(1);                       // use a 4-node quad mesh to support material list
  nsubs=1;                          // number of subcells
  nzones=M.NCells();                // number of cells
  nnodes=M.NNodes();                // number of nodes
  ndims=M.NDims();                  // number of dimensions
  lnodelist=N.nloc()*nzones;        // length of the nodelist
  nodelist=new int[lnodelist];      // allocate nodelist array
  xcoords=new double[nnodes];       // mesh vertex coordinates
  ycoords=new double[nnodes];       // mesh vertex coordinates
  matlist=new int[nzones];          // material numbers
  elnos=new int[nzones];            // element numbers
  var1=new double[nnodes];          // nodal variable

// store coordinates in correct format for silo

  for(long i=0;i<nnodes;i++){xcoords[i]=M.Coord(0,i);}
  for(long i=0;i<nnodes;i++){ycoords[i]=M.Coord(1,i);}

// set up an array of length ndims pointing to the coordinate arrays

  coords[0]=xcoords;
  coords[1]=ycoords;

// connectivities

  for(int i=0,j=0;i<nzones;i++){
    for(int k=0;k<N.nloc();k++,j++){
      if(k==0){nodelist[j]=M.Vertex(i,k);}
      if(k==1){nodelist[j]=M.Vertex(i,k);}
      if(k==2){nodelist[j]=M.Vertex(i,3);} // we need to flip the top 2 nodes around as the
      if(k==3){nodelist[j]=M.Vertex(i,2);} // nodelist goes anticlockwise around the element
    }
  }

// zone shapes

  for(int i=0;i<nshapetypes;i++){shapesize[i]=N.nloc();}
  for(int i=0;i<nshapetypes;i++){shapecounts[i]=nzones;}

// material numbers

  for(int i=0;i<nmat;i++){matnos[i]=i+1;}
  for(int i=0;i<nzones;i++){matlist[i]=m.at(i);}
  for(int i=0;i<nmat;i++){matname[i]="Air";}

// write out connectivity information

  dberr=DBPutZonelist(dbfile,"zonelist1",nzones,ndims,nodelist,lnodelist,origin,shapesize,shapecounts,nshapetypes);

// write out the mesh

  optlist=DBMakeOptlist(2);
  dberr=DBAddOption(optlist,DBOPT_DTIME,&time);
  dberr=DBAddOption(optlist,DBOPT_CYCLE,&step);
  dberr=DBPutUcdmesh(dbfile,"Generator",ndims,NULL,coords,nnodes,nzones,"zonelist1",NULL,DB_DOUBLE,optlist);

// write out the material

  optlist=DBMakeOptlist(1);
  dberr=DBAddOption(optlist,DBOPT_MATNAMES,matname);
  dberr=DBPutMaterial(dbfile,"Materials","Elements",nmat,matnos,matlist,matdims,ndims,NULL,NULL,NULL,NULL,mixlen,DB_INT,optlist);
  dberr=DBFreeOptlist(optlist);

// write the element numbers

//  for(int i=0;i<nzones;i++){elnos[i]=i;}
//  optlist=DBMakeOptlist(1);
//  dberr=DBAddOption(optlist,DBOPT_UNITS,(void*)"");
//  dberr=DBPutUcdvar1(dbfile,"elementID","Elements",elnos,nzones,NULL,0,DB_INT,DB_ZONECENT,optlist);
//  dberr=DBFreeOptlist(optlist);

// release the arrays ready for the next mesh

  delete[] nodelist;
  nodelist=NULL;

  delete[] xcoords;
  xcoords=NULL;

  delete[] ycoords;
  ycoords=NULL;

  delete[] matlist;
  matlist=NULL;

  delete[] elnos;
  elnos=NULL;

  delete[] var1;
  var1=NULL;

// draw a mesh called "Kinematics" to support the kinematic dataset
// to plot this we need to subdivide through the nodes as visit can
// not cope with nodes inside the element

  nsubs=S.order()*S.order();             // number of subcells
  nzones=M.NCells()*nsubs;               // number of cells
  nnodes=M.NNodes_CFEM();                // number of nodes
  ndims=M.NDims();                       // number of dimensions
  lnodelist=4*nzones;                    // length of the nodelist
  nodelist=new int[lnodelist];           // allocate nodelist array
  xcoords=new double[nnodes];            // mesh vertex coordinates
  ycoords=new double[nnodes];            // mesh vertex coordinates
  matlist=new int[nzones];               // material numbers
  elnos=new int[nzones];                 // element numbers
  nodnos=new long[nnodes];               // node numbers
  var1=new double[nnodes];               // nodal variable

// store coordinates in correct format for silo

  for(long i=0;i<nnodes;i++){xcoords[i]=x.at(0).at(i);}
  for(long i=0;i<nnodes;i++){ycoords[i]=x.at(1).at(i);}

// set up an array of length ndims pointing to the coordinate arrays

  coords[0]=xcoords;
  coords[1]=ycoords;

// connectivities divide the S.nloc() node element up into sub-zones

  for(int i=0,j=0;i<M.NCells();i++){
    for(int jloc=0;jloc<S.order();jloc++){
      for(int iloc=0;iloc<S.order();iloc++,j+=4){

// node numbers in the corners of the sub-zone

        int iloc0(jloc*(S.order()+1)+iloc);
        int iloc1(jloc*(S.order()+1)+iloc+1);
        int iloc2((jloc+1)*(S.order()+1)+iloc);
        int iloc3((jloc+1)*(S.order()+1)+iloc+1);

// store in the node list

        nodelist[j]=M.GlobalNode_CFEM(i,iloc0);
        nodelist[j+1]=M.GlobalNode_CFEM(i,iloc1);
        nodelist[j+2]=M.GlobalNode_CFEM(i,iloc3); // we need to flip the top 2 nodes around as the
        nodelist[j+3]=M.GlobalNode_CFEM(i,iloc2); // nodelist goes anticlockwise around the element

      }
    }

  }

// zone shapes

  for(int i=0;i<nshapetypes;i++){shapesize[i]=4;}
  for(int i=0;i<nshapetypes;i++){shapecounts[i]=nzones;}

// write out connectivity information

  dberr=DBPutZonelist(dbfile,"zonelist2",nzones,ndims,nodelist,lnodelist,origin,shapesize,shapecounts,nshapetypes);

// write out the mesh

  optlist=DBMakeOptlist(2);
  dberr=DBAddOption(optlist,DBOPT_DTIME,&time);
  dberr=DBAddOption(optlist,DBOPT_CYCLE,&step);
  dberr=DBPutUcdmesh(dbfile,"Kinematics",ndims,NULL,coords,nnodes,nzones,"zonelist2",NULL,DB_DOUBLE,optlist);
  dberr=DBPutPointmesh(dbfile,"Kinematic_Nodes",ndims,coords,nnodes,DB_DOUBLE,optlist);

// write the node numbers

  for(long i=0;i<nnodes;i++){nodnos[i]=i;}
  optlist = DBMakeOptlist(1);
  dberr=DBAddOption(optlist, DBOPT_UNITS, (void*)"");
  dberr=DBPutUcdvar1(dbfile,"nodeID","Kinematics",nodnos,nnodes,NULL,0,DB_LONG,DB_NODECENT,optlist);
  dberr=DBFreeOptlist(optlist);

// write the element numbers

  for(int i=0,k=0;i<M.NCells();i++){for(int iloc=0;iloc<nsubs;iloc++,k++){elnos[k]=i;}}
  optlist=DBMakeOptlist(1);
  dberr=DBAddOption(optlist,DBOPT_UNITS,(void*)"");
  dberr=DBPutUcdvar1(dbfile,"elementID","Kinematics",elnos,nzones,NULL,0,DB_INT,DB_ZONECENT,optlist);
  dberr=DBFreeOptlist(optlist);

// velocity field x-component

  for(long i=0;i<nnodes;i++){var1[i]=u.at(0).at(i);}
  optlist = DBMakeOptlist(1);
  dberr=DBAddOption(optlist, DBOPT_UNITS, (void*)"cm/s");
  dberr=DBPutUcdvar1(dbfile,"velocity_x","Kinematics",var1,nnodes,NULL,0,DB_DOUBLE,DB_NODECENT,optlist); 
  dberr=DBFreeOptlist(optlist);

// velocity field y-component

  for(long i=0;i<nnodes;i++){var1[i]=u.at(1).at(i);}
  optlist = DBMakeOptlist(1);
  dberr=DBAddOption(optlist, DBOPT_UNITS, (void*)"cm/s");
  dberr=DBPutUcdvar1(dbfile,"velocity_y","Kinematics",var1,nnodes,NULL,0,DB_DOUBLE,DB_NODECENT,optlist); 
  dberr=DBFreeOptlist(optlist);

// velocity field resultant

  for(long i=0;i<nnodes;i++){double uu(u.at(0).at(i)),vv(u.at(1).at(i));var1[i]=sqrt(uu*uu+vv*vv);}
  optlist = DBMakeOptlist(1);
  dberr=DBAddOption(optlist, DBOPT_UNITS, (void*)"cm/s");
  dberr=DBPutUcdvar1(dbfile,"velocity","Kinematics",var1,nnodes,NULL,0,DB_DOUBLE,DB_NODECENT,optlist); 
  dberr=DBFreeOptlist(optlist);

// release the arrays ready for the next mesh

  delete[] nodelist;
  nodelist=NULL;

  delete[] xcoords;
  xcoords=NULL;

  delete[] ycoords;
  ycoords=NULL;

  delete[] matlist;
  matlist=NULL;

  delete[] elnos;
  elnos=NULL;

  delete[] nodnos;
  nodnos=NULL;

  delete[] var1;
  var1=NULL;

// draw a mesh called "Thermodynamics" to support the thermodynamic dataset
// to plot this we need to subdivide through the nodes as visit can
// not cope with nodes inside the element

  nsubs=T.order()*T.order();             // number of subcells
  nzones=M.NCells()*nsubs;               // number of cells
  nnodes=M.NNodes_DFEM();                // number of nodes
  ndims=M.NDims();                       // number of dimensions
  lnodelist=4*nzones;                    // length of the nodelist
  nodelist=new int[lnodelist];           // allocate nodelist array
  xcoords=new double[nnodes];            // mesh vertex coordinates
  ycoords=new double[nnodes];            // mesh vertex coordinates
  matlist=new int[nzones];               // material numbers
  elnos=new int[nzones];                 // element numbers
  nodnos=new long[nnodes];               // node numbers
  var1=new double[nnodes];               // nodal variable

// store coordinates in correct format for silo

  for(long i=0;i<nnodes;i++){xcoords[i]=xt.at(0).at(i);}
  for(long i=0;i<nnodes;i++){ycoords[i]=xt.at(1).at(i);}

// set up an array of length ndims pointing to the coordinate arrays

  coords[0]=xcoords;
  coords[1]=ycoords;

// connectivities divide the S.nloc() node element up into sub-zones

  for(int i=0,j=0;i<M.NCells();i++){
    for(int jloc=0;jloc<T.order();jloc++){
      for(int iloc=0;iloc<T.order();iloc++,j+=4){

// node numbers in the corners of the sub-zone

        int iloc0(jloc*(T.order()+1)+iloc);
        int iloc1(jloc*(T.order()+1)+iloc+1);
        int iloc2((jloc+1)*(T.order()+1)+iloc);
        int iloc3((jloc+1)*(T.order()+1)+iloc+1);

// store in the node list

        nodelist[j]=M.GlobalNode_DFEM(i,iloc0);
        nodelist[j+1]=M.GlobalNode_DFEM(i,iloc1);
        nodelist[j+2]=M.GlobalNode_DFEM(i,iloc3); // we need to flip the top 2 nodes around as the
        nodelist[j+3]=M.GlobalNode_DFEM(i,iloc2); // nodelist goes anticlockwise around the element

      }
    }

  }

// zone shapes

  for(int i=0;i<nshapetypes;i++){shapesize[i]=4;}
  for(int i=0;i<nshapetypes;i++){shapecounts[i]=nzones;}

// write out connectivity information

  dberr=DBPutZonelist(dbfile,"zonelist3",nzones,ndims,nodelist,lnodelist,origin,shapesize,shapecounts,nshapetypes);

// write out the mesh

  optlist=DBMakeOptlist(2);
  dberr=DBAddOption(optlist,DBOPT_DTIME,&time);
  dberr=DBAddOption(optlist,DBOPT_CYCLE,&step);
  dberr=DBPutUcdmesh(dbfile,"Thermodynamics",ndims,NULL,coords,nnodes,nzones,"zonelist3",NULL,DB_DOUBLE,optlist);
  dberr=DBPutPointmesh(dbfile,"Thermodynamic_Nodes",ndims,coords,nnodes,DB_DOUBLE,optlist);








  delete[] nodelist;
  nodelist=NULL;

  delete[] xcoords;
  xcoords=NULL;

  delete[] ycoords;
  ycoords=NULL;

  delete[] matlist;
  matlist=NULL;

  delete[] elnos;
  elnos=NULL;

  delete[] nodnos;
  nodnos=NULL;

  delete[] var1;
  var1=NULL;

// draw a mesh called "Sample" to support sampled data

  nsubs=25;                                           // number of sub-cells
  int nsubx=sqrt(nsubs);                              // number of sub-cells in x direction
  int nsuby=sqrt(nsubs);                              // number of sub-cells in y direction
  nzones=M.NCells()*nsubs;                            // number of cells
  nnodes=NSampleNodes(M,nsubs,SampleNode);            // number of nodes
  ndims=M.NDims();                                    // number of dimensions
  lnodelist=4*nzones;                                 // length of the nodelist
  nodelist=new int[lnodelist];                        // allocate nodelist array
  xcoords=new double[nnodes];                         // mesh vertex coordinates
  ycoords=new double[nnodes];                         // mesh vertex coordinates
  matlist=new int[nzones];                            // material numbers
  elnos=new int[nzones];                              // element numbers
  nodnos=new long[nnodes];                            // node numbers
  var1=new double[nnodes];                            // nodal variable
  vout=new double[nzones];                            // sampled data

// store coordinates in correct format for silo

  for(int i=0;i<M.NCells();i++){

// subdivide the element

    for(int isuby=0,isub=0;isuby<nsuby+1;isuby++){
      for(int isubx=0;isubx<nsubx+1;isubx++,isub++){

        double dx(2.0/nsubx),dy(2.0/nsuby);

// aquire the parametric coordinates isubx,isuby of local node isub

        double xloc(-1.0+isubx*dx);
        double yloc(-1.0+isuby*dy);

// aquire the global vertices of sub-cell isub from a finite element method

        double xsum(0.0),ysum(0.0);
        for(int iloc=0;iloc<T.nloc();iloc++){
          xsum+=T.value(iloc,xloc,yloc)*xt.at(0).at(M.GlobalNode_DFEM(i,iloc));
          ysum+=T.value(iloc,xloc,yloc)*xt.at(1).at(M.GlobalNode_DFEM(i,iloc));
        }

// store the coordinates of sample point isub

        xcoords[SampleNode.at(i).at(isub)]=xsum;
        ycoords[SampleNode.at(i).at(isub)]=ysum;

      }
    }

  }

// set up an array of length ndims pointing to the coordinate arrays

  coords[0]=xcoords;
  coords[1]=ycoords;

// connectivities joining up the sub-cells

  for(int i=0,j=0;i<M.NCells();i++){

    for(int jloc=0;jloc<nsuby;jloc++){
      for(int iloc=0;iloc<nsubx;iloc++,j+=4){

// node numbers in the corners of the sub-zone

        int iloc0(jloc*(nsubx+1)+iloc);
        int iloc1(jloc*(nsubx+1)+iloc+1);
        int iloc2((jloc+1)*(nsubx+1)+iloc);
        int iloc3((jloc+1)*(nsubx+1)+iloc+1);

// store in the node list

        nodelist[j]=SampleNode.at(i).at(iloc0);
        nodelist[j+1]=SampleNode.at(i).at(iloc1);
        nodelist[j+2]=SampleNode.at(i).at(iloc3); // we need to flip the top 2 nodes around as the
        nodelist[j+3]=SampleNode.at(i).at(iloc2); // nodelist goes anticlockwise around the element

      }
    }

  }

// zone shapes

  for(int i=0;i<nshapetypes;i++){shapesize[i]=4;}
  for(int i=0;i<nshapetypes;i++){shapecounts[i]=nzones;}

// write out connectivity information

  dberr=DBPutZonelist(dbfile,"zonelist4",nzones,ndims,nodelist,lnodelist,origin,shapesize,shapecounts,nshapetypes);

// write out the mesh

  optlist=DBMakeOptlist(2);
  dberr=DBAddOption(optlist,DBOPT_DTIME,&time);
  dberr=DBAddOption(optlist,DBOPT_CYCLE,&step);
  dberr=DBPutUcdmesh(dbfile,"Sampling",ndims,NULL,coords,nnodes,nzones,"zonelist4",NULL,DB_DOUBLE,optlist);
  dberr=DBPutPointmesh(dbfile,"Sample_Points",ndims,coords,nnodes,DB_DOUBLE,optlist);

// write out the density

  vin=d;sample(M,nsubs,T,vin,vout);
  optlist = DBMakeOptlist(1);
  dberr=DBAddOption(optlist, DBOPT_UNITS, (void*)"g/cc");
  dberr=DBPutUcdvar1(dbfile,"density","Sampling",vout,nzones,NULL,0,DB_DOUBLE,DB_ZONECENT,optlist);
  dberr=DBFreeOptlist(optlist);














  delete[] nodelist;
  nodelist=NULL;

  delete[] xcoords;
  xcoords=NULL;

  delete[] ycoords;
  ycoords=NULL;

  delete[] matlist;
  matlist=NULL;

  delete[] elnos;
  elnos=NULL;

  delete[] nodnos;
  nodnos=NULL;

  delete[] var1;
  var1=NULL;

  delete[] vout;
  vout=NULL;

// close the database

  dberr=DBClose(dbfile);

  return;

}

// function to return the number of sample points and set up their global numbers

long NSampleNodes(Mesh const &M, int const &nsubs, vector<vector<long> > &SampleNode){

// first empty the vector so we can re-use this function if needed

  vector<vector<long> > vempty;
  SampleNode=vempty;

// set number of sub-divisions along each axis of the cell

  int p=sqrt(nsubs);

// first find number of cells in each direction

  int ncellsx(M.NCells(0));
  int ncellsy(M.NCells(1));

// set number of nodes to sample

  long mNSampleNodes((p*ncellsx+1)*(p*ncellsy+1));

// set up global node numbers for a mesh

  long k(0.0);
  for(int j=0;j<ncellsy;j++){
    for(int i=0;i<ncellsx;i++){
      vector<long> global_nodes;
      for(int jloc=0;jloc<p+1;jloc++){
        for(int iloc=0;iloc<p+1;iloc++){
          long k(j*p*(p*ncellsx+1)+jloc*(p*ncellsx+1)+(p*i+iloc));
          global_nodes.push_back(k);
        }
      }
      SampleNode.push_back(global_nodes);
    }
  }

  return mNSampleNodes;

}

// function to interpolate data on to sample points

void sample(Mesh const &M,int const &nsubs,Shape const &T, vector<double> const &vin,double vout[]){

// number of subdivisions in each direction

  int nsubx(sqrt(nsubs));
  int nsuby(sqrt(nsubs));

  for(int i=0;i<M.NCells();i++){

// subdivide the element

    for(int isuby=0,isub=0;isuby<nsuby;isuby++){
      for(int isubx=0;isubx<nsubx;isubx++,isub++){

        double dx(2.0/nsubx),dy(2.0/nsuby);

// aquire the parametric coordinates isubx,isuby of local node isub

        double xloc(-1.0+0.5*dx+isubx*dx);
        double yloc(-1.0+0.5*dy+isuby*dy);

// aquire the global vertices of sub-cell isub from a finite element method

        double varsum(0.0);
        for(int iloc=0;iloc<T.nloc();iloc++){
          varsum+=T.value(iloc,xloc,yloc)*vin.at(M.GlobalNode_DFEM(i,iloc));
        }

// store the coordinates of sample point isub

        vout[nsubs*i+isub]=varsum;

      }
    }

  }

  return;

}

void silo(VVD const &x,VD const &d,VD const &p,VD const &e,VD const &q,VD const &c,
          VVD const &u,VI const &m,int step,double time,Mesh const &M,VD const &g){

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

  Shape S(1);                          // declare a degree-1 shape function
  int nzones(M.NCells());              // number of zones
  int ndims(M.NDims());                // number of dimensions
  int nnodesk(M.NNodes());             // number of nodes on kinematic grid
  int origin(0);                       // first address in nodelist arrays
  int lnodelist=(S.nloc()*nzones);     // length of the node list
  int nodelist[lnodelist];             // the node list
  int nshapetypes(1);                  // number of different shape types on the mesh
  int shapesize[nshapetypes];          // number of nodes defining each shape
  int shapecounts[nshapetypes];        // number of zones of each shape type

// set up material data structure

  int nmat(M.NMaterials());             // number of materials
  int matdims[]={nzones,1};             // material dimensions
  char*matname[nmat];                   // material names
  int matnos[nmat];                     // materials numbers present
  int matlist[nzones];                  // material number in each zone
  int mixlen(0);                        // number of mixed cells
  int mix_next[mixlen];                 // indices into mixed data arrays
  int mix_mat[mixlen];                  // material numbers for mixed zones
  int mix_vf[mixlen];                   // volunme fractions

// output arrays

  int elnos[nzones];                    // element numbers
  long nknos[nnodesk];                  // node numbers
  double var1[nnodesk];                 // zone centred scalars
  double var2[nzones];                  // zone centred scalars

  cout<<"       silo(): Writing a silo graphics dump to file "<<filename<<endl;

// store coordinates in correct format for silo and repeat for each mesh

  double xcoordsk[nnodesk];for(int i=0;i<M.NNodes();i++){xcoordsk[i]=x.at(0).at(i);}
  double ycoordsk[nnodesk];for(int i=0;i<M.NNodes();i++){ycoordsk[i]=x.at(1).at(i);}
  double *coordsk[]={xcoordsk,ycoordsk};

// connectivities

  for(int i=0,j=0;i<nzones;i++){
    for(int k=0;k<S.nloc();k++,j++){
//      nodelist[j]=M.Vertex(i,k);
      if(k==0){nodelist[j]=M.Vertex(i,k);}
      if(k==1){nodelist[j]=M.Vertex(i,k);}
      if(k==2){nodelist[j]=M.Vertex(i,3);} // we need to flip the top 2 nodes around as the
      if(k==3){nodelist[j]=M.Vertex(i,2);} // nodelist goes anticlockwise around the element
    }
  }

// zone shapes

  for(int i=0;i<nshapetypes;i++){shapesize[i]=S.nloc();}
  for(int i=0;i<nshapetypes;i++){shapecounts[i]=nzones;}

// material numbers

  for(int i=0;i<nmat;i++){matnos[i]=i+1;}
  for(int i=0;i<nzones;i++){matlist[i]=m.at(i);}
  for(int i=0;i<nmat;i++){matname[i]="Air";}

// disengage deprecation signalling

  dberr=DBSetDeprecateWarnings(0);

// create the silo database - this opens it also

  dbfile=DBCreate(filename,DB_CLOBBER,DB_LOCAL,fileinfo,DB_HDF5);

// write out connectivity information. //

  dberr=DBPutZonelist(dbfile, "zonelist", nzones, ndims, nodelist, lnodelist, origin ,shapesize, shapecounts, nshapetypes);

// write out the meshes //

  optlist=DBMakeOptlist(2);
  dberr=DBAddOption(optlist,DBOPT_DTIME,&time);
  dberr=DBAddOption(optlist,DBOPT_CYCLE,&step);
  dberr=DBPutUcdmesh(dbfile,"Elements",ndims,NULL,coordsk,nnodesk,nzones,"zonelist",NULL,DB_DOUBLE,optlist);
  dberr=DBPutPointmesh(dbfile,"Kinematics",ndims,coordsk,nnodesk,DB_DOUBLE,optlist);
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

  for(long i=0;i<nzones;i++){var2[i]=d.at(i);}
  optlist = DBMakeOptlist(1);
  dberr=DBAddOption(optlist, DBOPT_UNITS, (void*)"g/cc");
  dberr=DBPutUcdvar1(dbfile,"density","Elements",var2,nzones,NULL,0,DB_DOUBLE,DB_ZONECENT,optlist);
  dberr=DBFreeOptlist(optlist);

  for(long i=0;i<nzones;i++){var2[i]=p.at(i);}
  optlist = DBMakeOptlist(1);
  dberr=DBAddOption(optlist, DBOPT_UNITS, (void*)"Mb");
  dberr=DBPutUcdvar1(dbfile,"pressure","Elements",var2,nzones,NULL,0,DB_DOUBLE,DB_ZONECENT,optlist);
  dberr=DBFreeOptlist(optlist);

  for(long i=0;i<nzones;i++){var2[i]=q.at(i);}
  optlist = DBMakeOptlist(1);
  dberr=DBAddOption(optlist, DBOPT_UNITS, (void*)"Mb");
  dberr=DBPutUcdvar1(dbfile,"bulk_q","Elements",var2,nzones,NULL,0,DB_DOUBLE,DB_ZONECENT,optlist);
  dberr=DBFreeOptlist(optlist);

  for(long i=0;i<nzones;i++){var2[i]=c.at(i);}
  optlist = DBMakeOptlist(1);
  dberr=DBAddOption(optlist, DBOPT_UNITS, (void*)"cm/s");
  dberr=DBPutUcdvar1(dbfile,"sound_speed","Elements",var2,nzones,NULL,0,DB_DOUBLE,DB_ZONECENT,optlist);
  dberr=DBFreeOptlist(optlist);

  for(int i=0;i<nzones;i++){var2[i]=e.at(i);}
  optlist = DBMakeOptlist(1);
  dberr=DBAddOption(optlist, DBOPT_UNITS, (void*)"Mbcc");
  dberr=DBPutUcdvar1(dbfile,"energy","Elements",var2,nzones,NULL,0,DB_DOUBLE,DB_ZONECENT,optlist);
  dberr=DBFreeOptlist(optlist);

  for(int i=0;i<nzones;i++){var2[i]=g.at(matlist[i]-1);}
  optlist = DBMakeOptlist(1);
  dberr=DBAddOption(optlist, DBOPT_UNITS, (void*)"");
  dberr=DBPutUcdvar1(dbfile,"gamma","Elements",var2,nzones,NULL,0,DB_DOUBLE,DB_ZONECENT,optlist);
  dberr=DBFreeOptlist(optlist);

// write vector data fields

  for(long i=0;i<nnodesk;i++){var1[i]=u.at(0).at(i);}
  optlist = DBMakeOptlist(1);
  dberr=DBAddOption(optlist, DBOPT_UNITS, (void*)"cm/s");
  dberr=DBPutUcdvar1(dbfile,"velocity_x","Elements",var1,nnodesk,NULL,0,DB_DOUBLE,DB_NODECENT,optlist); 
  dberr=DBFreeOptlist(optlist);

  for(long i=0;i<nnodesk;i++){var1[i]=u.at(1).at(i);}
  optlist = DBMakeOptlist(1);
  dberr=DBAddOption(optlist, DBOPT_UNITS, (void*)"cm/s");
  dberr=DBPutUcdvar1(dbfile,"velocity_y","Elements",var1,nnodesk,NULL,0,DB_DOUBLE,DB_NODECENT,optlist); 
  dberr=DBFreeOptlist(optlist);

// write velocity resultant

  for(long i=0;i<nnodesk;i++){var1[i]=sqrt(u.at(0).at(i)*u.at(0).at(i)+u.at(1).at(i)*u.at(1).at(i));}
  optlist = DBMakeOptlist(1);
  dberr=DBAddOption(optlist, DBOPT_UNITS, (void*)"cm/s");
  dberr=DBPutUcdvar1(dbfile,"velocity","Elements",var1,nnodesk,NULL,0,DB_DOUBLE,DB_NODECENT,optlist); 
  dberr=DBFreeOptlist(optlist);

// close the database

  dberr=DBClose(dbfile);

  return;

}

// return today's date

//std::string date(){

//  time_t now = time(0);

//  return ctime(&now);

//}
