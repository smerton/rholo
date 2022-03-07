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
#include "eos.h"     // eos lookups
#include <cmath>

// function signatures

std::string date();

using namespace std;

void silo(VVD const &x,VVD const &xt,VVD const &xinit,VD const &d,VD const &l,VD const &V,VD const &e,
          VVD const &u,VI const &m,int step,double time,Mesh const &M,VD const &g,Shape const &S,Shape const &T){

// local function signatures

  long NSampleNodes(Mesh const &M, int const &nsubs, vector<vector<long> > &SampleNode); // obtain number of sample points and set up their global numbers
  void sample(Mesh const &M,int const &nsubs,Shape const &T, vector<double> const &vin,double vout[]); // interpolate vin to obtain vout at the sample points
  void J(Mesh const &M,Shape const &S,VVD const &x,int const &i,int const &nsubs,vector<double> &detj); //compute a jacobian and return the determinant at the sample points
  void J(Mesh const &M,Shape const &S,VVD const &x,int const &i,vector<double> &detj); //compute a jacobian and return the determinant at the integration points
  void J(Mesh const &M,Shape const &S,VVD const &x0,VVD const &x,int const &i,int const &nsubs,VD &detJs); //compute a jacobian for the lagrangian motion and return the determinant at the sample points

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

  int  nmat(M.NMaterials());             // number of materials
  int  matdims[2];                       // material dimensions
  char *matname[nmat];                   // material names
  int  matnos[nmat];                     // materials numbers present
  int  *matlist;                         // material number in each zone
  int  mixlen(0);                        // number of mixed cells
  int  mix_next[mixlen];                 // indices into mixed data arrays
  int  mix_mat[mixlen];                  // material numbers for mixed zones
  int  mix_vf[mixlen];                   // volunme fractions

// reuse the following, several different meshes are required to support the different lengths of data
// some of these meshes subdivide the elements and so the following data can vary for the different meshes

  int nzones;                       // number of cells
  int nsubs;                        // number of sub-cells
  int nsubx;                        // number of sub-cells on x-axis
  int nsuby;                        // number of sub-cells on y-axis
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
  double *var1;                     // nodal/zonal variable to output against the mesh
  vector<vector<long> > SampleNode; // global numbers of the sample points
  vector<double> vin;               // vector to sample at the mesh sample points
  vector<double> detj0,detj,detjs;  // determinants at the sample points for writing out d and l as functions

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
  nsubs=1;                          // number of sub-cells
  nsubx=sqrt(nsubs);                // number of sub-cells on x-axis
  nsuby=sqrt(nsubs);                // number of sub-cells on y-axis
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

  for(int j=0,iel=0;j<M.NCells(1);j++){
    for(int i=0;i<M.NCells(0);i++,iel++){

// corner node numbers in the corners of the order p element

      int iloc0(0),iloc1(S.order()),iloc2(S.nloc()-1-S.order()),iloc3(S.nloc()-1);

// corner node numbers in the corners of the p1 element

      int jnod0(j*(M.NCells(0)+1)+i),jnod1(jnod0+1),jnod2(jnod0+(M.NCells(0)+1)),jnod3(jnod1+(M.NCells(0)+1));

// copy corner coordinates from kinematic mesh on to a p1 mesh

      xcoords[jnod0]=x.at(0).at(M.GlobalNode_CFEM(iel,iloc0));
      xcoords[jnod1]=x.at(0).at(M.GlobalNode_CFEM(iel,iloc1));
      xcoords[jnod2]=x.at(0).at(M.GlobalNode_CFEM(iel,iloc2));
      xcoords[jnod3]=x.at(0).at(M.GlobalNode_CFEM(iel,iloc3));

      ycoords[jnod0]=x.at(1).at(M.GlobalNode_CFEM(iel,iloc0));
      ycoords[jnod1]=x.at(1).at(M.GlobalNode_CFEM(iel,iloc1));
      ycoords[jnod2]=x.at(1).at(M.GlobalNode_CFEM(iel,iloc2));
      ycoords[jnod3]=x.at(1).at(M.GlobalNode_CFEM(iel,iloc3));

    }
  }

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
  for(int i=0,k=0;i<M.NCells();i++){for(int j=0;j<nsubs;j++,k++){matlist[k]=m.at(i);}}
  for(int i=0;i<nmat;i++){matname[i]="Air";}
  matdims[0]=nzones;
  matdims[1]=1;

// write out connectivity information

  dberr=DBPutZonelist(dbfile,"zonelist1",nzones,ndims,nodelist,lnodelist,origin,shapesize,shapecounts,nshapetypes);

// write out the mesh

  optlist=DBMakeOptlist(2);
  dberr=DBAddOption(optlist,DBOPT_DTIME,&time);
  dberr=DBAddOption(optlist,DBOPT_CYCLE,&step);
  dberr=DBPutUcdmesh(dbfile,"Generator",ndims,NULL,coords,nnodes,nzones,"zonelist1",NULL,DB_DOUBLE,optlist);
  dberr=DBPutPointmesh(dbfile,"Generator_Nodes",ndims,coords,nnodes,DB_DOUBLE,optlist);

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

  nsubs=S.order()*S.order();             // number of sub-cells
  nsubx=sqrt(nsubs);                     // number of sub-cells on x-axis
  nsuby=sqrt(nsubs);                     // number of sub-cells on y-axis
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

// material numbers

  for(int i=0;i<nmat;i++){matnos[i]=i+1;}
  for(int i=0,k=0;i<M.NCells();i++){for(int j=0;j<nsubs;j++,k++){matlist[k]=m.at(i);}}
  for(int i=0;i<nmat;i++){matname[i]="Air";}
  matdims[0]=nzones;
  matdims[1]=1;

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

// write out the materials

  optlist=DBMakeOptlist(1);
  dberr=DBAddOption(optlist,DBOPT_MATNAMES,matname);
  dberr=DBPutMaterial(dbfile,"Materials","Kinematics",nmat,matnos,matlist,matdims,ndims,NULL,NULL,NULL,NULL,mixlen,DB_INT,optlist);
  dberr=DBFreeOptlist(optlist);

// write the material IDs

  for(int i=0,k=0;i<M.NCells();i++){for(int iloc=0;iloc<nsubs;iloc++,k++){elnos[k]=m.at(i);}}
  optlist=DBMakeOptlist(1);
  dberr=DBAddOption(optlist,DBOPT_UNITS,(void*)"");
  dberr=DBPutUcdvar1(dbfile,"materialID","Kinematics",elnos,nzones,NULL,0,DB_INT,DB_ZONECENT,optlist);
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

  nsubs=T.order()*T.order();             // number of sub-cells
  nsubx=sqrt(nsubs);                     // number of sub-cells on x-axis
  nsuby=sqrt(nsubs);                     // number of sub-cells on y-axis
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

// material numbers

  for(int i=0;i<nmat;i++){matnos[i]=i+1;}
  for(int i=0,k=0;i<M.NCells();i++){for(int j=0;j<nsubs;j++,k++){matlist[k]=m.at(i);}}
  for(int i=0;i<nmat;i++){matname[i]="Air";}
  matdims[0]=nzones;
  matdims[1]=1;

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

  nsubs=1;                                            // number of sub-cells
  nsubx=sqrt(nsubs);                                  // number of sub-cells on x-axis
  nsuby=sqrt(nsubs);                                  // number of sub-cells on y-axis
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
  var1=new double[nzones];                            // nodal/zonal variable

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
        for(int iloc=0;iloc<S.nloc();iloc++){
          xsum+=S.value(iloc,xloc,yloc)*x.at(0).at(M.GlobalNode_CFEM(i,iloc));
          ysum+=S.value(iloc,xloc,yloc)*x.at(1).at(M.GlobalNode_CFEM(i,iloc));
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

// material numbers

  for(int i=0;i<nmat;i++){matnos[i]=i+1;}
  for(int i=0,k=0;i<M.NCells();i++){for(int j=0;j<nsubs;j++,k++){matlist[k]=m.at(i);}}
  for(int i=0;i<nmat;i++){matname[i]="Air";}
  matdims[0]=nzones;
  matdims[1]=1;

// write out connectivity information

  dberr=DBPutZonelist(dbfile,"zonelist4",nzones,ndims,nodelist,lnodelist,origin,shapesize,shapecounts,nshapetypes);

// write out the mesh

  optlist=DBMakeOptlist(2);
  dberr=DBAddOption(optlist,DBOPT_DTIME,&time);
  dberr=DBAddOption(optlist,DBOPT_CYCLE,&step);
  dberr=DBPutUcdmesh(dbfile,"Sampling",ndims,NULL,coords,nnodes,nzones,"zonelist4",NULL,DB_DOUBLE,optlist);
  dberr=DBPutPointmesh(dbfile,"Sample_Points",ndims,coords,nnodes,DB_DOUBLE,optlist);

// internal energy

  vin=e;sample(M,nsubs,T,vin,var1);
  optlist = DBMakeOptlist(1);
  dberr=DBAddOption(optlist, DBOPT_UNITS, (void*)"Mbcc");
  dberr=DBPutUcdvar1(dbfile,"energy","Sampling",var1,nzones,NULL,0,DB_DOUBLE,DB_ZONECENT,optlist);
  dberr=DBFreeOptlist(optlist);

// density

  for(int i=0,k=0;i<M.NCells();i++){J(M,S,xinit,i,nsubs,detj0);J(M,S,x,i,nsubs,detj);for(int isub=0;isub<nsubs;isub++,k++){var1[k]=d.at(i)*detj0.at(isub)/detj.at(isub);}}
  optlist = DBMakeOptlist(1);
  dberr=DBAddOption(optlist, DBOPT_UNITS, (void*)"g/cc");
  dberr=DBPutUcdvar1(dbfile,"density","Sampling",var1,nzones,NULL,0,DB_DOUBLE,DB_ZONECENT,optlist);
  dberr=DBFreeOptlist(optlist);

// pressure

  vin=e;sample(M,nsubs,T,vin,var1);
  for(int i=0,k=0;i<M.NCells();i++){J(M,S,xinit,i,nsubs,detj0);J(M,S,x,i,nsubs,detj);for(int isub=0;isub<nsubs;isub++,k++){double dk(d.at(i)*detj0.at(isub)/detj.at(isub)),ek(var1[k]),gk(g.at(m.at(i)-1)),pk(P(dk,ek,gk));var1[k]=pk;}}
  optlist = DBMakeOptlist(1);
  dberr=DBAddOption(optlist, DBOPT_UNITS, (void*)"Mb");
  dberr=DBPutUcdvar1(dbfile,"pressure","Sampling",var1,nzones,NULL,0,DB_DOUBLE,DB_ZONECENT,optlist);
  dberr=DBFreeOptlist(optlist);

// gamma

  for(int i=0,k=0;i<M.NCells();i++){for(int isub=0;isub<nsubs;isub++,k++){var1[k]=g.at(m.at(i)-1);}}
  optlist = DBMakeOptlist(1);
  dberr=DBAddOption(optlist, DBOPT_UNITS, (void*)"");
  dberr=DBPutUcdvar1(dbfile,"gamma","Sampling",var1,nzones,NULL,0,DB_DOUBLE,DB_ZONECENT,optlist);
  dberr=DBFreeOptlist(optlist);

// length

  for(int i=0,k=0;i<M.NCells();i++){J(M,S,xinit,x,i,nsubs,detjs);for(int isub=0;isub<nsubs;isub++,k++){var1[k]=l.at(i)*detjs.at(isub);}}
  optlist = DBMakeOptlist(1);
  dberr=DBAddOption(optlist, DBOPT_UNITS, (void*)"cm");
  dberr=DBPutUcdvar1(dbfile,"length","Sampling",var1,nzones,NULL,0,DB_DOUBLE,DB_ZONECENT,optlist);
  dberr=DBFreeOptlist(optlist);

// volume

  for(int i=0,k=0;i<M.NCells();i++){J(M,S,xinit,x,i,nsubs,detjs);for(int isub=0;isub<nsubs;isub++,k++){var1[k]=sqrt(V.at(i))/S.order();}}
  optlist = DBMakeOptlist(1);
  dberr=DBAddOption(optlist, DBOPT_UNITS, (void*)"cc");
  dberr=DBPutUcdvar1(dbfile,"volume","Sampling",var1,nzones,NULL,0,DB_DOUBLE,DB_ZONECENT,optlist);
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

// store the data at the sample point isub

        vout[nsubs*i+isub]=varsum;

      }
    }

  }

  return;

}

// compute a jacobian and return a vector containing the determinant at each sample point

void J(Mesh const &M,Shape const &N,VVD const &x, int const &i,int const &nsubs,vector<double> &detj){

// clear out the vector passed in so it is safe to push

  vector<double> vempty;
  detj=vempty;

// subdivide the element i

  int nsubx=sqrt(nsubs),nsuby=nsubx;

  for(int isuby=0,isub=0;isuby<nsuby;isuby++){
    for(int isubx=0;isubx<nsubx;isubx++,isub++){

      double dx(2.0/nsubx),dy(2.0/nsuby);

// aquire the parametric coordinates isubx,isuby of local node isub

      double xloc(-1.0+0.5*dx+isubx*dx);
      double yloc(-1.0+0.5*dy+isuby*dy);

// aquire the global vertices of sub-cell isub from a finite element method

      double dxdu(0.0),dxdv(0.0),dydu(0.0),dydv(0.0);
      for(int iloc=0;iloc<N.nloc();iloc++){
        dxdu+=N.dvalue(0,iloc,xloc,yloc)*x.at(0).at(M.GlobalNode_CFEM(i,iloc));
        dxdv+=N.dvalue(1,iloc,xloc,yloc)*x.at(0).at(M.GlobalNode_CFEM(i,iloc));
        dydu+=N.dvalue(0,iloc,xloc,yloc)*x.at(1).at(M.GlobalNode_CFEM(i,iloc));
        dydv+=N.dvalue(1,iloc,xloc,yloc)*x.at(1).at(M.GlobalNode_CFEM(i,iloc));
      }

// determinant at the sample point

      detj.push_back(dxdu*dydv-dxdv*dydu);

    }
  }

  return;

}

// compute a jacobian and return a vector containing the determinant at each integration point

void J(Mesh const &M,Shape const &N,VVD const &x, int const &i,vector<double> &detj){

// clear out the vector passed in so it is safe to push

  vector<double> vempty;
  detj=vempty;

// compute jacobian at each integration point

  for(int gi=0;gi<N.ngi();gi++){

    double dxdu(0.0),dxdv(0.0),dydu(0.0),dydv(0.0);
    for(int iloc=0;iloc<N.nloc();iloc++){
      dxdu+=N.dvalue(0,iloc,gi)*x.at(0).at(M.GlobalNode_CFEM(i,iloc));
      dxdv+=N.dvalue(1,iloc,gi)*x.at(0).at(M.GlobalNode_CFEM(i,iloc));
      dydu+=N.dvalue(0,iloc,gi)*x.at(1).at(M.GlobalNode_CFEM(i,iloc));
      dydv+=N.dvalue(1,iloc,gi)*x.at(1).at(M.GlobalNode_CFEM(i,iloc));
    }

// determinant at the sample point

    detj.push_back(dxdu*dydv-dxdv*dydu);

  }

  return;

}

// calculate a jacobian for the Lagrangian motion and return the determinant

void J(Mesh const &M,Shape const &S,VVD const &x0,VVD const &x,int const &i,int const &nsubs,VD &detJs){

// clear out the vector passed in so it is safe to push

  vector<double> vempty;
  detJs=vempty;

// subdivide the element i

  int nsubx=sqrt(nsubs),nsuby=nsubx;

  for(int isuby=0,isub=0;isuby<nsuby;isuby++){
    for(int isubx=0;isubx<nsubx;isubx++,isub++){

      double dx(2.0/nsubx),dy(2.0/nsuby);

// aquire the parametric coordinates isubx,isuby of local node isub

      double xloc(-1.0+0.5*dx+isubx*dx);
      double yloc(-1.0+0.5*dy+isuby*dy);

// compute the deriavtives

      double dx0du(0.0),dy0du(0.0),dx0dv(0.0),dy0dv(0.0);
      double dxdu(0.0),dydu(0.0),dxdv(0.0),dydv(0.0);

// derivatives of the physical coordinates at the quadrature points

      for(int iloc=0;iloc<S.nloc();iloc++){

        long gloc=(S.type()==CONTINUOUS)?M.GlobalNode_CFEM(i,iloc):M.GlobalNode_DFEM(i,iloc);

// at time 0

        dx0du+=x0.at(0).at(gloc)*S.dvalue(0,iloc,xloc,yloc); // dx0/du
        dx0dv+=x0.at(0).at(gloc)*S.dvalue(1,iloc,xloc,yloc); // dx0/dv
        dy0du+=x0.at(1).at(gloc)*S.dvalue(0,iloc,xloc,yloc); // dy0/du
        dy0dv+=x0.at(1).at(gloc)*S.dvalue(1,iloc,xloc,yloc); // dy0/dv

// at time t

        dxdu+=x.at(0).at(gloc)*S.dvalue(0,iloc,xloc,yloc); // dx/du
        dxdv+=x.at(0).at(gloc)*S.dvalue(1,iloc,xloc,yloc); // dx/dv
        dydu+=x.at(1).at(gloc)*S.dvalue(0,iloc,xloc,yloc); // dy/du
        dydv+=x.at(1).at(gloc)*S.dvalue(1,iloc,xloc,yloc); // dy/dv

      }

// define the jacobian for time-0, this maps to the isoparametric element from the time-0 element

      vector<vector<double> > J0{{dx0du,dx0dv},
                                 {dy0du,dy0dv}};

// determinant of J0

      double detJ0(dx0du*dy0dv-dx0dv*dy0du);

// the inverse of J0 is simply the adjugate divided by the determinant

      vector<vector<double> > invJ0{{ dy0dv/detJ0,-dx0dv/detJ0},
                                    {-dy0du/detJ0, dx0du/detJ0}};

// define the jacobian for time-t, this maps to the isoparametric element from the time-t element

      vector<vector<double> > J{{dxdu,dxdv},
                                {dydu,dydv}};

// determinant of J

      double detJ(dxdu*dydv-dxdv*dydu);

// define a jacobian of the lagrangian motion, this maps to the time-0 element from the time-t element

      vector<vector<double> > Js{{invJ0[0][0]*J[0][0]+invJ0[0][1]*J[1][0],invJ0[0][0]*J[0][1]+invJ0[0][1]*J[1][1]},
                                 {invJ0[1][0]*J[0][0]+invJ0[1][1]*J[1][0],invJ0[1][0]*J[0][1]+invJ0[1][1]*J[1][1]}};

// determinant of Js

      detJs.push_back(Js[0][0]*Js[1][1]-Js[0][1]*Js[1][0]);

    }
  }

  return;

}

// return today's date

//std::string date(){

//  time_t now = time(0);

//  return ctime(&now);

//}
