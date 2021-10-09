// silo graphics IO function

// Author S. R. Merton

#define FILENAME "silo/graphics.silo"
#define CODENAME "rholo"
#define VERSION "1.0"
#define TITLE "MyFirstTitle"
#define FILEINFO "This is a silo file created by rholo and it contains polyhedral meshes." // define info for reader app
#define MESHNAME "MyFirstMesh" // name of a test mesh
#define VD vector<double>      // vector of doubles

#include <iostream>
#include <vector>
#include <ctime>
#include <filesystem>
#include "silo.h"
#include "shape.h"

// function signatures

std::string date();

using namespace std;

void silo(VD*x,VD*d,VD*p,VD*e,VD*u,int step,double time,Shape*K){

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

  int nx((x->size()-1)/K->order());   // number of cells in x direction
  int ny(1);                          // number of cells in y direction
  int nnodes = 2*x->size();           // number of nodes
  int nzones = nx*ny;                 // number of zones
  int ndims = 2;                      // number of dimensions
  int origin = 0;                     // first address in nodelist array
  int lnodelist=(2*K->nloc()*nzones); // length of the node list
  int nodelist[lnodelist];            // the node list
  int nshapetypes(1);                 // number of different shape types on the mesh
  int shapesize[nshapetypes];         // number of nodes defining each shape
  int shapecounts[nshapetypes];       // number of zones of each shape type

  cout<<"Writing a silo graphics dump to file "<<filename<<endl;

// store coordinates in correct format for silo

  double xcoords[nnodes];for(int i=0;i<x->size();i++){xcoords[i]=x->at(i);};for(int i=0;i<x->size();i++){xcoords[x->size()+i]=x->at(i);}
  double ycoords[nnodes];for(int i=0;i<x->size();i++){ycoords[i]=0.0; };for(int i=0;i<x->size();i++){ycoords[x->size()+i]=0.25; }
  double *coords[]={xcoords,ycoords};

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

// zone shapes

  for(int i=0;i<nshapetypes;i++){shapesize[i]=2*K->nloc();}
  for(int i=0;i<nshapetypes;i++){shapecounts[i]=nzones;}

// create the silo database - this opens it also

  dbfile=DBCreate(filename,DB_CLOBBER,DB_LOCAL,fileinfo,DB_HDF5);

// write out connectivity information. //

  dberr=DBPutZonelist(dbfile, "zonelist", nzones, ndims, nodelist, lnodelist, origin ,shapesize, shapecounts, nshapetypes);

// write out the meshes //

  dberr=DBPutUcdmesh(dbfile, "Kinematics", ndims, NULL, coords, nnodes, nzones,"zonelist", NULL, DB_DOUBLE, NULL);
//  dberr=DBPutUcdmesh(dbfile, "Thermodynamics", ndims, NULL, coords, nnodes, nzones,"zonelist", NULL, DB_DOUBLE, NULL);
//  dberr=DBPutPointmesh(dbfile, "Nodes", ndims, coords, nnodes,DB_DOUBLE, NULL);
//  dberr=DBPutPointmesh(dbfile, "Gauss", ndims, coords, nnodes,DB_DOUBLE, NULL);

// close the database

  dberr=DBClose(dbfile);

  return;

}

// return today's date

std::string date(){

  time_t now = time(0);

  return ctime(&now);

}
