// TyphonIO graphics output function

// Author S. R. Merton

#include <iostream>
#include <vector>
#include "typhonio.h"

#define FILENAME "graphics.h5"
#define CODENAME "rholo"
#define VERSION "1.0"
#define DATE "06/10/21"
#define TITLE "My First One"
#define VD vector<double>       // vector of doubles

using namespace std;

void tio(VD*x0,VD*d0,VD*p,VD*e0,VD*u0,int step,double time){

  const char*filename(FILENAME);
  const char*codename(CODENAME);
  const char*version(VERSION);
  const char*date(DATE);
  const char*title(TITLE);

  TIO_File_t fileID;    // file handle
  TIO_Object_t stateID; // state handle
  TIO_Object_t meshID;  // mesh handle
  TIO_t err;            // error handle

// create the file

  if(step!=0){

// open the file

    cout<<"Opening graphics file FILENAME on step "<<step<<"..."<<endl;

    err=TIO_Open(FILENAME,&fileID,TIO_ACC_READWRITE,CODENAME,VERSION,DATE,TITLE,MPI_COMM_NULL,MPI_INFO_NULL,MPI_PROC_NULL);

    if(err!=TIO_SUCCESS){
      cout<<"tio(): Error opening graphics output file."<<endl;
      exit(1);
    }

  }else{

// create the file

    cout<<"Creating graphics file FILENAME..."<<endl;

    err=TIO_Create(FILENAME,&fileID,TIO_ACC_REPLACE,CODENAME,VERSION,DATE,TITLE,MPI_COMM_NULL,MPI_INFO_NULL,MPI_PROC_NULL);

    if(err!=TIO_SUCCESS){
      cout<<"tio(): Error creating graphics output file."<<endl;
      exit(1);
    }

  }

// create a new state for current time level

  string str1("step_"+to_string(step));
  const char*statename(str1.c_str());
  const char*units("s");

  err=TIO_Create_State(fileID,statename,&stateID,step,time,units);

  if(err!=TIO_SUCCESS){
    cout<<"tio(): Error creating new state."<<endl;
    exit(1);
  }

// output the mesh

  const char*meshname("Mesh");
  const char*group("Continuous");
  TIO_Dims_t ndims(TIO_2D); // no. dimensions
  int n1(4); // no. nodes
  int n2(1); // no. cells
  int n3(1); // no. shapes
  int n4(4); // no. connectivities
  int nchunks(1); // no. chunks
  const char*iunits("cm");
  const char*junits("cm");
  const char*kunits("cm");
  const char*ilabel("x");
  const char*jlabel("y");
  const char*klabel("z");

  err=TIO_Create_Mesh(fileID,stateID,meshname,&meshID,TIO_MESH_UNSTRUCT,TIO_COORD_CARTESIAN,TIO_FALSE,group, \
                      1,TIO_INT,TIO_DOUBLE,ndims,n1,n2,n3,n4,nchunks,iunits,junits,kunits,ilabel,jlabel,klabel);
  if(err!=TIO_SUCCESS){
    cout<<"tio(): Error creating mesh."<<endl;
    exit(1);
  }



// create the chunk

// debug - set up a dummy mesh
  int idx(0);
  int nghost_nodes(0);
  int nghost_cells(0);
  int nghost_shapes(0);
  int nghost_connectivity(0);
  int nmixcell(0);
  int nmixcomp(0);
  int nodeIDs[n1];for(int i=0;i<n1;i++){nodeIDs[i]=i;}
  int cellIDs[n2];for(int i=0;i<n2;i++){cellIDs[i]=i;}
  TIO_Shape_t shapes[n3];for(int i=0;i<n3;i++){shapes[i]=TIO_SHAPE_QUAD4;}
  int ncells[n3];for(int i=0;i<n3;i++){ncells[i]=n2;}
  int connectivity[]={0,1,2,3};
  double icoords[]={0.0,1.0,0.0,1.0};
  double jcoords[]={0.0,0.0,1.0,1.0};
  double kcoords[]={0.0,0.0,0.0,0.0};
// debug - set up a dummy mesh

  err=TIO_Set_Unstr_Chunk(fileID,meshID,idx,ndims,n1,n2,n3,n4,nghost_nodes,nghost_cells,nghost_shapes,nghost_connectivity,nmixcell,nmixcomp);
  if(err!=TIO_SUCCESS){
    cout<<"tio(): Error setting chunk."<<endl;
    exit(1);
  }

// write the chunk

  err=TIO_Write_UnstrMesh_Chunk(fileID,meshID,idx,TIO_XFER_INDEPENDENT,TIO_INT,TIO_DOUBLE,nodeIDs,cellIDs,\
                                shapes,ncells,connectivity,icoords,jcoords,kcoords);
  if(err!=TIO_SUCCESS){
    cout<<"tio(): Error writing chunk."<<endl;
    exit(1);
  }



// close the mesh

  err=TIO_Close_Mesh(fileID,meshID );

  if(err!=TIO_SUCCESS){
    cout<<"tio(): Error closing mesh."<<endl;
    exit(1);
  }


// close the state

  err=TIO_Close_State(fileID,stateID);

  if(err!=TIO_SUCCESS){
    cout<<"tio(): Error closing new state."<<endl;
    exit(1);
  }

// close the file

  err=TIO_Close(fileID);

  if(err!=TIO_SUCCESS){
    cout<<"tio(): Error closing graphics output file."<<endl;
    exit(1);
  }

  cout<<"Graphics file FILENAME closed."<<endl;

  return;

}
