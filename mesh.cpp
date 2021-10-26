// Function definitions for members of the mesh class

// Author S. R. Merton

#include "mesh.h"
#include <vector>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <algorithm>

using namespace std;

// enumerator to map the different input blocks

enum Block{invalid,dimension,elements,boundary,vertices,nodes};
Block resolveBlock(string line);

// Constructor to instantiate a new mesh from the data in the file meshfile

Mesh::Mesh(char* meshfile){

  cout<<"  Mesh::Mesh(): Loading a mesh from file "<<meshfile<<endl;

  long lines(0);
  string line;

// open the mesh file

  ifstream meshdata(meshfile); // open the mesh file

// read every line into the string

  while(getline(meshdata,line)){

// read input in each block

    switch(resolveBlock(line)){
      case dimension:{

        cout<<"  Mesh::Mesh(): parsing dimension input block..."<<endl;

// read in number of dimensions

        getline(meshdata,line);
        stringstream ss(line);
        ss>>mNDims;

        break;
      }

      case elements:{

        cout<<"  Mesh::Mesh(): parsing elements input block..."<<endl;

// read in number of cells

        getline(meshdata,line);
        stringstream ss(line);
        ss>>mNCells;

// read in data for each cell

        for(int i=0;i<mNCells;i++){
          int mat,type,v1,v2,v3,v4;
          vector<int> vectmp;
          getline(meshdata,line);
          stringstream ss(line);
          ss>>mat>>type>>v1>>v2>>v3>>v4;
          mMaterial.push_back(mat);
          mType.push_back(type);
          vectmp.push_back(v1);
          vectmp.push_back(v2);
          vectmp.push_back(v3);
          vectmp.push_back(v4);
          mVertex.push_back(vectmp);
        }

        mNMaterials=*max_element(mMaterial.begin(),mMaterial.end());

        break;
      }

      case boundary:{

        cout<<"  Mesh::Mesh(): parsing boundary input block..."<<endl;

// read in number of cells sides on the mesh edge

        getline(meshdata,line);
        stringstream ss(line);
        ss>>mNSides;

// read in boundary data for each cell side coincident with the mesh edge

        for(int i=0;i<mNSides;i++){
          int battr,btype,n1,n2;
          vector<int> vectmp;
          getline(meshdata,line);
          stringstream ss(line);
          ss>>battr>>btype>>n1>>n2;
          mSideAttr.push_back(battr);
          mSideType.push_back(btype);
          vectmp.push_back(n1);
          vectmp.push_back(n2);
          mSideNode.push_back(vectmp);
        }

        break;
      }

      case vertices:{

        cout<<"  Mesh::Mesh(): parsing vertices input block..."<<endl;

// read in number of mesh vertices (== nodes in the mesh as we can only read in quads)

        getline(meshdata,line);
        stringstream ss(line);
        ss>>mNNodes;

// read in serialisation of coordinate list

        int ncolumns;
        {getline(meshdata,line);
        stringstream ss(line);
        ss>>ncolumns;}
        vector<double> vectmp[ncolumns];

// read in coordiates

        for(int i=0;i<mNNodes;i++){
          double r;
          getline(meshdata,line);
          stringstream ss(line);
          for(int j=0;j<ncolumns;j++){
            ss>>r;
            vectmp[j].push_back(r);
          }
        }

// unpack into mesh object

        for(int j=0;j<ncolumns;j++){mCoord.push_back(vectmp[j]);}

        break;
      }

      case nodes:{

        cout<<"  Mesh::Mesh(): parsing nodes input block..."<<endl;

        getline(meshdata,line); // read past "nodes"
        getline(meshdata,line); // read past "FiniteElementSpace"
        getline(meshdata,line); // read past "FiniteElementCollection"
        getline(meshdata,line); // read past "VDim"
        getline(meshdata,line); // read past "Ordering"

// read in node coordinates

        for(int idim=0;idim<mNDims;idim++){
          double r;
          vector<double> vectmp;
          for(int i=0;i<mNNodes;i++){
            getline(meshdata,line);
            stringstream ss(line);
            ss>>r;
            vectmp.push_back(r);
          }
          mCoord.push_back(vectmp);
        }

        break;
      }

      default:{ // handle for unmapped input block
      break;

      }

    }

    lines+=1;

  }

// end of file read

  cout<<"  Mesh::Mesh(): Done. Mesh loaded, "<<lines<<" lines read."<<endl;
  cout<<endl;

// print out what we have found

  cout<<"  Mesh::Mesh(): Data found in mesh file:"<<endl;

  cout<<"    Number of dimensions: "<<NDims()<<endl;
  cout<<"    Number of materials: "<<NMaterials()<<endl;
  cout<<"    Number of cells: "<<NCells()<<endl;
  cout<<"    Number of nodes: "<<NNodes()<<endl;
  cout<<"    Number of cell sides coincident with a mesh edge: "<<NSides()<<endl;
  cout<<endl;

  cout<<"  Mesh::Mesh(): Cell Data:"<<endl;
  cout<<"    i m t vertices"<<endl;
  for(int i=0;i<NCells();i++){
    cout<<"    "<<i<<" "<<Material(i)<<" "<<Type(i)<<" ";
    for(int j=0;j<NVertices(i);j++){
      cout<<Vertex(i,j)<<" ";
    }
    cout<<endl;
  }

  cout<<endl;

  cout<<"  Mesh::Mesh(): Boundary Data:"<<endl;
  cout<<"    s a t nodes"<<endl;
  for(int i=0;i<NSides();i++){
    cout<<"    "<<i<<" "<<SideAttr(i)<<" "<<SideType(i)<<" ";
    for(int j=0;j<NSideNodes(i);j++){
      cout<<SideNode(i,j)<<" ";
    }
    cout<<endl;
  }

  cout<<endl;

  cout<<"  Mesh::Mesh(): Vertices:"<<endl;
  cout<<"    n coords"<<endl;
  for(int j=0;j<NNodes();j++){
    cout<<"    "<<j<<" ";
    for(int idim=0;idim<NDims();idim++){
      cout<<Coord(idim,j)<<" ";
    }
    cout<<endl;
  }

  return;

}

// function to split input data into blocks

 Block resolveBlock(std::string input) {
    if( input == "dimension" ) return dimension;
    if( input == "elements" ) return elements;
    if( input == "boundary" ) return boundary;
    if( input == "vertices" ) return vertices;
    if( input == "nodes" ) return nodes;
    return invalid;
 }

// function to convert a string to an int

constexpr unsigned int str2int(const char* str, int h = 0){return !str[h] ? 5381 : (str2int(str, h+1) * 33) ^ str[h];}

// member function to return the number of mesh dimensions

int Mesh::NDims(){return mNDims;}

// member function to return the number of nodes

int Mesh::NNodes(){return mNNodes;}

// member function to return the number of materials

int Mesh::NMaterials(){return mNMaterials;}

// member function to return the number of cells

int Mesh::NCells(){return mNCells;}

// member function to return the material in cell i

int Mesh::Material(int i){return mMaterial.at(i);}

// member function to return the geometric type of cell i

int Mesh::Type(int i){return mType.at(i);}

// member function to return the number of vertices of cell i

int Mesh::NVertices(int i){return mVertex.at(i).size();}

// member function to return the node number of vertex j of cell i

int Mesh::Vertex(int i, int j){return mVertex[i][j];}

// member function to return the number of element sides coincident with a mesh edge

int Mesh::NSides(){return mNSides;}

// member function to return the attribute of side i on edge of mesh

int Mesh::SideAttr(int i){return mSideAttr.at(i);}

// member function to return the type of side i on edge of mesh

int Mesh::SideType(int i){return mSideType.at(i);}

// member function to return the number of nodes on side i on edge of mesh

int Mesh::NSideNodes(int i){return mSideNode.at(i).size();}

// member function to return the node number of node j on side i on edge of mesh

int Mesh::SideNode(int i, int j){return mSideNode[i][j];}

// member function to return coordinate idim of node i

double Mesh::Coord(int idim, int i){return mCoord[idim][i];}

// Destructor function to release storage associated with a mesh class object

Mesh::~Mesh(){}
