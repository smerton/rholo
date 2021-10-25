// Function definitions for members of the mesh class

// Author S. R. Merton

#include "mesh.h"
#include <vector>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>

using namespace std;

enum Block{invalid,dimension,elements,boundary,vertices,nodes};
Block resolveBlock(string line);

// Constructor to instantiate a new mesh from the data in the file meshfile

Mesh::Mesh(char* meshfile){

  cout<<"  Mesh::Mesh(): Loading a mesh from file "<<meshfile<<endl;

  string line;

// open the mesh file

  ifstream meshdata(meshfile); // open the mesh file

// read every line into the string

  while(getline(meshdata,line)){

// read input in each block
//cout<<" block= "<<resolveBlock(line)<<endl;
    switch(resolveBlock(line)){
      case dimension:{
        cout<<"Mesh::Mesh(): parsing dimension input block..."<<endl;

// read in number of dimensions

        getline(meshdata,line);
        stringstream ss(line);
        ss>>mNDims;

        break;
      }

      case elements:{
        cout<<"Mesh::Mesh(): parsing elements input block..."<<endl;

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

        break;
      }

      case boundary:{
        cout<<"Mesh::Mesh(): parsing boundary input block..."<<endl;

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
        cout<<"Mesh::Mesh(): parsing vertices input block..."<<endl;

// read in number of mesh vertices (== nodes in the mesh as we can only read in quads)

        getline(meshdata,line);
        stringstream ss(line);
        ss>>mNNodes;

        break;
      }

      case nodes:{
        cout<<"Mesh::Mesh(): parsing nodes input block..."<<endl;

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

  }

// print out what we have found

  cout<<"mNDims= "<<mNDims<<endl;
  cout<<"mNCells= "<<mNCells<<endl;
  cout<<"mNNodes= "<<mNNodes<<endl;

  cout<<"Cell Data"<<endl;
  for(int i=0;i<mNCells;i++){
    cout<<i<<" "<<mMaterial.at(i)<<" "<<mType.at(i)<<" ";
    for(int ivrt=0;ivrt<mVertex.at(i).size();ivrt++){
      cout<<mVertex[i][ivrt]<<" ";
    }
    cout<<endl;
  }

  cout<<"Boundary Data"<<endl;
  for(int i=0;i<mNSides;i++){
    cout<<i<<" "<<mSideAttr.at(i)<<" "<<mSideType.at(i)<<" ";
    for(int j=0;j<mSideNode.at(i).size();j++){
      cout<<mSideNode[i][j]<<" ";
    }
    cout<<endl;
  }

  cout<<"Node Coordinates"<<endl;
  for(int j=0;j<mNNodes;j++){
    for(int idim=0;idim<mNDims;idim++){
      cout<<mCoord[idim][j]<<" ";
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

// Destructor function to release storage associated with a mesh class object

Mesh::~Mesh(){}
