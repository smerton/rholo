// Function definitions for members of the mesh class

// Author S. R. Merton

#include "mesh.h"
#include "shape.h"
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
          vectmp.push_back(v4);
          vectmp.push_back(v3);
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

// compute element volume

  Shape S(1); // load a shape function

// calculate a jacobian

  for(int i=0;i<NCells();i++){

    double ivol(0.0);

    for(int gi=0;gi<S.ngi();gi++){

      double dxdu(0.0),dydu(0.0),dxdv(0.0),dydv(0.0);

// derivatives of the physical coordinates

      for(int j=0;j<S.nloc();j++){
        dxdu+=Coord(0,Vertex(i,j))*S.dvalue(0,j,gi); // dx/du
        dydu+=Coord(1,Vertex(i,j))*S.dvalue(0,j,gi); // dy/du
        dxdv+=Coord(0,Vertex(i,j))*S.dvalue(1,j,gi); // dx/dv
        dydv+=Coord(1,Vertex(i,j))*S.dvalue(1,j,gi); // dy/dv
      }

// calculate the determinant

      double detJ(dxdu*dydv-dydu*dxdv);

// integrate the jacobian to obtain the cell volume

      ivol+=detJ*S.wgt(gi);

    }

// commit to volume address space in the mesh class

    mVolume.push_back(ivol);

  }

// connectivity traversal

  set_E2E();

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

// member function to set up the element->element connectivities

void Mesh::set_E2E(){

  for(int iel=0;iel<NCells();iel++){

// vertices of current element

    int n[5]={Vertex(iel,0),Vertex(iel,1),Vertex(iel,3),Vertex(iel,2),Vertex(iel,0)};

    vector<int> iel_neighbours={-1,-1,-1,-1};

// loop over faces of element iel

    for(int iface=0;iface<NVertices(iel);iface++){

// search for a neighbour on face iface of iel

      for(int ieln=0;ieln<NCells();ieln++){

        if(iel==ieln){continue;}

// vertices of element ieln

        int nn[5]={Vertex(ieln,0),Vertex(ieln,1),Vertex(ieln,3),Vertex(ieln,2),Vertex(ieln,0)};

// loop over faces of element ieln

        for(int ifacen=0;ifacen<NVertices(iel);ifacen++){

// check if faces iface and ifacen are coincident

          bool z1(n[iface]==nn[ifacen] || n[iface]==nn[ifacen+1]);
          bool z2(n[iface+1]==nn[ifacen] || n[iface+1]==nn[ifacen+1]);

          if(z1&&z2){iel_neighbours.at(iface)=ieln;} // ieln on face iface of iel

        }

      }

    }

    mE2E.push_back(iel_neighbours);

  }

// now connect boundary segments to the mesh

  for(int ib=0;ib<NSides();ib++){

// vertices of boundary segment

    int nb[2]={SideNode(ib,0),SideNode(ib,1)};

// search mesh for the adjacent element

    for(int iel=0;iel<NCells();iel++){

// vertices of element iel

      int n[5]={Vertex(iel,0),Vertex(iel,1),Vertex(iel,3),Vertex(iel,2),Vertex(iel,0)};

// search element for a face coincident with boundary segment

      for(int iface=0;iface<NVertices(iel);iface++){

// check if face iface lies on the domain boundary

        if(E2E(iel,iface)<0){

// check if face iface and boundary segment ib are coincident

          bool z1(n[iface]==nb[0] || n[iface]==nb[1]);
          bool z2(n[iface+1]==nb[0] || n[iface+1]==nb[1]);

// boundary segment ib on face iface of iel

          if(z1&&z2){
            mE2E.at(iel).at(iface)=NCells()+ib; // connect ib to face iface of iel
            vector<int> ib_neighbours={iel};// connect iel to boundary segment ib
            mE2E.push_back(ib_neighbours);
          }


        }

      }

    }

  }

  return;

}

// function to convert a string to an int

constexpr unsigned int str2int(const char* str, int h = 0){return !str[h] ? 5381 : (str2int(str, h+1) * 33) ^ str[h];}

// member function to return the number of mesh dimensions

int Mesh::NDims() const {return mNDims;}

// member function to return the number of nodes

int Mesh::NNodes() const {return mNNodes;}

// member function to return the number of materials

int Mesh::NMaterials() const {return mNMaterials;}

// member function to return the number of cells

int Mesh::NCells() const {return mNCells;}

// member function to return the material in cell i

int Mesh::Material(int i) const {return mMaterial.at(i);}

// member function to return the geometric type of cell i

int Mesh::Type(int i) const {return mType.at(i);}

// member function to return the number of vertices of cell i

int Mesh::NVertices(int i) const {return mVertex.at(i).size();}

// member function to return the node number of vertex j of cell i

int Mesh::Vertex(int i, int j) const {return mVertex[i][j];}

// member function to return the number of element sides coincident with a mesh edge

int Mesh::NSides() const {return mNSides;}

// member function to return the attribute of side i on edge of mesh

int Mesh::SideAttr(int i) const {return mSideAttr.at(i);}

// member function to return the type of side i on edge of mesh

int Mesh::SideType(int i) const {return mSideType.at(i);}

// member function to return the number of nodes on side i on edge of mesh

int Mesh::NSideNodes(int i) const {return mSideNode.at(i).size();}

// member function to return the node number of node j on side i on edge of mesh

int Mesh::SideNode(int i, int j) const {return mSideNode[i][j];}

// member function to return coordinate idim of node i

double Mesh::Coord(int idim, int i) const {return mCoord[idim][i];}

// member function to initialise a vector with the mesh coordinates

void Mesh::InitCoords(vector<vector<double> > &v){

// first resize the first stride

  v.resize(NDims());

// replace each dimension with the vector of coordinates

  for(int i=0;i<NDims();i++){
    v.at(i)=mCoord.at(i);
  }

  return;

}

// member function to return the element volume

double Mesh::Volume(int i) const {return mVolume.at(i);}

// member function to push a new boundary condition

void Mesh::bc_set(int iedge,int bc){mbc_edge.at(iedge)=bc;}

// member function to push a new boundary condition and assign a value to the mesh edge

void Mesh::bc_set(int iedge,int bc,double bcvalue){mbc_edge.at(iedge)=bc;mbc_value.at(iedge)=bcvalue;}

// member function to return the boundary condition on mesh edge iedge

int Mesh::bc_edge(int iedge) const {return mbc_edge[iedge];}

// member function to return the boundary value on mesh edge iedge

double Mesh::bc_value(int iedge) const {return mbc_value[iedge];}

// member function to return the number of boundary conditions

int Mesh::nbcs() const {return mbc_edge.size();}

// member function to return the element on face iface of element iel

int Mesh::E2E(int iel,int iface) const {return mE2E[iel][iface];}

// member function to return mesh boundary minimum in dimension idim

double Mesh::Min(int idim) const {return *min_element(mCoord.at(idim).begin(),mCoord.at(idim).end());}

// member function to return mesh boundary maximum in dimension idim

double Mesh::Max(int idim) const {return *max_element(mCoord.at(idim).begin(),mCoord.at(idim).end());}

// Destructor function to release storage associated with a mesh class object

Mesh::~Mesh(){}
