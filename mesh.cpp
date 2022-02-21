// Function definitions for members of the mesh class

// Author S. R. Merton

#define VVD vector<vector<double> > // laziness
#define VTOL 1.0e-10                // threshold for volume errors
#define ECUT 1.0e-8                 // cut-off on the energy field

#include "mesh.h"
#include "shape.h"
#include "eos.h"     // eos lookups
#include "matrix.h"
#include <vector>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <algorithm>
#include <iomanip>
#include <cmath>     // sqrt

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

// initialise ghost size

  mNGCells=0;
  mNGNodes=0;

// opposing faces

  int oface[4]={2,3,0,1};

// number of cells on each mesh edge

  int cells_on_edge[4]={0,0,0,0};

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

// find the number of elements along each edge of the mesh

  for(int iel=0;iel<NCells();iel++){
    for(int iface=0;iface<NVertices(iel);iface++){
      if(E2E(iel,iface)<0){
        cells_on_edge[iface]++;
      }
    }
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
            vector<int> ib_neighbours={-1,-1,-1,-1}; // vector address space for neighbours
            mE2E.push_back(ib_neighbours); // append vector address space for neighbours
            mE2E.at(NCells()+ib).at(oface[iface])=iel; // connect iel to boundary segment ib
            mNGCells++; // create space for a ghost cell
          }

        }

      }

    }

  }

// now construct ghost cells around perimeter of the mesh

  int inod(NNodes()-1);

  for(int ielg=NCells();ielg<NCells()+NGCells();ielg++){
    for(int ifaceg=0;ifaceg<4;ifaceg++){
      int iel(E2E(ielg,ifaceg));
      if(iel!=-1){

        inod++;

// ghost cell vertex numbering

        vector<int> vectmp1;
        vector<double> vectmp2[2];

// generate vertices for the ghost cells by reflecting the physical element

        switch(ifaceg){
          case(0):
            mCoord.at(0).push_back(Coord(0,Vertex(iel,1)));
            mCoord.at(1).push_back(Coord(1,Vertex(iel,3))+(Coord(1,Vertex(iel,3))-Coord(1,Vertex(iel,1))));
            vectmp1.push_back(Vertex(iel,2));   // ghost node 0 coincident with a physical node 2
            vectmp1.push_back(Vertex(iel,3));   // ghost node 1 coincident with a physical node 3
            vectmp1.push_back(inod+1);   // equates to vertex number of ghost node 2
            vectmp1.push_back(inod); // equates to vertex number of ghost node 3
            mNGNodes++;

// last ghost on edge so skip corner node

            if(ielg==(NCells()+cells_on_edge[0]+cells_on_edge[1]+cells_on_edge[2]-1)){
              mCoord.at(0).push_back(Coord(0,Vertex(iel,0)));
              mCoord.at(1).push_back(Coord(1,Vertex(iel,2))+(Coord(1,Vertex(iel,2))-Coord(1,Vertex(iel,0))));
              mNGNodes++;
              inod++;
            }

          break;
          case(1):
            mCoord.at(0).push_back(Coord(0,Vertex(iel,2))-(Coord(0,Vertex(iel,3))-Coord(0,Vertex(iel,2))));
            mCoord.at(1).push_back(Coord(1,Vertex(iel,3)));
            vectmp1.push_back(inod+1);   // equates to vertex number of ghost node 0
            vectmp1.push_back(Vertex(iel,0));   // ghost node 1 coincident with a physical node 0
            vectmp1.push_back(inod);   // equates to vertex number of ghost node 2
            vectmp1.push_back(Vertex(iel,2));   // ghost node 3 coincident with a physical node 2

            mNGNodes++;

// last ghost on edge so skip corner node

            if(ielg==(NCells()+cells_on_edge[0]+cells_on_edge[1]+cells_on_edge[2]+cells_on_edge[3]-1)){
              mCoord.at(0).push_back(Coord(0,Vertex(iel,0))-(Coord(0,Vertex(iel,1))-Coord(0,Vertex(iel,0))));
              mCoord.at(1).push_back(Coord(1,Vertex(iel,1)));
              mNGNodes++;
            }

          break;
          case(2):
            mCoord.at(0).push_back(Coord(0,Vertex(iel,2)));
            mCoord.at(1).push_back(Coord(1,Vertex(iel,0))-(Coord(1,Vertex(iel,2))-Coord(1,Vertex(iel,0))));
            vectmp1.push_back(inod);   // equates to vertex number of ghost node 0
            vectmp1.push_back(inod+1); // equates to vertex number of ghost node 1
            vectmp1.push_back(Vertex(iel,0));   // ghost node 2 coincident with a physical node 0
            vectmp1.push_back(Vertex(iel,1));   // ghost node 3 coincident with a physical node 1
            mNGNodes++;

// last ghost on edge so skip corner node

            if(ielg==(NCells()+cells_on_edge[0]-1)){
              mCoord.at(0).push_back(Coord(0,Vertex(iel,3)));
              mCoord.at(1).push_back(Coord(1,Vertex(iel,1))-(Coord(1,Vertex(iel,3))-Coord(1,Vertex(iel,1))));
              mNGNodes++;
              inod++;
            }

          break;
          case(3):
            mCoord.at(0).push_back(Coord(0,Vertex(iel,1))+(Coord(0,Vertex(iel,1))-Coord(0,Vertex(iel,0))));
            mCoord.at(1).push_back(Coord(1,Vertex(iel,0)));
            vectmp1.push_back(Vertex(iel,1));   // ghost node 0 coincident with a physical node 1
            vectmp1.push_back(inod);   // equates to vertex number of ghost node 1
            vectmp1.push_back(Vertex(iel,3));   // ghost node 2 coincident with a physical node 3
            vectmp1.push_back(inod+1); // equates to vertex number of ghost node 3
            mNGNodes++;

 // last ghost on edge so skip corner node

            if(ielg==(NCells()+cells_on_edge[0]+cells_on_edge[1]-1)){
              mCoord.at(0).push_back(Coord(0,Vertex(iel,3))+(Coord(0,Vertex(iel,3))-Coord(0,Vertex(iel,2))));
              mCoord.at(1).push_back(Coord(1,Vertex(iel,2)));
              mNGNodes++;
              inod++;
            }

          break;
        }

// store the ghost vertices

        mVertex.push_back(vectmp1);

// assign material number of the physical element iel to the ghost element ielg

        mMaterial.push_back(mMaterial.at(iel));

      }
    }
  }

// add a ghost cell to each corner of the mesh

  int v1(NNodes()),v2(NNodes()+1); // vertex addresses either side of each corner ghost, used to find coordinates of the corner
  vector<int> vectmp1(4),neighbours={-1,-1,-1,-1};

  for(int iface=0;iface<4;iface++){
    switch(iface){
      case(0):
        v1+=cells_on_edge[0];v2=v1+1; // vertices at end of face/start of next face
        mCoord.at(0).push_back(Coord(0,v2)); // set coordinate of corner ghost
        mCoord.at(1).push_back(Coord(1,v1)); // set coordinate of corner ghost
        vectmp1.at(0)=v1;
        vectmp1.at(1)=NNodes()+NGNodes();
        vectmp1.at(2)=cells_on_edge[0];
        vectmp1.at(3)=v2;
      break;
      case(1):
        v1+=cells_on_edge[1]+1;v2=v1+1; // vertices at end of face/start of next face
        mCoord.at(0).push_back(Coord(0,v1)); // set coordinate of corner ghost
        mCoord.at(1).push_back(Coord(1,v2)); // set coordinate of corner ghost
        vectmp1.at(0)=NNodes()-1;
        vectmp1.at(1)=v1;
        vectmp1.at(2)=v2;
        vectmp1.at(3)=NNodes()+NGNodes();
      break;
      case(2):
        v1+=cells_on_edge[2]+1;v2=v1+1; // vertices at end of face/start of next face
        mCoord.at(0).push_back(Coord(0,v2)); // set coordinate of corner ghost
        mCoord.at(1).push_back(Coord(1,v1)); // set coordinate of corner ghost
        vectmp1.at(0)=v2;
        vectmp1.at(1)=NNodes()-cells_on_edge[2]-1;
        vectmp1.at(2)=NNodes()+NGNodes();
        vectmp1.at(3)=v1;
      break;
      case(3):
        v1+=cells_on_edge[3]+1;v2=NNodes(); // vertices at end of face/start of next face
        mCoord.at(0).push_back(Coord(0,v1)); // set coordinate of corner ghost
        mCoord.at(1).push_back(Coord(1,v2)); // set coordinate of corner ghost
        vectmp1.at(0)=NNodes()+NGNodes();
        vectmp1.at(1)=v2;
        vectmp1.at(2)=v1;
        vectmp1.at(3)=0;
      break;
    }

// add the ghost cell

    mNGCells++;
    mNGNodes++;
    mVertex.push_back(vectmp1);
    mE2E.push_back(neighbours);

  }

// connect ghost cells

  for(int ielg=NCells();ielg<NCells()+NGCells();ielg++){

// vertices of current ghost

    int n[5]={Vertex(ielg,0),Vertex(ielg,1),Vertex(ielg,3),Vertex(ielg,2),Vertex(ielg,0)};

// loop over faces of ghost ielg

    for(int ifaceg=0;ifaceg<NVertices(ielg);ifaceg++){

// search for a neighbour ghost on face ifaceg of ielg

      for(int ieln=NCells();ieln<NCells()+NGCells();ieln++){


        if(ielg==ieln){continue;}

// vertices of element ieln

        int nn[5]={Vertex(ieln,0),Vertex(ieln,1),Vertex(ieln,3),Vertex(ieln,2),Vertex(ieln,0)};

// loop over faces of element ieln

        for(int ifacen=0;ifacen<NVertices(ieln);ifacen++){

// check if faces ifaceg and ifacen are coincident

          bool z1(n[ifaceg]==nn[ifacen] || n[ifaceg]==nn[ifacen+1]);
          bool z2(n[ifaceg+1]==nn[ifacen] || n[ifaceg+1]==nn[ifacen+1]);

          if(z1&&z2){mE2E.at(ielg).at(ifaceg)=ieln;} // ieln on face ifaceg of ielg

        }

      }

    }

  }

// add material to corner ghosts

  int iel1(cells_on_edge[0]-1);         // bottom right corner
  mMaterial.push_back(mMaterial.at(iel1));

  int iel2(NCells()-1);                 // top right corner
  mMaterial.push_back(mMaterial.at(iel2));

  int iel3(NCells()-cells_on_edge[2]);  // top left corner
  mMaterial.push_back(mMaterial.at(iel3));

  int iel4(0);                          // bottom left corner
  mMaterial.push_back(mMaterial.at(iel4));

  cout<<" mE2E.size() "<<mE2E.size()<<endl;
  cout<<" mVertex.size() "<<mVertex.size()<<endl;
  cout<<" NGCells() "<<NGCells()<<endl;
  cout<<" NCells() "<<NCells()<<endl;
  cout<<" NGNodes() "<<NGNodes()<<endl;
  cout<<" NNodes() "<<NNodes()<<endl;

  return;

}

// function to convert a string to an int

constexpr unsigned int str2int(const char* str, int h = 0){return !str[h] ? 5381 : (str2int(str, h+1) * 33) ^ str[h];}

// member function to return the number of mesh dimensions

int Mesh::NDims() const {return mNDims;}

// member function to return the number of nodes in the generator mesh

long Mesh::NNodes() const {return mNNodes;}

// member functino to return the number of nodes in a FEM mesh containing order p shapes

long Mesh::NNodes(int p) const {

// first we need the number of cells in each direction

    long ncellsx(0),ncellsy(0);
    for(int i=0;i<NSides();i++){
      if(SideAttr(i)==0){ncellsx+=1;}
      if(SideAttr(i)==1){ncellsy+=1;}
    }

// return number of nodes in a continuous finite element method

    return ((ncellsx*p+1)*(ncellsy*p+1));

}

// member function to insert order p shape of type t=CONTINUOUS or t=DISCONTINUOUS into the mesh and return the number of nodes

long Mesh::NNodes(int p,int t) {

// number of nodes to be determined - this depends on whether the FEM is continuous or discontinuous

  long nnodes(0);

  if(t==DISCONTINUOUS){

// number of nodes in a discontinuous finite element method

    mNNodes_DFEM=NCells()*(p+1)*(p+1);
    int nloc((p+1)*(p+1));

// set up global node numbers for a discontinuous finite element method

    for(int i=0;i<NCells();i++){
      vector<long> global_nodes;
      for(int iloc=0;iloc<nloc;iloc++){
        global_nodes.push_back(i*nloc+iloc);
      }
      mGlobalNode_DFEM.push_back(global_nodes);
    }

// set return value

    nnodes=mNNodes_DFEM;

  }else{

// first we need the number of cells in each direction

    int ncellsx(NCells(0)),ncellsy(NCells(1));

// number of nodes in a continuous finite element method

    mNNodes_CFEM=(ncellsx*p+1)*(ncellsy*p+1);

// set up global node numbers for a continuous finite element method

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
        mGlobalNode_CFEM.push_back(global_nodes);
      }
    }

// set return value

    nnodes=mNNodes_CFEM;

  }

  return nnodes;

}

// member function to return the number of global nodes on the continuous finite element mesh

long Mesh::NNodes_CFEM() const {return mNNodes_CFEM;}

// member function to return the number of global nodes on the discontinuous finite element mesh

long Mesh::NNodes_DFEM() const {return mNNodes_DFEM;}

// member function to return the number of ghost nodes

int Mesh::NGNodes() const {return mNGNodes;}

// member function to return the number of materials

int Mesh::NMaterials() const {return mNMaterials;}

// member function to return the number of cells

int Mesh::NCells() const {return mNCells;}

// member function to return the number of cells on axis idim

int Mesh::NCells(int idim) const {

  int ncells[NDims()];
  ncells[0]=0;ncells[1]=0;

  for(int i=0;i<NSides();i++){
    if(SideAttr(i)==0){ncells[0]+=1;}
    if(SideAttr(i)==1){ncells[1]+=1;}
  }

  return ncells[idim];

}

// member function to return the number of ghost cells

int Mesh::NGCells() const {return mNGCells;}

// member function to return the material in cell i

int Mesh::Material(int i) const {return mMaterial.at(i);}

// member function to return the geometric type of cell i

int Mesh::Type(int i) const {return mType.at(i);}

// member function to return the number of vertices of cell i

int Mesh::NVertices(int i) const {return mVertex.at(i).size();}

// member function to return the node number of vertex j of cell i

int Mesh::Vertex(int i, int j) const {return mVertex[i][j];}

// member function to return the global node number of local node j in cell i in a continuous finite element method

long Mesh::GlobalNode_CFEM(int i, int j) const {return mGlobalNode_CFEM[i][j];}

// member function to return the global node number of local node j in cell i in a discontinuous finite element method

long Mesh::GlobalNode_DFEM(int i, int j) const {return mGlobalNode_DFEM[i][j];}

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

// member function to initialise vector x to the mesh coordinates

void Mesh::InitCoords(vector<vector<double> > &v,int const p,int const t){

// here we are taking the coordinates from the generator which only knows the
// cell corners (effectively a 4 node p1 element) and using a finite element
// method to construct the coordinates of all the nodes in the target element
// which is order p and has (p+1)*(p+1) nodes.

// First step is to declare a shape function M for the order p target element
// Second step is to construct a p1 4 node element with shape N that has the
// same number of integration points as shape M

// Then the 4 coordinates from the generator are expanded onto the nodes of the
// target element using the prolongation operator in the shape class

// first resize the first stride

  v.resize(NDims());

// set the number of nodes

  for(int idim=0;idim<NDims();idim++){
    v.at(idim).resize((t==DISCONTINUOUS)?NNodes_DFEM():NNodes_CFEM());
  }

// declare an order p shape M whose nodes we need to populate with coordinates

  Shape M(p); // for M the type doesn't really matter here and the default number of integration points is OK

// declare p1 element to source the generator coordinates

  Shape N(1,sqrt(M.ngi()),CONTINUOUS); // N must have the same number of integration points as M, so declare with M.ngi()

// expand the coordinates of stencil N on to stencil M 

  for(int idim=0;idim<NDims();idim++){

    for(int i=0;i<NCells();i++){

      double xn[N.nloc()]; // cell corner coordinates
      double xm[M.nloc()]; // element node coordinates

// load generator coordinates from the mesh class on to the 4 node stencil N

      for(int iloc=0;iloc<N.nloc();iloc++){xn[iloc]=Coord(idim,Vertex(i,iloc));}

// apply prolongation operator to populate xm

      N.prolongate(xn,xm,M.order());

// commit new coordinates to return address space, each coordinate address is unique if discontinuous so we can append

      if(t==CONTINUOUS){
        for(int iloc=0;iloc<M.nloc();iloc++){v.at(idim).at(GlobalNode_CFEM(i,iloc))=xm[iloc];}
      }else{
        for(int iloc=0;iloc<M.nloc();iloc++){v.at(idim).at(GlobalNode_DFEM(i,iloc))=xm[iloc];}
      }

    }

  }

  return;

}

// member function to advect coordinate x with velocity u a distance u*dt

void Mesh::UpdateCoords(VVD &x, VVD const &u, double const dt) const{

// loop over dimension and advect nodes a distance u*dt

  for(int idim=0;idim<NDims();idim++){

    for(long i=0;i<x.at(idim).size();i++){

      x.at(idim).at(i)=x.at(idim).at(i)+u.at(idim).at(i)*dt;

    }

  }

  return;

}

// member function to map coordinates from an order p mesh to an order q mesh

void Mesh::MapCoords(VVD const &xp,VVD &xq,int const &p,int const &q) const{

// declare a shape for each element

  Shape P(p,p+1,CONTINUOUS);
  Shape Q(q,sqrt(P.ngi()),DISCONTINUOUS);

// loop over dimension and advect nodes a distance u*dt

  for(int idim=0;idim<NDims();idim++){

    for(long i=0;i<NCells();i++){

      double xptmp[P.nloc()]; // element P node coordinates
      double xqtmp[Q.nloc()]; // element Q node coordinates

// load coordinates from the mesh on stencil p

      for(int iloc=0;iloc<P.nloc();iloc++){xptmp[iloc]=xp.at(idim).at(GlobalNode_CFEM(i,iloc));}

// apply prolongation operator to populate xqtmp

      P.prolongate(xptmp,xqtmp,Q.order());

// commit new coordinates to return address space, each coordinate address is unique if discontinuous so we can append

      if(Q.type()==CONTINUOUS){
        for(int iloc=0;iloc<Q.nloc();iloc++){xq.at(idim).at(GlobalNode_CFEM(i,iloc))=xqtmp[iloc];}
      }else{
        for(int iloc=0;iloc<Q.nloc();iloc++){xq.at(idim).at(GlobalNode_DFEM(i,iloc))=xqtmp[iloc];}
      }

    }
  }

  return;

}

// initialise element length scale

void Mesh::InitLength(VD &l,int const &p,VD const &V) const{

// partition each cell volume based on p - this will normalise it to the element order

  for(int i=0;i<NCells();i++){
    l.at(i)=sqrt(V.at(i))/p;
  }

  return;

}

// update element length scale, just use the same definition as for the initial length until we understand how to normalise it properly for high-order

double Mesh::UpdateLength(int const &p,double const &V){return sqrt(V)/p;}

// update volume field

void Mesh::UpdateVolume(VD &V,VVD const &x, int const &p) const{

// declare a shape function

  Shape S(p);

// initialise a Jacobian

  vector<double> detJ(S.ngi());

// loop over the mesh and update volume in each cell

  for(int i=0;i<V.size();i++){

// loop over quadrature points and calculate the jacobian

    for(int gi=0;gi<S.ngi();gi++){

      double dxdu(0.0),dydu(0.0),dxdv(0.0),dydv(0.0);

// derivatives of the physical coordinates at the quadrature points

      for(int j=0;j<S.nloc();j++){
        dxdu+=x.at(0).at(GlobalNode_CFEM(i,j))*S.dvalue(0,j,gi); // dx/du
        dxdv+=x.at(0).at(GlobalNode_CFEM(i,j))*S.dvalue(1,j,gi); // dx/dv
        dydu+=x.at(1).at(GlobalNode_CFEM(i,j))*S.dvalue(0,j,gi); // dy/du
        dydv+=x.at(1).at(GlobalNode_CFEM(i,j))*S.dvalue(1,j,gi); // dy/dv
      }

// calculate the determinant at the quadrature point and commit to the vector

      detJ.at(gi)=dxdu*dydv-dxdv*dydu;

    }

// integrate the jacobian and reset the cell volume

    V.at(i)=0.0;
    for(int gi=0;gi<S.ngi();gi++){
      V.at(i)+=detJ[gi]*S.wgt(gi);
    }

    if(V.at(i)<VTOL){cout<<"Mesh::UpdateVolume(): -'ve volume detected in cell "<<i<<endl;exit(1);}

  }

  return;

}

// update density field

void Mesh::UpdateDensity(VD &d,VD const &V,VD const &m) const{

  for(int i=0;i<d.size();i++){
    d.at(i)=m.at(i)/V.at(i);
  }

  return;

}

// update internal energy field

void Mesh::UpdateEnergy(VD const &e0,VD &e1,VD const &p,VD const &q,VD const &V0,VD const &V1,VD const &m) const{

  for(int i=0;i<e0.size();i++){
    e1.at(i)=max(ECUT,e0.at(i)-((p.at(i)+q.at(i))*(V1.at(i)-V0.at(i)))/m.at(i));
  }

  return;

}

// load pressure field

void Mesh::UpdatePressure(VD &p,VD const &d,VD const &e,VD const &gamma,vector<int> const &mat){

  for(int i=0;i<p.size();i++){
    p.at(i)=P(d[i],e[i],gamma[mat[i]-1]);
    if(p[i]<0.0){cout<<"Mesh::UpdatePressure(): -'ve pressure detected in cell "<<i<<" e= "<<e[i]<<endl;exit(1);}
  }

  return;

}

// member function to return pressure field value at a point

double Mesh::UpdatePressure(double const&d,double const&e,double const&g){return P(d,e,g);}

// load new sound speeds

void Mesh::UpdateSoundSpeed(VD &c,VD const &g,vector<int> const &mat,VD const &p,VD const &d) const{

  for(int i=0;i<c.size();i++){
    c.at(i)=C(p.at(i),d.at(i),g[mat[i]-1]);
  }

  return;

}

// member function to return sound speed at a point

double Mesh::UpdateSoundSpeed(double const &g,double const &p,double const &d){return C(p,d,g);}

// member function to return artificial viscosity at a point


double Mesh::UpdateQ(double const&l,double const&d,double const&c,double const&cq,double const&cl,double const&divu){

return( (divu<0.0)?(d*l*divu*((cq*l*divu)-cl*c)):0.0);

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
