// Function definitions for members of the mesh class

// Author S. R. Merton

#include "mesh.h"
#include "math.h" // for pow() function
#include "shape.h"
#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>
#include <iomanip>
#include <string>

#define V1 vector<double>

void nodal_mass(Shape*S[],Mesh*M,V1*nodmass); // compute nodal masses

using namespace std;

Mesh::Mesh(char* meshfile){

// Constructor to instantiate a new mesh from data in the meshfile
//
// Arguments:
// meshfile is a pointer to a character containing the filename of the mesh file

  cout<<"Mesh::Mesh(): Loading a mesh from file "<<meshfile<<endl;

// load the mesh (google "c++ file iterator" for more interesting ways...)

  vector<string> mvec; string str;                    // vector to store file and a string to store a line
  ifstream meshdata(meshfile);                        // open the mesh file
  while(getline(meshdata, str)){mvec.push_back(str);} // read the whole file into mVec
  meshdata.close();                                   // close the mesh file

  cout<<"Mesh file "<<meshfile<<" accessed, "<<mvec.size()<<" lines loaded."<<endl;

// now unpack the mesh into the member data working backwards as we don't know how many comment lines are at the top

  long iaddr=mvec.size()-1;              // get address of the end of the data
  mNDims=stoi(mvec[iaddr]);

  iaddr=iaddr-NDims()-1;                 // move backwards to next data address
  for(int idim=0;idim<NDims();idim++){
    mDim.push_back(stol(mvec[iaddr+idim]));
  }

  iaddr=iaddr-2*NDims()-1;                // move backwards to next data address
  for(int ibc=0;ibc<2*NDims();ibc++){
    mBoundaryType.push_back(mvec[iaddr+ibc]);
  }

  iaddr=iaddr-2;                         // move backwards to next data address
  mNVertices=stol(mvec[iaddr]);

  iaddr=iaddr-2;                         // move backwards to next data address
  mNCells=stol(mvec[iaddr]);

  iaddr=iaddr-2;                         // move backwards to next data address
  mNMaterials=stoi(mvec[iaddr]);

  iaddr=iaddr-NMaterials()-1;            // move backwards to next data address
  for(int imat=0;imat<NMaterials();imat++){
    mMaterials.push_back(stoi(mvec[iaddr+imat]));
  }

  iaddr=iaddr-NMaterials()-1;            // move backwards to next data address
  for(int imat=0;imat<NMaterials();imat++){
    mVelocityZInit.push_back(stod(mvec[iaddr+imat]));
  }

  iaddr=iaddr-NMaterials()-1;            // move backwards to next data address
  for(int imat=0;imat<NMaterials();imat++){
    mVelocityYInit.push_back(stod(mvec[iaddr+imat]));
  }

  iaddr=iaddr-NMaterials()-1;            // move backwards to next data address
  for(int imat=0;imat<NMaterials();imat++){
    mVelocityXInit.push_back(stod(mvec[iaddr+imat]));
  }

  iaddr=iaddr-NMaterials()-1;            // move backwards to next data address
  for(int imat=0;imat<NMaterials();imat++){
    mPressureInit.push_back(stod(mvec[iaddr+imat]));
  }

  iaddr=iaddr-NMaterials()-1;            // move backwards to next data address
  for(int imat=0;imat<NMaterials();imat++){
    mDensityInit.push_back(stod(mvec[iaddr+imat]));
  }

  iaddr=iaddr-NMaterials()-1;            // move backwards to next data address
  for(int imat=0;imat<NMaterials();imat++){
    mEnergyInit.push_back(stod(mvec[iaddr+imat]));
  }

  iaddr=iaddr-2;                         // move backwards to next data address
  mNTypes=stoi(mvec[iaddr]);

  iaddr=iaddr-NTypes()-1;                // move backwards to next data address
  for(int itype=0;itype<NTypes();itype++){
    mTypes.push_back(stoi(mvec[iaddr+itype]));
  }

  iaddr=iaddr-NTypes()-1;                // move backwards to next data address
  for(int itype=0;itype<NTypes();itype++){
    mNames.push_back(mvec[iaddr+itype]);
  }

  iaddr=iaddr-2;                         // move backwards to next data address
  mNOrders=stoi(mvec[iaddr]);

  iaddr=iaddr-NOrders()-1;               // move backwards to next data address
  for(int iorder=0;iorder<NOrders();iorder++){
    int pinput=stoi(mvec[iaddr+iorder]);
    if(pinput==0)mUsePdV=true;           // turn on PdV energy field
    mTOrders.push_back(max(1,pinput));
    mOrders.push_back(1);
  }

  iaddr=iaddr-NOrders()-1;               // move backwards to next data address
  for(int iorder=0;iorder<NOrders();iorder++){
    mKOrders.push_back(stoi(mvec[iaddr+iorder]));
  }

  iaddr=iaddr-NCells()-1;                // move backwards to next data address
  for(long iel=0;iel<NCells();iel++){
    int pinput=stoi(mvec[iaddr+iel]);
    if(pinput==0)mUsePdV=true;           // turn on PdV energy field
    mTOrder.push_back(max(1,pinput));
    mOrder.push_back(1);
  }

  iaddr=iaddr-NCells()-1;                // move backwards to next data address
  for(long iel=0;iel<NCells();iel++){
    mKOrder.push_back(stoi(mvec[iaddr+iel]));
  }

  iaddr=iaddr-NCells()-1;                // move backwards to next data address
  for(long iel=0;iel<NCells();iel++){
    mType.push_back(stoi(mvec[iaddr+iel]));
  }

  iaddr=iaddr-NCells()-1;                // move backwards to next data address
  for(long iel=0;iel<NCells();iel++){
    mMaterial.push_back(stoi(mvec[iaddr+iel]));
  }

  iaddr=iaddr-NCells()-1;                // move backwards to next data address
  for(long iel=0;iel<NCells();iel++){
    mNLoc.push_back(stoi(mvec[iaddr+iel]));
  }

  for(auto&n:mNLoc){mNCorners+=n;}       // accumulate above vector to get length of the vertex list

  iaddr=iaddr-NCorners()-1;              // move backwards to next data address
  for(long icorner=0;icorner<NCorners();icorner++){
    mVertex.push_back(stol(mvec[iaddr+icorner]));
  }

  vector<double> vtxtmp[NDims()];
  for(int idim=0;idim<NDims();idim++){
    iaddr=iaddr-NVertices()-1;             // move backwards to next data address
    for(long ivert=0;ivert<NVertices();ivert++){
      vtxtmp[idim].push_back(stod(mvec[iaddr+ivert]));
    }
    mCoord.push_back(vtxtmp[idim]);
  }

// sanity checks

  for(long iel=0;iel<NCells();iel++){
    if(!UsePdV() && KOrder(iel)<=TOrder(iel)){
      cout<<"Error: Kinematics require higher element order than thermodynamics:"<<endl;
      cout<<"  Kinematics order= "<<KOrder(iel)<<endl;
      cout<<"  Thermodynamics order= "<<TOrder(iel)<<endl;
      exit(1);
    }
  }

// print out what we have found

  cout<<fixed<<setprecision(9); // set the output precision

  cout<<"Mesh::Mesh(): NDims= "<<NDims()<<endl;;

  for(int idim=0;idim<NDims();idim++){
    cout<<"Mesh::Mesh():  dimension "<<idim+1<<" = "<<Dim(idim)<<endl;
  }

  cout<<"Mesh::Mesh(): NVertices= "<<NVertices()<<endl;
  cout<<"Mesh::Mesh(): NCells= "<<NCells()<<endl;
  cout<<"Mesh::Mesh(): NMaterials= "<<NMaterials()<<endl;

  for(int imat=0;imat<NMaterials();imat++){
    cout<<"Mesh::Mesh():  material "<<imat+1<<" = "<<Materials(imat)<<endl;
  }

  cout<<"Mesh::Mesh(): NTypes= "<<NTypes()<<endl;

  for(int itype=0;itype<NTypes();itype++){
    cout<<"Mesh::Mesh():  polyhedral type "<<itype+1<<" = "<<Types(itype)<<" (name = "<<Names(itype)<<")"<<endl;
  }

  cout<<"Mesh::Mesh(): NOrders= "<<NOrders()<<endl;

  for(int iorder=0;iorder<NOrders();iorder++){
    cout<<"Mesh::Mesh():  kinematic polyhedral order "<<iorder+1<<" = "<<KOrders(iorder)<<endl;
  }

  for(int iorder=0;iorder<NOrders();iorder++){
    cout<<"Mesh::Mesh():  thermodynamics polyhedral order "<<iorder+1<<" = "<<TOrders(iorder)<<endl;
  }

  for(long iel=0;iel<NCells();iel++){
    cout<<"Mesh::Mesh():  element "<<iel+1<<" has material "<<Material(iel)<<endl;
  }

  cout<<"Mesh::Mesh(): NCorners= "<<NCorners()<<endl;

// insert requested element into each cell loaded in from the meshfile

  fill_dg();

  if(UsePdV()){
    cout<<"Mesh::Mesh(): PdV energy field active."<<endl;
  }else{
    cout<<"Mesh::Mesh(): Finite element energy field active, "<<TLoc(0)<<" nodes on the thermodynamic grid"<<endl;
  }

// generate connectivity

  set_conn();

// generate halo

  set_ghosts();

  cout<<"Mesh::Mesh():  "<<NNodes()<<" DG elements inserted in to the mesh !"<<endl;
  cout<<"Mesh::Mesh():  "<<NGhosts()<<" Ghost cells generated."<<endl;
  cout<<"Mesh::Mesh(): Mesh has been loaded."<<endl;

// insert dummy elements for silo

  insert_dummies();

// form node list along mesh boundaries

  boundary_nodes();

// form list of coincident nodes

  coincident_nodes();

  return;

}

// Destructor function to release storage associated with a mesh class object

Mesh::~Mesh(){}

// function to insert requested element into each cell read in from the dump

void Mesh::fill_dg(){

  cout<<"Mesh::fill_dg():    inserting requested DG element type..."<<endl;

  vector<double> vectmp[NDims()];  // vector temporary to store node coordinates for all elements
  vector<double> vectmpk[NDims()]; // vector temporary to store kinematic node coordinates for all elements
  vector<double> vectmpt[NDims()]; // vector temporary to store thermodynamic node coordinates for all elements
  mNNodes=0l;                      // count global nodes as they are inserted into the mesh
  mKNodes=0l;                      // count kinematic nodes as they are inserted into the mesh
  mTNodes=0l;                      // count thermodynamic nodes as they are inserted into the mesh

  for(long iel=0;iel<NCells();iel++){

    int p[3]={Order(iel),TOrder(iel),KOrder(iel)}; // polyhedral orders of the element we want to insert into mesh cell iel
    vector<long> vecnod;           // vector of global node numbers to insert for mesh cell iel
    vector<long> vecnodk;          // vector of global kinematic node numbers to insert for mesh cell iel
    vector<long> vecnodt;          // vector of global thermodynamic node numbers to insert for mesh cell iel
    vector<int> nfnodes;           // number of local nodes on each face
    vector<vector<int> > fnodes;   // local nodes on each face

    switch(Types(Type(iel)-1)){
      case 2001:{

        mMLoc.push_back((p[0]+1)*(p[0]+1)); // set the number of DG nodes to insert into the element
        mTLoc.push_back((p[1]+1)*(p[1]+1)); // set the number of thermodynamic nodes to insert into the element
        mKLoc.push_back((p[2]+1)*(p[2]+1)); // set the number of kinematic nodes to insert into the element
        vector<double> n[NLoc(iel)];    // shape function for the element from the mesh file
        vector<double> nk[NLoc(iel)];   // shape function for the element from the mesh file
        vector<double> nt[NLoc(iel)];   // shape function for the element from the mesh file
        vector<double> u;               // edge position of each node to insert into iel
        vector<double> uk;              // edge position of each node to insert into iel
        vector<double> ut;              // edge position of each node to insert into iel

// generate edge position of each DG node to insert into element iel

        for(int i=0;i<p[0]+1;i++){u.push_back(i*1.0d/p[0]);}

// evaluate a linear (4-node) shape function at each position where there is a DG node to insert

        for(int i=0;i<p[0]+1;i++){
          for(int j=0;j<p[0]+1;j++){
            n[0].push_back( (1.0d-u[j])*(1.0d-u[i]) );
            n[1].push_back( u[j]*(1.0d-u[i]) );
            n[2].push_back( u[j]*u[i] );
            n[3].push_back( (1.0d-u[j])*u[i] );
          }
        }

// generate edge position of each kinematic node to insert into element iel

        for(int i=0;i<p[2]+1;i++){uk.push_back(i*1.0d/(p[2]));}

// evaluate a linear (4-node) shape function at each position where there is a kinematic node to insert

        for(int i=0;i<p[2]+1;i++){
          for(int j=0;j<p[2]+1;j++){
            nk[0].push_back( (1.0d-uk[j])*(1.0d-uk[i]) );
            nk[1].push_back( uk[j]*(1.0d-uk[i]) );
            nk[2].push_back( uk[j]*uk[i] );
            nk[3].push_back( (1.0d-uk[j])*uk[i] );
          }
        }

// generate edge position of each thermodynamic node to insert into element iel

        for(int i=0;i<p[1]+1;i++){ut.push_back(i*1.0d/p[1]);}

// evaluate a linear (4-node) shape function at each position where there is a thermodynamic node to insert

        for(int i=0;i<p[1]+1;i++){
          for(int j=0;j<p[1]+1;j++){
            nt[0].push_back( (1.0d-ut[j])*(1.0d-ut[i]) );
            nt[1].push_back( ut[j]*(1.0d-ut[i]) );
            nt[2].push_back( ut[j]*ut[i] );
            nt[3].push_back( (1.0d-ut[j])*ut[i] );
          }
        }

// interpolate mesh vertices to find each DG node position to insert into element iel

        for(int idim=0;idim<NDims();idim++){
          for(int iloc=0;iloc<MLoc(iel);iloc++){
            double new_coord(0.0d);
            for(int j=0;j<NLoc(iel);j++){
              new_coord+=n[j][iloc]*Coord(idim,Vertex(iel,j));
            }
            vectmp[idim].push_back(new_coord);
          }
        }

// interpolate mesh vertices to find each kinematic node position to insert into element iel

        for(int idim=0;idim<NDims();idim++){
          for(int iloc=0;iloc<KLoc(iel);iloc++){
            double new_coord(0.0d);
            for(int j=0;j<NLoc(iel);j++){
              new_coord+=nk[j][iloc]*Coord(idim,Vertex(iel,j));
            }
            vectmpk[idim].push_back(new_coord);
          }
        }

// interpolate mesh vertices to find each thermodynamic node position to insert into element iel

        for(int idim=0;idim<NDims();idim++){
          for(int iloc=0;iloc<TLoc(iel);iloc++){
            double new_coord(0.0d);
            for(int j=0;j<NLoc(iel);j++){
              new_coord+=nt[j][iloc]*Coord(idim,Vertex(iel,j));
            }
            vectmpt[idim].push_back(new_coord);
          }
        }

// store the global address of the nodes inserted into element iel

        for(int iloc=0;iloc<MLoc(iel);iloc++){
          vecnod.push_back(mNNodes);
          mNNodes++;
        }

// store the global address of the kinematic nodes inserted into element iel

        for(int iloc=0;iloc<KLoc(iel);iloc++){
          vecnodk.push_back(mKNodes);
          mKNodes++;
        }

// store the global address of the thermodynamic nodes inserted into element iel

        for(int iloc=0;iloc<TLoc(iel);iloc++){
          vecnodt.push_back(mTNodes);
          mTNodes++;
        }

// set up face nodes on the surface of the thermodynamic element

        mNFaces.push_back(4);           // set the number of faces on the surface of this element

        for(int iface=0;iface<4;iface++){nfnodes.push_back(p[1]+1);}

        mNFLoc.push_back(nfnodes);

        vector<int> face1,face2,face3,face4;
        vector<vector<int> >faces;

        for(int i=0;i<p[1]+1;i++){
          face1.push_back(i);
          face2.push_back((i+1)*(p[1]+1)-1);
          face3.push_back((p[1]+1)*(p[1]+1)-i-1);
          face4.push_back((p[1]+1)*(p[1]+1)-(i+1)*(p[1]+1));
        }

        faces.push_back(face1);
        faces.push_back(face2);
        faces.push_back(face3);
        faces.push_back(face4);

        mFLoc.push_back(faces);

// set up face nodes on the surface of the kinematic element

        vector<int> face1k,face2k,face3k,face4k;
        vector<vector<int> >facesk;

        for(int i=0;i<p[2]+1;i++){
          face1k.push_back(i);
          face2k.push_back((i+1)*(p[2]+1)-1);
          face3k.push_back((p[2]+1)*(p[2]+1)-i-1);
          face4k.push_back((p[2]+1)*(p[2]+1)-(i+1)*(p[2]+1));
        }

        facesk.push_back(face1k);
        facesk.push_back(face2k);
        facesk.push_back(face3k);
        facesk.push_back(face4k);

        mFLock.push_back(facesk);

// append global node numbers generated in element iel

        mNode.push_back(vecnod);
        mKNode.push_back(vecnodk);
        mTNode.push_back(vecnodt);

        break;}
      default:
        cout<<"Mesh::fill_dg():  Type "<<Types(Type(iel)-1)<<" not supported - aborting."<<endl;
        exit(1);
    }

  }

// insert the new coordinates into the mesh

  for(int idim=0;idim<NDims();idim++){mDGCoord.push_back(vectmp[idim]);}
  for(int idim=0;idim<NDims();idim++){mKCoord0.push_back(vectmpk[idim]);}
  for(int idim=0;idim<NDims();idim++){mKCoord.push_back(vectmpk[idim]);}
  for(int idim=0;idim<NDims();idim++){mTCoord0.push_back(vectmpt[idim]);}
  for(int idim=0;idim<NDims();idim++){mTCoord.push_back(vectmpt[idim]);}

  return;

}

// perform connectivity traversal to collect adjacency information between neighbouring elements

void Mesh::set_conn(){

  cout<<"Mesh::set_conn():  performing adjacency search..."<<endl;

// sweep over elements

  for(long iel=0;iel<NCells();iel++){

// set up an array to store the neighbour face that each face of element iel connects to

    int faces[NFaces(iel)]; // neighbour faces
    for(int iface=0;iface<NFaces(iel);iface++){faces[iface]=-1;} // add ghost cell neighbour face if no neighbour face is found

// set up an array to store the element number of each neighbour that element iel has

    long neighbours[NFaces(iel)]; // neighbour elements
    for(int iface=0;iface<NFaces(iel);iface++){neighbours[iface]=-1l;} // add ghost cell neighbour if no neighbour element is found

// vector to store a vector of neighbour nodes for each face of iel

    vector<vector<int> > n2nvec;

// visit each face on the surface of element iel

    for(int iface=0;iface<NFaces(iel);iface++){

// set up an array to store the nodes on neighbour face that each node on face iface connects to

      int nodes[NFLoc(iel,iface)]; // nodes connecting to neighbour face

// visit all other elements on the mesh and check for adjacency on iface of iel

      for(long ieln=0;ieln<NCells();ieln++){

        if(iel==ieln)continue;

        for(int ifacen=0;ifacen<NFaces(ieln);ifacen++){

// check for coincident nodes on faces iface and ifacen

          if(NFLoc(ieln,ifacen)==NFLoc(iel,iface)){ // Need to update this for different NFLoc() on ieln...

// sum the connections and store face locations where there is mapping 
// note that ieln and iel are only adjacent if ALL nodes on a face are coincident

            int nconns=0;int flocs[NFLoc(iel,iface)];
            for(int iloc=0;iloc<NFLoc(iel,iface);iloc++){
              for(int jloc=0;jloc<NFLoc(ieln,ifacen);jloc++){

                int check=0;
                long inod=Node(iel,FLoc(iel,iface,iloc));
                long inodn=Node(ieln,FLoc(ieln,ifacen,jloc));
                for(int idim=0;idim<NDims();idim++){
                  check+=(DGCoord(idim,inod)==DGCoord(idim,inodn));
                }

// make a connection between inod and inodn if all coordinate components are coincident

                if(check==NDims()){flocs[iloc]=jloc;nconns++;}

              }
            }

// store the adjacency information across face iface of element iel

            if(nconns==NFLoc(iel,iface)){
              faces[iface]=ifacen;neighbours[iface]=ieln;
              for(int iloc=0;iloc<NFLoc(iel,iface);iloc++){nodes[iloc]=flocs[iloc];}
            }

          } // if(NFLoc()...)

        } // ifacen

      } // ieln

// pack neighbour node numbers on face iface into nodal connectivity vector

      vector<int> n2nvec_iface;
      for(int iloc=0;iloc<NFLoc(iel,iface);iloc++){n2nvec_iface.push_back(nodes[iloc]);}
      n2nvec.push_back(n2nvec_iface);

    } // iface

// pack faces of neighbours into face connectivity vector

    vector<int> f2fvec;

    for(int iface=0;iface<NFaces(iel);iface++){f2fvec.push_back(faces[iface]);}

    mF2F.push_back(f2fvec);

// pack element number of each neighbour into element connectivity vector

    vector<long> e2evec;

    for(int iface=0;iface<NFaces(iel);iface++){e2evec.push_back(neighbours[iface]);}

    mE2E.push_back(e2evec);

// pack neighbour node numbers on all faces into nodal connectivity vector

    mN2N.push_back(n2nvec);

  } // iel

  return;

}

// set up halo address space and associated connectivity

void Mesh::set_ghosts(){

// set up halo address space

  long mNGNodes=NNodes()-1; // don't modify mNNodes here as this will mess things up for silo
  int faces[6]={0,1,2,3,0,1}; // opposing faces for quads use iface+2 to access
//  int faces[6]={0,1,2,3,4,5,0,1,2}; // opposing faces for hexes use iface+3 to access

  long mGKNodes=mKNodes; // global kinematic node number in ghost cells

  for(long iel=0;iel<NCells();iel++){
    for(int iface=0;iface<NFaces(iel);iface++){
      long ieln=E2E(iel,iface); // element number of neighbour
      if(ieln<0){

        vector<long> vecnodgk; // all global kinematic node numbers in a ghost cell

// assign connectivity information to the ghost cell

        mNGhosts++; // increase ghost count
        long ielg=mNGhosts+NCells()-1; // assign a ghost cell number offset from the physical cell number
        long ifaceg=faces[iface+2]; // face number of the ghost cell ielg  connecting to face iface of iel
        mE2E[iel].at(iface)=ielg; // modify initial data with the ghost cell address
        mF2F[iel].at(iface)=ifaceg; // modify initial face data with the ghost cell face adjacent to physical face
        for(int isloc=0;isloc<NFLoc(iel,iface);isloc++){
          mN2N[iel][iface].at(isloc)=isloc; // assume a 1-1 mapping between surface nodes on the ghost face
        }

// insert global nodes in to the ghost element ielg

        vector<long> vecnod;
        mMLoc.push_back(MLoc(iel)); // the number of nodes matches the physical cell
        for(int iloc=0;iloc<MLoc(ielg);iloc++){
          mNGNodes++;
          vecnod.push_back(mNGNodes);
        }
        mNode.push_back(vecnod);

// set up face nodes on the surface of the ghost element

        mNFaces.push_back(4);           // set the number of faces on the surface of this ghost element

        vector<int> nfnodes;
        int p(Order(iel));            // polyhedral order of the physical element
        for(int iface=0;iface<4;iface++){nfnodes.push_back(p+1);}

        mNFLoc.push_back(nfnodes);

        vector<int> face1,face2,face3,face4;
        vector<vector<int> >ghost_faces;

        for(int i=0;i<p+1;i++){
//          face1.push_back(i); // had to tweak the numbering for ghosts as there are no coordinates to match up the nodes
          face1.push_back(p-i); // had to tweak the numbering for ghosts as there are no coordinates to match up the nodes
//          face2.push_back((i+1)*(p+1)-1);
          face2.push_back((p+1)*(p+1)-1-i*(p+1));
//          face3.push_back((p+1)*(p+1)-i-1);
          face3.push_back((p+1)*(p+1)-(p+1)+i);
//          face4.push_back((p+1)*(p+1)-(i+1)*(p+1));
          face4.push_back(i*(p+1));
        }

        ghost_faces.push_back(face1);
        ghost_faces.push_back(face2);
        ghost_faces.push_back(face3);
        ghost_faces.push_back(face4);

        mFLoc.push_back(ghost_faces);

// store the global address of the kinematic nodes inserted into this ghost element

        for(int iloc=0;iloc<KLoc(iel);iloc++){
          vecnodgk.push_back(mGKNodes);
          mGKNodes++;
        }

// append global node numbers generated in ghost element

        mKNode.push_back(vecnodgk);

      }

    }
  }

  cout<<"Number of ghosts in the level 1 halo = "<<NGhosts()<<endl;

// test

  for(long iel=0;iel<NCells();iel++){
    for(int iface=0;iface<NFaces(iel);iface++){
      int ifacen=F2F(iel,iface); // face of neighbour
      long ieln=E2E(iel,iface); // element number of neighbour
      if(ieln>NCells()-1){
        cout<<"*element "<<iel<<" face "<<iface<<" -> element "<<ieln<<" face "<<ifacen;
      }else{
        cout<<" element "<<iel<<" face "<<iface<<" -> element "<<ieln<<" face "<<ifacen;
      }
      cout<<" : node connections to neighbour= ";
      for(int iloc=0;iloc<NFLoc(iel,iface);iloc++){
        int jloc=N2N(iel,iface,iloc);
        cout<<Node(iel,FLoc(iel,iface,iloc))<<"->"<<Node(ieln,FLoc(ieln,ifacen,jloc))<<",";
//        cout<<Node(iel,FLoc(iel,iface,iloc))<<"->"<<FLoc(ieln,ifacen,jloc)<<",";
      }
      cout<<endl;
    }
  }

  return;

}

// traverse mesh boundary and store the nodes that are found there, this is useful for setting boundary conditions

void Mesh::boundary_nodes(){

  cout<<"Mesh::boundary_nodes(): Generating node list along mesh boundary..."<<endl;

// vector to store nodes on each face

  vector<long> nodes_on_face[NFaces(0)];

// traverse mesh and store the kinematic nodes that are on each boundary

  for(long iel=0;iel<NCells();iel++){

    for(int iface=0;iface<NFaces(iel);iface++){

// look up neighbour element on this face

      long ieln=E2E(iel,iface);

// check if face lies on a mesh edge

      if(ieln>NCells()-1){

// loop over nodes on face and fetch global number

        for(int iloc=0;iloc<KOrder(iel)+1;iloc++){
          long inod=KNode(iel,FLock(iel,iface,iloc));
          nodes_on_face[iface].push_back(inod);
        }

      }

    }

  }

// commit to mesh class address space

  for(int iface=0;iface<NFaces(0);iface++){
    mBoundaryNode.push_back(nodes_on_face[iface]);
  }

// output what has been found

  for(int iface=0;iface<NFaces(0);iface++){
    cout<<"Mesh::boundary_nodes(): Mesh edge "<<iface<<" contains "<<NBoundaryNodes(iface)<<" nodes:"<<endl;
    for(long inod=0;inod<NBoundaryNodes(iface);inod++){
      cout<<"Mesh::boundary_nodes():   node "<<BoundaryNode(iface,inod)<<" found on mesh face "<<iface<<endl;
    }
  }

  cout<<"Mesh::boundary_nodes(): Done."<<endl;

  return;

}

// member function to generate a list of coincident nodes
// This is useful for condensing the velocity field for
// mesh movement schemes e.g. mass weighted averaging.

void Mesh::coincident_nodes(){

// vector to store coincident nodes

  vector<long> tmp[KNodes()];

// opposite face

  int const oface[4]={2,3,0,1};

// traverse the mesh and visit each kniematic node on each face

  for(long iel=0;iel<NCells();iel++){

// polyhedral order of the kinematic grid

    int p(KOrder(iel));

    for(int iface=0;iface<NFaces(iel);iface++){

// neighbour element on this face

      long ieln=E2E(iel,iface);

// avoid coincident nodes that are in ghost cells

      if(ieln>NCells()-1) continue;

      for(int iloc=0;iloc<p+1;iloc++){

// global number of local node on face iface of element iel

        long iglo=KNode(iel,FLock(iel,iface,iloc));

// global number of coincident node on adjacent neighbour face

        long jglo=KNode(ieln,FLock(ieln,oface[iface],p-iloc));

// store coincident node numbers

        tmp[iglo].push_back(jglo);

      }

    }

  }

// add connections that were too far away to be found on first pass

  for(long inod=0;inod<KNodes();inod++){
    for(int i=0;i<tmp[inod].size();i++){
      long jnod=tmp[inod][i]; // 15,20
      for(int j=0;j<tmp[jnod].size();j++){
        long knod=tmp[jnod][j]; // 8,27 then 8,27

// add if not on existing list

        bool zfound(false);
        for(int k=0;k<tmp[inod].size();k++){
          long lnod=tmp[inod][k];
          if(lnod==knod){zfound=true;}
        }
//        if(!zfound && knod!=inod){tmp[inod].push_back(knod);}
        if(!zfound){tmp[inod].push_back(knod);} // bit more convenient ??

      }
    }
  }

// commit to mesh class

  for(long inod=0;inod<KNodes();inod++){
    if(tmp[inod].size()==0){tmp[inod].push_back(inod);}// include parent node
    mCoincidentNode.push_back(tmp[inod]);
  }

// debug
//  for(long inod=0;inod<KNodes();inod++){
//    cout<<" Node "<<inod<<" coincident with";
//    for(int j=0;j<NCoincidentNodes(inod);j++){
//      cout<<" "<<CoincidentNode(inod,j);
//    }
//    cout<<endl;
//  }
//
//  exit(1);
// debug

  return;

}

// member function to swap pointers to start-of-step values

void Mesh::end_of_step(){

  setVolume0(&mVolume);
  setVelocity0(&mVelocity);
  setEnergy0(&mEnergy);
  setCCEnergy0(&mCCEnergy);

//  mKCoord0.swap(mKCoord); // no member function to do this
  mKCoord0.assign(mKCoord.begin(),mKCoord.end()); // no member function to do this

//  mTCoord0.swap(mTCoord); // no member function to do this
  mTCoord0.assign(mTCoord.begin(),mTCoord.end()); // no member function to do this

  return;

}

// member function to initialise the pressure field for the start of the calculation

void Mesh::pressure_init(){

// reset pressure field in physical cells

  for(long iel=0;iel<NCells();iel++){
    int mat=Material(iel);
    mPressure.push_back(PressureInit(mat));
  }

  cout<<"Mesh::pressure_init(): "<<mPressure.size()<<" pressures set..."<<endl;

  return;

}

// member function to initialise the density field for the start of the calculation

void Mesh::density_init(){

// reset density field in physical cells

  for(long iel=0;iel<NCells();iel++){
    int mat=Material(iel);
    mDensity.push_back(DensityInit(mat));
  }

  cout<<"Mesh::density_init(): "<<mDensity.size()<<" densities set..."<<endl;

  return;

}

// member function to initialise the energy field based on the initial pressure

void Mesh::energy_init(){

// get masses for starting a nodal energy field
// commented out as J() uses kinematic grid which gives wrong nodal volumes for the thermodynamics
//  Shape S1(Shape::QUAD2D,Order(0),Order(0)+2);
//  Shape*S[2]={&S1,&S1}; // nodal_mass calls J which requires 2nd address in shape array
//  V1 nodmass; // nodal masses
//  nodal_mass(S,this,&nodmass); // masses should all be 0.015625000 for p1 (vol=0.5,rho=0.125)
//  for(long iel=0;iel<NCells();iel++){
//    int mat=Material(iel);
//    double e(PressureInit(mat)/(0.4*DensityInit(mat))); // invert the eos to get e from p
//    for(int iloc=0;iloc<TLoc(iel);iloc++){
//      double e(PressureInit(mat)/(0.4*DensityInit(mat))); // invert the eos to get e from p
//      mEnergy0.push_back(e);
//      mEnergy.push_back(e);
//    }
//  }

// reset energy field in physical cells

  for(long iel=0;iel<NCells();iel++){
    int mat=Material(iel);
    double e(PressureInit(mat)/(0.4*DensityInit(mat))); // invert the eos to get e from p
    for(int iloc=0;iloc<TLoc(iel);iloc++){
      mEnergy0.push_back(e);
      mEnergy.push_back(e);
    }
    mCCEnergy0.push_back(e);
    mCCEnergy.push_back(e);
  }

  cout<<"Mesh::energy_init(): "<<mEnergy.size()<<" internal energies set..."<<endl;

  return;

}

// member function to initialise the velocity field for the start of the calculation

void Mesh::velocity_init(){

// vector temporary with dimension in coordinate number

  vector<double> vtmp[3];

// reset velocty field in physical cells

  for(long iel=0;iel<NCells();iel++){
    int mat=Material(iel);
    for(int iloc=0;iloc<KLoc(iel);iloc++){
      vtmp[0].push_back(VelocityXInit(mat));
      vtmp[1].push_back(VelocityYInit(mat));
      vtmp[2].push_back(VelocityZInit(mat));
    }
  }

  cout<<"Mesh::velocity_init(): "<<vtmp[0].size()<<" velocities set..."<<endl;

// reset velocity field in ghost cells, use the data of the parent cell

  for(long iel=0;iel<NCells();iel++){
    int mat=Material(iel);
    for(int iface=0;iface<NFaces(iel);iface++){
      long ieln=E2E(iel,iface);
      if(ieln>NCells()-1){
//        for(int iloc=0;iloc<KLoc(ieln);iloc++){ // BUG 22/02/21: mKLoc[] not set in ghosts cells !!
        for(int iloc=0;iloc<KLoc(iel);iloc++){
          vtmp[0].push_back(VelocityXInit(mat));
          vtmp[1].push_back(VelocityYInit(mat));
          vtmp[2].push_back(VelocityZInit(mat));
        }
      }
    }
  }

// construct each velocity component from the temporary

  for(int idim=0;idim<NDims();idim++){
    mVelocity0.push_back(vtmp[idim]);
    mVelocity.push_back(vtmp[idim]);
  }

  cout<<"Mesh::velocity_init(): "<<mVelocity0[0].size()<<" velocities set including level 1 halo."<<endl;

  return;

}

// member function to impose boundary constraint on the initial velocity field,
// this is to prevent nodes on the mesh boundary moving perpendicular to the mesh.

void Mesh::velocity_bcs(){

// traverse edges of the mesh

  for(int iface=0;iface<NFaces(0);iface++){
    if(BoundaryType(iface).compare("REFLECTIVE")==0){
      cout<<"Mesh::velocity_bcs(): Mesh facet "<<iface<<" reflects, constraint imposed on perpendicular component of velocity field."<<endl;
      switch(iface){
        case 0:
          for(long ibnod=0;ibnod<NBoundaryNodes(iface);ibnod++){
            mVelocity0[1].at(BoundaryNode(iface,ibnod))=0.0;
            mVelocity[1].at(BoundaryNode(iface,ibnod))=0.0;
          }
          break;
        case 1:
          for(long ibnod=0;ibnod<NBoundaryNodes(iface);ibnod++){
            mVelocity0[0].at(BoundaryNode(iface,ibnod))=0.0;
            mVelocity[0].at(BoundaryNode(iface,ibnod))=0.0;
          }
          break;
        case 2:
          for(long ibnod=0;ibnod<NBoundaryNodes(iface);ibnod++){
            mVelocity0[1].at(BoundaryNode(iface,ibnod))=0.0;
            mVelocity[1].at(BoundaryNode(iface,ibnod))=0.0;
          }
          break;
        case 3:
          for(long ibnod=0;ibnod<NBoundaryNodes(iface);ibnod++){
            mVelocity0[0].at(BoundaryNode(iface,ibnod))=0.0;
            mVelocity[0].at(BoundaryNode(iface,ibnod))=0.0;
          }
          break;
      }
    }else{
      cout<<"Mesh::velocity_bcs(): Mesh facet "<<iface<<" transmits, no constraint imposed."<<endl;
    }
  }

  return;

}

// member function to update halo pressure, mass, density, energy, volume and velocity fields

void Mesh::update_halo(){

//  cout<<"Mesh::update_halo(): updating halo pressure,mass,density,energy, volume and velocity fields..."<<endl;

// buffer to store data to export to the halo address space

  vector<double> buffer[15]; // 15 member fields to update

// sweep physical elements and their faces to reach elements in the halo

  for(int iel=0;iel<NCells();iel++){
    for(int iface=0;iface<NFaces(iel);iface++){
      long ieln=E2E(iel,iface); // neighbour element on face iface

      if(ieln>NCells()-1){

// ieln is in the halo, action depends on boundary condition imposed on iface

        if(BoundaryType(iface).compare("TRANSMISSIVE")==0){

// impose transmissive constraint, allow mesh to expand - commit data to halo buffer

          buffer[0].push_back(Volume0(iel));
          buffer[1].push_back(Volume(iel));
//          buffer[2].push_back(0.0); // vacuum breaks the pvrs Riemann solver
          buffer[2].push_back(Pressure(iel));
//          buffer[3].push_back(0.0); // void breaks the pvrs Riemann solver
          buffer[3].push_back(Mass(iel));
//          buffer[4].push_back(0.0); // void breaks the pvrs Riemann solver
          buffer[4].push_back(Density(iel));

          for(int iloc=0;iloc<TLoc(iel);iloc++){
//            buffer[5].push_back(0.0);
            buffer[5].push_back(Energy0(TNode(iel,iloc)));
//            buffer[6].push_back(0.0);
            buffer[6].push_back(Energy(TNode(iel,iloc)));
          }

//          buffer[7].push_back(0.0);
          buffer[7].push_back(CCEnergy0(iel));
//          buffer[8].push_back(0.0);
          buffer[8].push_back(CCEnergy(iel));

          for(int iloc=0;iloc<KLoc(iel);iloc++){
//            buffer[9].push_back(0.0);
            buffer[9].push_back(Velocity0(0,KNode(iel,iloc)));
//            buffer[10].push_back(0.0);
            buffer[10].push_back(Velocity(0,KNode(iel,iloc)));
//            buffer[11].push_back(0.0);
            buffer[11].push_back(Velocity0(1,KNode(iel,iloc)));
//            buffer[12].push_back(0.0);
            buffer[12].push_back(Velocity(1,KNode(iel,iloc)));
//            buffer[13].push_back(0.0);
            buffer[13].push_back(Velocity0(2,KNode(iel,iloc)));
//            buffer[14].push_back(0.0);
            buffer[14].push_back(Velocity(2,KNode(iel,iloc)));
          }

        }else{

// impose reflective constraint, enforce zero gradient across iface - commit data to halo buffer

          buffer[0].push_back(Volume0(iel));
          buffer[1].push_back(Volume(iel));
          buffer[2].push_back(Pressure(iel));
          buffer[3].push_back(Mass(iel));
          buffer[4].push_back(Density(iel));

          for(int iloc=0;iloc<TLoc(iel);iloc++){
            buffer[5].push_back(Energy0(TNode(iel,Reflect(iloc,iface,TOrder(iel)))));
            buffer[6].push_back(Energy(TNode(iel,Reflect(iloc,iface,TOrder(iel)))));
          }

          buffer[7].push_back(CCEnergy0(iel));
          buffer[8].push_back(CCEnergy(iel));

          for(int iloc=0;iloc<KLoc(iel);iloc++){
            buffer[9].push_back(Velocity0(0,KNode(iel,Reflect(iloc,iface,KOrder(iel)))));
            buffer[10].push_back(Velocity(0,KNode(iel,Reflect(iloc,iface,KOrder(iel)))));
            buffer[11].push_back(Velocity0(1,KNode(iel,Reflect(iloc,iface,KOrder(iel)))));
            buffer[12].push_back(Velocity(1,KNode(iel,Reflect(iloc,iface,KOrder(iel)))));
            buffer[13].push_back(Velocity0(2,KNode(iel,Reflect(iloc,iface,KOrder(iel)))));
            buffer[14].push_back(Velocity(2,KNode(iel,Reflect(iloc,iface,KOrder(iel)))));
          }

        }

      }

    }
  }

// replace class member data across the halo address space with the data in the buffer
//
//  copy_n(buffer[0].begin(),buffer[0].size(),&mVolume0[NCells()]);
//  copy_n(buffer[1].begin(),buffer[1].size(),&mVolume[NCells()]);
//  copy_n(buffer[2].begin(),buffer[2].size(),&mPressure[NCells()]);
//  copy_n(buffer[3].begin(),buffer[3].size(),&mMass[NCells()]);
//  copy_n(buffer[4].begin(),buffer[4].size(),&mDensity[NCells()]);
//  copy_n(buffer[5].begin(),buffer[5].size(),&mEnergy0[NCells()*TLoc(0)]);
//  copy_n(buffer[6].begin(),buffer[6].size(),&mEnergy[NCells()*TLoc(0)]);
//  copy_n(buffer[7].begin(),buffer[7].size(),&mCCEnergy0[NCells()]);
//  copy_n(buffer[8].begin(),buffer[8].size(),&mCCEnergy[NCells()]);
//  copy_n(buffer[9].begin(),buffer[9].size(),&mVelocity0[0][NCells()*KLoc(0)]);
//  copy_n(buffer[10].begin(),buffer[10].size(),&mVelocity[0][NCells()*KLoc(0)]);
//  copy_n(buffer[11].begin(),buffer[11].size(),&mVelocity0[1][NCells()*KLoc(0)]);
//  copy_n(buffer[12].begin(),buffer[12].size(),&mVelocity[1][NCells()*KLoc(0)]);
//  copy_n(buffer[13].begin(),buffer[13].size(),&mVelocity0[2][NCells()*KLoc(0)]);
//  copy_n(buffer[14].begin(),buffer[14].size(),&mVelocity[2][NCells()*KLoc(0)]);

// insert halo buffer into class member address space

  mVolume0.insert(end(mVolume0),begin(buffer[0]),end(buffer[0]));
  mVolume.insert(end(mVolume),begin(buffer[1]),end(buffer[1]));
  mPressure.insert(end(mPressure),begin(buffer[2]),end(buffer[2]));
  mMass.insert(end(mMass),begin(buffer[3]),end(buffer[3]));
  mDensity.insert(end(mDensity),begin(buffer[4]),end(buffer[4]));
  mEnergy0.insert(end(mEnergy0),begin(buffer[5]),end(buffer[5]));
  mEnergy.insert(end(mEnergy),begin(buffer[6]),end(buffer[6]));
  mCCEnergy.insert(end(mCCEnergy),begin(buffer[7]),end(buffer[7]));
  mCCEnergy0.insert(end(mCCEnergy0),begin(buffer[8]),end(buffer[8]));

  mVelocity0[0].insert(end(mVelocity0[0]),begin(buffer[9]),end(buffer[9]));
  mVelocity[0].insert(end(mVelocity[0]),begin(buffer[10]),end(buffer[10]));
  mVelocity0[1].insert(end(mVelocity0[1]),begin(buffer[11]),end(buffer[11]));
  mVelocity[1].insert(end(mVelocity[1]),begin(buffer[12]),end(buffer[12]));
  mVelocity0[2].insert(end(mVelocity0[2]),begin(buffer[13]),end(buffer[13]));
  mVelocity[2].insert(end(mVelocity[2]),begin(buffer[14]),end(buffer[14]));

  return;

}

// member function to create dummy mesh for silo
// this inserts dummy elements of zero volume inbetween each finite element
// and splits each physical element up through the nodes so a zonelist can be generated
// that is plottable in silo

void Mesh::insert_dummies(){

// this has only been done for 2D, stop if input mesh is 3D

  if(Dim(NDims()-1)>1){
    cout<<"Mesh::insert_dummies():  Dummy mesh not set up yet !"<<endl;
    exit(1);
  }

  cout<<"Mesh::insert_dummies():  Generating dummy mesh for graphics output..."<<endl;

// add a dummy cells inbetween each element

  int ncellsx=2*(Dim(0)-1)-1;
  int ncellsy=2*(Dim(1)-1)-1;
  long iel=0;
//  vector<vector<long> > zonelist[4]; // 0=phys, 1=right dummy, 2=top dummy, 3=corner dummy

// visit all cells including dummies

  for(int iely=0;iely<ncellsy;iely++){

    int ieldx=0;

    for(int ielx=0;ielx<ncellsx;ielx++){

// determine if cell ielx,iely is a physical cell or a dummy

      if(ielx%2==0 && iely%2==0){

// physical cell - sweep over subcells and draw the zonelist

        int p[3]={Order(iel),TOrder(iel),KOrder(iel)}; // polyhedral order of physical element iel

        vector<long> ielzonelist; // nodelist for current element

// add the global node of each corner of each subcell

        for(int j=0;j<p[0];j++){
          for(int i=0;i<p[0];i++){

// nodes at each end of each face of the subcell

//            ielzonelist.push_back(Node(iel,j*(p[0]+1)+i));       // bottom left corner of subcell 
//            ielzonelist.push_back(Node(iel,j*(p[0]+1)+i+1));     // bottom right corner of subcell 
//            ielzonelist.push_back(Node(iel,j*(p[0]+1)+i+1));     // bottom right corner of subcell 
//            ielzonelist.push_back(Node(iel,(j+1)*(p[0]+1)+i+1)); // top right corner of subcell 
//            ielzonelist.push_back(Node(iel,(j+1)*(p[0]+1)+i+1)); // top right corner of subcell
//            ielzonelist.push_back(Node(iel,(j+1)*(p[0]+1)+i));   // top left corner of subcell
//            ielzonelist.push_back(Node(iel,(j+1)*(p[0]+1)+i));   // top left corner of subcell
//            ielzonelist.push_back(Node(iel,j*(p[0]+1)+i));       // bottom left corner of subcell

// nodes at each corner of the subcell

            ielzonelist.push_back(Node(iel,j*(p[0]+1)+i));       // bottom left corner of subcell
            ielzonelist.push_back(Node(iel,j*(p[0]+1)+i+1));     // bottom right corner of subcell
            ielzonelist.push_back(Node(iel,(j+1)*(p[0]+1)+i+1)); // top right corner of subcell 
            ielzonelist.push_back(Node(iel,(j+1)*(p[0]+1)+i));   // top left corner of subcell

// append to the nodelist

            mNodeList.push_back(Node(iel,j*(p[0]+1)+i));       // bottom left corner of subcell
            mNodeList.push_back(Node(iel,j*(p[0]+1)+i+1));     // bottom right corner of subcell
            mNodeList.push_back(Node(iel,(j+1)*(p[0]+1)+i+1)); // top right corner of subcell
            mNodeList.push_back(Node(iel,(j+1)*(p[0]+1)+i));   // top left corner of subcell

          }
        }

// append the zonelists for all subcells

//        zonelist[0].push_back(ielzonelist);

// add the kinematic node of each corner of each subcell

        for(int j=0;j<p[2];j++){
          for(int i=0;i<p[2];i++){

// append to the nodelist

            mKNodeList.push_back(KNode(iel,j*(p[2]+1)+i));       // bottom left corner of subcell
            mKNodeList.push_back(KNode(iel,j*(p[2]+1)+i+1));     // bottom right corner of subcell
            mKNodeList.push_back(KNode(iel,(j+1)*(p[2]+1)+i+1)); // top right corner of subcell
            mKNodeList.push_back(KNode(iel,(j+1)*(p[2]+1)+i));   // top left corner of subcell

// append element number in the subcell

            mElement.push_back(iel);

          }
        }

// add the thermodynamic node of each corner of each subcell

        for(int j=0;j<p[1];j++){
          for(int i=0;i<p[1];i++){

// append to the nodelist

            mTNodeList.push_back(TNode(iel,j*(p[1]+1)+i));       // bottom left corner of subcell
            mTNodeList.push_back(TNode(iel,j*(p[1]+1)+i+1));     // bottom right corner of subcell
            mTNodeList.push_back(TNode(iel,(j+1)*(p[1]+1)+i+1)); // top right corner of subcell
            mTNodeList.push_back(TNode(iel,(j+1)*(p[1]+1)+i));   // top left corner of subcell

          }
        }

        iel++;

      }else if(ielx%2!=0 && iely%2==0){

// dummy cell to the right of a physical cell - sweep over subcells and draw the zonelist
 
        int iell=iel-1;           // physical element to the left (iel just visited)
        int ielr=iell+1;          // physical element to the left (iel just visited)
        int p[3]={Order(iell),TOrder(iell),KOrder(iell)}; // polyhedral order of physical element iell
        vector<long> ielzonelist; // nodelist for current element

// add the global node of each corner of each subcell

        for(int j=0;j<p[0];j++){

// nodes at each end of each face of the subcell

//          ielzonelist.push_back(Node(iell,j*(p[0]+1)+p[0]));
//          ielzonelist.push_back(Node(ielr,j*(p[0]+1)));
//          ielzonelist.push_back(Node(ielr,j*(p[0]+1)));
//          ielzonelist.push_back(Node(ielr,(j+1)*(p[0]+1)));
//          ielzonelist.push_back(Node(ielr,(j+1)*(p[0]+1)));
//          ielzonelist.push_back(Node(iell,(j+1)*(p[0]+1)+p[0]));
//          ielzonelist.push_back(Node(iell,(j+1)*(p[0]+1)+p[0]));
//          ielzonelist.push_back(Node(iell,j*(p[0]+1)+p[0]));

// nodes at each corner of the subcell

          ielzonelist.push_back(Node(iell,j*(p[0]+1)+p[0]));
          ielzonelist.push_back(Node(ielr,j*(p[0]+1)));
          ielzonelist.push_back(Node(ielr,(j+1)*(p[0]+1)));
          ielzonelist.push_back(Node(iell,(j+1)*(p[0]+1)+p[0]));

// append to the nodelist

          mNodeList.push_back(Node(iell,j*(p[0]+1)+p[0]));
          mNodeList.push_back(Node(ielr,j*(p[0]+1)));
          mNodeList.push_back(Node(ielr,(j+1)*(p[0]+1)));
          mNodeList.push_back(Node(iell,(j+1)*(p[0]+1)+p[0]));

        }

// append the zonelists for all subcells

//        zonelist[1].push_back(ielzonelist);

// add the kinematic node of each corner of each subcell

        for(int j=0;j<p[2];j++){

// append to the nodelist

          mKNodeList.push_back(KNode(iell,j*(p[2]+1)+p[2]));
          mKNodeList.push_back(KNode(ielr,j*(p[2]+1)));
          mKNodeList.push_back(KNode(ielr,(j+1)*(p[2]+1)));
          mKNodeList.push_back(KNode(iell,(j+1)*(p[2]+1)+p[2]));

// append negative element number to flag dummy

          mElement.push_back(-iel);

        }

// add the thermodynamic node of each corner of each subcell

        for(int j=0;j<p[1];j++){

// append to the nodelist

          mTNodeList.push_back(TNode(iell,j*(p[1]+1)+p[1]));
          mTNodeList.push_back(TNode(ielr,j*(p[1]+1)));
          mTNodeList.push_back(TNode(ielr,(j+1)*(p[1]+1)));
          mTNodeList.push_back(TNode(iell,(j+1)*(p[1]+1)+p[1]));

        }

      }else if(ielx%2==0 && iely%2!=0){

// dummy cell to the top of a physical cell - sweep over subcells and draw the zonelist
 
        int ielt=iel+ieldx;       // physical element above the dummy
        int ielb=ielt-(Dim(0)-1); // physical element below the dummy
        int p[3]={Order(ielb),TOrder(ielb),KOrder(ielb)}; // polyhedral order of physical element ielb
        vector<long> ielzonelist; // nodelist for current element

// add the global node of each corner of each subcell

        for(int i=0;i<p[0];i++){

// nodes at each end of each face of the subcell

//          ielzonelist.push_back(Node(ielb,p[0]*(p[0]+1)+i));
//          ielzonelist.push_back(Node(ielb,p[0]*(p[0]+1)+i+1));
//          ielzonelist.push_back(Node(ielb,p[0]*(p[0]+1)+i+1));
//          ielzonelist.push_back(Node(ielt,i+1));
//          ielzonelist.push_back(Node(ielt,i+1));
//          ielzonelist.push_back(Node(ielt,i));
//          ielzonelist.push_back(Node(ielt,i));
//          ielzonelist.push_back(Node(ielb,p[0]*(p[0]+1)+i));

// nodes at each corner of the subcell

          ielzonelist.push_back(Node(ielb,p[0]*(p[0]+1)+i));
          ielzonelist.push_back(Node(ielb,p[0]*(p[0]+1)+i+1));
          ielzonelist.push_back(Node(ielt,i+1));
          ielzonelist.push_back(Node(ielt,i));

// append to the nodelist

          mNodeList.push_back(Node(ielb,p[0]*(p[0]+1)+i));
          mNodeList.push_back(Node(ielb,p[0]*(p[0]+1)+i+1));
          mNodeList.push_back(Node(ielt,i+1));
          mNodeList.push_back(Node(ielt,i));

        }

// append the zonelists for all subcells

//        zonelist[2].push_back(ielzonelist);

// add the kinematic node of each corner of each subcell

        for(int i=0;i<p[2];i++){

// append to the nodelist

          mKNodeList.push_back(KNode(ielb,p[2]*(p[2]+1)+i));
          mKNodeList.push_back(KNode(ielb,p[2]*(p[2]+1)+i+1));
          mKNodeList.push_back(KNode(ielt,i+1));
          mKNodeList.push_back(KNode(ielt,i));

// append negative element number to flag dummy

          mElement.push_back(-iel);

        }

// add the thermodynamic node of each corner of each subcell

        for(int i=0;i<p[1];i++){

// append to the nodelist

          mTNodeList.push_back(TNode(ielb,p[1]*(p[1]+1)+i));
          mTNodeList.push_back(TNode(ielb,p[1]*(p[1]+1)+i+1));
          mTNodeList.push_back(TNode(ielt,i+1));
          mTNodeList.push_back(TNode(ielt,i));

        }

        ieldx++;

      }else if(ielx%2!=0 && iely%2!=0){

// dummy cell on the corner of 4 physical cells

        int ieltr=iel+ieldx;         // physical element to the top right
        int ieltl=ieltr-1;           // physical element to the top left
        int ielbl=ieltl-(Dim(0)-1); // physical element to the bottom left
        int ielbr=ieltr-(Dim(0)-1); // physical element to the bottom right
        int p[3]={Order(ielbl),TOrder(ielbl),KOrder(ielbl)}; // polyhedral order of physical element ielbl
        vector<long> ielzonelist; // nodelist for current element

// add the global node of each corner

// nodes at each end of each face of the subcell

//        ielzonelist.push_back(Node(ielbl,(p[0]+1)*(p[0]+1)-1));
//        ielzonelist.push_back(Node(ielbr,(p[0]+1)*(p[0]+1)-(p[0]+1)));
//        ielzonelist.push_back(Node(ielbr,(p[0]+1)*(p[0]+1)-(p[0]+1)));
//        ielzonelist.push_back(Node(ieltr,0));
//        ielzonelist.push_back(Node(ieltr,0));
//        ielzonelist.push_back(Node(ieltl,p[0]));
//        ielzonelist.push_back(Node(ieltl,p[0]));
//        ielzonelist.push_back(Node(ielbl,(p[0]+1)*(p[0]+1)-1));

// nodes at each corner of the subcell

        ielzonelist.push_back(Node(ielbl,(p[0]+1)*(p[0]+1)-1));
        ielzonelist.push_back(Node(ielbr,(p[0]+1)*(p[0]+1)-(p[0]+1)));
        ielzonelist.push_back(Node(ieltr,0));
        ielzonelist.push_back(Node(ieltl,p[0]));


// append to the nodelist

        mNodeList.push_back(Node(ielbl,(p[0]+1)*(p[0]+1)-1));
        mNodeList.push_back(Node(ielbr,(p[0]+1)*(p[0]+1)-(p[0]+1)));
        mNodeList.push_back(Node(ieltr,0));
        mNodeList.push_back(Node(ieltl,p[0]));

// append the zonelists for all subcells

//        zonelist[3].push_back(ielzonelist);

// add the kinematic node of each corner

        mKNodeList.push_back(KNode(ielbl,(p[2]+1)*(p[2]+1)-1));
        mKNodeList.push_back(KNode(ielbr,(p[2]+1)*(p[2]+1)-(p[2]+1)));
        mKNodeList.push_back(KNode(ieltr,0));
        mKNodeList.push_back(KNode(ieltl,p[2]));

// append negative element number to flag dummy

        mElement.push_back(-iel);

// add the thermodynamic node of each corner

        mTNodeList.push_back(TNode(ielbl,(p[1]+1)*(p[1]+1)-1));
        mTNodeList.push_back(TNode(ielbr,(p[1]+1)*(p[1]+1)-(p[1]+1)));
        mTNodeList.push_back(TNode(ieltr,0));
        mTNodeList.push_back(TNode(ieltl,p[1]));

      }

    }
  }

  cout<<"Mesh::insert_dummies():  done."<<endl;

  return;

}

// accessor function to return number of mesh dimensions

int Mesh::NDims(){return mNDims;}

// accessor function to return the length of mesh dimension idim

long Mesh::Dim(int idim){return mDim[idim];}

// accessor function to return number of mesh vertices

long Mesh::NVertices(){return mNVertices;}

// accessor function to return number of cells on the mesh

long Mesh::NCells(){return mNCells;}

// accessor function to return number of ghosts on the mesh

long Mesh::NGhosts(){return mNGhosts;}

// accessor function to return number of materials

int Mesh::NMaterials(){return mNMaterials;}

// accessor function to return material number imat

int Mesh::Materials(int imat){return mMaterials[imat];}

// accessor function to return the number of polyhedral types

int Mesh::NTypes(){return mNTypes;}

// accessor function to return polyhedral type itype

int Mesh::Types(int itype){return mTypes[itype];}

// accessor function to return the name of polyhedral type itype

string Mesh::Names(int itype){return mNames[itype];}

// accessor function to return the number of polyhedral orders

int Mesh::NOrders(){return mNOrders;}

// accessor function to return polyhedral order iorder

int Mesh::Orders(int iorder){return mOrders[iorder];}

// accessor function to return polyhedral order iorder for kinematics

int Mesh::KOrders(int iorder){return mKOrders[iorder];}

// accessor function to return polyhedral order iorder for thermodynamics

int Mesh::TOrders(int iorder){return mTOrders[iorder];}

// accessor function to return polyhedral order of element iel

int Mesh::Order(long iel){return mOrder[iel];}

// accessor function to return polyhedral order of element iel to use for kinematics

int Mesh::KOrder(long iel){return mKOrder[iel];}

// accessor function to return polyhedral order of element iel to use for thermodynamics

int Mesh::TOrder(long iel){return mTOrder[iel];}

// accessor function to return polyhedral type of element iel

int Mesh::Type(long iel){return mType[iel];}

// accessor function to return material number in element iel

int Mesh::Material(long iel){return mMaterial[iel];}

// accessor function to return the element corner sum

long Mesh::NCorners(){return mNCorners;}

// accessor function to return the vertex in corner icorner of element iel

long Mesh::Vertex(long iel,int icorner){return mVertex[iel*NLoc(iel)+icorner];}

// accessor function to return coordinate idim of mesh vertex ivert

double Mesh::Coord(int idim,long ivert){return mCoord[idim][ivert];}

// accessor function to return coordinate idim of global node inod

double Mesh::DGCoord(int idim,long inod){return mDGCoord[idim][inod];}

// accessor function to return coordinate idim of kinematic node inod

double Mesh::KCoord0(int idim,long inod){return mKCoord0[idim][inod];}

// accessor function to return half-step/full-step coordinate idim of kinematic node inod

double Mesh::KCoord(int idim,long inod){return mKCoord[idim][inod];}

// accessor function to return coordinate idim of thermodynamic node inod

double Mesh::TCoord0(int idim,long inod){return mTCoord0[idim][inod];}

// accessor function to return half-step/full-step coordinate idim of thermodynamic node inod

double Mesh::TCoord(int idim,long inod){return mTCoord[idim][inod];}

// accessor function to return the number of local nodes in meshfile element iel

int Mesh::NLoc(long iel){return mNLoc[iel];}

// accessor function to return the number of local finite element nodes inserted in to mesh ile element iel

int Mesh::MLoc(long iel){return mMLoc[iel];}

// accessor function to return the number of local kinematic nodes inserted in to mesh ile element iel

int Mesh::KLoc(long iel){return mKLoc[iel];}

// accessor function to return the number of local thermodynamic nodes inserted in to mesh ile element iel

int Mesh::TLoc(long iel){return mTLoc[iel];}

// accessor function to return the number of global nodes in the mesh

long Mesh::NNodes(){return mNNodes;}

// accessor function to return the number of kinematic nodes in the mesh

long Mesh::KNodes(){return mKNodes;}

// accessor function to return the number of thermodynamic nodes in the mesh

long Mesh::TNodes(){return mTNodes;}

// accessor function to return the density field in element iel

double Mesh::Density(long iel){return mDensity[iel];}

// accessor function to return the pressure field in element iel

double Mesh::Pressure(long iel){return mPressure[iel];}

// accessor function to return component idim of the start-of-step velocity field at global node inod

double Mesh::Velocity0(int idim,long inod){return mVelocity0[idim][inod];}

// accessor function to return component idim of the velocity field at global node inod

double Mesh::Velocity(int idim,long inod){return mVelocity[idim][inod];}

// accessor function to return component idim of the mass-weighted velocity field at global node inod

double Mesh::MWVelocity(int idim,long inod){return mMWVelocity[idim][inod];}

// accessor function to return the start-of-step energy field at global node inod

double Mesh::Energy0(long inod){return mEnergy0[inod];}

// accessor function to return the half-step energy field at global node inod

double Mesh::Energy(long inod){return mEnergy[inod];}

// accessor function to return the start-of-step cell-centred energy field in element iel

double Mesh::CCEnergy0(long iel){return mCCEnergy0[iel];}

// accessor function to return the half-step cell-centred energy field in element iel

double Mesh::CCEnergy(long iel){return mCCEnergy[iel];}

// accessor function to return the global node number of local node iloc in element iel

long Mesh::Node(long iel,int iloc){return mNode[iel][iloc];}

// accessor function to return the kinematic node number of local node iloc in element iel

long Mesh::KNode(long iel,int iloc){return mKNode[iel][iloc];}

// accessor function to return the thermodynamic node number of local node iloc in element iel

long Mesh::TNode(long iel,int iloc){return mTNode[iel][iloc];}

// accessor function to return the element to element connectivity on face iface of element iel

long Mesh::E2E(long iel,int iface){return mE2E[iel][iface];}

// accessor function to return the node to node connectivity on node iloc of face iface of element iel

int Mesh::N2N(long iel,int iface,int iloc){return mN2N[iel][iface][iloc];}

// accessor function to return the face to face connectivity on face iface of element iel

int Mesh::F2F(long iel,int iface){return mF2F[iel][iface];}

// accessor function to return the number of faces on the surface of element iel

int Mesh::NFaces(long iel){return mNFaces[iel];}

// accessor functin to return the number of local nodes on face iface of element iel

int Mesh::NFLoc(long iel,int iface){return mNFLoc[iel][iface];}

// accessor function to return the local node number of node ifloc on face iface of thermodynamic element iel

int Mesh::FLoc(long iel,int iface,int ifloc){return mFLoc[iel][iface][ifloc];}

// accessor function to return the local node number of node ifloc on face iface of kinematic element iel

int Mesh::FLock(long iel,int iface,int ifloc){return mFLock[iel][iface][ifloc];}

// accessor function to return element ilist on the silo nodelist

long Mesh::NodeList(long ilist){return mNodeList[ilist];}

// accessor function to return element ilist on the silo kinematic nodelist

long Mesh::KNodeList(long ilist){return mKNodeList[ilist];}

// accessor function to return element ilist on the silo thermodynamic nodelist

long Mesh::TNodeList(long ilist){return mTNodeList[ilist];}

// accessor function to return length of the nodelist

long Mesh::lNodeList(){return mNodeList.size();}

// accessor function to return length of the kinematic nodelist

long Mesh::lKNodeList(){return mKNodeList.size();}

// accessor function to return length of the thermodynamic nodelist

long Mesh::lTNodeList(){return mTNodeList.size();}

// member function to return the initial energy field in material imat

double Mesh::EnergyInit(int imat){return mEnergyInit[imat-1];} // offset by 1 as first address in vector is 0

// member function to return the initial pressure field in material imat

double Mesh::PressureInit(int imat){return mPressureInit[imat-1];} // offset by 1 as first address in vector is 0

// member function to return the initial density field in material imat

double Mesh::DensityInit(int imat){return mDensityInit[imat-1];} // offset by 1 as first address in vector is 0

// member function to return the initial velocity in the x-direction in material imat

double Mesh::VelocityXInit(int imat){return mVelocityXInit[imat-1];} // offset by 1 as first address in vector is 0

// member function to return the initial velocity in the y-direction in material imat

double Mesh::VelocityYInit(int imat){return mVelocityYInit[imat-1];} // offset by 1 as first address in vector is 0

// member function to return the initial velocity in the z-direction in material imat

double Mesh::VelocityZInit(int imat){return mVelocityZInit[imat-1];} // offset by 1 as first address in vector is 0

// member function to return the start-of-step volume of element iel

double Mesh::Volume0(long iel){return mVolume0[iel];}

// member function to return the volume of element iel

double Mesh::Volume(long iel){return mVolume[iel];}

// member function to return the mass of element iel

double Mesh::Mass(long iel){return mMass[iel];}

// member function to update the start-of-step volume field to the new volumes passed in

//void Mesh::setVolume0(vector<double>*new_volumes){mVolume0.swap((*new_volumes));}
void Mesh::setVolume0(vector<double>*new_volumes){mVolume0.assign((*new_volumes).begin(),(*new_volumes).end());}

// member function to update the volume field to the new volumes passed in

//void Mesh::setVolume(vector<double>*new_volumes){mVolume.swap((*new_volumes));}
void Mesh::setVolume(vector<double>*new_volumes){mVolume.assign((*new_volumes).begin(),(*new_volumes).end());}

// member function to update all element pressures to the new pressures passed in

void Mesh::setPressure(vector<double>*new_pressures){mPressure.swap((*new_pressures));}

// member function to update all element densities to the new densities passed in

//void Mesh::setDensity(vector<double>*new_densities){mDensity.swap((*new_densities));}
void Mesh::setDensity(vector<double>*new_densities){mDensity.assign((*new_densities).begin(),(*new_densities).end());}

// member function to update the mass field (this should only be called once as mass is constant across the entire calculation)

void Mesh::setMass(vector<double>*new_mass){mMass.swap((*new_mass));}

// member function to update start-of-step energy field to the new energies passed in

//void Mesh::setEnergy0(vector<double> *new_energies){mEnergy0.swap((*new_energies));}
void Mesh::setEnergy0(vector<double> *new_energies){mEnergy0.assign((*new_energies).begin(),(*new_energies).end());}

// member function to update half-step energy field to the new energies passed in

//void Mesh::setEnergy(vector<double> *new_energies){mEnergy.swap((*new_energies));}
void Mesh::setEnergy(vector<double> *new_energies){mEnergy.assign((*new_energies).begin(),(*new_energies).end());}

// member function to update start-of-step cell-centred energy field to the new energies passed in

//void Mesh::setCCEnergy0(vector<double> *new_energies){mCCEnergy0.swap((*new_energies));}
void Mesh::setCCEnergy0(vector<double> *new_energies){mCCEnergy0.assign((*new_energies).begin(),(*new_energies).end());}

// member function to update half-step cell-centred energy field to the new energies passed in

//void Mesh::setCCEnergy(vector<double> *new_energies){mCCEnergy.swap((*new_energies));}
void Mesh::setCCEnergy(vector<double> *new_energies){mCCEnergy.assign((*new_energies).begin(),(*new_energies).end());}

// member function to move kinematic node inod in direction idim to new position pos

void Mesh::setKCoord0(int idim,long inod,double pos){mKCoord0[idim][inod]=pos;}

// member function to move half-step/full-step kinematic node inod in direction idim to new position pos

void Mesh::setKCoord(int idim,long inod,double pos){mKCoord[idim][inod]=pos;}

// member function to move thermodynamic node inod in direction idim to new position pos

void Mesh::setTCoord0(int idim,long inod,double pos){mTCoord0[idim][inod]=pos;}

// member function to move half-step/full-step thermodynamic node inod in direction idim to new position pos

void Mesh::setTCoord(int idim,long inod,double pos){mTCoord[idim][inod]=pos;}

// member function to update the start-of-step velocity field to the new vector passed in

//void Mesh::setVelocity0(vector<vector<double> >*new_velocities){mVelocity0.swap((*new_velocities));}
void Mesh::setVelocity0(vector<vector<double> >*new_velocities){mVelocity0.assign((*new_velocities).begin(),(*new_velocities).end());}

// member function to update the full-step velocity field to the new vector passed in

void Mesh::setVelocity(vector<vector<double> >*new_velocities){mVelocity.swap((*new_velocities));}

// member function to update the mass-weighted velocity field to the new vector passed in

void Mesh::setMWVelocity(vector<vector<double> >*new_velocities){mMWVelocity.assign((*new_velocities).begin(),(*new_velocities).end());}

// member function to return the size of the start-of-step velocity field in dimension idim

long Mesh::NVelocities0(int idim){return (mVelocity0.size()==0) ? 0 : mVelocity0[idim].size();}

// member function to return the size of the full-step velocity field in dimension idim

long Mesh::NVelocities(int idim){return (mVelocity.size()==0) ? 0 : mVelocity[idim].size();}

// member function to return number of  nodes found along mesh boundary iface

long Mesh::NBoundaryNodes(int iface){return mBoundaryNode[iface].size();}

// member function to return global number of node inod found along mesh boundary iface

long Mesh::BoundaryNode(int iface,long inod){return mBoundaryNode[iface][inod];}

// member function to return the number of nodes that are coincident with global node inod

int Mesh::NCoincidentNodes(long inod){return mCoincidentNode[inod].size();}

// member function to return global number of node i that is coincident with global node inod

long Mesh::CoincidentNode(long inod,int i){return mCoincidentNode[inod][i];}

// member function to return the boundary condition imposed on mesh face iface

string Mesh::BoundaryType(int iface){return mBoundaryType[iface];}

// accessor function to reflect local node iloc across face iface of an order p element

int Mesh::Reflect(int iloc,int iface,int p){

// extract row and column address of the local node

  int i(iloc%(p+1)),j(iloc/(p+1)),jloc(iloc);

// choose a face and reflect across it

  switch(iface){
    case 0:
      jloc=(p-j)*(p+1)+i;
      break;
    case 1:
      jloc=(j+1)*(p+1)-i-1;
      break;
    case 2:
      jloc=(p-j)*(p+1)+i;
      break;
    case 3:
      jloc=(j+1)*(p+1)-i-1;
      break;
  }

  return jloc;

}

// member function to switch between PdV and finite element energy fields

bool Mesh::UsePdV(){return mUsePdV;}
