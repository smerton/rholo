// Function definition to implement boundary conditions

// Author S. R. Merton

#include <iostream>
#include <vector>
#include "globals.h" // defines
#include "matrix.h"  // matrix class
#include "mesh.h"    // mesh class
#include "shape.h"   // shape class
#include "bcs.h"     // boundary condition signatures

using namespace std;

// insert boundary conditions in to the mass matrix via row elimination

void bc_insert(Matrix &A,Mesh const &M,Shape const &S,VD const &d,VD const &detJ,VVD &u0,VVD &u1,VD &b0,VD &b1,long const &nknodes){

// loop over boundary elements and choose what type of boundary needs to be applied

  cout<<"There are "<<M.NSides()<<" element sides on the mesh boundary:"<<endl;

// initialise boundary vectors

  fill(b0.begin(),b0.end(),0.0); // value on the boundary
  fill(b1.begin(),b1.end(),0.0); // eliminated row

// find cell sides coincident with the edges of the mesh

  for(long i=0;i<M.NCells();i++){
    for(int iside=0;iside<M.NVertices(i);iside++){

      if(M.E2E(i,iside)<M.NCells()){continue;}

// side iside is on a mesh edge, impose boundary conditions

      string bcname;

      int idim=(iside==0||iside==2)?1:0; // direction perpendicular to mesh boundary

// acquire local node numbers on side iside

      int local_node[S.order()+1];

      switch(iside){

        case(0): // bottom edge of mesh

          for(int iloc=0;iloc<S.order()+1;iloc++){local_node[iloc]=iloc;}

          break;

        case(1): // right edge of mesh

          for(int iloc=0;iloc<S.order()+1;iloc++){local_node[iloc]=iloc*(S.order()+1)+S.order();}

          break;

        case(2): // top edge of mesh

          for(int iloc=0;iloc<S.order()+1;iloc++){local_node[iloc]=S.nloc()-S.order()-1+iloc;}

          break;

        case(3): // left edge of mesh

          for(int iloc=0;iloc<S.order()+1;iloc++){local_node[iloc]=iloc*(S.order()+1);}

          break;

      }

// select boundary condition type on side iside

      switch(M.bc_edge(iside)){

        case(VACUUM):

// do nothing so mesh expands into the void

          bcname="vacuum";

          break;

        case(REFLECTIVE):

// set v.n=0 on domain boundary and impose a constraint on the acceleration field

          bcname="reflective";

          break;

        case(VELOCITY):

// set v.n=<value> on domain boundary and impose a constraint on the acceleration field

          bcname="velocity";

// set v.n=<value> on domain boundary

          for(int isloc=0;isloc<S.order()+1;isloc++){
            long k(M.GlobalNode_CFEM(i,local_node[isloc]));
            u0.at(idim).at(k)=M.bc_value(iside);
            u1.at(idim).at(k)=M.bc_value(iside);
          }

// eliminate k'th solution as we are imposing a condition on it

          for(int isloc=0;isloc<S.order()+1;isloc++){

            long k(M.GlobalNode_CFEM(i,local_node[isloc])); // boundary node
            double bval(1.0e-200); // boundary value

// collect known information - is this double counting ??

            for(long irow=0;irow<nknodes;irow++){b1.at(idim*nknodes+irow)+=A.read(idim*nknodes+irow,idim*nknodes+k)*bval;}

// store boundary value

            b0.at(idim*NROWS+k)=bval;

          }

// modify mass matrix to restore symmetry, this is to try and avoid costly changes to the solution strategy

          for(int isloc=0;isloc<S.order()+1;isloc++){
            long k(M.GlobalNode_CFEM(i,local_node[isloc])); // boundary node
            for(long l=0;l<M.NDims()*nknodes;l++){
              if(l!=k){
                A.write(l,idim*nknodes+k,0.0);
                A.write(idim*nknodes+k,l,0.0);
              }
            }
            A.write(idim*nknodes+k,idim*nknodes+k,1.0);
          }

          break;

        case(FLUID):

// add in additional boundary fluid masses to give zero pressure gradient across the edges of the domain

          bcname="fluid";

          break;

        case(ACCELERATION):

// impose a.n=<value> on domain boundary via row elimination of the mass matrix

          bcname="acceleration";

          break;

        default:

          bcname="undefined";

          cout<<"bc_insert(): "<<bcname<<" boundary conditions not coded, stopping."<<endl;

          exit(1);

      }

      cout<<"Edge "<<iside<<" boundary type : "<<bcname<<endl;
    }
  }

  return;

}
