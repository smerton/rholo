// Signature of functions to insert boundary conditions

// Author S. R. Merton

#define VACUUM 1              // vacuum boundary
#define REFLECTIVE 2          // reflective boundary
#define VELOCITY 3            // velocity v.n applied to boundary
#define FLUID 4               // fluid on the boundary
#define ACCELERATION 5        // velocity a.n applied to boundary

// insert boundary conditions in to the mass matrix via row elimination

void bc_insert(Matrix &A,Mesh const &M,Shape const &S,VD const &d,VD const &detJ,
               VVD &u0,VVD &u1,VD &b0,VD &b1,long const &nnodes);
