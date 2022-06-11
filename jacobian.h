// Signature for the jacobian of the lagrangian motion

// Author S. R. Merton

void jacobian(int const &i,VVD const &x,Mesh const &M,Shape const &S,VD &detJ,VVVD &detDJ);            // calculate a jacobian and determinant
void jacobian(int const &i,VVD const &x0,VVD const &x,Mesh const &M,Shape const &S,VD &detJ,VVVD &Js); // calculate a jacobian for the Lagrangian motion
void jacobian(int const &i,VVD const &x,Mesh const &M,Shape const &S,VD &detJ);                        // calculate a jacobian at the local nodes
