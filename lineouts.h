// Signatures for lineout functions

// Author S. R. Merton

// top level function to choose between 1D and 2D lines

void lineouts(Mesh const &M,Shape const &S,Shape const &T,VD const &dinit,VD const &e,VVD const &xinit,VVD const &x,
              VVD const &xt,VVD const &u,int const &test_problem,vector<int> const &mat,VD const &g);

// function to code for 1D lineouts

void lineouts_1d(Mesh const &M,Shape const &S,Shape const &T,VD const &dinit,VD const &e,VVD const &xinit,VVD const &x,
              VVD const &xt,VVD const &u,int const &test_problem,vector<int> const &mat,VD const &g);

// function to code for 2D lineouts

void lineouts_2d(Mesh const &M,Shape const &S,Shape const &T,VD const &dinit,VD const &e,VVD const &xinit,VVD const &x,
              VVD const &xt,VVD const &u,int const &test_problem,vector<int> const &mat,VD const &g);
