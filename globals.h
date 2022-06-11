// Global parameters

// Author S. R. Merton

#define VD vector<double>              // vector of doubles
#define VVD vector<VD>                 // vector of VD
#define VVVD vector<VVD>               // vector of VVD
#define VI vector<int>                 // vector of ints

#define VTOL 1.0e-10                   // threshold for volume errors
#define ECUT 1.0e-8                    // cut-off on the energy field

#define NROWS nknodes                  // number of rows in the global matrix
#define NCOLS nknodes                  // number of columns in the global matrix
#define NGI S.ngi()*n                  // number of integration points on the mesh
#define GPNT i*S.ngi()+gi              // global address of integration point gi
#define ROW M.GlobalNode_CFEM(i,iloc)  // row address in global matrix
#define COL M.GlobalNode_CFEM(i,jloc)  // column address in global matrix

#define LS_PSEUDO_1D 1                 // use a pseudo-1D length scale, best for all 1D tests on 2D meshes
#define LS_LOCAL 2                     // use sqrt(volume(t))/p in the length scale definition
#define LS_AVERAGE 3                   // use sqrt(average volume)/p in the length scale definition where average volume is locally smooth
#define LS_DIRECTIONAL 4               // use a directional length scale definition for improved symmetry

#define INCLUDE_GHOSTS 1               // used in ghost updates to update ghost data as well as physical data
#define EXCLUDE_GHOSTS 0               // used in ghost updates to update physical data only
