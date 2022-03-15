// Export the signature of the mesh class

#define VD vector<double>           // laziness
#define VVD vector<vector<double> > // laziness

#define INCLUDE_GHOSTS 1 // used in ghost updates to update ghost data as well as physical data
#define EXCLUDE_GHOSTS 0 // used in ghost updates to update physical data only

#define LS_PSEUDO_1D 1       // use a pseudo-1D length scale, best for all 1D tests on 2D meshes
#define LS_LOCAL 2           // use sqrt(volume(t))/p in the length scale definition
#define LS_AVERAGE 3         // use sqrt(average volume)/p in the length scale definition where average volume is locally smooth
#define LS_DIRECTIONAL 4     // use a directional length scale definition for improved symmetry

#include <vector>

using namespace std;

class Mesh{

  public:

  Mesh(char* meshfile);  // constructor for a new mesh
  ~Mesh();               // destructor function for storage release

// accessor functions to member data

  int NDims() const; // returns the number of mesh dimensions

  long NNodes() const; // returns the number of nodes on the mesh read in from the generator
  long NNodes(int p) const; // returns the number of nodes in a FEM mesh containing order p shapes
  long NNodes(int p,int t); // returns the number of nodes in a FEM mesh containing order p shapes of type t=CONTINUOUS or t=DISCONTINUOUS and sets up global node numbers
  long NNodes_CFEM() const; // returns the number of global nodes on the continuous finite element mesh
  long NNodes_DFEM() const; // returns the number of global nodes on the discontinuous finite element mesh
  int NGNodes() const; // returns the number of ghost nodes on the mesh edge 
  int NMaterials() const; // returns the number of materials on the mesh
  int NCells() const; // returns the number of cells on the mesh
  int NCells(int idim) const; // returns the number of cells on axis idim
  int NGCells() const; // returns the number of ghost cells on the mesh edge
  int Material(int i) const; // returns material number in cell i
  int Type(int i) const; // returns the geometric type of cell i
  int NVertices(int i) const; // returns the number of vertices of cell i
  int Vertex(int i,int j) const; // returns the node number of vertex j of cell i
  int NSides() const; // returns the number of cell sides coincident with a mesh edge
  int SideAttr(int i) const; // returns the attribute of side i on edge of mesh
  int SideType(int i) const; // returns the type of side i on edge of mesh
  int NSideNodes(int i) const; // returns the numnber of nodes on side i on edge of mesh
  int SideNode(int i,int j) const; // returns the node number of node j on side i on edge of mesh
  double Coord(int idim,int i) const; // returns coordinate idim of node i
  void InitCoords(VVD &v,int const p,int const t); // initialise the mesh coordinates onto the order p stencil of type t
  void InitLength(VD &l,int const &p,VD const &V,int const &ls_type) const; // initialise the element length scale
  double Volume(int i) const; // returns the volume of the element
  void bc_set(int iedge,int bc); // push new boundary condition bc on to mesh edge iedge
  void bc_set(int iedge,int bc,double bcvalue); // push new boundary condition bc and value bcvalue on to mesh edge iedge
  int bc_edge(int iedge) const; // returns the boundary condition on edge iedge of the mesh 
  double bc_value(int iedge) const; // returns the boundary value on edge iedge of the mesh 
  int nbcs() const; // returns the number of boundary conditions that have been set
  int E2E(int iel,int iface) const; // returns the element on face iface of element iel
  double Min(int idim) const; // returns mesh boundary minimum in dimension idim
  double Max(int idim) const; // returns mesh boundary maximum in dimension idim
  void UpdateCoords(VVD &x,VVD const &u,double const dt) const; // advect coordinate x a distance u*dt with velocity u
  void MapCoords(VVD const &xp,VVD &xq,int const &p,int const &q) const; // map coordinates from an order p mesh to an order q mesh
  double UpdateLength(int const &p,double const &V,double const &l0,double const &detJs,int const &length_scale_type); // update length scale for an element of polyhedral order p, volume V and initial length l0
  void UpdateVolume(VD &V,VVD const &x,int const &p) const; // update volume field V given coordinate x and polyhedral element order p
  void UpdateDensity(VD &d,VD const &V,VD const &m) const; // update denisty field d given a volume field V and a mass field m
  void UpdateEnergy(VD const &e0,VD &e1,VD const &p,VD const &q,VD const &V0,VD const &V1,VD const &m) const ; // update mesh energy field
  void UpdatePressure(VD &p,VD const &d,VD const &e,VD const &gamma,vector<int> const &mat); // load pressure field
  double UpdatePressure(double const &d,double const &e,double const&g); // returns pressure field value at a point
  void UpdateSoundSpeed(VD &c,VD const &g,vector<int> const &mat,VD const &p,VD const &d) const; // load new sound speeds
  double UpdateSoundSpeed(double const &g,double const &p,double const &d); // returns sound speed at a point
  double UpdateQ(double const&l,double const&d,double const&c,double const&cq,double const&cl,double const&divu); // returns artificial viscosity at a point
  long GlobalNode_CFEM(int const i,int const j) const; // global node number of local node j in element i in a continuous finite element method
  long GlobalNode_DFEM(int const i,int const j) const; // global node number of local node j in element i in a discontinuous finite element method

  private:

// member data

  int mNDims; // number of mesh dimensions
  long mNNodes; // number of nodes on the mesh read in from the generator
  long mNNodes_CFEM; // number of nodes on the mesh in a continuous finite element method after polyhedral order has been set (may differ from generator)
  long mNNodes_DFEM; // number of nodes on the mesh in a discontinuous finite element method after polyhedral order has been set (may differ from generator)
  long mNGNodes; // number of ghost nodes on the mesh edge
  int mNMaterials; // number of materials on the mesh
  int mNCells; // number of cells on the mesh
  int mNGCells; // number of ghost cells on the mesh edge
  int mNSides; // number of cell sides coinciding with the mesh edge
  vector<int> mMaterial; // material in each cell
  vector<int> mType; // polyhedral type of each cell
  vector<vector<int> > mVertex; // vertex number in each cell on the generator mesh
  vector<vector<long> > mGlobalNode_CFEM; // global node numbers in a continuous finite element method
  vector<vector<long> > mGlobalNode_DFEM; // global node numbers in a discontinuous finite element method
  vector<int> mSideAttr; // cell side boundary attribute
  vector<int> mSideType; // cell side boundary type
  vector<vector<int> > mSideNode; // node numbers on each cell side coincident with the mesh edge
  VVD mCoord; // coordinates of each node
  VD mVolume; // element volume
  vector<int> mbc_edge=vector<int>(4); // boundary condition on each edge of the mesh
  VD mbc_value=vector<double>(4); // boundary value on each edge of mesh
  vector<vector<int> > mE2E; // element->element connectivities

// member function signatures

  constexpr unsigned int str2int(const char* s); // to convert string to integer
  void set_E2E(); // a function to set up connectivities

};
