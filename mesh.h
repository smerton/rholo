// Export the signature of the mesh class

#define VD vector<double>           // laziness
#define VVD vector<vector<double> > // laziness

#define INCLUDE_GHOSTS 1 // used in ghost updates to update ghost data as well as physical data
#define EXCLUDE_GHOSTS 0 // used in ghost updates to update physical data only

#include <vector>

using namespace std;

class Mesh{

  public:

  Mesh(char* meshfile);  // constructor for a new mesh
  ~Mesh();               // destructor function for storage release

// accessor functions to member data

  int NDims() const; // returns the number of mesh dimensions
  int NNodes() const; // returns the number of nodes on the mesh
  int NGNodes() const; // returns the number of ghost nodes on the mesh edge 
  int NMaterials() const; // returns the number of materials on the mesh
  int NCells() const; // returns the number of cells on the mesh
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
  void InitCoords(VVD &v,int const flag); // initialise the mesh coordinates and include ghosts with flag=INCLUDE_GHOSTS
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
  void UpdateLength(VD &l,int const &p,VVD const &x) const; // update length of each element of polyhedral order p
  void UpdateVolume(VD &V,VVD const &x,int const &p) const; // update volume field V given coordinate x and polyhedral element order p
  void UpdateDensity(VD &d,VD const &V,VD const &m) const; // update denisty field d given a volume field V and a mass field m
  void UpdateEnergy(VD const &e0,VD &e1,VD const &p,VD const &q,VD const &V0,VD const &V1,VD const &m) const ; // update mesh energy field
  void UpdatePressure(VD &p,VD const &d,VD const &e,VD const &gamma,vector<int> const &mat); // load pressure field
  void UpdateSoundSpeed(VD &c,VD const &g,vector<int> const &mat,VD const &p,VD const &d) const; // load new sound speeds

  private:

// member data

  int mNDims; // number of mesh dimensions
  int mNNodes; // number of nodes on the mesh
  int mNGNodes; // number of ghost nodes on the mesh edge
  int mNMaterials; // number of materials on the mesh
  int mNCells; // number of cells on the mesh
  int mNGCells; // number of ghost cells on the mesh edge
  int mNSides; // number of cell sides coinciding with the mesh edge
  vector<int> mMaterial; // material in each cell
  vector<int> mType; // polyhedral type of each cell
  vector<vector<int> > mVertex; // node number of each vertex in each cell
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
