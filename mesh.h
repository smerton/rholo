// Export the signature of the mesh class

#include <vector>

using namespace std;

class Mesh{

  public:

  Mesh(char* meshfile);  // constructor for a new mesh
  ~Mesh();               // destructor function for storage release

// accessor functions to member data

  int NDims() const; // returns the number of mesh dimensions
  int NNodes() const; // returns the number of nodes on the mesh
  int NMaterials() const; // returns the number of materials on the mesh
  int NCells() const; // returns the number of cells on the mesh
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
  void InitCoords(vector<vector<double> > &v); // initialise the mesh coordinates
  double Volume(int i) const; // returns the volume of the element

  private:

// member data

  int mNDims; // number of mesh dimensions
  int mNNodes; // number of nodes on the mesh
  int mNMaterials; // number of materials on the mesh
  int mNCells; // number of cells on the mesh
  int mNSides; // number of cell sides coinciding with the mesh edge
  vector<int> mMaterial; // material in each cell
  vector<int> mType; // polyhedral type of each cell
  vector<vector<int> > mVertex; // node number of each vertex in each cell
  vector<int> mSideAttr; // cell side boundary attribute
  vector<int> mSideType; // cell side boundary type
  vector<vector<int> > mSideNode; // node numbers on each cell side coincident with the mesh edge
  vector<vector<double> > mCoord; // coordinates of each node
  vector<double > mVolume; // element volume

// member function signatures

  constexpr unsigned int str2int(const char* s); // to convert string to integer

};
