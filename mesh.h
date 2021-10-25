// Export the signature of the mesh class

#include <vector>

using namespace std;

class Mesh{

  public:

  Mesh(char* meshfile);  // constructor for a new mesh
  ~Mesh();               // destructor function for storage release

// accessor functions to member data

  int NDims(); // returns the number of mesh dimensions
  int NNodes(); // returns the number of nodes on the mesh
  int NMaterials(); // returns the number of materials on the mesh
  int NCells(); // returns the number of cells on the mesh
  int Material(int iel); // returns material number in cell iel
  double Coord(int idim,int j); // returns coordinate idim of node j

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

// member function signatures

  constexpr unsigned int str2int(const char* s); // to convert string to integer

};
