// signature for the input class

// Author S. R. Merton

#include <vector>
#include <string>

using namespace std;

class Input{

  public:

    Input(char* inputfile); // constructor for a new input
    ~Input(); // destructor for storage release

// accessor functions to member data

    string Title(); // returns problem title
    int NDescriptions(); // returns number of lines in the problem description
    string Description(int i); // returns line i of the problem description
    int Nodes(int idim); // returns the number of nodes along direction idim
    int NCells();  // returns the number of cells in the mesh
    double XMin(); // returns the mesh extent
    double XMax(); // returns the mesh extent
    double YMin(); // returns the mesh extent
    double YMax(); // returns the mesh extent
    double ZMin(); // returns the mesh extent
    double ZMax(); // returns the mesh extent
    string Boundary(int iface); // returns the boundary condition of face iface
    int Pk(); // returns the polyhedral order of the kinematic element
    int Pt(); // returns the polyhedral order of the thermodynamic element
    int NMaterials(); // returns the number of materials
    int Material(int imat); // returns material number of material imat
    int Ambient(); // returns ambient material number
    double Density(int imat); // returns density of material imat
    double Pressure(int imat); // returns pressure in material imat
    double Velocity(int idim,int imat); // returns velocity in direction idim of material imat
    double Range(int idim,int imat); // returns coordinate idim of material imat
    int NDims(); // number of mesh dimensions
    bool zNoh(); // returns true if mesh distortion is set for Noh problem
    bool zSedov(); // returns true if mesh distortion is set for Sedov problem
    bool zSaltzmann(); // returns true if mesh distortion is set for Saltzmann problem

  private:

// member data

    string mKeyword; // keyword to register an input line
    string mTitle;   // title of the problem
    string mbuffer;  // allows multiple lines with a common keyword
    vector<string> mDescription; // description of the problem
    int mNodes[3];   // number of nodes in each direction
    int mNCells;   // number of cells in the mesh
    double mXMin;    // mesh extent
    double mXMax;
    double mYMin;
    double mYMax;
    double mZMin;
    double mZMax;
    string mBoundaries[6]; // boundary conditions
    int mPk;               // polyhedral order for kinematics
    int mPt;               // polyhedral order for thermodynamics
    vector<int> mMaterial; // material number
    int mAmbient=1;        // ambient material
    vector<double> mDensity;     // density of material
    vector<double> mPressure;    // pressure of material
    vector<double> mVelocity[3]; // velocity of material
    vector<double> mRange[6];    // spatial boundary of material
    int mNDims; // number of mesh dimensions
    bool mzNoh; // adjust mesh for the Noh problem
    bool mzSedov; // adjust mesh for the Sedov problem
    bool mzSaltzmann; // adjust mesh for the Saltzmann problem

};
