// Export the signature of the mesh class
// prefix m denotes member data
// prefix set denotes a function to write to member data so that the class member can remain private

#include <vector>
#include <string>

using namespace std;

class Mesh{

  public:

    Mesh(char* meshfile);  // constructor for a new mesh
    ~Mesh();               // destructor function for storage release
    void insert_dummies(); // function to insert dummies for a ucd silo mesh

// accessor functions to member data

    int NDims(); // returns the number of mesh dimensions
    long Dim(int idim); // returns the length of mesh dimension idim
    long NVertices(); // returns the number of mesh vertices
    long NCells(); // returns the number of mesh cells
    long NGhosts(); // returns the number of ghosts
    int NMaterials(); // returns the number of materials on the mesh
    int Materials(int imat); // returns the mesh material number of material imat
    int NTypes(); // returns the number of polyhedral element types on the mesh
    int Types(int itype); // returns the polyhedral type of type itype
    string Names(int itype); // returns the name of polyhedral type itype
    int NOrders(); // returns the number of polyhedral orders on the mesh
    int Orders(int iorder); // returns the polyhedral order of order iorder
    int KOrders(int iorder); // returns the polyhedral order of order iorder for the kinematics
    int TOrders(int iorder); // returns the polyhedral order of order iorder for the thermodynamics
    int Order(long iel); // returns the polyhedral order of element iel
    int KOrder(long iel); // returns the polyhedral order of element iel to use for kinematics
    int TOrder(long iel); // returns the polyhedral order of element iel to use for thermodynamics
    int Type(long iel); // returns the polyhedral type of element iel
    int Material(long iel); // returns the material number in element iel
    long NCorners(); // returns the element corner sum
    long Vertex(long iel,int icorner); // returns the vertex of corner icorner in element iel
    double Coord(int idim,long ivert); // returns the coordinate idim of vertex ivert
    double DGCoord(int idim,long inod); // returns the coordinate idim of global node inod
    double KCoord0(int idim,long inod); // returns the coordinate idim of kinematic node inod at the start of a time step
    double KCoord(int idim,long inod); // returns the coordinate idim of kinematic node inod at the half-step
    double TCoord0(int idim,long inod); // returns the coordinate idim of thermodynamic node inod at the start of a time step
    double TCoord(int idim,long inod); // returns the coordinate idim of thermodynamic node inod at the half-step
    int NLoc(long iel); // returns the number of local nodes in element iel in the meshfile
    int MLoc(long iel); // returns the number of local nodes to insert into the finite element
    int KLoc(long iel); // returns the number of local kinematic nodes to insert into the finite element
    int TLoc(long iel); // returns the number of local thermodyanmic nodes to insert into the finite element
    long NNodes(); // returns the number of finite element nodes in the mesh
    long KNodes(); // returns the number of kinematic nodes in the mesh
    long TNodes(); // returns the number of thermodynamic nodes in the mesh
    double Density(long iel); // returns the density in element iel
    double Pressure(long iel); // returns the pressure in element iel
    double Velocity0(int idim,long inod); // returns component idim of the start-of-step velocity at global node inod
    double Velocity(int idim,long inod); // returns component idim of the full-step velocity at global node inod
    double MWVelocity(int idim,long inod); // returns component idim of the mass-weighted velocity at global node inod
    double Energy0(long inod); // returns the start-of-step energy field at global node inod
    double Energy(long inod); // returns the half-step energy field at global node inod
    double CCEnergy0(long iel); // returns the start-of-step cell-centred energy field in element iel
    double CCEnergy(long iel); // returns the half-step cell-centred energy field in element iel
    long Node(long iel,int iloc); // returns the global node of local node iloc in element iel
    long KNode(long iel,int iloc); // returns the kinematic node of local node iloc in element iel
    long TNode(long iel,int iloc); // returns the thermodynamic node of local node iloc in element iel
    long E2E(long iel,int iface); // returns the neighbour element on face iface of element iel
    int N2N(long iel,int iface,int ifloc); // returns local face node on neighbour face connecting face node ifloc on face iface of iel
    int F2F(long iel,int iface); // returns face of the neighbour on face iface of iel
    int NFaces(long iel); // returns number of faces on the surface of inserted  element iel
    int NFLoc(long iel,int iface); // returns number of local nodes on face iface of inserted element iel
    int FLoc(long iel,int iface,int ifloc); // returns local node number of node ifloc on face iface of thermodynamic element iel
    int FLock(long iel,int iface,int ifloc); // returns local node number of node ifloc on face iface of kinematic element iel
    long NodeList(long ilist); // returns entry ilist on the silo nodelist
    long KNodeList(long ilist); // returns entry ilist on the silo kinematic nodelist
    long TNodeList(long ilist); // returns entry ilist on the silo thermodynamic nodelist
    long lNodeList(); // returns the length of the silo nodelist
    long lKNodeList(); // returns the length of the silo kinematic nodelist
    long lTNodeList(); // returns the length of the silo thermodynamic nodelist
    double EnergyInit(int imat); // returns the initial energy in material imat
    double PressureInit(int imat); // returns the initial pressure in material imat
    double DensityInit(int imat); // returns the initial density in material imat
    double VelocityXInit(int imat); // returns the initial velocity in the x-direction in material imat
    double VelocityYInit(int imat); // returns the initial velocity in the y-direction in material imat
    double VelocityZInit(int imat); // returns the initial velocity in the z-direction in material imat
    double Volume0(long iel); // returns the start-of-step volume in element iel
    double Volume(long iel); // returns the volume in element iel
    double Mass(long iel); // returns the mass of element iel
    void setVolume0(vector<double> *new_volumes); // updates the start-of-step element volumes to the new vector passed in
    void setVolume(vector<double> *new_volumes); // updates the element volumes to the new vector passed in
    void setPressure(vector<double> *new_pressures); // updates the element pressures to the new vector passed in
    void setDensity(vector<double> *new_densities); // updates the element densities to the new vector passed in
    void setMass(vector<double> *new_mass); // updates the element mass to the new vector passed in
    void setEnergy0(vector<double> *new_energies); // updates the start-of-step energy field to the new vector passed in
    void setEnergy(vector<double> *new_energies); // updates the half-step energy field to the new vector passed in
    void setCCEnergy0(vector<double> *new_energies); // updates the start-of-step cell-centred energy field to the new vector passed in
    void setCCEnergy(vector<double> *new_energies); // updates the half-step cell-centred energy field to the new vector passed in
    void setKCoord0(int idim,long inod,double pos); // move kinematic node inod to position pos along direction idim
    void setKCoord(int idim,long inod,double pos); // move kinematic node inod to position pos along direction idim
    void setTCoord0(int idim,long inod,double pos); // move thermodynamic node inod to position pos along direction idim
    void setTCoord(int idim,long inod,double pos); // move thermodynamic node inod to position pos along direction idim
    void setVelocity0(vector<vector<double> > *new_velocities); // update the start-of-step velocity field to the new vector passed in
    void setVelocity(vector<vector<double> > *new_velocities); // update the full-step velocity field to the new vector passed in
    void setMWVelocity(vector<vector<double> > *new_velocities); // update the mass-weighted velocity field to the new vector passed in
    long NVelocities0(int idim); // returns the size of the start-of-step velocity field
    long NVelocities(int idim); // returns the size of the full-step velocity field
    long NBoundaryNodes(int iface); // returns number of nodes found along mesh boundary containing element face iface
    long BoundaryNode(int iface,long inod); // global node number of node inod on mesh boundary containing element face iface
    int NCoincidentNodes(long inod); // returns the numebr of nodes that are coincident with global node inod
    long CoincidentNode(long inod,int i); // returns the global number of node i that is coincident with global node inod
    string BoundaryType(int iface); // returns boundary condition to apply to face iface of the mesh (REFLECTIVE || TRANSMISSIVE)
    int Reflect(int iloc,int iface,int p); // reflect local node iloc across face iface of an order p element
    bool UsePdV(); // swicthes between cell-centred PdV and nodal finite element energy fields

// nodelists for silo

    vector<long> mElement; // element number on the silo mesh

// signatures of member functions to initialise the mesh data fields to input values

    void pressure_init(); // initial pressure field
    void density_init(); // initial density field
    void energy_init(); // initial energy field
    void velocity_init(); // initial velocity field
    void velocity_bcs(); // initial velocity boundary conditions
    void update_halo(); // update the ghost halo pressure, mass, density, energy and volume fields

// signature of function to swap start-of-step values

    void end_of_step();

  private:

// member data

    int mNDims; // number of mesh dimensions
    vector<long> mDim; // length of each mesh dimension
    long mNVertices; // number of mesh vertices
    long mNCells; // number of mesh cells
    long mNGhosts=0; // number of ghost cells
    int mNMaterials; // number of materials
    vector<int> mMaterials; // materials in the mesh
    int mNTypes; // number of polyhedral types in the mesh
    vector<int> mTypes; // polyhedral types on the mesh
    vector<string> mNames; // name of each polyhedral type on the mesh
    int mNOrders; // number of polyhedral orders present
    vector<int> mOrders; // polyhedral orders present
    vector<int> mKOrders; // polyhedral orders present for kinematics
    vector<int> mTOrders; // polyhedral orders present for thermodynamics
    vector<int> mOrder; // polyhedral order of each element
    vector<int> mKOrder; // polyhedral order to use for kinematics 
    vector<int> mTOrder; // polyhedral order to use for thermodynamics
    vector<int> mType; // polyhedral type of each element
    vector<int> mMaterial; // material number in each element
    long mNCorners=0l; // element corner sum to define the length of the vertex list
    vector<long> mVertex; // mesh vertex of each element corner
    vector<vector<double> > mCoord; // vertex coordinates
    vector<vector<double> > mDGCoord; // finite element node coordinates
    vector<vector<double> > mKCoord0; // finite element node coordinates on kinematic grid at the start of a time step
    vector<vector<double> > mKCoord; // finite element node coordinates on kinematic grid at the half-step
    vector<vector<double> > mTCoord0; // finite element node coordinates on thermodynamic grid at the start of a time step
    vector<vector<double> > mTCoord; // finite element node coordinates on thermodynamic grid at the half-step
    vector<int> mNLoc; // number of local nodes in each meshfile element
    vector<int> mMLoc; // number of finite element nodes in each element
    vector<int> mKLoc; // number of finite element nodes in each element on the kinematic grid
    vector<int> mTLoc; // number of finite element nodes in each element on the thermodynamic grid
    long mNNodes; // number of finite element nodes in the mesh
    long mKNodes; // number of finite element nodes in the mesh on the kinematic grid
    long mTNodes; // number of finite element nodes in the mesh on the thermodynamic grid
    vector<double> mDensity; // element density
    vector<double> mPressure; // element pressure
    vector<vector<double> > mVelocity0; // start-of-step nodal velocity components
    vector<vector<double> > mVelocity; // full-step nodal velocity components
    vector<vector<double> > mMWVelocity; // mass-weighted nodal velocity components
    vector<double> mEnergy0; // nodal energy field at the start of a time step
    vector<double> mEnergy; // nodal energy field at the half-step
    vector<double> mCCEnergy0; // cell-centred energy field at the start of a time step
    vector<double> mCCEnergy; // cell-centred energy field at the half-step
    vector<vector<long> > mNode; // global node numbers
    vector<vector<long> > mKNode; // kinematic node numbers
    vector<vector<long> > mTNode; // thermodynamic node numbers
    vector<vector<long> > mE2E; // element-element connectivities
    vector<vector<vector<int> > > mN2N; // node-node connectivities
    vector<vector<int> > mF2F; // face-face connectivities
    vector<int> mNFaces; // number of faces on each inserted element
    vector<vector<int> > mNFLoc; // number of local nodes on each face of each inserted element
    vector<vector<vector<int> > > mFLoc; // local node numbers on each face of each thermodynamic element
    vector<vector<vector<int> > > mFLock; // local node numbers on each face of each kinematic element
    vector<long> mNodeList; // silo nodelist array for dummy element mesh
    vector<long> mKNodeList; // silo nodelist array for dummy element mesh
    vector<long> mTNodeList; // silo nodelist array for dummy element mesh
    vector<double> mEnergyInit; // initial energy in each material
    vector<double> mPressureInit; // initial pressure in each material
    vector<double> mDensityInit; // initial density in each material
    vector<double> mVelocityXInit; // initial velocity in x direction in each material
    vector<double> mVelocityYInit; // initial velocity in y direction in each material
    vector<double> mVelocityZInit; // initial velocity in z direction in each material
    vector<double> mVolume0; // element volume at start of a time step
    vector<double> mVolume; // element volume
    vector<double> mMass; // element mass
    vector<vector<long> > mBoundaryNode; // global nodes along mesh boundaries
    vector<vector<long> > mCoincidentNode; // list of coincident kinematic nodes
    vector<string> mBoundaryType; // boundary condition to apply to mesh edges
    bool mUsePdV=false; // switches between PdV and FE energy fields

// signature of member function to insert requested element in to each meshfile cell

    void fill_dg();

// signature of member function to perform connectivity traversal

    void set_conn();

// signature of member function to generate halo of ghost cells

    void set_ghosts();

// signature of function to generate list of nodes along the mesh boundaries

    void boundary_nodes();

// signature of function to generate list of coincident nodes

    void coincident_nodes();

};
