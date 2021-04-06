// Export the signature of the shape class

// Author S. R. Merton

#include <vector>

using namespace std;

class Shape{

  public:

    Shape(int n); // constructor function for a new shape of order n
    ~Shape(); // destructor to release class storage

// accessor functions to member data

    int order(); // returns the polyhedral order of the shape
    int nloc(); // returns number of local nodes
    int sloc(); // returns number of nodes on the surface
    int ngi(); // number of Gauss integration points
    int nfaces(); // number of element faces
    int reflect(int iloc,int iface); // refelct iloc across face iface

    double wgt(int gi); // quadrature weight of integration point gi
    double value(int i,int gi); // shape i value at Gauss point gi
    double value(int i,double x); // shape i value at local coordinate x
    double dvalue(int i,int gi); // derivative i value at Gauss point gi
    double dvalue(int i,double); // derivative i value at local coordinate x

  private:

// member data

    int morder; // polyhedral order of the shape
    int mnloc; // number of local nodes
    int msloc; // number of surface nodes
    int mngi; // number of Gauss points
    int mnfaces; // number of faces

    vector<vector<int> > mreflect; // reflection across a face

    vector<double> mwgt; // quadrature weight
    vector<vector<double> > mvalue; // shape value
    vector<vector<double> > mdvalue; // derivative value

};
