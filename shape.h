// Export the signature of the shape class

// Author S. R. Merton

#define CONTINUOUS 1
#define DISCONTINUOUS 2

#include <vector>

using namespace std;

class Shape{

  public:

    Shape(int n);                                                  // constructor function for a new shape of order n in isoparametric coordinates
    Shape(int n,int ngi,int t);                                    // new isoparametric shape of order n, type t=DISCONTINUOUS or t=CONTINUOUS using ngi integration points per axis
    Shape(int n,vector<vector<double> > x);                        // constructor function for a new shape of order n in global coordinates x

    ~Shape();                                                      // destructor to release class storage

// accessor functions to member data

    int order() const;                                             // returns the polyhedral order of the shape
    int type() const;                                              // returns the type of the shape
    int ndims() const;                                             // returns the number of dimensions of the shape
    int nloc() const;                                              // returns number of local nodes
    int sloc() const;                                              // returns number of nodes on the surface
    int ngi() const;                                               // number of Gauss integration points
    int nfaces() const;                                            // number of element faces
    int reflect(int iloc,int iface) const;                         // refelct iloc across face iface
    int pos(int idim,int iloc) const;                              // node number in dimension idim

    double wgt(int gi) const;                                      // quadrature weight of integration point gi
    double value(int i,int gi) const;                              // shape i value at Gauss point gi
    double dvalue(int idim,int i,int gi) const;                    // shape i derivative idim at Gauss point gi
    double value(int i,double u,double v) const;                   // shape i value at local coordinates u,v
    double value(int i,vector<double> x) const;                    // shape i value at global coordinate x
    double dvalue(int idim,int i,double u,double v) const;         // shape i derivative idim at local coordinates u,v
    double dvalue(int idim,int i,vector<double> x) const;          // shape i derivative idim at global coordinate x
    double coeff(int i,int j) const;                               // coefficient j in the polynomial expansion of shape i
    double integrate(int i) const;                                 // integrate shape i in global coordinates
    double integrate(int idim,int i) const;                        // integrate derivative idim of shape i in global coordinates

// accessor function to prolongation operator

    void prolongate(double*u,double*v,int p) const;                // prolongation operator to map vector u[] to an order p element

  private:

// member data

    int morder;                                         // polyhedral order of the shape
    int mtype;                                          // type of shape can be DISCONTINUOUS or CONTINUOUS
    int mndims;                                         // number of dimensions of the shape
    int mnloc;                                          // number of local nodes
    int msloc;                                          // number of surface nodes
    int mngi;                                           // number of Gauss points
    int mnfaces;                                        // number of faces
    vector<int> mpos[3];                                // node number in each dimension
    vector<vector<int> > mreflect;                      // reflection across a face
    vector<double> mwgt;                                // quadrature weights
    vector<vector<double> > mvalue;                     // shape values
    vector<vector<double > > mdvalue[3];                // derivative values
    vector<vector<double> > mcoeff;                     // coefficients of the polynomial
    vector<vector<double> > mr;                         // coordinates of the shape function nodes

};
