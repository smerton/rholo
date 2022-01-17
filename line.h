// Export the signature of the line class

// Author S. R. Merton

#include <vector>

using namespace std;

class Line{

  public:

    Line(vector<double> r1,vector<double> r2); // constructor for a new line from start point r1 to end point r2

    ~Line();                                   // destructor to release class storage

// accessor functions to member data

  double start(int idim) const;       // returns coordinate idim of the start
  double end(int idim) const;         // returns coordinate idim of the end
  double m(int idim) const;           // returns coordinate idim of the gradient
  double m() const;                   // returns the gradient
  double length() const;              // returns the length of the line
  int nsegments() const;              // returns the number of segments
  double coord(int idim,int i) const; // returns segment i end point coordinate idim

// a function to divide a line into n segments

  void divide(int n);

  private:

// member data

  vector<double> mstart;          // coordinates of start point
  vector<double> mend;            // coordinates of end point
  vector<double> mm;              // gradient of the line
  double mlength;                 // length of the line
  int mnsegments;                 // number of segments
  vector<vector<double> > mcoord; // segment end point coordinates

};
