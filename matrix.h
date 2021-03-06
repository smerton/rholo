// Export the signature of the matrix class

using namespace std;

class Matrix{

  public:

  Matrix(int n);                           // constructor function
  ~Matrix();                               // destructor functions

  Matrix transpose();                      // returns transpose
  Matrix adjoint();                        // returns adjoint
  Matrix product(Matrix*B);                // returns the product with B
  double det();                            // returns determinant
  void solve(double*x,double*b);           // solves linear system Ax=b given pointers to x and b

  Matrix inverse();                        // returns a Matrix object containing the inverse
  void copy(double**A);                    // copy to a matrix class object
  double read(int i,int j);                // reads element i,j of the member object
  void write(int i,int j,double dat);      // writes data dat to element i,j of the member object

  int NRows();                             // returns the number of rows
  int NCols();                             // returns the number of columns

  private:

  double** mMat;                           // the matrix object
  int mN;                                  // number of rows(=number of columns)

};
