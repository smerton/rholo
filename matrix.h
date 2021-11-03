// Export the signature of the matrix class

using namespace std;

class Matrix{

  public:

  Matrix(int n);                           // constructor function
  ~Matrix();                               // destructor functions

  Matrix transpose();                      // returns transpose
  Matrix adjoint();                        // returns adjoint
  Matrix adjugate();                       // returns adjugate
  int active=0;                            // tests whether the current object is active
  void product(Matrix *A,Matrix *B);       // stores the matrix product AB
  double det();                            // returns determinant
  void solve(double*x,double*b);           // solves linear system Ax=b given pointers to x and b
  void inverse(Matrix *A);                 // generates the inverse for the matrix A
//  void inverse2(double* A, int N);       // generates the inverse for the matrix A using lapack
  void inverse2(Matrix *A);                // generates the inverse for the matrix A using lapack
  void copy(double**A);                    // copy to a matrix class object
  double read(int i,int j);                // reads element i,j of the member object
  void write(int i,int j,double dat);      // writes data dat to element i,j of the member object
  void add(int i,int j,double dat);        // adds data dat to element i,j of the member object

  int NRows();                             // returns the number of rows
  int NCols();                             // returns the number of columns

  private:

  double** mMat;                           // the matrix object
  int mN;                                  // number of rows(=number of columns)

};
