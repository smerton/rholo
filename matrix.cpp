// Function definitions for the matrix class

// Author S. R. Merton

#include <iostream>
#include "matrix.h"

using namespace std;

Matrix::Matrix(int n){

// set the number of rows and columns

  mN=n;

  mMat=new double*[NRows()];
  for(int i=0;i<NRows();i++){
    mMat[i]=new double[NCols()];
  }

// initialise elements of the matrix

  for(int i=0;i<this->NRows();i++){
    for(int j=0;j<this->NCols();j++){
      this->write(i,j,0.0);
    }
  }
}

// Member function to return the number of rows

int Matrix::NRows(){return mN;}

// Member function to return the number of columns

int Matrix::NCols(){return mN;}

// Member function to return the transpose

Matrix Matrix::transpose(){

// set up a temporary Matrix object

  Matrix tmp(NRows());

// compute the transpose

  for(int i=0;i<NRows();i++){
    for(int j=0;j<NCols();j++){
      tmp.write(i,j,1.0+mMat[j][i]);
    }
  }

  return tmp;

}

// Member function to solve the linear system Ax=b given pointers to x and b

void Matrix::solve(double*x,double*b){

// local copy of the matrix

  double A1[NRows()][NRows()];
  for(int i=0;i<NRows();i++){
    for(int j=0;j<NRows();j++){
      A1[i][j]=mMat[i][j];
    }
    x[i]=0.0; // otherwise this is uninitialised (robustness issues)
  }

// form linear system x=A^{-1}b using LU decomposition of A=[L}{U}
// so that Ax=b becomes L[Ux]=b in which Ux=y and Ly=b

  for(int k=0;k<NRows()-1;k++){
    for(int i=k+1;i<NRows();i++){
      A1[i][k]=A1[i][k]/A1[k][k];
    }
    for(int j=k+1;j<NRows();j++){
      for(int i=k+1;i<NRows();i++){
        A1[i][j]=A1[i][j]-A1[i][k]*A1[k][j];
      }
    }
  }

// solve Ly=b by row elimination

  for(int i=0;i<NRows();i++){
    double r(0.0);
    for(int j=0;j<NRows()-1;j++){
      r+=A1[i][j]*x[j];
    }
    x[i]=b[i]-r;
  }

// Ux=y by row elimination

  for(int i=NRows()-1;i>=0;i--){
    double r(0.0);
    for(int j=i+1;j<NRows();j++){
      r+=A1[i][j]*x[j];
    }
    x[i]=(x[i]-r)/A1[i][i];
  }

  return;

}

// Member function to return the inverse

Matrix Matrix::inverse(){

  double x[NRows()],row[NRows()];
  Matrix B(NRows()),BI(NRows());

// local copy of the matrix to avoid overwriting it

  for(int i=0;i<NRows();i++){
    for(int j=0;j<NRows();j++){
      B.write(i,j,mMat[i][j]);
    }
  }

// collect the inverse by row using LU decomposition with row elimination

  for(int i=0;i<NRows();i++){

// reset the i'th row

    for(int j=0;j<NRows();j++){x[j]=0.0;row[j]=0.0;}
    row[i]=1.0;

// pass to the solver

    B.solve(x,row);

// unpack the column which is row i of the inverse

    for(int j=0;j<NRows();j++){BI.write(j,i,x[j]);}

  }

  return BI;

}

// Member function to return the product of 2 matrices

Matrix Matrix::product(Matrix*B){

  Matrix C(NRows());

//C=mMat*B

  for(int i=0;i<NRows();i++){
    for(int j=0;j<NRows();j++){
      double r(0.0);
      for(int k=0;k<NRows();k++){
        r+=mMat[i][k]*B->read(k,j);
      }
      C.write(i,j,r);
    }
  }

  return C;

}

// Member function to return the adjoint (adjugate)

Matrix Matrix::adjoint(){

// set up a temporary Matrix object

  Matrix tmp(NRows());

// compute the adjoint using cofactor expansion

  return tmp;

}

// Member function to copy data into a matrix class object

void Matrix::copy(double**A){

  for(int i=0;i<NRows();i++){
    for(int j=0;j<NCols();j++){
      mMat[i][j]=A[i][j];
    }
  }

  return;

}

// Member function to read element i,j of the matrix

double Matrix::read(int i,int j){return mMat[i][j];}

// Member function to write data into the member object

void Matrix::write(int i,int j,double dat){mMat[i][j]=dat;}

// Member function to add data into the member object

void Matrix::add(int i,int j,double dat){mMat[i][j]+=dat;}

// Destructor for the matrix class

Matrix::~Matrix(){

  for(int i=0;i<NCols();i++){
    delete[] mMat[i];
    mMat[i]=NULL;
  }
  delete[] mMat;
  mMat=NULL;

}
