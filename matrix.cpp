// Function definitions for the matrix class

// Author S. R. Merton

#include <iostream>
#include "matrix.h"

// interface to LAPACK

extern "C"{

// LU decomposition of a general matrix

void dgetrf_(int* M,int* N,double* A,int* lda,int* IPIV,int *INFO);

// generate inverse of a matrix given its LU decomposition

void dgetri_(int* N,double* A,int* lda,int* IPIV,double* WORK,int* lwork,int* INFO);

}

using namespace std;

//void Matrix::inverse2(double* A, int N){
//
//// emit the inverse via lapack
//
//  int *IPIV=new int[N];
//  int LWORK=N*N;
//  double *WORK=new double[LWORK];
//  int INFO;
//
//  dgetrf_(&N,&N,A,&N,IPIV,&INFO);
//  dgetri_(&N,A,&N,IPIV,WORK,&LWORK,&INFO);
//
//  delete[] IPIV;
//  delete[] WORK;
//
//}

void Matrix::inverse2(Matrix *A){

// emit the inverse via lapack driver routines dgetrf() & dgetri()

  int N(A->NRows());
  int *IPIV=new int[N];
  int LWORK=N*N;
  double *WORK=new double[LWORK];
  int INFO;
  double *B=new double[N*N];

// flatten

//  cout<<"  Matrix::inverse2(): flattening..."<<endl;

  for(int i=0,k=0;i<N;i++){
    for(int j=0;j<N;j++,k++){
      B[k]=A->read(i,j);
    }
  }

//  cout<<"  Matrix::inverse2(): performing LU factorisation..."<<endl;
  dgetrf_(&N,&N,B,&N,IPIV,&INFO);

//  cout<<"  Matrix::inverse2(): finding inverse using LU factorisation..."<<endl;
  dgetri_(&N,B,&N,IPIV,WORK,&LWORK,&INFO);

// unflatten

//  cout<<"  Matrix::inverse2(): unflattening..."<<endl;
  for(int i=0,k=0;i<N;i++){
    for(int j=0;j<N;j++,k++){
      mMat[i][j]=B[k];
    }
  }

  delete[] IPIV;
  delete[] WORK;
  delete[] B;

}

// constructor for a new square matrix object
// A[n][n] this constructor would be suitable
// for a type of matrix whose inverse may be
// required

Matrix::Matrix(int n){

// set the number of rows and columns

  mNCols=n;
  mNRows=n;

  mMat=new double*[NRows()];
  for(long i=0;i<NRows();i++){
    mMat[i]=new double[NCols()];
  }

// initialise elements of the matrix

  for(int i=0;i<this->NRows();i++){
    for(int j=0;j<this->NCols();j++){
      this->write(i,j,0.0);
    }
  }

// mark object as active

  active=1;

}

// constructor for a new non-square matrix object
// A[nrows][ncols] this constructor would be
// suitable for a type of matrix whose inverse
// is not required and needs to be more general

Matrix::Matrix(int nrows,int ncols){

// set the number of rows and columns

  mNRows=nrows;
  mNCols=ncols;

  mMat=new double*[NRows()];
  for(long i=0;i<NRows();i++){
    mMat[i]=new double[NCols()];
  }

// initialise elements of the matrix

  for(int i=0;i<this->NRows();i++){
    for(int j=0;j<this->NCols();j++){
      this->write(i,j,0.0);
    }
  }

// mark object as active

  active=1;

}

// Member function to return the number of rows

int Matrix::NRows(){return mNRows;}

// Member function to return the number of columns

int Matrix::NCols(){return mNCols;}

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

//  double A1[NRows()][NRows()]; // on stack

  double** A1=new double*[NRows()]; // on heap
  for(long i=0;i<NRows();i++){
    A1[i]=new double[NRows()];
  }

  for(long i=0;i<NRows();i++){
    for(long j=0;j<NRows();j++){
      A1[i][j]=mMat[i][j];
    }
    x[i]=0.0; // otherwise this is uninitialised (robustness issues)
  }

// form linear system x=A^{-1}b using LU decomposition of A=[L}{U}
// so that Ax=b becomes L[Ux]=b in which Ux=y and Ly=b

  for(long k=0;k<NRows()-1;k++){
    for(long i=k+1;i<NRows();i++){
      A1[i][k]=A1[i][k]/A1[k][k];
    }
    for(long j=k+1;j<NRows();j++){
      for(long i=k+1;i<NRows();i++){
        A1[i][j]=A1[i][j]-A1[i][k]*A1[k][j];
      }
    }
  }

// solve Ly=b by row elimination

  for(long i=0;i<NRows();i++){
    double r(0.0);
    for(long j=0;j<NRows()-1;j++){
      r+=A1[i][j]*x[j];
    }
    x[i]=b[i]-r;
  }

// Ux=y by row elimination

  for(long i=NRows()-1;i>=0;i--){
    double r(0.0);
    for(long j=i+1;j<NRows();j++){
      r+=A1[i][j]*x[j];
    }
    x[i]=(x[i]-r)/A1[i][i];
  }

// release heap storage

  for(long i=0;i<NRows();i++){
    delete[] A1[i];
    A1[i]=NULL;
  }
  delete[] A1;
  A1=NULL;

  return;

}

// Member function to form the inverse

void Matrix::inverse(Matrix *A){

  double x[NRows()],row[NRows()];

// collect the inverse by row using LU decomposition with row elimination

  for(int i=0;i<NRows();i++){

// reset the i'th row

    for(int j=0;j<NRows();j++){x[j]=0.0;row[j]=0.0;}
    row[i]=1.0;

// pass to the solver

    A->solve(x,row);

// unpack the column which is row i of the inverse

    for(int j=0;j<NRows();j++){mMat[j][i]=x[j];}

  }

  return;

}

// Member function to form the product of 2 matrices

void Matrix::product(Matrix *A,Matrix *B){

  for(int i=0;i<A->NRows();i++){
    for(int j=0;j<B->NCols();j++){
      mMat[i][j]=0.0;
      for(int k=0;k<A->NCols();k++){
        mMat[i][j]+=A->read(i,k)*B->read(k,j);
      }
    }
  }

  return;

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

// mark object as inactive

  active=0;

// release storage

  for(long i=0;i<this->NRows();i++){
    delete[] mMat[i];
    mMat[i]=NULL;
  }
  delete[] mMat;
  mMat=NULL;

}
