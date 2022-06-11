// Function definitions to code for some utilities

// Author S. R. Merton

#include <vector>

using namespace std;

// type safe function to return the sign of the argument

//template <typename T> int sgn(T val) {return(T(0)<val)-(val<T(0));} // -1,0 or 1
template <typename T> int sgn(T val) {return( (val>=T(0))?T(1):T(-1));} // -1 or 1

// function to empty a vector

void vempty(vector<double>&v){vector<double> e;v.swap(e);return;}
