// Function definitions for eos lookups

// Author S. R. Merton

#include <iostream>
#include "eos.h"
#include <cmath>

using namespace std;

// return pressure given the energy

double P(double d,double e,double g){return (g-1.0)*d*e;}

// invert the eos to return energy given the pressure

double E(double d,double p,double g){return p/((g-1.0)*d);}

// sound speed

double C(double p,double d,double g){return sqrt(g*p/d);}
