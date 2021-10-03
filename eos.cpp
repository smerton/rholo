// Function definitions for eos lookups

// Author S. R. Merton

#include <iostream>
#include "eos.h"

using namespace std;

// return pressure given the energy

double P(double d,double e){return (GAMMA-1.0)*d*e;}

// invert the eos to return energy given the pressure

double E(double d,double p){return p/((GAMMA-1.0)*d);}
