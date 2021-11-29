// Signature for eos lookups to return pressure, energy and sound speed

// Author S. R. Merton

#define GAMMA (5.0/3.0)          // ratio of specific heats, usually 1.4 or 5/3

using namespace std;

double P(double d,double e); // pressure as a function of energy
double E(double d,double p); // invert the eos to get energy as a function of pressure
double C(double p,double d); // sound speed
