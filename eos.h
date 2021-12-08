// Signature for eos lookups to return pressure, energy and sound speed

// Author S. R. Merton

using namespace std;

double P(double d,double e,double g); // pressure as a function of energy
double E(double d,double p,double g); // invert the eos to get energy as a function of pressure
double C(double p,double d,double g); // sound speed
