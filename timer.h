// Export the signature of the timer class

#define TIMER_MAIN 0          // timer for main
#define TIMER_ASSEMBLY 1      // timer for acceleration field assembly
#define TIMER_INVERSE 2       // timer for matrix inverter
#define TIMER_ENERGY 3        // timer for energy field assembly
#define TIMER_FORCE 4         // timer for force assembly
#define TIMER_KSOLVE 5        // timer for kinematic solve
#define TIMER_GRAPHICS 6      // timer for graphics output
#define TIMER_CFL 7           // timer for CFL calculation
#define TIMER_MOTION 8        // timer for Lagrangian motion
#define TIMER_ECHECK 9        // timer for energy conservation check
#define TIMER_OUTPUT 10       // timer for output
#define TIMER_VISCOSITY 11    // timer for viscous terms
#define TIMER_LOCATE 12       // timer to locate untimed stuff

using namespace std;

class Timer{

  public:

  Timer(int n);              // constructor function for n timers
  ~Timer();                  // destructor function

  void Init();               // reset the time
  void Start(int i);         // start timer i
  void Stop(int i);          // stop timer i
  double Span(int i);        // returns the time span for timer i
  double Total();            // returns the sum of all time spans

  private:

  std::chrono::time_point<std::chrono::high_resolution_clock> *mTime;     // current time
  std::chrono::time_point<std::chrono::high_resolution_clock> *mLastTime; // current time
  std::chrono::duration<double, std::milli> *mduration;                   // duration

  int mntimers; // number of timers
  double *mSpan; // pointer to the time spans
  double mTotal; // sum of all time spans

};
