// Export the signature of the timer class

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
