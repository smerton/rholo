// Function definitions for the timer class

// Author S. R. Merton

#include <iostream>
#include <chrono>
#include "timer.h"

using namespace std;
using namespace std::chrono;

Timer::Timer(int n){

// set number of timers

  mntimers=n;

  cout<<mntimers<<" timers have been requested"<<endl;

  mSpan=new double[mntimers];

  mTime=new std::chrono::time_point<std::chrono::high_resolution_clock>[mntimers];

  mLastTime=new std::chrono::time_point<std::chrono::high_resolution_clock>[mntimers];

  mduration=new std::chrono::duration<double, std::milli>[mntimers];

}

// Member function to reset all timers

void Timer::Init(){


//  mTime=std::chrono::high_resolution_clock::now();
//  mLastTime=std::chrono::high_resolution_clock::now();

//  mduration=mTime-mLastTime;

  for(int i=0;i<mntimers;i++){

    mTime[i]=std::chrono::high_resolution_clock::now();
    mLastTime[i]=std::chrono::high_resolution_clock::now();

    mduration[i]=mTime[i]-mLastTime[i];

    mSpan[i]=0.0;
  }

  mTotal=0.0;

  return;

}

// Member function to start a timer

void Timer::Start(int i){

  mLastTime[i]=std::chrono::high_resolution_clock::now();

  return;

}

// Member function to stop a timer

void Timer::Stop(int i){

  mTime[i]=std::chrono::high_resolution_clock::now();
  mduration[i]=mTime[i]-mLastTime[i];
  mSpan[i]+=mduration[i].count()*1.0e-03;

  return;

}

// Member function to sum the time spans to find the total

double Timer::Total(){

  mTotal=0.0;

  for(int i=0;i<mntimers;i++){
    mTotal+=mSpan[i];
  }

  return mTotal;

}

// Accessor function to return the span of a timer

double Timer::Span(int i){return mSpan[i];}

// destructor for the timer class

Timer::~Timer(){

  delete[] mSpan;

}
