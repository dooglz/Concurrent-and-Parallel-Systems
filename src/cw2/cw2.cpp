#include <stdio.h>
#include <stdlib.h>
#include <thread>
#include <cvmarkersobj.h>

#include <vector>
#include <algorithm>
#include <iostream>
#include <chrono>
#include <string>
#include "Renderer.h"
#include "nbody.h"
#include "timer.h"
#include <fstream>

using namespace std;
using namespace Concurrency::diagnostic;

uint16_t simStalls = 0;
uint16_t visStalls = 0;
mutex simThreadLock;
bool run;

long long simTimes[100];
long PARTICLESIZE = 28096;
long CORECOUNT = 7;
long SAMPLESIZE = 50000;

std::ofstream csv;
template <typename T> T average(T t[], int n) {
  T s = t[n - 1];
  for (int i = 0; i < (n - 1); i++)
    s += t[i];
  return s / n;
}


long bb = 0;
void renderThread() {
	marker_series series;
  string aa = to_string(CORECOUNT) + "_" + to_string(PARTICLESIZE) + ".csv";
  cout << aa << std::endl;
  csv = std::ofstream(aa, std::ofstream::out);
  sim::Init();
  uint16_t avgc = 0;

  while (true) {
    {
      lock_guard<std::mutex> lock(simThreadLock);
      if (!run) {
        return;
      }
    }

    sim::Tick();
	series.write_flag(_T("Tick"));
    // auto cc = sim::Tick();

    // csv << cc << std::endl;
    // cout <<cc << std::endl;
     bb++;
    if (bb > SAMPLESIZE) {
      //	csv.close();
      //	exit(0);
    }

    /*
simTimes[avgc++] = cc;
if (avgc == 100) {
  avgc = 0;
  csv << average(simTimes, 100) << std::endl;
  bb++;
  if (bb > 5) {
    csv.close();
    exit(0);
  }
}
    */
  }
}

int main(int argc, char *argv[]) {
 // const auto n = std::chrono::steady_clock::now();
  if (argc > 1) {
    CORECOUNT = (unsigned int)std::stoi(__argv[1]);
  }
  if (argc > 2) {
    PARTICLESIZE = (unsigned int)std::stoi(__argv[2]);
  }
  if (argc > 3) {
    SAMPLESIZE = (unsigned int)std::stoi(__argv[3]);
  }


  vis::Init();
  // sim::Init();
  std::thread simThread(renderThread);
  run = true;
  uint16_t a = 0;
  while (!vis::ShouldQuit()) {
    // sim::Tick();
    vis::Start();
    /*
a++;
if (a % 500 == 0) {
  std::cout << "Stalls: " << simStalls << " : " << visStalls
            << " Avg SimTime: " << Timer::format(average(simTimes, 1000))
            << std::endl;
  a = 0;
}*/
   // const auto n2 = std::chrono::steady_clock::now();
   // if (std::chrono::duration_cast<std::chrono::seconds>(n2 - n).count() > 10) {
  //    break;
  //  }
	if (bb > 3){
	//	break;
	}
  }
  {
    std::lock_guard<std::mutex> lock(simThreadLock);
    run = false;
  }
  simThread.join();
  vis::Stop();
  vis::Shutdown();

  return 0;
}
