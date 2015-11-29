#include <stdio.h>
#include <stdlib.h>
#include <thread>

#include <vector>
#include <algorithm>
#include <iostream>
#include <chrono>
#include <string>
#include "Renderer.h"
#include "nbody.h"
#include "timer.h"

using namespace std;

uint16_t simStalls = 0;
uint16_t visStalls = 0;
mutex simThreadLock;
bool run;

long long simTimes[1000];


template <typename T> T average(T t[], int n) {
  T s = t[n - 1];
  for (int i = 0; i < (n - 1); i++)
    s += t[i];
  return s / n;
}

void renderThread() {
  sim::Init();
  uint16_t avgc = 0;

  while (true) {
    {
      lock_guard<std::mutex> lock(simThreadLock);
      if (!run) {
        return;
      }
    }
  simTimes[avgc++] = sim::Tick();
  if (avgc == 1000){ avgc = 0; }
  }
}

int main(void) {
  vis::Init();
  sim::Init();
  //std::thread simThread(renderThread);
  run = true;
  uint16_t a = 0;
  while (!vis::ShouldQuit()) {
    sim::Tick();
    vis::Start();
    a++;
    if (a % 500 == 0) {
    std::cout << "Stalls: " << simStalls << " : " << visStalls << " Avg SimTime: " << Timer::format(average(simTimes, 1000)) << std::endl;
      a = 0;
    }
  }
  {
    std::lock_guard<std::mutex> lock(simThreadLock);
    run = false;
  }
 // simThread.join();
  vis::Stop();
  vis::Shutdown();

  return 0;
}
