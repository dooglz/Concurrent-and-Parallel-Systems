#include <stdio.h>
#include <stdlib.h>
#include <thread>

#include <vector>
#include <algorithm>
#include <iostream>
#include "Renderer.h"
#include "nbody.h"

size_t simStalls = 0;
size_t visStalls = 0;
std::mutex simThreadLock;
bool run;

void renderThread() {
  sim::Init();

  while (true) {
    {
      std::lock_guard<std::mutex> lock(simThreadLock);
      if (!run) {
        return;
      }
    }
    sim::Tick();
  }
}

int main(void) {
  vis::Init();
  std::thread simThread(renderThread);
  run = true;
  size_t a = 0;
  while (!vis::ShouldQuit()) {
    vis::Start();
    a++;
    if (a % 1000 == 0) {
      std::cout << simStalls << " : " << visStalls << std::endl;
      a = 0;
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
