#include <stdio.h>
#include <stdlib.h>

#include <vector>
#include <algorithm>
#include "Renderer.h"


int main(void)
{
  vis::Init();
  while (!vis::ShouldQuit()){
    vis::Start();
  }
  vis::Stop();
  vis::Shutdown();

  return 0;
}
