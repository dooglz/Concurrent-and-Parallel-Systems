#pragma once
#include <glm/glm.hpp>

namespace vis{

  struct RenderParticle{
    glm::vec3 pos;
    unsigned char r, g, b; // Color
  };


  void Init();
  void Start();
  void Stop();
  void Shutdown();
  bool ShouldQuit();
}