#pragma once
#include <mutex>
#include <glm/glm.hpp>

#define RENDERBUFFERSIZE 3
#define PARTICLESIZE 1000


namespace vis{

  enum rbState
  {
    Rendering,UPDATING,READYTORENDER,READYTOUPDATE
  };

  struct RenderParticle{
    glm::vec3 pos;
    unsigned char r, g, b; // Color
  };

  extern RenderParticle* renderBuffer;
  extern rbState* renderBufferStates;
  extern std::mutex renderBufferStateLocks;

  void Init();
  void Start();
  void Stop();
  void Shutdown();
  bool ShouldQuit();
}