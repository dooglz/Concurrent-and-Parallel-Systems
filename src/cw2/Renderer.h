#pragma once
#include <mutex>
#include <glm/glm.hpp>

#define RENDERBUFFERSIZE 3
#define PARTICLESIZE 100000

extern size_t simStalls;
extern size_t visStalls;

namespace vis {

enum rbState { RENDERING, UPDATING, READYTORENDER, READYTOUPDATE };

struct RenderParticle {
  glm::vec3 pos;
  unsigned char r, g, b; // Color
};

extern RenderParticle *renderBuffer;
extern rbState *renderBufferStates;
extern std::mutex renderBufferStateLocks;
extern glm::vec3 CameraPosition;

void Init();
void Start();
void Stop();
void Shutdown();
bool ShouldQuit();
}