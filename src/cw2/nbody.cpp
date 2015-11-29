#include "nbody.h"
#include "Renderer.h"
#include <algorithm>
#include <glm/gtx/norm.hpp>
#include <chrono>
#include <thread>
#include <iostream>
#include <omp.h>

using namespace glm;
#define SOFTENING 1e-9f
//#define FRICITON 0.999998f
#define FRICITON 1.0f
struct Body {
  glm::vec3 pos, speed;
  unsigned char r, g, b; // Color
};

Body *bodies;

void sim::Init() {
  bodies = new Body[PARTICLESIZE];
#pragma omp parallel for
  for (int i = 0; i < PARTICLESIZE; i++) {
    bodies[i].r = rand() % 256;
    bodies[i].g = rand() % 256;
    bodies[i].b = rand() % 256;
    bodies[i].pos = vec3((rand() % 1000) - 500, (rand() % 1000) - 500,
                         (rand() % 1000) - 500);
  }
  omp_set_num_threads(8);
}
long long sim::Tick() {
  const auto n = std::chrono::steady_clock::now();
  const float delta = 0.1;

// Simulate all particles
#pragma omp parallel for
  for (int i = 0; i < PARTICLESIZE; i++) {
    vec3 newVelo(0, 0, 0);
    for (int j = 0; j < PARTICLESIZE; j++) {
      vec3 r = bodies[j].pos - bodies[i].pos;
      float distSqr = dot(r, r) + SOFTENING;
      if (distSqr > 0.1f) {
        float invDist = 1.0f / sqrtf(distSqr);
        float invDist3 = invDist * invDist * invDist;

        newVelo += r * invDist3;
      }
    }

    bodies[i].speed += delta * newVelo * FRICITON;
    bodies[i].speed -= 0.000001f * bodies[i].pos;
    bodies[i].pos += bodies[i].speed;

    bodies[i].pos.x = clamp(bodies[i].pos.x, -1000.0f, 1000.0f);
    bodies[i].pos.y = clamp(bodies[i].pos.y, -1000.0f, 1000.0f);
    bodies[i].pos.z = clamp(bodies[i].pos.z, -1000.0f, 1000.0f);
  }

  const auto n2 = std::chrono::steady_clock::now();
  // aquire lock
  int offset = -1;
  while (offset == -1) {
    { // scope for lock
      std::lock_guard<std::mutex> lock(vis::renderBufferStateLocks);
      for (size_t i = 0; i < RENDERBUFFERSIZE; i++) {
        if (vis::renderBufferStates[i] == vis::READYTOUPDATE) {
          offset = i;
          vis::renderBufferStates[i] = vis::UPDATING;
          break;
        }
      }
    }
    if (offset == -1) {
      simStalls++;
      std::this_thread::sleep_for(std::chrono::milliseconds(1));
    }
  }

  const size_t addrOffset = offset * PARTICLESIZE;
#pragma omp parallel for
  for (int i = 0; i < PARTICLESIZE; i++) {
    Body &p = bodies[i];
    vis::renderBuffer[addrOffset + i] = {p.pos, p.r, p.g, p.b};
  }
  { // scope for lock
    std::lock_guard<std::mutex> lock(vis::renderBufferStateLocks);
    vis::renderBufferStates[offset] = vis::READYTORENDER;
  }

  return std::chrono::duration_cast<std::chrono::nanoseconds>(n2 - n).count();
}