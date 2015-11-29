#include "nbody.h"
#include "Renderer.h"
#include <algorithm>
#include <glm/gtx/norm.hpp>
#include <chrono>
#include <thread>
#include <iostream>

struct Particle {
  glm::vec3 pos, speed;
  unsigned char r, g, b; // Color
  float life; // Remaining life of the particle. if <0 : dead and unused.
  float cameradistance; // *Squared* distance to the camera. if dead : -1.0f

  bool operator<(const Particle &that) const {
    // Sort in reverse order : far particles drawn first.
    return this->cameradistance > that.cameradistance;
  }
};
Particle *ParticlesContainer;

void SortParticles() {
  std::sort(&ParticlesContainer[0], &ParticlesContainer[PARTICLESIZE]);
}

void sim::Init() {
  ParticlesContainer = new Particle[PARTICLESIZE];
  for (int i = 0; i < PARTICLESIZE; i++) {
    ParticlesContainer[i].life = -1.0f;
    ParticlesContainer[i].cameradistance = -1.0f;
    ParticlesContainer[i].r = rand() % 256;
    ParticlesContainer[i].g = rand() % 256;
    ParticlesContainer[i].b = rand() % 256;
  }
}
long long sim::Tick() {
	const auto n = std::chrono::steady_clock::now();
  const float delta = 0.01;

  // Simulate all particles
  for (int i = 0; i < PARTICLESIZE; i++) {

    Particle &p = ParticlesContainer[i]; // shortcut

    if (p.life > 0.0f) {

      // Decrease life
      p.life -= delta;
      if (p.life > 0.0f) {
        // Simulate simple physics : gravity only, no collisions
        p.speed += glm::vec3(0.0f, -9.81f, 0.0f) * delta * 0.5f;
        p.pos += p.speed * delta;
        // p.cameradistance = glm::length2(p.pos - vis::CameraPosition);
      } else {
        p.cameradistance = -1.0f;
      }
    } else {
      // reset
      p.life = (rand() % 100) * 0.05f;
      p.pos = glm::vec3(0, 0, -20.0f);
      glm::vec3 maindir = glm::vec3(0.0f, 10.0f, 0.0f);
      glm::vec3 randomdir = glm::vec3((rand() % 2000 - 1000.0f) / 1000.0f,
                                      (rand() % 2000 - 1000.0f) / 1000.0f,
                                      (rand() % 2000 - 1000.0f) / 1000.0f);
      p.speed = maindir + randomdir * 1.5f;
    }
  }
  const auto n2 = std::chrono::steady_clock::now();
  // SortParticles();
  // copy

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
  for (int i = 0; i < PARTICLESIZE; i++) {
    Particle &p = ParticlesContainer[i];
    vis::renderBuffer[addrOffset + i] = {p.pos, p.r, p.g, p.b};
  }
  { // scope for lock
    std::lock_guard<std::mutex> lock(vis::renderBufferStateLocks);
    vis::renderBufferStates[offset] = vis::READYTORENDER;
  }

  return std::chrono::duration_cast<std::chrono::nanoseconds>(n2 - n).count();
}