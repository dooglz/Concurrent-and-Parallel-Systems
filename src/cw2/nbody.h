#pragma once
#include <glm/glm.hpp>

struct Body {
  glm::vec3 pos, speed;
  unsigned char r, g, b; // Color
  unsigned int code;
  Body(float x, float y, float z);
  Body(glm::vec3 pos);
  Body();
};

namespace sim {
long long  Tick();
void Init();
}