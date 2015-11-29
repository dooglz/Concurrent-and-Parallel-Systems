#include <stdio.h>
#include <stdlib.h>
#include <vector>
#include <algorithm>

#include <GL/glew.h>

#define glfw3dll
#include <GLFW/glfw3.h>

#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtx/norm.hpp>
#include <thread>
#include <iostream>

#include "shader.hpp"
#include "Renderer.h"

using namespace vis;
using namespace glm;

GLFWwindow *window;

GLfloat *g_particule_position_size_data;
GLubyte *g_particule_color_data;
// The VBO containing the positions and sizes of the particles
GLuint particles_color_buffer;
// The VBO containing the colors of the particles
GLuint particles_position_buffer;
GLuint billboard_vertex_buffer;
GLuint programID;

GLuint CameraRight_worldspace_ID;
GLuint CameraUp_worldspace_ID;
GLuint ViewProjMatrixID;
GLuint VertexArrayID;

RenderParticle *vis::renderBuffer;
rbState *vis::renderBufferStates;
std::mutex vis::renderBufferStateLocks;

const int MaxParticles = PARTICLESIZE;

double lastTime;
glm::mat4 ViewProjectionMatrix;
glm::vec3 vis::CameraPosition;
glm::mat4 ViewMatrix;
glm::mat4 ProjectionMatrix;

void vis::Init() {
  // Initialise GLFW
  if (!glfwInit()) {
    fprintf(stderr, "Failed to initialize GLFW\n");
    getchar();
    return;
  }

  glfwWindowHint(GLFW_SAMPLES, 4);
  glfwWindowHint(GLFW_RESIZABLE, GL_FALSE);
  glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
  glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
  glfwWindowHint(GLFW_OPENGL_FORWARD_COMPAT,
                 GL_TRUE); // To make MacOS happy; should not be needed
  glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);

  // vsync

  // Open a window and create its OpenGL context
  window = glfwCreateWindow(1024, 1024, "Tutorial 18 - Particules", NULL, NULL);
  if (window == NULL) {
    fprintf(stderr, "Failed to open GLFW window. If you have an Intel GPU, "
                    "they are not 3.3 compatible. Try the 2.1 version of the "
                    "tutorials.\n");
    getchar();
    glfwTerminate();
    return;
  }
  glfwMakeContextCurrent(window);

  // Initialize GLEW
  glewExperimental = true; // Needed for core profile
  if (glewInit() != GLEW_OK) {
    fprintf(stderr, "Failed to initialize GLEW\n");
    getchar();
    glfwTerminate();
    return;
  }

  // Ensure we can capture the escape key being pressed below
  glfwSetInputMode(window, GLFW_STICKY_KEYS, GL_TRUE);
  // Hide the mouse and enable unlimited mouvement
  glfwSetInputMode(window, GLFW_CURSOR, GLFW_CURSOR_NORMAL);

  // Set the mouse at the center of the screen
  glfwPollEvents();

  //glfwSwapInterval(1);

  glClearColor(0.0f, 0.0f, 0.0f, 0.0f);

  // Enable depth test
  glEnable(GL_DEPTH_TEST);
  // Accept fragment if it closer to the camera than the former one
  glDepthFunc(GL_LESS);

  VertexArrayID;
  glGenVertexArrays(1, &VertexArrayID);
  glBindVertexArray(VertexArrayID);

  // Create and compile our GLSL program from the shaders
  programID = LoadShaders("Particle.vertexshader", "Particle.fragmentshader");

  // Vertex shader
  CameraRight_worldspace_ID =
      glGetUniformLocation(programID, "CameraRight_worldspace");
  CameraUp_worldspace_ID =
      glGetUniformLocation(programID, "CameraUp_worldspace");
  ViewProjMatrixID = glGetUniformLocation(programID, "VP");

  g_particule_position_size_data = new GLfloat[MaxParticles * 3];
  g_particule_color_data = new GLubyte[MaxParticles * 3];

  // The VBO containing the 4 vertices of the particles.
  // Thanks to instancing, they will be shared by all particles.
  static const GLfloat g_vertex_buffer_data[] = {
      -0.5f, -0.5f, 0.0f, 0.5f, -0.5f, 0.0f,
      -0.5f, 0.5f,  0.0f, 0.5f, 0.5f,  0.0f,
  };
  billboard_vertex_buffer;
  glGenBuffers(1, &billboard_vertex_buffer);
  glBindBuffer(GL_ARRAY_BUFFER, billboard_vertex_buffer);
  glBufferData(GL_ARRAY_BUFFER, sizeof(g_vertex_buffer_data),
               g_vertex_buffer_data, GL_STATIC_DRAW);

  // The VBO containing the positions and sizes of the particles
  particles_position_buffer;
  glGenBuffers(1, &particles_position_buffer);
  glBindBuffer(GL_ARRAY_BUFFER, particles_position_buffer);
  // Initialize with empty (NULL) buffer : it will be updated later, each frame.
  glBufferData(GL_ARRAY_BUFFER, MaxParticles * 3 * sizeof(GLfloat), NULL,
               GL_STREAM_DRAW);

  // The VBO containing the colors of the particles
  particles_color_buffer;
  glGenBuffers(1, &particles_color_buffer);
  glBindBuffer(GL_ARRAY_BUFFER, particles_color_buffer);
  // Initialize with empty (NULL) buffer : it will be updated later, each frame.
  glBufferData(GL_ARRAY_BUFFER, MaxParticles * 3 * sizeof(GLubyte), NULL,
               GL_STREAM_DRAW);

  lastTime = glfwGetTime();

  ProjectionMatrix =
      // glm::perspective(glm::radians(95.0f), 1.0f, 0.1f, 1000.0f);
      glm::ortho(-1000.0, 1000.0, -1000.0, 1000.0, 0.1, 1000.0);
  ViewMatrix = glm::lookAt(glm::vec3(0, 0, 1000.0), glm::vec3(0, 0, -1),
                           glm::vec3(0, 1, 0));

  CameraPosition = vec3(glm::inverse(ViewMatrix)[3]);

  ViewProjectionMatrix = ProjectionMatrix * ViewMatrix;

  renderBuffer = new RenderParticle[RENDERBUFFERSIZE * PARTICLESIZE];
  renderBufferStates = new rbState[RENDERBUFFERSIZE];
  for (size_t i = 0; i < RENDERBUFFERSIZE; i++) {
    renderBufferStates[i] = READYTOUPDATE;
  }
}

void vis::Start() {

  int offset = -1;
  while (offset == -1) {
    { // scope for lock
      std::lock_guard<std::mutex> lock(vis::renderBufferStateLocks);
      for (size_t i = 0; i < RENDERBUFFERSIZE; i++) {
        if (vis::renderBufferStates[i] == vis::READYTORENDER) {
          offset = i;
          vis::renderBufferStates[i] = vis::RENDERING;
          break;
        }
      }
    }
    if (offset == -1) {
      visStalls++;
      std::this_thread::sleep_for(std::chrono::milliseconds(2));
    }
  }
  const size_t addrOffset = offset * PARTICLESIZE;

  // Simulate all particles
  int ParticlesCount = 0;
  for (int i = 0; i < MaxParticles; i++) {
    RenderParticle &p = renderBuffer[addrOffset + i]; // shortcut

    // Fill the GPU buffer
    g_particule_position_size_data[3 * ParticlesCount + 0] = p.pos.x;
    g_particule_position_size_data[3 * ParticlesCount + 1] = p.pos.y;
    g_particule_position_size_data[3 * ParticlesCount + 2] = p.pos.z;

    g_particule_color_data[3 * ParticlesCount + 0] = p.r;
    g_particule_color_data[3 * ParticlesCount + 1] = p.g;
    g_particule_color_data[3 * ParticlesCount + 2] = p.b;
    ParticlesCount++;
  }

  // Clear the screen
  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

  // printf("%d ",ParticlesCount);

  // Update the buffers that OpenGL uses for rendering.
  // There are much more sophisticated means to stream data from the CPU to the
  // GPU,
  // but this is outside the scope of this tutorial.
  // http://www.opengl.org/wiki/Buffer_Object_Streaming

  glBindBuffer(GL_ARRAY_BUFFER, particles_position_buffer);
  glBufferData(GL_ARRAY_BUFFER, MaxParticles * 3 * sizeof(GLfloat), NULL,
               GL_STREAM_DRAW); // Buffer orphaning, a common way to improve
                                // streaming perf. See above link for details.
  glBufferSubData(GL_ARRAY_BUFFER, 0, ParticlesCount * sizeof(GLfloat) * 3,
                  g_particule_position_size_data);

  glBindBuffer(GL_ARRAY_BUFFER, particles_color_buffer);
  glBufferData(GL_ARRAY_BUFFER, MaxParticles * 3 * sizeof(GLubyte), NULL,
               GL_STREAM_DRAW); // Buffer orphaning, a common way to improve
                                // streaming perf. See above link for details.
  glBufferSubData(GL_ARRAY_BUFFER, 0, ParticlesCount * sizeof(GLubyte) * 3,
                  g_particule_color_data);

  glEnable(GL_BLEND);
  glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

  // Use our shader
  glUseProgram(programID);

  // Same as the billboards tutorial
  glUniform3f(CameraRight_worldspace_ID, ViewMatrix[0][0], ViewMatrix[1][0],
              ViewMatrix[2][0]);
  glUniform3f(CameraUp_worldspace_ID, ViewMatrix[0][1], ViewMatrix[1][1],
              ViewMatrix[2][1]);

  glUniformMatrix4fv(ViewProjMatrixID, 1, GL_FALSE,
                     &ViewProjectionMatrix[0][0]);

  // 1rst attribute buffer : vertices
  glEnableVertexAttribArray(0);
  glBindBuffer(GL_ARRAY_BUFFER, billboard_vertex_buffer);
  glVertexAttribPointer(0, // attribute. No particular reason for 0, but must
                           // match the layout in the shader.
                        3, // size
                        GL_FLOAT, // type
                        GL_FALSE, // normalized?
                        0,        // stride
                        (void *)0 // array buffer offset
                        );

  // 2nd attribute buffer : positions of particles' centers
  glEnableVertexAttribArray(1);
  glBindBuffer(GL_ARRAY_BUFFER, particles_position_buffer);
  glVertexAttribPointer(1, // attribute. No particular reason for 1, but must
                           // match the layout in the shader.
                        3, // size : x + y + z + size => 4
                        GL_FLOAT, // type
                        GL_FALSE, // normalized?
                        0,        // stride
                        (void *)0 // array buffer offset
                        );

  // 3rd attribute buffer : particles' colors
  glEnableVertexAttribArray(2);
  glBindBuffer(GL_ARRAY_BUFFER, particles_color_buffer);
  glVertexAttribPointer(2, // attribute. No particular reason for 1, but must
                           // match the layout in the shader.
                        3, // size : r + g + b => 3
                        GL_UNSIGNED_BYTE, // type
                        GL_TRUE,  // normalized?    *** YES, this means that the
                                  // unsigned char[4] will be accessible with a
                                  // vec4 (floats) in the shader ***
                        0,        // stride
                        (void *)0 // array buffer offset
                        );

  // These functions are specific to glDrawArrays*Instanced*.
  // The first parameter is the attribute buffer we're talking about.
  // The second parameter is the "rate at which generic vertex attributes
  // advance when rendering multiple instances"
  // http://www.opengl.org/sdk/docs/man/xhtml/glVertexAttribDivisor.xml
  glVertexAttribDivisor(
      0, 0); // particles vertices : always reuse the same 4 vertices -> 0
  glVertexAttribDivisor(1, 1); // positions : one per quad (its center) -> 1
  glVertexAttribDivisor(2, 1); // color : one per quad -> 1

  // Draw the particules !
  // This draws many times a small triangle_strip (which looks like a quad).
  // This is equivalent to :
  // for(i in ParticlesCount) : glDrawArrays(GL_TRIANGLE_STRIP, 0, 4),
  // but faster.
  glDrawArraysInstanced(GL_TRIANGLE_STRIP, 0, 4, ParticlesCount);

  glDisableVertexAttribArray(0);
  glDisableVertexAttribArray(1);
  glDisableVertexAttribArray(2);

  // Swap buffers
  glfwSwapBuffers(window);
  glfwPollEvents();

  { // scope for lock
    std::lock_guard<std::mutex> lock(vis::renderBufferStateLocks);
    vis::renderBufferStates[offset] = vis::READYTOUPDATE;
  }
}
void vis::Stop() {}
void vis::Shutdown() {
  delete[] g_particule_position_size_data;

  // Cleanup VBO and shader
  glDeleteBuffers(1, &particles_color_buffer);
  glDeleteBuffers(1, &particles_position_buffer);
  glDeleteBuffers(1, &billboard_vertex_buffer);
  glDeleteProgram(programID);
  glDeleteVertexArrays(1, &VertexArrayID);

  // Close OpenGL window and terminate GLFW
  glfwTerminate();

  delete[] renderBuffer;
  delete[] renderBufferStates;
}

bool vis::ShouldQuit() {
  return ((glfwGetKey(window, GLFW_KEY_ESCAPE) == GLFW_PRESS) ||
          (glfwWindowShouldClose(window) == 1));
}