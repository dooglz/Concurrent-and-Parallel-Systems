#pragma once
#include "nbody.h"
#include <string>
#include <glm/glm.hpp>
/*
using namespace glm;

// Paul Nettle

// (http://www.flipcode.com/askmid/)

class Octree;
typedef bool (*callback)(const Octree &o, void *data);

// -----------------------------------------------------------------------------
// This defines a cubic bounding volume (center, radius)
// -----------------------------------------------------------------------------

struct Bounds {
  glm::vec3 center; // Center of a cubic bounding volume
  float radius;     // Radius of a cubic bounding volume
  Bounds(vec3 p, float r){
    center = p;
    radius = r;
  }
  Bounds(float x, float y, float z, float r) :Bounds(vec3(x,y,z),r){}
  Bounds() :Bounds(0.0,0.0,0.0,0.0){}
};


class Octree {
public:
  Octree();
  virtual ~Octree();

  // Accessors
  Body** Bodys() const { return _Bodys; }
  const uint32_t BodyCount() const { return _BodyCount; }

  // Implementation
  const void build(Body **Bodys, const unsigned int count,
    const unsigned int maximumDepth,
    const Bounds &bounds,
    const unsigned int currentDepth = 0);

  const void Add(Body* b);

protected:
  bool _activeChilderen[8];
  Octree *_childeren[8];
  uint32_t _BodyCount;
  Body **_Bodys;
  Bounds _bounds;
};


Octree::Octree()
  : _BodyCount(0), _Bodys(0), _bounds(0.0, 0.0, 0.0, 0.0) {
  memset(_childeren, 0, sizeof(_childeren));
}

*/
/*
class Octree {
public:
  Octree();
  virtual ~Octree();

  // Accessors
  inline const Body *const *Bodys() const { return _Bodys; }
  inline const unsigned int BodyCount() const { return _BodyCount; }

  // Implementation
  virtual const bool build(Body **Bodys, const unsigned int count,
                           const unsigned int maximumDepth,
                           const Bounds &bounds,
                           const unsigned int currentDepth = 0);
  virtual const bool traverse(callback proc, void *data) const;

protected:
  Octree *_child[8];
  unsigned int _BodyCount;
  Body **_Bodys;
  Body _center;
  double _radius;
};

// -----------------------------------------------------------------------------
// Construction -- Just "nullify" the class
// -----------------------------------------------------------------------------

Octree::Octree()
    : _BodyCount(0), _Bodys(0), _center(0.0, 0.0, 0.0), _radius(0.0) {
  memset(_child, 0, sizeof(_child));
}

Octree::~Octree() { delete[] _Bodys; }

// -----------------------------------------------------------------------------
// Build the octree
// -----------------------------------------------------------------------------

const bool Octree::build(Body **Bodys, const unsigned int count,
                         const unsigned int maximumDepth, const Bounds &bounds,
                         const unsigned int currentDepth) {
  // You know you're a leaf when...
  //
  // 1. The number of Bodys is <= 1
  // 2. We've recursed too deep into the tree
  //    (currentDepth >= maximumDepth)
  //
  //    NOTE: We specifically use ">=" for the depth comparison so that we
  //          can set the maximumDepth depth to 0 if we want a tree with
  //          no depth.

  if (count <= 1 || currentDepth >= maximumDepth) {
    // Just store the Bodys in the node, making it a leaf

    _BodyCount = count;
    _Bodys = new Body *[count];
    memcpy(_Bodys, Bodys, sizeof(Body *) * count);
    return true;
  }

  // We'll need this (see comment further down in this source)

  unsigned int childBodyCounts[8];
  memset(childBodyCounts, 0, 8 * sizeof(unsigned int));

  // Classify each Body to a child node

  for (unsigned int i = 0; i < count; i++) {
    // Current Body

    Body &p = *Bodys[i];

    // Here, we need to know which child each Body belongs to. To
    // do this, we build an index into the _child[] array using the
    // relative position of the Body to the center of the current
    // node

    p.code = 0;
    if (p.pos.x > bounds.center.x)
      p.code |= 1;
    if (p.pos.y > bounds.center.y)
      p.code |= 2;
    if (p.pos.z > bounds.center.z)
      p.code |= 4;

    // We'll need to keep track of how many Bodys get stuck in each child so
    // we'll just keep track of it here, since we have the information handy.

    childBodyCounts[p.code]++;
  }

  // Recursively call build() for each of the 8 children

  for (size_t i = 0; i < 8; i++) {
    // Don't bother going any further if there aren't any Bodys for this child

    if (!childBodyCounts[i])
      continue;

    // Allocate the child

    _child[i] = new Octree;

    // Allocate a list of Bodys that were coded JUST for this child only

    Body **newList = new Body *[childBodyCounts[i]];

    // Go through the input list of Bodys and copy over the Bodysthat were coded
    // for this child

    Body **ptr = newList;

    for (unsigned int j = 0; j < count; j++) {
      if (Bodys[j]->code == i) {
        *ptr = Bodys[j];
        ptr++;
      }
    }

    // At this Body, we have a list of Bodys that will belong to
    // this child node.
    // Generate a new bounding volume -- We do this with a touch of
    // trickery...
    //
    // We use a table of offsets. These offsets determine where a
    // node is, relative to it's parent. So, for example, if want to
    // generate the bottom-left-rear (-x, -y, -z) child for a node,
    // we use (-1, -1, -1).
    //
    // However, since the radius of a child is always half of its
    // parent's, we use a table of 0.5, rather than 1.0.
    //
    // These values are stored the following table. Note that this
    // won't compile because it assumes Bodys are structs, but you
    // get the idea.

    glm::vec3 boundsOffsetTable[8] = {{-0.5, -0.5, -0.5},
                                      {+0.5, -0.5, -0.5},
                                      {-0.5, +0.5, -0.5},
                                      {+0.5, +0.5, -0.5},
                                      {-0.5, -0.5, +0.5},
                                      {+0.5, -0.5, +0.5},
                                      {-0.5, +0.5, +0.5},
                                      {+0.5, +0.5, +0.5}};

    // Calculate our offset from the center of the parent's node to
    // the center of the child's node

    glm::vec3 offset = boundsOffsetTable[i] * bounds.radius;

    // Create a new Bounds, with the center offset and half the radius

    Bounds newBounds;
    newBounds.radius = bounds.radius * 0.5;
    newBounds.center = bounds.center + offset;

    // Recurse

    _child[i]->build(newList, childBodyCounts[i], maximumDepth, newBounds, currentDepth + 1);

    // Clean up

    delete[] newList;
  }

  return true;
}

// -----------------------------------------------------------------------------
// Generic tree traversal
// -----------------------------------------------------------------------------

const bool Octree::traverse(callback proc, void *data) const {
  // Call the callback for this node (if the callback returns false, then
  // stop traversing.

  if (!proc(*this, data))
    return false;

  // If I'm a node, recursively traverse my children

  if (!_BodyCount) {
    for (unsigned int i = 0; i < 8; i++) {
      // We store incomplete trees (i.e. we're not guaranteed that a node has
      // all 8 children)
      if (!_child[i])
        continue;

      if (!_child[i]->traverse(proc, data))
        return false;
    }
  }

  return true;
}
*/