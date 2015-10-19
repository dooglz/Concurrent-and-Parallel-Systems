#include "sequential.h"
#include "Timer.h"
#include <iostream>

// Defines the data size to work on
const unsigned int SIZE = 1000;
// Used to validate the result.  This is related to the data size
const double CHECK_VALUE = 12.0;

int main(int argc, char **argv) {
  SystemInfo.Print();
  seq::start(10);

  return 0;
}