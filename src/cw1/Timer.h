#pragma once
#include <stdint.h>
#include <vector>
#include <mutex>
#include <thread>
#include <chrono>
using namespace std;

static const char Spinner(const unsigned int t) {
  char spinners[] = {
      '|', '/', '-', '\\',
  };
  return (spinners[t % 4]);
}

struct SysInfo {
  uint16_t cpu_cores;
  uint16_t cpu_logical;
  uint16_t cpu_hardware_concurrency;
  string cpu_vendor;
  string cpu_name;
  bool cpu_hyperThreaded;
  const void Print() const;
  const string toString() const;
  SysInfo();
};
ostream &operator<<(std::ostream &os, const SysInfo &obj);
static const SysInfo SystemInfo;

struct Timer {
  chrono::steady_clock::time_point start;
  chrono::steady_clock::time_point end;
  Timer() { Start(); }
  void Start() { start = chrono::steady_clock::now(); }
  void Stop() { end = chrono::steady_clock::now(); }
  const chrono::steady_clock::duration Duration() { return end - start; }
  unsigned long long Duration_NS() {
    return chrono::duration_cast<chrono::nanoseconds>(Duration()).count();
  };
};

struct ResultFile {
  string name;
  vector<string> attributes;
  vector<string> headdings;
  vector<vector<unsigned long long>> times;
  vector<unsigned long long> averages;
  vector<float> averagePercentages;
  const void CalcAvg();
  const void PrintToCSV(const string &filename);
};