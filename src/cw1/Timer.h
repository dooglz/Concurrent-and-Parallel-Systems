#pragma once
#include <stdint.h>
#include <vector>
#include <mutex>
#include <thread>
#include <chrono>

using namespace std;

#if defined(__OPTIMIZE__)
#define OPTMODE true
#else
#define OPTMODE false
#endif

#if defined(_DEBUG)
#define DBGMODE true
#else
#define DBGMODE false
#endif

struct SysInfo {
private:
  const uint16_t get_cpu_cores() const;
  const uint16_t get_cpu_logical() const;
  const string get_cpu_vendor() const;
  const string get_cpu_name() const;
  const bool get_cpu_hyperThreaded() const;

public:
  const bool prog_optimisation = OPTMODE;
  const bool prog_debugMode = DBGMODE;
  const uint16_t cpu_cores;
  const uint16_t cpu_logical;
  const uint16_t cpu_hardware_concurrency;
  const string cpu_vendor;
  const string cpu_name;
  const bool cpu_hyperThreaded;
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
  const chrono::steady_clock::duration Duration() const{ return end - start; }
  const string toString() const;
  unsigned long long Duration_NS() const {
    return chrono::duration_cast<chrono::nanoseconds>(Duration()).count();
  };
  const static string format(const unsigned long long NS);
};
ostream &operator<<(std::ostream &os, const Timer &obj);

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