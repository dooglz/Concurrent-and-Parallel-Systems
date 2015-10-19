#include "Timer.h"
#include <fstream>
#include <iostream>
#include <algorithm>
#include <string>
#include <utility>
#include <regex>
#include <cstring>
using namespace std;
const void ResultFile::PrintToCSV(const string &filename) {
  time_t rawtime;
  time(&rawtime);
  string safefilename = filename + "_" + ctime(&rawtime) + ".csv";
  std::replace(safefilename.begin(), safefilename.end(), ' ', '_');
  std::replace(safefilename.begin(), safefilename.end(), ':', '-');
  safefilename.erase(
    std::remove(safefilename.begin(), safefilename.end(), '\n'),
    safefilename.end());
  safefilename.erase(
    std::remove(safefilename.begin(), safefilename.end(), '\r'),
    safefilename.end());
  ofstream data(safefilename, ofstream::out);

  data << "name," << name.c_str() << endl;
  data << "date," << ctime(&rawtime);
  data << SystemInfo << endl;
  for (auto a : attributes) {
    data << a.c_str() << endl;
  }
  if (headdings.size() > 0) {
    data << "Run,";
    for (auto h : headdings) {
      data << h.c_str() << ",";
    }
    data << endl;
  }
  if (averages.size() > 0) {
    data << "Averages,";
    for (auto h : averages) {
      data << h << ",";
    }
    data << endl;
  }
  if (averagePercentages.size() > 0) {
    data << "Average %,";
    for (auto h : averagePercentages) {
      data << h << ",";
    }
    data << endl;
  }
  for (size_t i = 0; i < times.size(); ++i) {
    data << i << ",";
    for (auto tt : times[i]) {
      data << tt << ",";
    }
    data << endl;
  }

  data.close();
  cout << "Printed to: " << safefilename << endl;
}


const void ResultFile::CalcAvg() {
  averages.clear();
  unsigned long long totalTime = 0;
  for (size_t i = 0; i < times[0].size(); i++) {
    long double total = 0;
    uint32_t count = 0;
    for (size_t j = 0; j < times.size(); j++) {
      total += times[j][i];
      count++;
    }
    unsigned long long avg = (unsigned long long)round(total / count);
    totalTime += avg;
    averages.push_back((unsigned long long)round(total / count));
  }
  averagePercentages.clear();
  for (auto a : averages) {
    averagePercentages.push_back((float)a * 100 / (float)totalTime);
  }
}

const void SysInfo::Print() const {
  cout << "Optimised Code:\t" << prog_optimisation << endl;
  cout << "Debug Mode:\t" << prog_debugMode << endl;
  cout << "Cpu Vendor:\t" << cpu_vendor << endl;
  cout << "Cpu Name:\t" << cpu_name << endl;
  cout << "Cpu Loigcal:\t" << cpu_logical << endl;
  cout << "Cpu Cores:\t" << cpu_cores << endl;
  cout << "Cpu HWC:\t" << cpu_hardware_concurrency << endl;
  cout << "Cpu HYPTH:\t" << (cpu_hyperThreaded ? "true" : "false") << endl;
}
const string SysInfo::toString() const {
  return ("Optimised Code," +
    (prog_optimisation ? std::string("true") : std::string("false")) +
    "\nDebug Mode," +
    (prog_debugMode ? std::string("true") : std::string("false")) +
    "\nCpu Vendor," + cpu_vendor + "\nCpu Name," + cpu_name +
    "\nCpu Loigcal," + to_string(cpu_logical) + "\nCpu Cores," +
    to_string(cpu_cores) + "\nCpu HYPTH," +
    (cpu_hyperThreaded ? "true" : "false") + "\nCpu HWC," +
    to_string(cpu_hardware_concurrency));
}

ostream &operator<<(std::ostream &os, const SysInfo &obj) {
  return os << obj.toString();
}

void cpuID(unsigned i, unsigned regs[4]) {
#ifdef _WIN32
  __cpuid((int *)regs, (int)i);

#else
  asm volatile("cpuid"
    : "=a"(regs[0]), "=b"(regs[1]), "=c"(regs[2]), "=d"(regs[3])
    : "a"(i), "c"(0));
  // ECX is set to zero for CPUID function 4
#endif
}

// https://stackoverflow.com/questions/150355/

const uint16_t SysInfo::get_cpu_cores()const {
  uint16_t cpu_cores = get_cpu_logical();
  unsigned int regs[4];
  if (get_cpu_vendor() == "GenuineIntel") {
    // Get DCP cache info
    cpuID(4, regs);
    cpu_cores = ((regs[0] >> 26) & 0x1f) + 1; // EAX[31:26] + 1

  }
  else if (get_cpu_vendor() == "AuthenticAMD") {
    // Get NC: Number of CPU cores - 1
    cpuID(0x80000008, regs);
    cpu_cores = ((unsigned)(regs[2] & 0xff)) + 1; // ECX[7:0] + 1
  }
  return cpu_cores;
}

const uint16_t SysInfo::get_cpu_logical()const {
  unsigned int regs[4];
  cpuID(1, regs);
  return (regs[1] >> 16) & 0xff; // EBX[23:16]
}

const string SysInfo::get_cpu_vendor()const {
  // Get vendor
  char vendor[12];
  unsigned int regs[4];
  cpuID(0, regs);
  ((unsigned *)vendor)[0] = regs[1]; // EBX
  ((unsigned *)vendor)[1] = regs[3]; // EDX
  ((unsigned *)vendor)[2] = regs[2]; // ECX
  return string(vendor, 12);
}

const string SysInfo::get_cpu_name()const {
  unsigned int regs[4];
  char CPUBrandString[64];
  cpuID(0x80000000, regs);
  unsigned int nExIds = regs[0];
  memset(CPUBrandString, 0, sizeof(CPUBrandString));
  for (unsigned int i = 0x80000000; i <= nExIds; ++i) {
    cpuID(i, regs);
    if (i == 0x80000002)
      memcpy(CPUBrandString, regs, sizeof(regs));
    else if (i == 0x80000003)
      memcpy(CPUBrandString + 16, regs, sizeof(regs));
    else if (i == 0x80000004)
      memcpy(CPUBrandString + 32, regs, sizeof(regs));
  }
  string cpu_name = string(CPUBrandString);
  // removing leading, trailing and extra spaces
  cpu_name = regex_replace(cpu_name, regex("^ +| +$|( ) +"), string("$1"));
  return cpu_name;

}
const bool SysInfo::get_cpu_hyperThreaded()const {
  unsigned int regs[4];
  cpuID(1, regs);
  unsigned cpuFeatures = regs[3]; // EDX
  return  cpuFeatures & (1 << 28) && cpu_cores < cpu_logical;
}

SysInfo::SysInfo()
  : cpu_cores(get_cpu_cores()), cpu_logical(get_cpu_logical()),
  cpu_hardware_concurrency(std::thread::hardware_concurrency()),
  cpu_name(get_cpu_name()), cpu_vendor(get_cpu_vendor()),
  cpu_hyperThreaded(get_cpu_hyperThreaded()) {
}