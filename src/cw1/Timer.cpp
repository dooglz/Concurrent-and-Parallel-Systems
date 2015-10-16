#include "Timer.h"
#include <fstream>
#include <iostream>
#include <algorithm>
#include <string>
#include <utility>

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
  data << "date," << ctime(&rawtime) << endl;
  for (auto a : attributes) {
    data << a.c_str() << endl;
  }
  if (headdings.size() > 0) {
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
  for (size_t i =0; i< times.size(); ++i) {
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
  unsigned long long totalTime =0;
  for (size_t i = 0; i < times[0].size(); i++)
  {
    long double total = 0;
    uint32_t count = 0;
    for (size_t j = 0; j < times.size(); j++)
    {
      total += times[j][i];
      count++;
    }
    unsigned long long avg = round(total / count);
    totalTime += avg;
    averages.push_back(round(total / count));
  }
  averagePercentages.clear();
  for (auto a : averages)
  {
    averagePercentages.push_back((float)a*100 / (float)totalTime);
  }
}