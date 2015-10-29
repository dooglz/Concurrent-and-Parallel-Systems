#include "sequential.h"
#include "sequentialOMP.h"
#include "parallel.h"
#include "Timer.h"
#include <iostream>
#include "ezOptionParser.hpp"

// Defines the data size to work on
const unsigned int SIZE = 1000;
// Used to validate the result.  This is related to the data size
const double CHECK_VALUE = 12.0;


void Usage(ez::ezOptionParser &opt) {
	std::string usage;
	opt.getUsage(usage);
	std::cout << usage;
};

int main(int argc, const char *argv[]) {
  SystemInfo.Print();


  ez::ezOptionParser opt;

  opt.overview = "Demo of automatic usage message creation.";
  opt.syntax = "usage [OPTIONS]";
  opt.example = "usage -h\n\n";
  opt.footer = "Sam Serrels 2015\n";

  opt.add("", 0, 0, 0, "Display usage", "-h", "-help", "--help", "--usage");
  opt.add("0", 0, 1, 0, "Thread Limit", "-t", "-threads", "--threads");
  opt.add("0", 0, 1, 0, "Iterations", "-i", "-iterations", "--iterations");
  opt.add("", 0, 0, 0, "Enable Simd128", "-simd128");
  opt.add("", 0, 0, 0, "Enable Simd256", "-simd256");
 

  opt.parse(argc, argv);

  if (opt.isSet("-h")) {
	  Usage(opt);
	  return 0;
  }
  std::vector<std::string> badOptions;

  if (!opt.gotRequired(badOptions)) {
	  for (size_t i = 0; i < badOptions.size(); ++i) {
		  std::cerr << "ERROR: Missing required option " << badOptions[i] << ".\n\n";
	  }
	  Usage(opt);
	  return 1;
  }

  if (!opt.gotExpected(badOptions)) {
	  for (size_t i = 0; i < badOptions.size(); ++i) {
		  std::cerr << "ERROR: Got unexpected number of arguments for option " << badOptions[i]
			  << ".\n\n";
	  }
	  Usage(opt);
	  return 1;
  }


  int cores = 1;
  if (opt.isSet("-t")){
	  opt.get("-t")->getInt(cores);
  }
  int iterations = 100;
  if (opt.isSet("-i")){
	  opt.get("-i")->getInt(iterations);
  }
  seqOMP::start(iterations, cores, opt.isSet("-simd128"), opt.isSet("-simd256"));

  return 0;
}