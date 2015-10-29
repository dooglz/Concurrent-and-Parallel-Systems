#pragma once
extern const unsigned int SIZE;
extern const double CHECK_VALUE;

namespace seqOMP {
	int start(const unsigned int runs, const unsigned int threadCount, const bool simd128, const bool simd256);
}
