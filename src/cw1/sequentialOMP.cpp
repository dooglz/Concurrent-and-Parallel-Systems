#include <iostream>
#include <string>
#include <assert.h>
#include <omp.h>
#include "sequentialOMP.h"
#include "Timer.h"
#include <algorithm>

#define MAX_THREADS 2
#define PAR_DAXPY 0
#define SIMD_DAXPY256 0
#define SIMD_DAXPY128 0
#define SIMD_L_Loop 1
using namespace std;
namespace seqOMP {

void (*cdaxpy)(unsigned int, const double, double *, double *, unsigned int);
int(*cgaussian)(double **, int, int *);
// Fills A with random doubles, ppopulates b with row sums, returns largest val
double fillArray(double **a, int n, double *b) {
  double largestValue = 0.0;
  int init = 1325;

  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < n; ++j) {
      init = 3125 * init % 65536; // cheap and nasty random generator
      a[j][i] = (static_cast<double>(init) - 32768.0) / 16384.0;
      largestValue = (a[j][i] > largestValue) ? a[j][i] : largestValue;
    }
  }
  // fill b with 0
  for (int i = 0; i < n; ++i) {
    b[i] = 0.0;
  }
  // add every element of each row of A to each row of B
  for (int j = 0; j < n; ++j) {
    for (int i = 0; i < n; ++i) {
      b[i] += a[j][i];
    }
  }

  return largestValue;
}

/* Purpose : Find largest component of double vector dx
n    : number of elements in input vector
dx   : double vector with n+1 elements, dx[0] is not used
dx_off : offset in reading dx
*/
int indexOfLargestElement(int n, double *dx, int dx_off) {
  double dmax, dtemp;
  int itemp = 0;

  if (n < 1) {
    itemp = -1;
  } else if (n == 1) {
    itemp = 0;
  } else {
    itemp = 0;
    dmax = abs(dx[0 + dx_off]);
    for (int i = 0; i < n; ++i) {
      dtemp = abs(dx[i + dx_off]);
      if (dtemp > dmax) {
        itemp = i;
        dmax = dtemp;
      }
    }
  }

  return itemp;
}

// Scales a vector by a constant
void scaleVecByConstant(int n, double da, double *dx, int dx_off, int incx) {
  if (n > 0) {
    if (incx != 1) {
      int nincx = n * incx;
      for (int i = 0; i < nincx; i += incx)
        dx[i + dx_off] *= da;
    } else {
      for (int i = 0; i < n; ++i)
        dx[i + dx_off] *= da;
    }
  }
}

/* Constant times a vector plus a vector
Purpose : To compute dy = da * dx + dy
--- Input ---
n       : number of elements in input vector(s)
scaler  : double scalar multiplier
dx      : double vector with n+1 elements
dy      : double vector with n+1 element
--- Output ---
dy = da * dx + dy, unchanged if n <= 0
*/
void daxpyS128(unsigned int n, const double scaler, double *dx, double *dy,
               unsigned int offset) {
  if ((n <= 0) || (scaler == 0)) {
    return;
  }
  double *const y = &dy[offset];
  double *const x = &dx[offset];
  const __m128d scalers = _mm_set1_pd(scaler);
  const int remainder = n % 2;
  const int nm1 = n - 1;

  for (int i = 0; i < nm1; i += 2) {
    // load X
    const __m128d xs = _mm_loadu_pd(&x[i]);
    // load y
    __m128d ys = _mm_loadu_pd(&y[i]);
    // mutliply X by scalers, add to Y
    ys = _mm_add_pd(ys, _mm_mul_pd(xs, scalers));
    // load back into y
    _mm_storeu_pd(&y[i], ys);
  }

  if (remainder != 0) {
    y[n - 1] += scaler * x[n - 1];
  }
}
void daxpyS256(unsigned int n, const double scaler, double *dx, double *dy,
               unsigned int offset) {
  if ((n <= 0) || (scaler == 0)) {
    return;
  }
  double *const y = &dy[offset];
  double *const x = &dx[offset];

  const __m256d scalers = _mm256_set1_pd(scaler);
  const int remainder = n % 4;
  const int nm1 = n - 3;

  for (int i = 0; i < nm1; i += 4) {
    // load X
    const __m256d xs = _mm256_loadu_pd(&x[i]);
    // load y
    __m256d ys = _mm256_loadu_pd(&y[i]);
    // mutliply X by scalers, add to Y
    ys = _mm256_add_pd(ys, _mm256_mul_pd(xs, scalers));
    // load back into y
    _mm256_storeu_pd(&y[i], ys);
  }

  for (int i = n - remainder; i < n; ++i) {
    y[i] += scaler * x[i];
  }
}
void daxpyPar(unsigned int n, const double scaler, double *dx, double *dy,
              unsigned int offset) {
  if ((n <= 0) || (scaler == 0)) {
    return;
  }
  double *const y = &dy[offset];
  double *const x = &dx[offset];
#pragma omp parallel
  {
    const int thread_id = omp_get_thread_num();
    const int thread_count = omp_get_num_threads();
    const int per_thread = n / thread_count;

    for (int i = (thread_id * per_thread); i < per_thread; ++i) {
      y[i] += scaler * x[i];
    }
  }
}

void daxpy(unsigned int n, const double scaler, double *dx, double *dy,
           unsigned int offset) {

  double *const y = &dy[offset];
  double *const x = &dx[offset];
  for (int i = 0; i < n; ++i) {
    y[i] += scaler * x[i];
  }
}

// Performs Gaussian elimination with partial pivoting
int gaussian_eliminate(double **a, int n, int *ipivot) {
  // Pointers to columns being worked on
  double *col_k;
  int nm1 = n - 1;
  int info = 0;

  if (nm1 >= 0) {
    int kp1, l;
    for (int k = 0; k < nm1; ++k) {
      // Set pointer for col_k to relevant column in a
      col_k = &a[k][0];
      kp1 = k + 1;

      // Find pivot index
      l = indexOfLargestElement(n - k, col_k, k) + k;
      ipivot[k] = l;

      // Zero pivot means that this column is already triangularized
      if (col_k[l] != 0) {
        double t;
        // Check if we need to interchange
        if (l != k) {
          t = col_k[l];
          col_k[l] = col_k[k];
          col_k[k] = t;
        }

        // Compute multipliers
        t = -1.0 / col_k[k];
        //Multiply collum by t
        scaleVecByConstant(n - kp1, t, col_k, kp1, 1);

        // Row elimination with column indexing
        for (int j = kp1; j < n; ++j) {
          // Set pointer for col_j to relevant column in a
          double *col_j = &a[j][0];

          double t = col_j[l];
          if (l != k) {
            col_j[l] = col_j[k];
            col_j[k] = t;
          }
		  cdaxpy(n - kp1, t, col_k, col_j, kp1);
        }

      } else
        info = k;
    }
  }

  ipivot[n - 1] = n - 1;
  if (a[n - 1][n - 1] == 0) {
    info = n - 1;
  }

  return info;
}
// Performs Gaussian elimination with partial pivoting
int gaussian_eliminatePAR(double **a, int n, int *ipivot) {
  // Pointers to columns being worked on
  double *col_k;
  int nm1 = n - 1;
  int info = 0;

  if (nm1 >= 0) {
    int kp1, l;
    for (int k = 0; k < nm1; ++k) {
      // Set pointer for col_k to relevant column in a
      col_k = &a[k][0];
      kp1 = k + 1;

      // Find pivot index
      l = indexOfLargestElement(n - k, col_k, k) + k;
      ipivot[k] = l;

      // Zero pivot means that this column is already triangularized
      if (col_k[l] != 0) {
        double t;
        // Check if we need to interchange
        if (l != k) {
          t = col_k[l];
          col_k[l] = col_k[k];
          col_k[k] = t;
        }

        // Compute multipliers
        t = -1.0 / col_k[k];
        scaleVecByConstant(n - kp1, t, col_k, kp1, 1);

// Row elimination with column indexing
#pragma omp parallel for
        for (int j = kp1; j < n; ++j) {
          // Set pointer for col_j to relevant column in a
          double *col_j = &a[j][0];

          double t = col_j[l];
          if (l != k) {
            col_j[l] = col_j[k];
            col_j[k] = t;
          }
		  cdaxpy(n - kp1, t, col_k, col_j, kp1);
        }

      } else
        info = k;
    }
  }

  ipivot[n - 1] = n - 1;
  if (a[n - 1][n - 1] == 0) {
    info = n - 1;
  }

  return info;
}

/*
// Performs a dot product calculation of two vectors
double ddot(int n, double *dx, int dx_off, double *dy, int dy_off) {
  double temp = 0.0;
  if (n > 0) {
    for (int i = 0; i < n; ++i) {
      temp += dx[i + dx_off] * dy[i + dy_off];
    }
  }
  return temp;
}
*/
// Solves the system a * x = b using the factors computed in dgeco or
// gaussian_eliminate
void dgesl(double **a, int n, int *ipivot, double *b) {
  int k, nm1;
  nm1 = n - 1;

  // Solve a * x = b.  First solve l * y = b
  if (nm1 >= 1) {
    for (k = 0; k < nm1; ++k) {

      int l = ipivot[k];
      double t = b[l];

      if (l != k) {

        b[l] = b[k];
        b[k] = t;
      }

	  cdaxpy(n - (k + 1), t, &a[k][0], b, (k + 1));
    }
  }

  // Now solve u * x = y
  for (int kb = 0; kb < n; ++kb) {
    k = n - (kb + 1);
    b[k] /= a[k][k];
    double t = -b[k];
	cdaxpy(k, t, &a[k][0], b, 0);
  }
}

// Multiply matrix m times vector x and add the result to vector y
void dmxpy(int n1, double *y, int n2, double *x, double **m) {
  for (int j = 0; j < n2; ++j) {
    for (int i = 0; i < n1; ++i) {
      y[i] += x[j] * m[j][i];
    }
  }
}

// Runs the benchmark
void run(double **a, double *b, int n, int *ipivot) {}

// Validates the result
void validate(double **a, double *b, double *x, int n) {
  // copy b into x
  for (int i = 0; i < n; ++i) {
    x[i] = b[i];
  }

  // reset A and B arrays to orignal rand values
  double biggestA = fillArray(a, n, b);

  for (int i = 0; i < n; ++i) {
    b[i] = -b[i];
  }

  // multipy a*x, add to b
  dmxpy(n, b, n, x, a);

  double biggestB = 0.0;
  double biggestX = 0.0;
  for (int i = 0; i < n; ++i) {
    biggestB = (biggestB > abs(b[i])) ? biggestB : abs(b[i]);
    biggestX = (biggestX > abs(x[i])) ? biggestX : abs(x[i]);
  }

  double residn =
      biggestB / (n * biggestA * biggestX * (2.2204460492503131e-016));
  assert(residn < CHECK_VALUE);
  /*
  if (residn > CHECK_VALUE) {
  assert(false);
  cout << "Validation failed!" << endl;
  cout << "Computed Norm Res = " << residn << endl;
  cout << "Reference Norm Res = " << CHECK_VALUE << endl;
  } else {
  cout << "Calculations are correct!" << endl;
  cout << "Computed Norm Res = " << residn << endl;
  cout << "Reference Norm Res = " << CHECK_VALUE << endl;
  }
  */
}

int start(const unsigned int runs, const unsigned int threadCount,
          const bool simd128, const bool simd256) {

  ResultFile r;
  r.name = "Sequential LinPack OMP" + to_string(runs);

  if (simd256) {
    cdaxpy = &daxpyS256;
    r.name += "Simd256";
  } else if (simd128) {
    cdaxpy = &daxpyS128;
    r.name += "Simd128";
  } else {
    cdaxpy = &daxpy;
  }
  if (threadCount > 1){
	  unsigned int t = min((unsigned int)omp_get_max_threads(), threadCount);
	  cgaussian = &gaussian_eliminatePAR;
	  omp_set_num_threads(t);
	  r.name += "ParLloopT" + to_string(t);
  }else{
	  cgaussian = &gaussian_eliminate;
  }
  cout << r.name<< endl;
 r.headdings = {"Allocate Memory", "Create Input Numbers",
                   " gaussian_eliminate", "Solve", "Validate"};
  Timer time_total;
  for (size_t i = 0; i < runs; i++) {
    cout << i << endl;
    // Allocate data on the heap
    Timer time_allocate;
    double **a = new double *[SIZE];
    for (size_t i = 0; i < SIZE; ++i) {
      // a[i] = new double[SIZE];
      a[i] = (double *)_aligned_malloc(SIZE * sizeof(double), sizeof(double));
    }

    double *b =
        (double *)_aligned_malloc(SIZE * sizeof(double), sizeof(double));

    double *x =
        (double *)_aligned_malloc(SIZE * sizeof(double), sizeof(double));

    int *ipivot = (int *)_aligned_malloc(SIZE * sizeof(int), sizeof(int));
    // double *b = new double[SIZE];
    // double *x = new double[SIZE];
    // int *ipivot = new int[SIZE];
    time_allocate.Stop();

    // Main application
    Timer time_genRnd;
    auto aa = fillArray(a, SIZE, b);
    time_genRnd.Stop();

    Timer time_gauss;
	cgaussian(a, SIZE, ipivot);
    time_gauss.Stop();

    Timer time_dgesl;
    dgesl(a, SIZE, ipivot, b);
    time_dgesl.Stop();

    Timer time_validate;
    validate(a, b, x, SIZE);
    time_validate.Stop();
    r.times.push_back({time_allocate.Duration_NS(), time_genRnd.Duration_NS(),
                       time_gauss.Duration_NS(), time_dgesl.Duration_NS(),
                       time_validate.Duration_NS()});
    // Free the memory
    for (size_t i = 0; i < SIZE; ++i) {
      _aligned_free(a[i]);
    }
    _aligned_free(b);
    _aligned_free(x);
    _aligned_free(ipivot);
    // delete[] b;
    // delete[] x;
    // delete[] ipivot;
  }
  r.CalcAvg();
  r.PrintToCSV(r.name);
  time_total.Stop();
  cout << "Total Time: " << Timer::format(time_total.Duration_NS());
  return 0;
}
}