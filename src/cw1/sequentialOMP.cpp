#include <iostream>
#include <string>
#include <assert.h>
#include <omp.h>
#include "sequentialOMP.h"
#include "Timer.h"
#include <algorithm>

#define MAX_THREADS 8
#define PAR_DAXPY 0
#define SIMD_DAXPY 0
#define SIMD_L_Loop 0
#define PAR_DAXPY_CACHE 1
#define cache_line_size 8

using namespace std;

namespace seqOMP {

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

void daxpy(const unsigned int n, const double scaler, double *dx, double *dy, const unsigned int offset) {
  if ((n <= 0) || (scaler == 0)) { return; }

#if PAR_DAXPY

// omp_set_num_threads(min(n, MAX_THREADS));
#if PAR_DAXPY_CACHE
    if (offset % cache_line_size != 0) {
      // if we start looping through this, we won't be on the cache boundaries!
      // So let's caclulate the odd values in on thread first, and change the
      // offset.
      int rem = offset % cache_line_size;
      int loops_to_do = cache_line_size - rem;
      for (int i = 0; (i < loops_to_do) && (i < n); ++i) {
        int a = i + offset;
        dy[i + offset] += scaler * dx[i + offset];
      }
      offset += loops_to_do;
      n -= loops_to_do;
      if (n < 1) {
        return;
      }
    }
#endif

#pragma omp parallel
    {
      const int thread_id = omp_get_thread_num();
      const int thread_count = omp_get_num_threads();
      const int stride = thread_count * cache_line_size;
// id(thread_id == 0) {cout << " n: " << }

#if PAR_DAXPY_CACHE

      for (int i = (thread_id * cache_line_size); i < n; i += stride) {
        for (int j = 0; (j < cache_line_size) && (j + i < n); ++j) {
          dy[j + i + offset] += scaler * dx[j + i + offset];
        }
      }
#else
      for (int i = thread_id; i < n; i = i + thread_count) {
        dy[i + offset] += scaler * dx[i + offset];
      }
#endif
    }
#elif SIMD_DAXPY
//register double           x0, x1, x2, x3, y0, y1, y2, y3;

register int              i =0;
const int                 
inc2 = 2 * offset,
inc3 = 3 * offset,
inc4 = 4 * offset;

  if (((n >> 2) << 2) != 0)
  {
    const long n1 = (n & -4) - inc4;
    
    while (i < n1)
    {
      /*
      x0 = (*dx); y0 = (*dy);     
      x1 = dx[offset]; y1 = dy[offset];
      x2 = dx[inc2]; y2 = dy[inc2]; 
      x3 = dx[inc3]; y3 = dy[inc3];
      */

      //simd multiply

      __m128d dd = _mm_set1_pd(scaler);
      __m128d x0x1 = _mm_set_pd((*dx), dx[offset]);
      __m128d x2x3 = _mm_set_pd(dx[inc2], dx[inc3]);
      __m128d y0y1 = _mm_set_pd((*dy), dy[offset]);
      __m128d y2y3 = _mm_set_pd(dy[inc2], dy[inc3]);

      x0x1 = _mm_mul_pd(x0x1, dd);
      x2x3 = _mm_mul_pd(x2x3, dd);
      y0y1 = _mm_add_pd(y0y1, x0x1);
      y2y3 = _mm_add_pd(y2y3, x2x3);

      *dy = _mm_cvtsd_f64(y0y1);
      dy[offset] = _mm_cvtsd_f64(_mm_unpackhi_pd(y0y1, y0y1));
      dy[inc2] = _mm_cvtsd_f64(y2y3);
      dy[inc3] = _mm_cvtsd_f64(_mm_unpackhi_pd(y2y3, y2y3));

      //v.m2 = _mm_mul_pd(v.m2, dd);
/*
      
      *dy = y0 + alpha * x0; 
      dy[offset] = y1 + alpha * x1;
      dy[inc2] = y2 + alpha * x2; 
      dy[inc3] = y3 + alpha * x3;
      */
      dx += inc4;
      dy += inc4;
      i += inc4;
    }
  }

  while (i < n)
  {
    register double x0 = (*dx);
    register double y0 = (*dy);

    *dy = y0 + scaler * x0;

    dx += offset;
    dy += offset;
    i+= offset;
   
  }

#else
    for (size_t i = 0; i < n; ++i) {
      dy[i + offset] += scaler * dx[i + offset];
    }
#endif
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
        scaleVecByConstant(n - kp1, t, col_k, kp1, 1);

// Row elimination with column indexing
#if SIMD_L_Loop
#pragma omp parallel for
#endif
        for (int j = kp1; j < n; ++j) {
          // Set pointer for col_j to relevant column in a
          double *col_j = &a[j][0];

          double t = col_j[l];
          if (l != k) {
            col_j[l] = col_j[k];
            col_j[k] = t;
          }
          daxpy(n - kp1, t, col_k, col_j, kp1);
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

      daxpy(n - (k + 1), t, &a[k][0], b, (k + 1));
    }
  }

  // Now solve u * x = y
  for (int kb = 0; kb < n; ++kb) {
    k = n - (kb + 1);
    b[k] /= a[k][k];
    double t = -b[k];
    daxpy(k, t, &a[k][0], b, 0);
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

int start(const unsigned int runs) {
  omp_set_num_threads(MAX_THREADS);
  ResultFile r;
  r.name = "Sequential LinPack OMP" + to_string(runs);
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
    gaussian_eliminate(a, SIZE, ipivot);
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
  // r.CalcAvg();
  // r.PrintToCSV(r.name);
  time_total.Stop();
  cout << "Total Time: " << Timer::format(time_total.Duration_NS());
  return 0;
}
}