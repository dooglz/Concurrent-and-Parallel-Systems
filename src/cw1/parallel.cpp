#include <iostream>
#include <assert.h>
#include <string>
#include "parallel.h"
#include "Timer.h"

using namespace std;

namespace par {

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

  double fillArray2(double *a, int n, double *b) {
    double largestValue = 0.0;
    int init = 1325;

    for (int i = 0; i < n; ++i) {
      for (int j = 0; j < n; ++j) {
        init = 3125 * init % 65536; // cheap and nasty random generator
        a[j*n + i] = (static_cast<double>(init) - 32768.0) / 16384.0;
        largestValue = (a[j*n + i] > largestValue) ? a[j*n + i] : largestValue;
      }
    }
    // fill b with 0
    for (int i = 0; i < n; ++i) {
      b[i] = 0.0;
    }
    // add every element of each row of A to each row of B
    for (int j = 0; j < n; ++j) {
      for (int i = 0; i < n; ++i) {
        b[i] += a[j*n + i];
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
    }
    else if (n == 1) {
      itemp = 0;
    }
    else {
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
      }
      else {
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

  void daxpy(int n, double scaler, double *dx, double *dy, int offset) {
    if ((n > 0) && (scaler != 0)) {
      for (int i = 0; i < n; ++i) {
        dy[i + offset] += scaler * dx[i + offset];
      }
    }
  }

  /*
  void msaxpy(int nr, int nc, double dx[], int offset, double scaler[], double dy[])
  {
    int i, j;
# pragma omp parallel \
  shared ( dx, nc, nr, scaler, dy ) \
  private ( i, j ) 

# pragma omp for 

    for (j = 0; j < nc; j++)
    {
      for (i = 0; i < nr; i++)
      {
        dy[i + j*offset] += dx[j*offset] * scaler[i];
      }
    }
    return;
  }
*/

  void msaxpy(int nr, double * a, int n, double * x, double * y)
  {
    int i, j;
#pragma omp parallel shared(nr,a,x,y) private(i,j) 
    /* #pragma omp for schedule(dynamic,1) nowait */
#pragma omp for 
    for (j = 0; j<nr; j++) {
      for (i = 0; i<nr; i++) {
        y[i, j] += a[j*n] * x[i];
      }
    }

  }

  // Performs Gaussian elimination with partial pivoting
  int gaussian_eliminate(double **a, int n, int *ipivot) {
    // Pointers to columns being worked on
    double *col_k, *col_j;
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
          for (int j = kp1; j < n; ++j) {
            // Set pointer for col_j to relevant column in a
            col_j = &a[j][0];

            t = col_j[l];
            if (l != k) {
              col_j[l] = col_j[k];
              col_j[k] = t;
            }
            daxpy(n - kp1, t, col_k, col_j, kp1);
          }
        }
        else
          info = k;
      }
    }

    ipivot[n - 1] = n - 1;
    if (a[n - 1][n - 1] == 0) {
      info = n - 1;
    }

    return info;
  }

  void sswap(int n, double * px, int ix, double * py, int iy)
  {
    int i, ipx, ipy;
    double t;
    ipx = 0; 
    ipy = 0;
    for (i = 0; i<n; i++) {
      t = *(px + ipx);
      *(px + ipx) = *(py + ipy);
      *(py + ipy) = t;
      ipx += ix;
      ipy += iy;
    }
  }

  void sscal(int n,double sa,double*  x, int ix)
  {
    int i, ipx;

    if (n <= 0) return;
    if (ix != 1) {
      for (i = 0; i<n; i++) {
        ipx = i*ix;
        x[ipx] *= sa;
      }
    }
    else {
      for (i = 0; i<n; i++) {
        x[i] *= sa;
      }
    }
  }

  void saxpy(int n, double a, double x[], double y[])

  {
    int i,ix,iy,m;

    if (n <= 0 || a == 0.0) {
      return;
    }
    m = n % 4;

    for (i = 0; i < m; i++) {
      y[i] = y[i] + a * x[i];
    }

    for (i = m; i < n; i = i + 4) {
      y[i] = y[i] + a * x[i];
      y[i + 1] = y[i + 1] + a * x[i + 1];
      y[i + 2] = y[i + 2] + a * x[i + 2];
      y[i + 3] = y[i + 3] + a * x[i + 3];
    }

    return;
  }

#define al(i,j) a[i+(n*j)]
  void sgesl(double *a, int n, int *ipvt, double *b)
  {
    /*
    C version of Cleve Moler's SGESL - with Fortran ordering
    and call by address convention
    */
    double t;
    int k, k1, kb, l, nk, nm1, nv;
    nv = n;
    if (nv > 0) {
      nm1 = nv - 1;
      for (k = 0; k<nm1; k++) {
        l = ipvt[k];
        t = b[l];
        if (l != k) {
          b[l] = b[k];
          b[k] = t;
        }
        nk = nv - k - 1;
        k1 = k + 1;
        saxpy(nk, t, &al(k1, k), &b[k1]);
      }
    }
    for (kb = 0; kb<nv; kb++) {
      k = nv - kb - 1;
      b[k] = b[k] / al(k, k);
      t = -b[k];
      saxpy(k, t, &al(0, k), b);
    }
  }



#define am(p,q) (a+p+q*n)
  // Performs Gaussian elimination with partial pivoting
  int gaussian_eliminate3(double *a, int n, int *ipivot) {
    int info = 0;
    int nm1 = n - 1;
    for (int k = 0; k<nm1; k++) {
      int kp1 = k + 1;
     // int l = isamax(n - k, am(k, k), 1) + k;
      int l = indexOfLargestElement(n - k, am(k, k), k) + k;
      ipivot[k] = l;
      if ((*am(l, k)) == 0.0) {
        info = k;
        return 0;
      }
      double t;
      if (l != k) {
        t = *am(l, k);
        *am(l, k) = *am(k, k);
        *am(k, k) = t;
      }
      t = -1.0 / (*am(k, k));

      sscal(n - k - 1, t, am(kp1, k), 1);
      if (l != k) {
        sswap(n - k - 1, am(l, kp1), n, am(k, kp1), n);
      }
      msaxpy(n - k - 1, am(k, kp1), n, am(kp1, k), am(kp1, kp1));
    }

    return 0;
  }


  // Solves the system a * x = b using the factors computed in dgeco or
  // gaussian_eliminate
  void dgesl(double **a, int n, int *ipivot, double *b) {
    double t;
    int k, l, nm1, kp1;

    nm1 = n - 1;

    // Solve a * x = b.  First solve l * y = b
    if (nm1 >= 1) {
      for (k = 0; k < nm1; ++k) {
        l = ipivot[k];
        t = b[l];
        if (l != k) {
          b[l] = b[k];
          b[k] = t;
        }
        kp1 = k + 1;
        daxpy(n - kp1, t, &a[k][0], b, kp1);
      }
    }

    // Now solve u * x = y
    for (int kb = 0; kb < n; ++kb) {
      k = n - (kb + 1);
      b[k] /= a[k][k];
      t = -b[k];
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
  void dmxpy2(int n, double *y,  double *x, double *m) {
    for (int j = 0; j < n; ++j) {
      for (int i = 0; i < n; ++i) {
        y[i] += x[j] * m[j*n+i];
      }
    }
  }

  // Runs the benchmark
  void run(double **a, double *b, int n, int *ipivot) {}


  void validate2(double *a, double *b, double *x, int n) {
    // copy b into x
    for (int i = 0; i < n; ++i) {
      x[i] = b[i];
    }

    // reset A and B arrays to orignal rand values
    double biggestA = fillArray2(a, n, b);

    for (int i = 0; i < n; ++i) {
      b[i] = -b[i];
    }

    // multipy a*x, add to b
    dmxpy2(n, b, x, a);

    double biggestB = 0.0;
    double biggestX = 0.0;
    for (int i = 0; i < n; ++i) {
      biggestB = (biggestB > abs(b[i])) ? biggestB : abs(b[i]);
      biggestX = (biggestX > abs(x[i])) ? biggestX : abs(x[i]);
    }

    double residn =
      biggestB / (n * biggestA * biggestX * (2.2204460492503131e-016));
    assert(residn < CHECK_VALUE);
  }

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
//    assert(residn < CHECK_VALUE);
    
    if (residn > CHECK_VALUE) {
   
    cout << "Validation failed!" << endl;
    cout << "Computed Norm Res = " << residn << endl;
    cout << "Reference Norm Res = " << CHECK_VALUE << endl;
    assert(false);
    } else {
    cout << "Calculations are correct!" << endl;
    cout << "Computed Norm Res = " << residn << endl;
    cout << "Reference Norm Res = " << CHECK_VALUE << endl;
    }
    
  }

  int start(const unsigned int runs) {
    ResultFile r;
    r.name = "Sequential LinPack" + to_string(runs);
    r.headdings = { "Allocate Memory", "Create Input Numbers",
      " gaussian_eliminate", "Solve", "Validate" };
    Timer time_total;
    for (size_t i = 0; i < runs; i++) {
      cout << i << endl;
      // Allocate data on the heap
      Timer time_allocate;
      /*
      double **a = new double *[SIZE];
      for (size_t i = 0; i < SIZE; ++i) {
        a[i] = new double[SIZE];
      }
      */
	  double *a = new double[NSIZE*NSIZE];
	  double *b = new double[NSIZE];
	  double *x = new double[NSIZE];
	  int *ipivot = new int[NSIZE];
      time_allocate.Stop();

      // Main application
      Timer time_genRnd;
	  auto aa = fillArray2(a, NSIZE, b);
      time_genRnd.Stop();

      Timer time_gauss;
	  gaussian_eliminate3(a, NSIZE, ipivot);
      time_gauss.Stop();

      Timer time_dgesl;
	  sgesl(a, NSIZE, ipivot, b);
      time_dgesl.Stop();

      Timer time_validate;
	  validate2(a, b, x, NSIZE);
      time_validate.Stop();
      r.times.push_back({ time_allocate.Duration_NS(), time_genRnd.Duration_NS(),
        time_gauss.Duration_NS(), time_dgesl.Duration_NS(),
        time_validate.Duration_NS() });
      // Free the memory
	  for (size_t i = 0; i < NSIZE; ++i) {
      //  delete[] a[i];
      }
      delete[] a;
      delete[] b;
      delete[] x;
      delete[] ipivot;
    }
   // r.CalcAvg();
   // r.PrintToCSV(r.name);
    time_total.Stop();
    cout << "Total Time: " << Timer::format(time_total.Duration_NS());
    return 0;
  }
}