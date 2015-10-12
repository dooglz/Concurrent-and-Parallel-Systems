#include <iostream>
#include "sequential.h"

using namespace std;

namespace seq {
double matgen(double **a, int lda, int n, double *b) {
  double norma = 0.0;
  int init = 1325;

  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < n; ++j) {
      init = 3125 * init % 65536;
      a[j][i] = (static_cast<double>(init) - 32768.0) / 16384.0;
      norma = (a[j][i] > norma) ? a[j][i] : norma;
    }
  }
  for (int i = 0; i < n; ++i)
    b[i] = 0.0;
  for (int j = 0; j < n; ++j) {
    for (int i = 0; i < n; ++i)
      b[i] += a[j][i];
  }

  return norma;
}

/* Purpose : Find largest component of double vector dx
--- Input ---
n    : number of elements in input vector
dx   : double vector with n+1 elements, dx[0] is not used
dx_off : offset in reading dx
stride : storage spacing between elements of dx
*/
int indexOfLargestElement(int n, double *dx, int dx_off, int stride) {
  double dmax, dtemp;
  int ix, itemp = 0;

  if (n < 1)
    itemp = -1;
  else if (n == 1)
    itemp = 0;
  else if (stride != 1) {
    dmax = abs(dx[0 + dx_off]);
    ix = 1 + stride;
    for (int i = 1; i < n; ++i) {
      dtemp = abs(dx[ix + dx_off]);
      if (dtemp > dmax) {
        itemp = i;
        dmax = dtemp;
      }
      ix += stride;
    }
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
  int nincx;

  if (n > 0) {
    if (incx != 1) {
      nincx = n * incx;
      for (int i = 0; i < nincx; i += incx)
        dx[i + dx_off] *= da;
    } else {
      for (int i = 0; i < n; ++i)
        dx[i + dx_off] *= da;
    }
  }
}

// Constant times a vector plus a vector
void daxpy(int n, double da, double *dx, int dx_off, int incx, double *dy,
           int dy_off, int incy) {
  int ix, iy;

  if ((n > 0) && (da != 0)) {
    if (incx != 1 || incy != 1) {
      ix = 0;
      iy = 0;
      if (incx < 0)
        ix = (-n + 1) * incx;
      if (incy < 0)
        iy = (-n + 1) * incy;
      for (int i = 0; i < n; ++i) {
        dy[iy + dy_off] += da * dx[ix + dx_off];
        ix += incx;
        iy += incy;
      }
    } else {
      for (int i = 0; i < n; ++i)
        dy[i + dy_off] += da * dx[i + dx_off];
    }
  }
}

// Performs Gaussian elimination with partial pivoting
int gaussian_eliminate(double **a, int lda, int n, int *ipvt) {
  // Pointers to columns being worked on
  double *col_k, *col_j;

  double t;
  int kp1, l;
  int nm1 = n - 1;
  int info = 0;

  if (nm1 >= 0) {
    for (int k = 0; k < nm1; ++k) {
      // Set pointer for col_k to relevant column in a
      col_k = &a[k][0];
      kp1 = k + 1;

      // Find pivot index
      l = indexOfLargestElement(n - k, col_k, k, 1) + k;
      ipvt[k] = l;

      // Zero pivot means that this column is already triangularized
      if (col_k[l] != 0) {
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
          daxpy(n - kp1, t, col_k, kp1, 1, col_j, kp1, 1);
        }
      } else
        info = k;
    }
  }

  ipvt[n - 1] = n - 1;
  if (a[n - 1][n - 1] == 0)
    info = n - 1;

  return info;
}

// Performs a dot product calculation of two vectors
double ddot(int n, double *dx, int dx_off, int incx, double *dy, int dy_off,
            int incy) {
  double temp = 0.0;
  int ix, iy;

  if (n > 0) {
    if (incx != 1 || incy != 1) {
      ix = 0;
      iy = 0;
      if (incx < 0)
        ix = (-n + 1) * incx;
      if (incy < 0)
        iy = (-n + 1) * incy;
      for (int i = 0; i < n; ++i) {
        temp += dx[ix + dx_off] * dy[iy + dy_off];
        ix += incx;
        iy += incy;
      }
    } else
      for (int i = 0; i < n; ++i)
        temp += dx[i + dx_off] * dy[i + dy_off];
  }

  return temp;
}

// Solves the system a * x = b using the factors computed in dgeco or
// gaussian_eliminate
void dgesl(double **a, int lda, int n, int *ipvt, double *b, int job) {
  double t;
  int k, l, nm1, kp1;

  nm1 = n - 1;

  if (job == 0) {
    // Solve a * x = b.  First solve l * y = b
    if (nm1 >= 1) {
      for (k = 0; k < nm1; ++k) {
        l = ipvt[k];
        t = b[l];
        if (l != k) {
          b[l] = b[k];
          b[k] = t;
        }
        kp1 = k + 1;
        daxpy(n - kp1, t, &a[k][0], kp1, 1, b, kp1, 1);
      }
    }

    // Now solve u * x = y
    for (int kb = 0; kb < n; ++kb) {
      k = n - (kb + 1);
      b[k] /= a[k][k];
      t = -b[k];
      daxpy(k, t, &a[k][0], 0, 1, b, 0, 1);
    }
  } else {
    // Solve trans(a) * x = b.  First solve trans(u) * y = b
    for (k = 0; k < n; ++k) {
      t = ddot(k, &a[k][0], 0, 1, b, 0, 1);
      b[k] = (b[k] - t) / a[k][k];
    }
    // Solve trans(l) * x = y
    if (nm1 >= 1) {
      for (int kb = 1; kb < nm1; ++kb) {
        k = n - (kb + 1);
        kp1 = k + 1;
        b[k] += ddot(n - kp1, &a[k][0], kp1, 1, b, kp1, 1);
        l = ipvt[k];
        if (l != k) {
          t = b[l];
          b[l] = b[k];
          b[k] = t;
        }
      }
    }
  }
}

// Multiply matrix m times vector x and add the result to vector y
void dmxpy(int n1, double *y, int n2, int ldm, double *x, double **m) {
  for (int j = 0; j < n2; ++j)
    for (int i = 0; i < n1; ++i)
      y[i] += x[j] * m[j][i];
}

// Estimates roundoff in quantities of size x
double epslon(double x) {
  double eps = 0.0;
  double a = 4.0 / 3.0;
  double b, c;

  while (eps == 0) {
    b = a - 1.0;
    c = b + b + b;
    eps = abs(c - 1.0);
  }

  return eps * abs(x);
}

// Initialises the system
void initialise(double **a, double *b, double &ops, double &norma, double lda) {
  long long nl = static_cast<long long>(SIZE);
  ops = (2.0 * static_cast<double>((nl * nl * nl))) / 3.0 +
        2.0 * static_cast<double>((nl * nl));
  norma = matgen(a, lda, SIZE, b);
}

// Runs the benchmark
void run(double **a, double *b, int &info, double lda, int n, int *ipvt) {
  info = gaussian_eliminate(a, lda, n, ipvt);
  dgesl(a, lda, n, ipvt, b, 0);
}

// Validates the result
void validate(double **a, double *b, double *x, double &norma, double &normx,
              double &resid, double lda, int n) {
  double eps, residn;
  double ref[] = {6.0, 12.0, 20.0};

  for (int i = 0; i < n; ++i)
    x[i] = b[i];

  norma = matgen(a, lda, n, b);

  for (int i = 0; i < n; ++i)
    b[i] = -b[i];

  dmxpy(n, b, n, lda, x, a);
  resid = 0.0;
  normx = 0.0;
  for (int i = 0; i < n; ++i) {
    resid = (resid > abs(b[i])) ? resid : abs(b[i]);
    normx = (normx > abs(x[i])) ? normx : abs(x[i]);
  }

  eps = epslon(1.0);
  residn = resid / (n * norma * normx * eps);
  if (residn > CHECK_VALUE) {
    cout << "Validation failed!" << endl;
    cout << "Computed Norm Res = " << residn << endl;
    cout << "Reference Norm Res = " << CHECK_VALUE << endl;
  } else {
    cout << "Calculations are correct!" << endl;
    cout << "Computed Norm Res = " << residn << endl;
    cout << "Reference Norm Res = " << CHECK_VALUE << endl;
  }
}

int start() {
  // Allocate data on the heap
  double **a = new double *[SIZE];
  for (int i = 0; i < SIZE; ++i)
    a[i] = new double[SIZE];
  double *b = new double[SIZE];
  double *x = new double[SIZE];
  int *ipvt = new int[SIZE];

  double ldaa = static_cast<double>(SIZE);
  double lda = ldaa + 1;
  double ops, norma, normx;
  double resid;
  int info;

  // Main application
  initialise(a, b, ops, norma, lda);
  run(a, b, info, lda, SIZE, ipvt);
  validate(a, b, x, norma, normx, resid, lda, SIZE);

  // Free the memory
  for (int i = 0; i < SIZE; ++i)
    delete[] a[i];
  delete[] a;
  delete[] b;
  delete[] x;
  delete[] ipvt;

  return 0;
}
}