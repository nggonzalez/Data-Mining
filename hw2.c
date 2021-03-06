#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <limits.h>

double norm(int length, double *sum);
void initColToValue(double *x, int length, double value);
void multiplyMatrixByColumn(double *r, double *a, double *x, int length);
double calculateAlpha(double *r, double *a, int length);
void subtractCols(double *sum, double *y, double *r, int length);
double getValueInMatrix(int n, int row, int col, double* matrix);
void setValueInMatrix(int n, int row, int col, double newVal, double* matrix);
void cloneArray(double *destination, double *starting, int length);
void getNextXSoln(double *x, double alpha, double *r, int length);
void scaleMatrix(int length, double *matrix, double scalar);
void transposeA(double *at, double *a, int length);
void dumb_solve(double *a, double *y, int n, double eps, int numit, double *x, int *niter, double *discreps);


void dumb_solve(double *a, double *y, int n, double eps, int numit, double *x, int *niter, double *discreps) {
  *niter = 0;
  initColToValue(x, n, 0);
  initColToValue(discreps, numit, 0);

  double currentDiscreps = INT_MAX;
  double *ax = malloc(sizeof(double) * n);
  double *temp = malloc(sizeof(double) * n); // col of size n
  double *r = malloc(sizeof(double) * n); // col of size n
  double *sum = malloc(sizeof(double) * n);
  double alpha;
  int i;
  double *xs = malloc(sizeof(double) * n);
  double *rs = malloc(sizeof(double) * n);
  double *at = malloc(sizeof(double) * n * n);
  transposeA(at, a, n);
  double startDiscreps;

  while(*niter < numit && currentDiscreps > eps) {
    initColToValue(r, n, 0); // r = 0s
    initColToValue(sum, n, 0); // sum = 0s

    multiplyMatrixByColumn(ax, a, x, n); // ax = Ax
    subtractCols(temp, y, ax, n); // r = y - ax
    // printf("Y-ax:\n");
    // for(i = 0; i < n; i++)
    //   printf("%.5f\t", temp[i]);
    // printf("\n\n");

    multiplyMatrixByColumn(r, at, temp, n); // r = 2at(y-ax)


    // printf("Residuals:\n");
    // for(i = 0; i < n; i++)
    //   printf("%.5f\t", r[i]);
    // printf("\n\n");

    currentDiscreps = pow(norm(n, temp), 2);
    discreps[(*niter)++] = currentDiscreps;
    // printf("currentDiscreps: %.8f\n", currentDiscreps);

    // alpha = calculateAlpha(r, a, n);
    alpha = INT_MAX;
    // if(currentDiscreps > discreps[*niter - 1]) {
    //   alpha = .5*alpha;
    // } else {
    //   alpha = 1.2*alpha;
    // }
    // printf("Alpha %f\n", alpha);

    startDiscreps = currentDiscreps;
    cloneArray(xs, x, n);
    cloneArray(rs, r, n);
    while(startDiscreps <= currentDiscreps) {
      getNextXSoln(x, alpha, r, n);
      multiplyMatrixByColumn(ax, a, x, n); // ax = Ax
      subtractCols(temp, y, ax, n); // temp = y - ax
      currentDiscreps = pow(norm(n, temp), 2);
      alpha *= .5;
      cloneArray(x, xs, n);
      cloneArray(r, rs, n);
      // printf("currentDiscreps: %.8f\n", currentDiscreps);
    }
    getNextXSoln(x, alpha, r, n);
    // currentDiscreps = startDiscreps;
  }

  free(at);
  free(ax);
  free(temp);
  free(xs);
  free(rs);
  free(sum);
  free(r);
}


void transposeA(double *at, double *a, int length) {
  int row, col;
  double swap;
  for(row = 0; row < length; row++) {
    for (col = 0; col < length; col++) {
      swap = 2 * getValueInMatrix(length, row, col, a);
      setValueInMatrix(length, col, row, swap, at);
    }
  }

  // printf("AT\n");
  // for(row = 0; row < length; row++) {
  //   for (col = 0; col < length; col++) {
  //     printf("%f\t", getValueInMatrix(length, row, col, at));
  //   }
  //   printf("\n");
  // }
  // printf("\n\n");
}


void cloneArray(double *destination, double *starting, int length) {
  int i;
  for(i = 0; i < length; i++) {
    destination[i] = starting[i];
  }

  // for(i = 0; i < length; i++) {
  //   printf("%f\t" , destination[i]);
  // }
  // printf("\n\n");
}


void scaleMatrix(int length, double *matrix, double scalar) {
  int i, j;
  for(i = 0; i < length; i++) {
    for(j = 0; j < length; j++) {
      setValueInMatrix(length, i, j, scalar * getValueInMatrix(length, i, j, matrix), matrix);
    }
  }
}


void getNextXSoln(double *x, double alpha, double *r, int length) {
  int i;
  for(i = 0; i < length; i++) {
    x[i] = x[i] + (alpha * r[i]);
    // printf("x%d: %f\t", i, x[i]);
  }
  // printf("\n\n");
}


double calculateAlpha(double *r, double *a, int length) {
  double alpha = 0, divisor = 0;
  int i;
  for(i = 0; i < length; i++) {
    alpha += r[i] * r[i];
    // printf("Alpha %f\n", alpha);
  }

  double *temp = malloc(sizeof(double) * length);
  int row;
  double tempSum = 0;
  for(i = 0; i < length; i++) {
    for(row = 0; row < length; row++) {
      tempSum += r[row] * getValueInMatrix(length, row, i, a);
      // printf("tempSum: %f\n", tempSum);
    }

    temp[i] = tempSum;
    tempSum = 0;
  }

  for(i = 0; i < length; i++) {
    divisor += temp[i] * r[i];
  }

  free(temp);

  // printf("divisor: %f\n", divisor);

  if(divisor != 0) {
    alpha = alpha / divisor;
  }

  return alpha;
}


void subtractCols(double *sum, double *first, double *last, int length) {
  int i;
  for(i = 0; i < length; i++) {
    sum[i] = first[i] - last[i];
    // printf("Y-ax %.50f\n", sum[i]);
  }
  // printf("\n");
}


void initColToValue(double *x, int length, double value) {
  int i;
  for(i = 0; i < length; i++)
    x[i] = value;
}


void setValueInMatrix(int n, int row, int col, double newVal, double* matrix) {
  matrix[row * n + col] = newVal;
}


double getValueInMatrix(int n, int row, int col, double* matrix) {
  return matrix[row * n + col];
}


double norm(int length, double *sum) {
  double norm = 0;
  int k;
  for(k = 0; k < length; k++) {
    norm += pow(sum[k], 2);
  }

  return sqrt(norm);
}


void multiplyMatrixByColumn(double *r, double *a, double *x, int length) {
  int row, col;
  double sum;
  for(row = 0; row < length; row++) {
    sum = 0;
    for(col = 0; col < length; col++) {
      sum += getValueInMatrix(length, row, col, a) * x[col];
    }
    r[row] = sum;
  }
}

// will this multiply the matrix correctly
// a is n * n
// x is a col of size n


int main(int argc, char*argv[]) {
  // 5x + 4y –  z = 0
  //     10y – 3z = 11
  //            z = 3

  // int n = 10;
  // int numit = 1000;
  // int niter = 0;
  // double eps = 0.000001;

  // double* a = malloc(sizeof(double) * n * n);
  // double* y = malloc(sizeof(double) * n);
  // double* x = malloc(sizeof(double) * n);
  // double* discreps = malloc(sizeof(double) * numit);

  // n = 1;
  // a[0] = 2;
  // y[0] = 22;

  // n = 2;
  // a[0] = 1;
  // a[1] = 1;
  // a[2] = 2;
  // a[3] = 1;

  // y[0] = 3;
  // y[1] = 4;

  // n=3;
  // a[0] = 5;
  // a[1] = 4;
  // a[2] = -1;
  // a[3] = 0;
  // a[4] = 10;
  // a[5] = -3;
  // a[6] = 0;
  // a[7] = 0;
  // a[8] = 1;

  // y[0] = 0;
  // y[1] = 11;
  // y[2] = 3;

  // n = 3;
  // int row, col;
  // for(row = 0; row < n; row++) {
  //   for (col = 0; col < n; col++)
  //   {
  //     if(row == col) {
  //       setValueInMatrix(n, row, col, 1.0/ (double)((row+1)*(col+1)), a);
  //     } else {
  //       setValueInMatrix(n, row, col, 0, a);
  //     }
  //   }
  // }
  // initColToValue(y, n, 1.0);

  // n = 2;
  // a[0] = 1;
  // a[1] = 1;
  // a[2] = 2;
  // a[3] = -2;

  // y[0] = 3;
  // y[1] = -4;

  // n = 3;
  // a[0] = 1;
  // a[1] = -2;
  // a[2] = 3;
  // a[3] = 2;
  // a[4] = 1;
  // a[5] = 1;
  // a[6] = -3;
  // a[7] = 2;
  // a[8] = -2;

  // y[0] = 7;
  // y[1] = 4;
  // y[2] = -10;

  // dumb_solve(a, y, n, eps, numit, x, &niter, discreps);

  // printf("A:\n");
  // for(row = 0; row < n; row++) {
  //   for (col = 0; col < n; col++)
  //   {
  //     printf("%.20f ", getValueInMatrix(n, row, col, a));
  //   }
  //   printf("\n");
  // }
  // printf("\n");

  // printf("X:\n");
  // for(row = 0; row < n; row++) {
  //   printf("%.5f ", x[row]);
  // }
  // printf("\n");

  // printf("Y:\n");
  // for(row = 0; row < n; row++) {
  //   printf("%.5f ", y[row]);
  // }
  // printf("\n\n");

  // printf("Numit: %d,\tNiter: %d\n\n", numit, niter);

  // printf("Discreps:\n");
  // for(row = 0; row < niter; row++) {
  //   printf("%.5f ", discreps[row]);
  // }

  return 1;
}
