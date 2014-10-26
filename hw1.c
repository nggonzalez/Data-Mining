#include <stdio.h>
#include <stdlib.h>
#include <math.h>

void jacobi(double *a, int n, double *s, double *u, double *v);
double norm(int n, int col, double *u);
double calculateUSquared(int size, int col, double* u);
double calculateGamma(int size, int row, int col, double* u);
void bubbleSort(int size, double *u, double *s, double *v);
void setValueInMatrix(int n, int row, int col, double newVal, double* matrix);
double getValueInMatrix(int n, int row, int col, double* matrix);
double sign(double num);
// void print(double *u, int n);

int main(int argc, char*argv[]) {
  // int n = 3;
  // double* a = malloc(sizeof(double) * n * n);
  // double* u = malloc(sizeof(double) * n * n);
  // double* v = malloc(sizeof(double) * n * n);
  // double* s = malloc(sizeof(double) * n);

  //   a[0] = 1;
  //   a[1] = 2;
  //   a[2] = 3;
  //   a[3] = 2;
  //   a[4] = -3;
  //   a[5] = 4;
  //   a[6] = 3;
  //   a[7] = 4;
  //   a[8] = 5;


  // jacobi(a, n, s, u, v);
  
  // printf("U:\n");
  // print(u, n);
  // printf("V:\n");
  // print(v, n);
  // printf("S:\n");
  // for(int i = 0; i < n; i++) {
  //   printf("%.5f ", s[i]);
  // }
  // printf("\n");

  return 1;
}

void jacobi(double *a, int n, double *s, double *u, double *v) {
  int i;
  for(i = 0; i < n * n; i++) {
    u[i] = a[i];
    v[i] = 0;
  }

  for(i = 0; i < n; i++) {
    setValueInMatrix(n, i, i, 1.0, v);
  }

  double threshold = .000000001;
  double convergence = 1 + threshold;
  while(convergence > threshold) {
    convergence = 0;
    
    int col, row;
    for(col = 0; col < n; col++) {
      for(row = 0; row < n; row++) {
        if(row >= col) {
          continue;
        }
        double alpha = calculateUSquared(n, row, u);
        double beta = calculateUSquared(n, col, u);
        double gamma = calculateGamma(n, row, col, u);

        //printf("Alpha: %f, Beta: %f, Gamma: %f\n", alpha, beta, gamma);
        convergence = fmax(convergence, fabs(gamma)/sqrt(alpha * beta));
        // printf("convergence %f, threshold %f\n", convergence, threshold);

        double zeta = (beta - alpha) / (2.0 * gamma);
        double tangent = sign(zeta) / (fabs(zeta) + sqrt(1.0 + (zeta * zeta)));
        double cosine = 1.0 / sqrt(1.0 + tangent * tangent);
        double sine = cosine * tangent;
        // printf("Zeta: %f, Tangent: %f, Cosine: %f Sine: %f\n", zeta, tangent, cosine, sine);


        int k;
        for(k = 0; k < n; k++) {
          tangent = getValueInMatrix(n, k, row, u);
          //printf("Tangent %f\t", tangent);

          double Ukj = getValueInMatrix(n, k, col, u);
          double newVal = cosine*tangent - sine * Ukj;
          
          //printf("newVal Uki: %f\t", newVal);
          
          setValueInMatrix(n, k, row, newVal, u);
          
          newVal = sine*tangent + cosine * Ukj;
          //printf("newVal Ukj: %f\n", newVal);
          
          setValueInMatrix(n, k, col, newVal, u);
        }

        for(k = 0; k < n; k++) {
          tangent = getValueInMatrix(n, k, row, v);
          // printf("Tangent %f\t", tangent);
          double Vkj = getValueInMatrix(n, k, col, v);
          double newVal = cosine*tangent - sine * Vkj;
          // printf("newVal Vki: %f\t", newVal);

          setValueInMatrix(n, k, row, newVal, v);
          
          newVal = sine*tangent + cosine * Vkj;
          // printf("newVal Vkj: %f\n", newVal);
          setValueInMatrix(n, k, col, newVal, v);
        }

      }
    }
  }

  int col;
  for(col = 0; col < n; col++) {
    s[col] = norm(n, col, u);

    // printf("Norm %f\n", s[col]);

    int k;
    for(k = 0; k < n; k++) {
      double Ukj = getValueInMatrix(n, k, col, u);
      setValueInMatrix(n, k, col, Ukj / s[col], u);
    }
  }

  bubbleSort(n, u, s, v);
}

double sign(double num) {
  if(num > 0) {
    return 1;
  } else if (num < 0) {
    return -1;
  } else {
    return 0;
  }
}

// void print(double *u, int n) {
//     int i, j;
//     for (i = 0; i < n; i++) {
//         for (j = 0; j < n; j++) {
//             printf("%f ", u[i*n + j]);
//         }
//         printf("\n");
//     }
//     printf("\n");
// }

double norm(int n, int col, double *u) {
  double sum = 0;
  int k;
  for(k = 0; k < n; k++) {
    sum += pow(getValueInMatrix(n, k, col, u), 2);
  }

  return sqrt(sum);  
}

double calculateUSquared(int n, int col, double* u) {
  double sum = 0;
  int row;
  for(row = 0; row < n; row++) {
    double temp = getValueInMatrix(n, row, col, u);
    sum += temp*temp;
  }

  return sum;
}


double calculateGamma(int n, int row, int col, double* u) {
  double sum = 0;
  int k;
  for(k = 0; k < n; k++) {
    sum += u[k * n + row] * u[k * n + col];
  }

  return sum;
}


void bubbleSort(int size, double *u, double *s, double *v) {
  int i, j, k;
  for(i = 0; i < size - 1; i++) {
    for(j = 0; j < size - i - 1; j++) {
      if(s[j] < s[j + 1]) {
        double swap = s[j];
        s[j] = s[j+1];
        s[j+1] = swap;

        for(k = 0; k < size; k++) {
          double swapUkj = getValueInMatrix(size, k, j, u);
          double swapUkj1 = getValueInMatrix(size, k, j+1, u);
          setValueInMatrix(size, k, j, swapUkj1, u);
          setValueInMatrix(size, k, j+1, swapUkj, u);
        }

        for(k = 0; k < size; k++) {
          double swapVkj = getValueInMatrix(size, k, j, v);
          double swapVkj1 = getValueInMatrix(size, k, j+1, v);
          setValueInMatrix(size, k, j, swapVkj1, v);
          setValueInMatrix(size, k, j+1, swapVkj, v);
        }
      }
    }
  }
}


void setValueInMatrix(int n, int row, int col, double newVal, double* matrix) {
  matrix[row * n + col] = newVal;
}


double getValueInMatrix(int n, int row, int col, double* matrix) {
  return matrix[row * n + col];
}
