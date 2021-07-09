#include <stdlib.h>
#include <complex.h>
#include <string.h>
#include <stdio.h>

#define PI 3.1415926535897932384626433832795L

int powerOfTwo(int n);

void inplaceFFT(int length, double complex *x);

void inplaceInverseFFT(int length, double complex *x);


void oamfftConvolution(int lenH, int *h, int lenX, int *x, int *y) {
  if (lenH > lenX) {
    oamfftConvolution(lenX, x, lenH, h, y);
    return;
  }
  int N = 8 * powerOfTwo(lenH),
      step_size = N - (lenH - 1);
       
  double complex *ch = calloc(N, sizeof(double complex));
  double complex *cx = calloc(N, sizeof(double complex));
  for (int i=0; i < lenH; i++) ch[i] = h[i];
  inplaceFFT(N, ch);
  int position = 0;
  memset(y, 0, (lenX + lenH - 1) * sizeof(int));
  
  while (position + step_size <= lenX) {
    memset(cx, 0, N * sizeof(double complex));
    for (int i = 0; i < step_size; ++i) {
      cx[i] = x[position + i];
    }
    inplaceFFT(N, cx);
    for (int i = 0; i < N; i++) {
      cx[i] *= ch[i];
    }
    inplaceInverseFFT(N, cx);
    for (int i = 0; i < N; ++i) {
      y[position + i] += (int)(cx[i]+0.5);
    }
    position += step_size;
  }
  if (position != lenX) {
    step_size = lenX - position;
    memset(cx, 0, N * sizeof(double complex));
    for (int i = 0; i < step_size; ++i) {
        cx[i] = x[position + i];
      }
    }
    inplaceFFT(N, cx);
    for (int i = 0; i < N; i++) {
      cx[i] *= ch[i];
    }
    inplaceInverseFFT(N, cx);
    for (int i = 0; i < lenX + lenH - position - 1; ++i) {
      y[position + i] += (int)(cx[i]+0.5);
  }
  free(ch);
  free(cx);
}
