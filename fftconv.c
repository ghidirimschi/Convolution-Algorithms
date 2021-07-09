#include <stdlib.h>
#include <complex.h>

#define PI 3.1415926535897932384626433832795L

int powerOfTwo(int n) {
  int p2;
  for (p2=1; p2<n; p2*=2);
  return p2;
}

static void inplaceCooleyTukeyFFT(int length, double complex *a, double complex omega, double complex *wsp) {
  // wsp is a workspace array. This is used to avoid abundant allocas and frees of memory.
  if (length < 2) return;
  // split array a in even and odd part
  length /= 2;
  double complex *even = wsp;
  double complex *odd = wsp + length;
  int idx = 0;  
  for (int i=0; i < length; i++) {
    even[i] = a[idx++];
    odd[i] = a[idx++];    
  }
  // recursive calls, with w^2 as root of unity for the halved length
  // note that in the next recursion level, a[] plays the role of wsp.
  inplaceCooleyTukeyFFT(length, even, omega*omega, a);
  inplaceCooleyTukeyFFT(length, odd, omega*omega, a);
  // merge phase
  double complex h, x = 1;
  for (int i=0; i < length; i++) {
    h = x*odd[i];
    a[i] = even[i] + h;
    a[i+length] = even[i] - h;
    x = x*omega;
  }
}

void inplaceFFT(int length, double complex *x) {
  // note that length is expected to be a power of two!
  double complex *wsp = malloc(length*sizeof(double complex)); // workspace array
  double complex omega = cexp(-2.0*PI*I/length);     // note that I is the c99 representation of the complex number 0 + i1
  inplaceCooleyTukeyFFT(length, x, omega, wsp);
  free(wsp);
}

void inplaceInverseFFT(int length, double complex *x) {
  // note that length is expected to be a power of two!
  double complex *wsp = malloc(length*sizeof(double complex)); // workspace array
  double complex omega = cexp(2.0*PI*I/length);     // note that I is the c99 representation of the complex number 0 + i1
  inplaceCooleyTukeyFFT(length, x, omega, wsp);
  free(wsp);
  for (int i=0; i < length; i++) x[i] /= length;
}

void fftConvolution(int lenH, int *h, int lenX, int *x, int *y) {
  int lenY = lenX + lenH - 1;
  int len2 = powerOfTwo(lenY);
  double complex *cx = calloc(2*len2, sizeof(double complex));
  double complex *ch = cx + len2;
  for (int i=0; i < lenX; i++) cx[i] = x[i];
  for (int i=0; i < lenH; i++) ch[i] = h[i];
  inplaceFFT(len2, ch);
  inplaceFFT(len2, cx);
  for (int i=0; i < len2; i++) cx[i] *= ch[i];
  inplaceInverseFFT(len2, cx);
  for (int i=0; i < lenY; i++) y[i] = (int)(cx[i]+0.5); // round to integers
  free(cx);
}
