#include <stdlib.h>
#include <string.h>

int **D;
#define D(a,b) ((b) < (a) ? 0 : D[a][(b) - (a)])


int **allocD(int N) {
  int **D = malloc(N * sizeof(int*));
  for (int i = 0; i < N; ++i) {
    D[i] = malloc((N - i) * sizeof(int));
  }
  return D;
}

int E_dak(int p, int q) {
  return D(p, q) - D(p + 1, q);
}

int E_deuk(int p, int q) {
  return D(p, q) - D(p, q - 1);
}

void sietseConvolutionAux(int lenH, int *h, int lenX, int *x, int *y) {
  int lenY = lenH + lenX - 1;
  D = allocD(lenH);
  int aSum, bSum;
  
  for (int i = 0; i < lenH; ++i) {
    aSum = 0;
    bSum = 0;
    for (int j = i; j < lenX; ++j) {
      aSum += h[j];
      bSum += x[j];
      D[i][j - i] = aSum * bSum;
    }
  } 
  
  int *C = malloc(lenY * sizeof(int));
  for (int k = 0; k < lenH - 1; ++k) {
    C[k] = 0;
    for (int j = 0; j < k / 2 + 1; ++j) {
      C[k] += E_dak(j, k - j);
    }
  }
  for (int k = lenH; k < lenY; ++k) {
    C[k] = 0;
    for (int j = k - lenH + 1; j <= (k+1)/2; ++j) {
      C[k] += E_deuk(j, k - j);
    }
  }
  y[0] = C[0];  
  y[lenY - 1] = C[lenY - 1];  
  y[lenH - 1] = D(0, lenH - 1) - C[lenH  - 2] - C[lenH];
  
  for (int i = 1; i < lenH - 1; ++i) {
    y[i] = C[i] - C[i - 1];    
  }
  
  for (int i = lenH; i < lenY - 1; ++i) {
    y[i] = C[i] - C[i + 1];
  }
  
  y[lenY - 1] = C[lenY - 1];
  
  for (int i = 0; i < lenH; ++i) {
    free(D[i]);
  }
  free(D);
  free(C);
}

void sietseConvolution(int lenH, int *h, int lenX, int *x, int *y) {
  if (lenH > lenX) {
    sietseConvolution(lenX, x, lenH, h, y);
    return;
  }
  if (lenH == lenX) {
    sietseConvolutionAux(lenH, h, lenX, x, y);
    return;
  }
  int *th = calloc(lenX, sizeof(int));
  memcpy(th, h, lenH * sizeof(int));
  int *ty = malloc(((lenX << 1) - 1) * sizeof(int));
  sietseConvolutionAux(lenX, th, lenX, x, ty);
  for (int i = 0; i < lenX + lenH - 1; ++i) {
    y[i] = ty[i];
  }
  
  free(th);
  free(ty);
}
