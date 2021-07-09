#include <stdlib.h>
#include <string.h>

void duplConvolution(int lenH, int *h, int lenX, int *x, int *y);

void derivativeConvolution(int lenH, int *h, int lenX, int *x, int *y) {
  int lenY = lenH + lenX - 1;
  
  int *dH = malloc((lenH + 1) * sizeof(int));
  dH[0] = h[0];
  for(int i = 1; i < lenH; ++i) {
    dH[i] = h[i] - h[i - 1];
  }
  dH[lenH] = -h[lenH - 1];
  
  int *dY = malloc((lenY + 1) * sizeof(int));
  duplConvolution(lenH + 1, dH, lenX, x, dY);
  
  y[0] = dY[0];
  for (int i = 1; i < lenY; ++i) {
    y[i] = y[i - 1] + dY[i];
  }
  
  free(dH);
  free(dY);
} 
