#include <stdlib.h>
#include <stdio.h>

void karatsubaConvolutionAux(int lenH, int *h, int lenX, int *x, int *y, int *wsp) {
  if (lenH == 1) {
    for (int i = 0; i < lenX; ++i) {
      y[i] = h[0] * x[i];
    }
    return;
  }
  int split = lenH/2;
  int *x0 = x + split;  
  int *h0 = h + split; 
  karatsubaConvolutionAux(split, h, split, x, y, wsp);
  y[2 * split - 1] = 0; 
  karatsubaConvolutionAux(lenH - split, h0, lenX - split, x0, y + 2 * split, wsp);
  for (int i = 0; i < split; ++i) {
    x0[i] += x[i];
    h0[i] += h[i];
  }
  int z2Len = lenH - split + lenX - split - 1;
 
  int *z1 = wsp;
 
  karatsubaConvolutionAux(lenH - split, h0, lenX - split, x0, z1, z1 + z2Len);
  
  for (int i = 0; i < split; ++i) {
    x0[i] -= x[i];
    h0[i] -= h[i];
  }
  for (int i = 0; i < 2 * split - 1; ++i) {
    z1[i] -= y[i];
  }
  for (int i = 0; i < z2Len; ++i) {
    z1[i] -= y[2 * split + i];
    y[split + i] += z1[i];
  }
}   

void karatsubaConvolution(int lenH, int *h, int lenX, int *x, int *y) {
  if (lenH > lenX) {
    karatsubaConvolution(lenX, x, lenH, h, y);
    return;
  }
  int wspSize = (32 - __builtin_clz(lenH)) * (lenX - 1);
  int *wsp = malloc(wspSize * sizeof(int));
  karatsubaConvolutionAux(lenH, h, lenX, x, y, wsp);
  free(wsp);
}
