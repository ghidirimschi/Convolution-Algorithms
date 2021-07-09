#include <stdlib.h>
#include <string.h>
#include "LibList.h"

#define MIN_SAMPLE -100
#define MAX_SAMPLE  100


#if (MAX_SAMPLE - MIN_SAMPLE < 10000) //use histogram sort for a limited dynamic range
void duplConvolution(int lenH, int *h, int lenX, int *x, int *y) {
  List *hist = calloc(MAX_SAMPLE - MIN_SAMPLE + 1, sizeof(List));
  memset(y, 0, (lenH + lenX - 1) * sizeof(int));
  int idx;
  for (int i = 0; i < lenH; ++i) {
    idx = h[i] - MIN_SAMPLE;
    hist[idx] = addItem(i, hist[idx]);
  }
  int *scaledX = malloc(lenX * sizeof(int));
  int shift;
  for (int i = 0; i < MAX_SAMPLE - MIN_SAMPLE + 1; ++i) {
    if (isEmptyList(hist[i])) {
      continue;
    }
    shift = i + MIN_SAMPLE;
    for (int j = 0; j < lenX; ++j) {
      scaledX[j] = shift * x[j];
    }
    while (!isEmptyList(hist[i])) {
      shift = firstItem(hist[i]);
      hist[i] = removeFirstNode(hist[i]);
      for (int j = 0; j < lenX; ++j) {
        y[shift + j] += scaledX[j];
      }
    }
  }
  free(scaledX);
  free(hist);
}
#else

typedef struct ElPos {
  int el;
  int pos;
} ElPos;

int cmpElPos(const void *a, const void *b) {
  int aEl = (*(ElPos*)a).el, bEl = (*(ElPos*)b).el;
  if (aEl == -bEl) {
    return bEl - aEl;
  }
  return abs(aEl) - abs(bEl);
}

void duplConvolution(int lenH, int *h, int lenX, int *x, int *y) {
  
  ElPos *sortedH = malloc(lenH * sizeof(ElPos));
  for (int i = 0; i < lenH; ++i) {
    sortedH[i].el = h[i];
    sortedH[i].pos = i;
  }
  qsort(sortedH, lenH, sizeof(ElPos), cmpElPos);
  int *scaledX = malloc(lenX * sizeof(int));
  memset(y, 0, (lenH + lenX - 1) * sizeof(int));
  
  int i = 0,j;
  
  int shift, scale;
  
  while (!sortedH[i].el) {
    ++i;
  }
  
  while (sortedH[i].el == 1) {
    shift = sortedH[i].pos;
    for (j = 0; j < lenX; ++j) {
      y[j + shift] += x[j];
    }
    ++i;
  } 
  
  while (sortedH[i].el == -1) {
    shift = sortedH[i].pos;
    for (j = 0; j < lenX; ++j) {
      y[j + shift] -= x[j];
    }
    ++i;
  } 
  
  
  while (i < lenH) {
    shift = sortedH[i].pos;
    scale = sortedH[i].el;
    for (j = 0; j < lenX; ++j) {
      scaledX[j] = scale * x[j];
      y[j + shift] += scaledX[j];
    }  
    ++i;
    while (i < lenH && sortedH[i].el == scale) {
      shift = sortedH[i].pos;
      for (j = 0; j < lenX; ++j) {
        y[j + shift] += scaledX[j];
      } 
      ++i;
    }
    while (i < lenH && sortedH[i].el == -scale) {
      shift = sortedH[i].pos;
      for (j = 0; j < lenX; ++j) {
        y[j + shift] -= scaledX[j];
      } 
      ++i;
    }
  }
  free(sortedH);
  free(scaledX);
}
#endif
