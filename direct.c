#include <string.h>

void inputSideConvolution(int lenH, int *h, int lenX, int *x, int *y) {
  if (lenH > lenX) {
    inputSideConvolution(lenX, x, lenH, h, y);
    return;
  }
  // from here: lenH <= lenX
  int lenY = lenH + lenX - 1;
  memset(y, 0, lenY*sizeof(int)); // sets y[i]=0 for all i
  for(int i=0; i < lenX; i++) {
    for(int j=0; j < lenH; j++) {
      y[i+j] += x[i]*h[j];
    }
  }
}

void outputSideConvolution(int lenH, int *h, int lenX, int *x, int *y) {
  // from here: lenH <= lenX
  int lwb, upb, lenY = lenH + lenX - 1;
  for(int i = 0; i < lenY; i++) {
    // lwb and upb are chosen such that: 0<=j && j<lenH && 0<=i-j && i-j<lenX
    // which is equivalent with: MAX(0,i+1-lenX) <= j < MIN(i+1,lenH)
    lwb = (i < lenX ? 0: i-lenX+1);
    upb = (i < lenH ? i+1 : lenH);
    y[i] = 0;
    for (int j = lwb; j < upb; j++) {
      y[i] += h[j]*x[i-j];
    }
  }
}
