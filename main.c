#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <limits.h>
#include "timer.h"

#define CHECK_ALGORITHMS 1
#define LIMIT_CACHING   1


int *generateArray(int lowerLimit, int upperLimit, int size, double prob, int period) {
  int *array = malloc(size * sizeof(int));
  double rProb;
  for (int i = 0; i < period; ++i) {
    array[i] = (rand() % (upperLimit - lowerLimit + 1)) + lowerLimit;
  }
  if (!period) {
    array[0] = (rand() % (upperLimit - lowerLimit + 1)) + lowerLimit;
  }
  for(int i = !period ? 1 : period; i < size; i++) {
    rProb = (double)rand() / RAND_MAX;
    if(rProb <= prob) {
      if (period) {
        for (int j = 0; j < period && i+j < size; ++j) {
          array[i+j] = array[j];
        }
        i += period - 1;
      } else {
        array[i] = array[i - 1];
      }
    } else {
      array[i] = (rand() % (upperLimit - lowerLimit + 1)) + lowerLimit; 
    }
  }
  return array;
}

int *readSignal(int *len) {
  int *x;
  char c;
  scanf("%d:", len);
  x = calloc(*len, sizeof(int));
  do c = getchar(); while (c != '[');
  int temp;
  if (*len > 0) {
    scanf("%d", &temp);
    x[0] = temp;
    for (int i=1; i < *len; i++) {
      scanf(",%d", &temp);
      x[i] = temp;
    }
  }
  do c = getchar(); while (c != ']');
  return x;
}

void printSignal(int len, int *x) {
  if (len > 0) {
    printf("[%d", x[0]);
    for (int i = 1; i < len; i++) printf(",%d", x[i]);
  }
  printf("]");
  return;
}

/************* Prototypes convolution algorithms **************/
void inputSideConvolution(int lenH, int *h, int lenX, int *x, int *y);
void outputSideConvolution(int lenH, int *h, int lenX, int *x, int *y);
void fftConvolution(int lenH, int *h, int lenX, int *x, int *y);
void karatsubaConvolution(int lenH, int *h, int lenX, int *x, int *y);
void duplConvolution(int lenH, int *h, int lenX, int *x, int *y);
void sietseConvolution(int lenH, int *h, int lenX, int *x, int *y);
void lrsnoConvolution(int lenH, int *h, int lenX, int *x, int *y);
void lrsnoConvolutionSort(int lenH, int *h, int lenX, int *x, int *y);
void oamfftConvolution(int lenH, int *h, int lenX, int *x, int *y);
void derivativeConvolution(int lenH, int *h, int lenX, int *x, int *y);
void nttConvolutionBarret(int lenH, int *h, int lenX, int *x, int *y);
void nttConvolutionSoftware(int lenH, int *h, int lenX, int *x, int *y);
void nttcrtConv(int lenH, int *h, int lenX, int *x, int *y);
void nttcrtBarrettConv(int lenH, int *h, int lenX, int *x, int *y);
/************* Timing routine ********************************/
  
double timeConvolutionAlgorithm(void (algorithm)(int,int*,int,int*,int*),
				int lenH, int *h, int lenX, int *x, int *y) {
  Timer timer;
  double fastest = 999999999;
  const int size = 20*1024*1024;
  char *c = (char *)malloc(size);
  for (int iteration=0; iteration<50; iteration++) {
    #if LIMIT_CACHING
    int i = 1;
    for (int j = 0; j < size; j++)
      c[j] = j * (++i*2);
  
    #endif
    startTimer(&timer);
    algorithm(lenH, h, lenX, x, y);
    stopTimer(&timer);
    fastest = (seconds(timer) < fastest ? seconds(timer) : fastest);
  }
  free(c);
  return fastest;
}

void clearHashMap();

int main(int argc, char *argv[]) {
  void (*algorithm[])(int,int*,int,int*,int*) = {
     fftConvolution, nttcrtConv, nttConvolutionBarret, nttcrtBarrettConv, nttConvolutionSoftware, oamfftConvolution, lrsnoConvolution, derivativeConvolution, duplConvolution, inputSideConvolution,  outputSideConvolution, karatsubaConvolution, sietseConvolution
  };

  char name[][32] = {
    "fftConvolution", "nttcrtConv", "nttConvolutionBarret", "nttcrtBarrettConv", "nttConvolutionSoftware", "oamfftConvolution", "lrsnoConvolution",  "derivativeConvolution", "duplConvolution", "inputSideConvolution", "outputSideConvolution", "karatsubaConvolution", "sietseConvolution"
  };

  if ((argc < 3) || (argc > 5)) {
    fprintf(stderr, "Usage: %s <lenX> <lenH> [probability from [0..1]] [period of h]\n", argv[0]);
    return -1;
  }
  
  int lenX = atoi(argv[1]);
  int lenH = atoi(argv[2]);
  int lenY = lenX + lenH - 1;

  double probability = (argc < 4 ? 0.0 : atof(argv[3]));
  int period = (argc < 5 ? 0 : atoi(argv[4]));
  
  int *h = generateArray(0, 100, lenH, probability, period);
  int *x = generateArray(0, 100, lenX, probability, 0);
  int *y = malloc(lenY*sizeof(int));
  
#if CHECK_ALGORITHMS
  int *reference = malloc(lenY*sizeof(int));
  outputSideConvolution(lenH, h, lenX, x, reference);

  for (int alg=0; alg < sizeof(name)/32; alg++) {
    algorithm[alg](lenH, h, lenX, x, y);
    for (int i=0; i < lenY; i++) {
      if (y[i] != reference[i]) {
	printf("Fatal error: Algorithm %s produced wrong output (y[%d]=%d, while it should have been %d)\n", name[alg], i, y[i], reference[i]);
	exit(EXIT_FAILURE);
      }
    }
  }
  free(reference);
#endif
  
  for (int i=0; i < sizeof(name)/32; i++) {
    double elapsedTime = timeConvolutionAlgorithm(algorithm[i], lenH, h, lenX, x, y);
    printf("%32s:    %lf miliseconds.\n", name[i], elapsedTime * 1000.0);
    //printf("%lf\t", elapsedTime * 1000.0);

  }
  printf("\n");

  free(y);
  free(x);
  free(h);
  clearHashMap();
  return 0;
}
