#include <stdlib.h>
#include <stdio.h>
#include <stdint.h>



static int omegas65537[][2] = {
	{1, 1},
	{65536, 65536},
	{65281, 256},
	{4096, 65521},
	{64, 64513},
	{65529, 8192},
	{8224, 64509},
	{13987, 39504},
	{282, 64375},
	{15028, 35843},
	{19139, 23398},
	{61869, 29606},
	{54449, 3505},
	{6561, 58355},
	{81, 8091},
	{9, 7282},
	{3, 21846}
};

static int nInv65537[] = {
	1,
	32769,
	49153,
	57345,
	61441,
	63489,
	64513,
	65025,
	65281,
	65409,
	65473,
	65505,
	65521,
	65529,
	65533,
  65535,
	65536
};


static int omegas114689[][2] = {
	{1, 1},
	{114688, 114688},
	{56022, 58667},
	{108799, 9327},
	{28708, 75390},
	{1800, 24913},
	{11579, 24792},
	{26123, 107322},
	{31809, 87413},
	{73478, 4347},
	{110050, 19531},
	{50231, 11955},
	{50625, 8022},
	{225, 84615},
	{15, 7646}
};

static int nInv114689[] = {
	1,
	57345,
	86017,
	100353,
	107521,
	111105,
	112897,
	113793,
	114241,
	114465,
	114577,
	114633,
	114661,
	114675,
	114682
};


static int add(int a, int b, int M) {
  a += b;
  return (a > M) ? (a - M) : a;
}


static int multiply(int a, int b, int M) {
  return (int) (((long) a * b) % M);
}

static int msb(int x) {
  int r = 0;
  while (x >>= 1) {
    r++;
  }
  return r;
}

static int powerOfTwo(int n) {
  int p2;
  for (p2=1; p2<n; p2*=2);
  return p2;
}


static void inplaceCooleyTukeyNTT(int length, int *a, int omega, int *wsp, int M) {
  // wsp is a workspace array. This is used to avoid abundant allocas and frees of memory.
  if (length < 2) return;
  // split array a in even and odd part
  length /= 2;
  int *even = wsp;
  int *odd = wsp + length;
  int idx = 0;  
  for (int i=0; i < length; i++) {
    even[i] = a[idx++];
    odd[i] = a[idx++];    
  }
  // recursive calls, with w^2 as root of unity for the halved length
  // note that in the next recursion level, a[] plays the role of wsp.
  inplaceCooleyTukeyNTT(length, even, multiply(omega, omega, M), a, M);
  inplaceCooleyTukeyNTT(length, odd, multiply(omega, omega, M), a, M);
  // merge phase
  int h, x = 1;
  for (int i=0; i < length; i++) {
    h = multiply(x, odd[i], M);
    a[i] = add(even[i], h, M);
    a[i+length] = add(even[i], M - h, M);
    x = multiply(x, omega, M);
  }
}

static void inplaceNTT(int length, int *x, int omega, int M) {
  // note that length is expected to be a power of two!
  int *wsp = malloc(length * sizeof(int)); // workspace array
  inplaceCooleyTukeyNTT(length, x, omega, wsp, M);
  free(wsp);
}

static void inplaceInverseNTT(int length, int *x, int omega, int invN, int M) {
  int *wsp = malloc(length*sizeof(int)); // workspace array 
  inplaceCooleyTukeyNTT(length, x, omega, wsp, M);
  free(wsp);
  for (int i=0; i < length; i++) x[i] = multiply(x[i], invN, M);
}


static void nttConvolution(int lenH, int *h, int lenX, int *x, int *y, int M, int omegas[][2], int nInv[]) {
  int lenY = lenX + lenH - 1;
  int len2 = powerOfTwo(lenY);
  
  int log = msb(len2);
  int omega = omegas[log][0];
  
  int *cx = calloc(2*len2, sizeof(int));
  int *ch = cx + len2;
  for (int i=0; i < lenX; i++) cx[i] = x[i];
  for (int i=0; i < lenH; i++) ch[i] = h[i];
  inplaceNTT(len2, ch, omega, M);
  inplaceNTT(len2, cx, omega, M);
  for (int i=0; i < len2; i++) {
    cx[i] = multiply(cx[i], ch[i], M);
  }
  inplaceInverseNTT(len2, cx, omegas[log][1], nInv[log], M);
  for (int i=0; i < lenY; i++) y[i] = (int)(cx[i]); // round to integers
  free(cx);
}


static void chinese_remainder(int *n, int **a, int nrMods, long *y, int yLen)
{
	long prod = 7516372993, p[2] = {5010762410, 2505610584};
  
	for (int i = 0; i < nrMods; ++i) {
    for (int j = 0; j < yLen; ++j) {
      y[j] += ((long) a[i][j] * p[i]) % prod;
      if (y[j] > prod) {
        y[j] -= prod;
      }
      
    }
	}
}


void nttcrtConv(int lenH, int *h, int lenX, int *x, int *y) {
  int lenY = lenH + lenX - 1;
  int const nrMods = 2;
  int mods[] = {65537, 114689};
  
  int *nInv[] = {nInv65537, nInv114689};
  int (*omegas[])[] = {omegas65537, omegas114689};
  
  int **ym = malloc(nrMods * sizeof(int *));  
  for (int i = 0; i < nrMods; ++i) {
    ym[i] = malloc(lenY * sizeof(int));
    nttConvolution(lenH, h, lenX, x, ym[i], mods[i], omegas[i], nInv[i]);
  }
  
  long *yL = calloc(lenY, sizeof(long));
  
  chinese_remainder(mods, ym, nrMods, yL, lenY);
  
  for (int i = 0; i < lenY; ++i) {
    y[i] = (int) yL[i];
  }
  
  for (int i = 0; i < nrMods; ++i) {
    free(ym[i]);
  }
  free(ym);
  free(yL);
}
