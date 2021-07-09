#include <stdlib.h>
#include <stdint.h>

static const int64_t omegas[][2] = {
  {1, 1},
  {9223372036737335296, 9223372036737335296},
	{5821608562908302673, 3401763473829032624},
	{6183202344127710422, 4637966906727796075},
	{8232954128065750752, 5263173472767262291},
	{2506374561187543290, 2022096071381403462},
	{5175679149588297482, 2345448649146739134},
	{4529680260184893260, 7628857276699979704},
	{4693975314133565746, 7164255507103637110},
	{250208342857496214, 8429686828039810538},
	{4469250536659403476, 7714104868670012898},
	{1890336861624766990, 3124841710923827579},
	{935815637929945223, 4634860353951294624},
	{1236403532507684466, 2872752946200244613},
	{2949988598460004147, 2751275275530404409},
	{3207215361743254657, 5067346905827978431},
	{3793264300966149197, 514542188484315801},
	{5841790927866575890, 7348949947088264256},
	{2514491450482062217, 1108796242057739048},
	{3993053986736158084, 7387127451862881604},
	{6123203078656286151, 3990408202324679787},
	{7818225862990107593, 7306313761062697635},
	{4224198908362084779, 3276901089038766973},
	{4556310295883022695, 6050139748279733103},
  {2419180138865645092, 6543991442235915766}
};

static const int64_t invN[] = {
  1,
  4611686018368667649,
  6917529027553001473,
  8070450532145168385,
  8646911284441251841,
  8935141660589293569,
  9079256848663314433,
  9151314442700324865,
  9187343239718830081,
  9205357638228082689,
  9214364837482708993,
  9218868437110022145,
  9221120236923678721,
  9222246136830507009,
  9222809086783921153,
  9223090561760628225,
  9223231299248981761,
  9223301667993158529,
  9223336852365246913,
  9223354444551291105,
  9223363240644313201,
  9223367638690824249,
  9223369837714079773,
  9223370937225707535,
  9223371486981521416,
  4611685743490760708
};



static const int64_t M = 9223372036737335297;

static int64_t add(int64_t a, int64_t b) {
  a -= M;
  a += b;
  a += (a >> 63) & M;
  return a;
}

/*
int64_t multiply(int64_t a, int64_t b) {
  return ((__int128) a * b) % M;
}*/


static int64_t multiplyHighUnsigned(int64_t x, int64_t y) {
    int64_t x_high = (uint64_t) x >> 32;
    int64_t y_high = (uint64_t) y >> 32;
    int64_t x_low = x & 0xFFFFFFFFL;
    int64_t y_low = y & 0xFFFFFFFFL;
        
    int64_t z2 = x_low * y_low;
    int64_t t = x_high * y_low + ((uint64_t) z2 >> 32);
    int64_t z1 = t & 0xFFFFFFFFL;
    int64_t z0 = (uint64_t) t >> 32;
    z1 += x_low * y_high;
    return x_high * y_high + z0 + ((uint64_t) z1 >> 32);
}

static int64_t multiply(int64_t a, int64_t b) {
  int64_t xh = multiplyHighUnsigned(a, b); // high word of product
  int64_t xl = a * b; // low word of product
  int64_t BARR_R = -M;
  int64_t xrh = multiplyHighUnsigned(xh, BARR_R); // high word of xr
  int64_t xrm = multiplyHighUnsigned(xl, BARR_R); // middle word of xr, first part
  int64_t add = xh * BARR_R; // second part of middle word
  xrm += add; // add them
  if((add ^ (1L << 63)) > ((xrm ^ (1L << 63)))) {
    xrh++; // carry, see note 1
  }

  int64_t t = xl - ((xrh << 2) | ((uint64_t) xrm >> 62)) * M - M; // see note 2
  t += (t >> 63) & M;
  return t;
}


static int powerOfTwo(int n) {
  int p2;
  for (p2=1; p2<n; p2*=2);
  return p2;
}

static int msb(int x) {
  int r = 0;
  while (x >>= 1) {
    r++;
  }
  return r;
}

static void inplaceCooleyTukeyNTT(int length, int64_t *a, int64_t omega, int64_t *wsp) {
  // wsp is a workspace array. This is used to avoid abundant allocas and frees of memory.
  if (length < 2) return;
  // split array a in even and odd part
  length /= 2;
  int64_t *even = wsp;
  int64_t *odd = wsp + length;
  int idx = 0;  
  for (int i=0; i < length; i++) {
    even[i] = a[idx++];
    odd[i] = a[idx++];    
  }
  // recursive calls, with w^2 as root of unity for the halved length
  // note that in the next recursion level, a[] plays the role of wsp.
  inplaceCooleyTukeyNTT(length, even, multiply(omega, omega), a);
  inplaceCooleyTukeyNTT(length, odd, multiply(omega, omega), a);
  // merge phase
  int64_t h, x = 1;
  for (int i=0; i < length; i++) {
    h = multiply(x, odd[i]);
    a[i] = add(even[i], h);
    a[i+length] = add(even[i], M - h);
    x = multiply(x, omega);
  }
}

static void inplaceNTT(int length, int64_t *x) {
  // note that length is expected to be a power of two!
  int64_t *wsp = malloc(length * sizeof(int64_t)); // workspace array
  int64_t omega = omegas[msb(length)][0];
  inplaceCooleyTukeyNTT(length, x, omega, wsp);
  free(wsp);
}



static void inplaceInverseNTT(int length, int64_t *x) {
  // note that length is expected to be a power of two!
  int64_t *wsp = malloc(length*sizeof(int64_t)); // workspace array
  int64_t omega = omegas[msb(length)][1];     
  inplaceCooleyTukeyNTT(length, x, omega, wsp);
  free(wsp);
  for (int i=0; i < length; i++) x[i] = multiply(x[i], invN[msb(length)]);
}

void nttConvolutionBarret(int lenH, int *h, int lenX, int *x, int *y) {
  int lenY = lenX + lenH - 1;
  int len2 = powerOfTwo(lenY);
  int64_t *cx = calloc(2*len2, sizeof(int64_t));
  int64_t *ch = cx + len2;
  for (int i=0; i < lenX; i++) cx[i] = x[i];
  for (int i=0; i < lenH; i++) ch[i] = h[i];
  inplaceNTT(len2, ch);
  inplaceNTT(len2, cx);
  for (int i=0; i < len2; i++) {
    cx[i] = multiply(cx[i], ch[i]);
  }
  inplaceInverseNTT(len2, cx);
  for (int i=0; i < lenY; i++) y[i] = (int)(cx[i]); // round to integers
  free(cx);
}
