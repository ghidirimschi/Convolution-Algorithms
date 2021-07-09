#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define MIN_REPEAT 8

static const unsigned PRIME_BASE = 257;
static const unsigned PRIME_MOD = 1000000007;

typedef long long ll;

static long powMod(ll base, long exp)
{
    long t = 1L;
    while (exp > 0) 
    {
        if (exp % 2 != 0) {
            t = (t * base) % PRIME_MOD;
        }
        base = (base * base) % PRIME_MOD;
        exp /= 2;
    }
    return t % PRIME_MOD;
}
  
typedef struct {
  ll hashVal;
  int hashPos;
} HashStruct;


static ll hash(int *s, int k) {
  ll ret = 0;
  for (int i = 0; i < k; i++)
  {
      ret = ((ret * PRIME_BASE) % PRIME_MOD + s[i]) % PRIME_MOD;
  }
  if (ret < 0) {
    ret += PRIME_MOD;
  }
  return ret;
}

static int hashCmp(const void * a, const void * b) {
   HashStruct *aH = (HashStruct*)a, *bH = (HashStruct*)b;
   if (aH->hashVal == bH->hashVal) {
     return aH->hashPos - bH->hashPos;
   }
   return aH->hashVal - bH->hashVal;
}

static ll absll(ll a) {
  return a < 0 ? -a : a;
}

static int substrings(int k, int *x, int lenX, int *other) {
  HashStruct *hashes = malloc((lenX - k + 1) * sizeof(HashStruct));
  ll powK = powMod(PRIME_BASE, k); 
  hashes[0].hashVal = hash(x, k);
  hashes[0].hashPos = 0;
  ll tmpHash;
  for (int i = 0; i < lenX - k; ++i) {
    tmpHash = (hashes[i].hashVal * PRIME_BASE + x[k + i]) % PRIME_MOD;
    tmpHash -= (powK * x[i]) % PRIME_MOD;
    tmpHash %= PRIME_MOD;
    if (tmpHash < 0) {
      tmpHash += PRIME_MOD; 
    }
    hashes[i+1].hashVal = tmpHash;
    hashes[i+1].hashPos = i+1;
  }
  qsort(hashes, lenX - k + 1, sizeof(HashStruct), hashCmp);
  tmpHash = hashes[0].hashVal;
  int tmpPos = hashes[0].hashPos;
  for (int i = 1; i <= lenX - k; ++i) {
    if (hashes[i].hashVal == tmpHash) {
      if (absll(hashes[i].hashPos - tmpPos) >= k) {
        *other = hashes[i].hashPos;
        free(hashes);
        return tmpPos;
      }
    } else {
      tmpHash = hashes[i].hashVal;
      tmpPos = hashes[i].hashPos;
    }
  }
  free(hashes);
  return -1;  
}


static int findLongestSubstring(int *x, int lenX, int *s1, int *s2) {
  int l = 0, h = lenX, m;
  int ts1, ts2;
  int k = -1; 
  
  while (l <= h && h >= MIN_REPEAT) {
    m = l + (h-l)/2;
    ts1 = substrings(m, x, lenX, &ts2);
    if (ts1 != -1) {
      *s1 = ts1;
      *s2 = ts2;
      l = m + 1;
      k = m;
    } else {
      h = m - 1;
    }
  }
  return k;
} 

void inputSideConvolution(int lenH, int *h, int lenX, int *x, int *y);
void fftConvolution(int lenH, int *h, int lenX, int *x, int *y);



static void pickConvolution(int lenH, int *h, int lenX, int *x, int *y) {
  if (lenH < 1000 || lenX < 1000) {
    inputSideConvolution(lenH, h, lenX, x, y);
    return;
  }
 
  fftConvolution(lenH, h, lenX, x, y);
  
} 


void lrnopConvolutionSort(int lenH, int *h, int lenX, int *x, int *y) {
  int s1, s2, lenY = lenX + lenH - 1;
  int k = findLongestSubstring(h, lenH, &s1, &s2);
  memset(y, 0, lenY * sizeof(int));
  if (k == -1) {
    pickConvolution(lenH, h, lenX, x, y);
    return;
  }
  int *tY = malloc(lenY * sizeof(int));
  if (s1 != 0) {
    pickConvolution(s1, h, lenX, x, tY);
    for (int i = 0; i < s1 + lenX - 1; ++i) {
      y[i] += tY[i];
    }
  }
  
  lrnopConvolutionSort(k, h + s1, lenX, x, tY);
  for (int i = 0; i < lenX + k - 1; ++i) {
    y[s1 + i] += tY[i];
    y[s2 + i] += tY[i];
  }
  
  if (s2 > s1 + k) {
    pickConvolution(s2 - s1 - k, h + s1 + k, lenX, x, tY);
    for (int i = 0; i < s2 - s1 - k + lenX - 1; ++i) {
      y[s1 + k + i] += tY[i];
    }
  }
  
  if (s2 + k < lenH) {
    pickConvolution(lenH - s2 - k, h + s2 + k, lenX, x, tY);
    for (int i = 0; i < lenH - s2 - k + lenX - 1; ++i) {
      y[s2 + k + i] += tY[i]; 
    }
  }
  free(tY);
}
