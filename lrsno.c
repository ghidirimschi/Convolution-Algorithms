#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "LibList.h"

const int PRIME_BASE = 13;
const int PRIME_MOD = 11099107;

typedef long long ll;

#define MIN_REPEAT  1
List hm[11099107];
List occupied;

ll hash(int *s, int k) {
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

int absll(int a) {
  return a < 0 ? -a : a;
}



int insertHash(int hash, int pos, int *str, int k) {
  List l = hm[hash];
  while (l != NULL) {
    if ((absll(pos - l->item) >= k) && !memcmp(str + pos, str + l->item, k * sizeof(int))) {
      return l->item;
    }
    l = l->next;
  }
  hm[hash] = addItem(pos, hm[hash]);
  occupied = addItem(hash, occupied);
  return -1;
}

void clearHashMap() {
  List occTmp  = occupied;
  while (occTmp != NULL) {
    freeList(hm[occTmp->item]);
    hm[occTmp->item] = NULL;
    occTmp = occTmp->next;
  }
  freeList(occupied);
  occupied = NULL;
}

int powMod(long base, int exp)
{
    int t = 1;
    while (exp > 0) 
    {
        if (exp % 2 != 0) {
            t = (t * base) % PRIME_MOD;
        }
        base = (base * base) % PRIME_MOD;
        exp >>= 1;
    }
    return t % PRIME_MOD;
}

int substrings(int k, int *x, int lenX, int *other) {
  int powK = powMod(PRIME_BASE, k); 
  int prevHash = hash(x, k);
  
  clearHashMap();
  insertHash(prevHash, 0, x, k);
  int tmpHash;
  int rtn;
  for (int i = 0; i < lenX - k; ++i) {
    tmpHash = (prevHash * PRIME_BASE + x[k + i]) % PRIME_MOD;
    tmpHash -= (powK * x[i]); 
    tmpHash %= PRIME_MOD;
    if (tmpHash < 0) {
      tmpHash += PRIME_MOD; 
    }
    rtn = insertHash(tmpHash, i + 1, x, k);
    if (rtn != -1) {
      *other = rtn;
      return i + 1;
    }
    prevHash = tmpHash;
  }
  return -1;  
}

int findLongestSubstring(int *x, int lenX, int *s1, int *s2) {
  int l = 0, h = lenX, m;
  int ts1, ts2;
  int k = -1; 
  
  while (l <= h && h >= MIN_REPEAT) {
    m = l + (h-l)/2;
    ts1 = substrings(m, x, lenX, &ts2);
    if (ts1 != -1) {
      *s1 = ts1 < ts2 ? ts1 : ts2;
      *s2 = ts2 > ts1 ? ts2 : ts1;
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
void oamfftConvolution(int lenH, int *h, int lenX, int *x, int *y);


void pickConvolution(int lenH, int *h, int lenX, int *x, int *y) {
  if (lenH < 1000 || lenX < 1000) {
    inputSideConvolution(lenH, h, lenX, x, y);
    return;
  }
 
  fftConvolution(lenH, h, lenX, x, y);
  
} 


void lrsnoConvolution(int lenH, int *h, int lenX, int *x, int *y) {
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
  
  
  
  lrsnoConvolution(k, h + s1, lenX, x, tY);
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
