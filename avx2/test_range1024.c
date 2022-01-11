#include <stdint.h>
#include "params.h"
#include "nttrange.h"

int main(void) {
  int32_t i;
  uint32_t bounds0[NTT_N], bounds1[NTT_N];

  for(i=0;i<KEM_N;i++) {
    bounds0[i] = 2047;
    bounds1[i] = 1;
  }
  for(i=KEM_N;i<NTT_N;i++) {
    bounds0[i] = 0;
    bounds1[i] = 0;
  }

  range_mul(bounds0,bounds1);
  return 0;
}
