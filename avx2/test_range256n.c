#include <stdint.h>
#include "params.h"
#include "nttrange.h"

int main(void) {
  int32_t i;
  uint32_t bounds0[NTT_N], bounds1[NTT_N];

  for(i=0;i<NTT_N;i++) {
    bounds0[i] = 4096;
    bounds1[i] = 5;
  }

  range_mul(bounds0,bounds1);
  return 0;
}
