/* Taken from the official AVX2 optimized implementation of LAC */

#include <stdint.h>
#include <immintrin.h>
#include "lacmul.h"
#include "../params.h"
#include "../poly.h"

#define DIM_N KEM_N
#define Q KEM_Q
#define BIG_Q (1024*4*Q)

void lac_polysmall_mul(uint8_t *c, const uint8_t *a, const int8_t *b)
{
	int i,j;
	uint8_t v[DIM_N+DIM_N],*v_p;
	__m256i tmp0, tmp1, tmp2, tmp_one;
	int32_t sum_all;
	int32_t r[8];

	tmp_one=_mm256_set_epi16(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1);
	//construct vector of a
	for(i=0;i<DIM_N;i++)
	{
		v[i]=a[DIM_N-1-i];
		v[i+DIM_N]=Q-v[i];
	}

	//compute c[i]
	for(i=0;i<DIM_N;i++)
	{
		v_p=(v+DIM_N-i-1);
		tmp2 =_mm256_setzero_si256();
		for(j=0;j<DIM_N;j+=128)
		{
			tmp0 = _mm256_loadu_si256((__m256i *)(v_p+j));
			tmp1 = _mm256_loadu_si256((__m256i *)(b+j));
			tmp0 = _mm256_maddubs_epi16(tmp0, tmp1);
			tmp2 = _mm256_add_epi16(tmp2, tmp0);

			tmp0 = _mm256_loadu_si256((__m256i *)(v_p+j+32));
			tmp1 = _mm256_loadu_si256((__m256i *)(b+j+32));
			tmp0 = _mm256_maddubs_epi16(tmp0, tmp1);
			tmp2 = _mm256_add_epi16(tmp2, tmp0);

			tmp0 = _mm256_loadu_si256((__m256i *)(v_p+j+64));
			tmp1 = _mm256_loadu_si256((__m256i *)(b+j+64));
			tmp0 = _mm256_maddubs_epi16(tmp0, tmp1);
			tmp2 = _mm256_add_epi16(tmp2, tmp0);

			tmp0 = _mm256_loadu_si256((__m256i *)(v_p+j+96));
			tmp1 = _mm256_loadu_si256((__m256i *)(b+j+96));
			tmp0 = _mm256_maddubs_epi16(tmp0, tmp1);
			tmp2 = _mm256_add_epi16(tmp2, tmp0);
		}

		tmp2 = _mm256_madd_epi16(tmp2,tmp_one);
		_mm256_storeu_si256((__m256i*)r,tmp2);
		sum_all=r[0]+r[1]+r[2]+r[3]+r[4]+r[5]+r[6]+r[7];

		c[i]=(sum_all+BIG_Q)%Q;
	}
}

void lac_poly_mul(poly *r, const poly *a, const poly *b) {
  unsigned int i;
  uint8_t as[KEM_N];
  int8_t bs[KEM_N];

  for(i=0;i<KEM_N;i++) {
    as[i] = a->coeffs[i];
    bs[i] = b->coeffs[i];
  }

  lac_polysmall_mul(as,as,bs);

  for(i=0;i<KEM_N;i++)
    r->coeffs[i] = as[i];
}
