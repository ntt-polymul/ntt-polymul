void transpose_n1(__m256i *M)
{
	register __m256i r0, r1, r2, r3, r4, r5, r6, r7, r8, r9, r10, r11;
	register __m256i temp, temp0, temp1, temp2;

	r0 = _mm256_unpacklo_epi16(M[0], M[1]);
	r1 = _mm256_unpacklo_epi16(M[2], M[3]);
	r2 = _mm256_unpacklo_epi16(M[4], M[5]);
	r3 = _mm256_unpacklo_epi16(M[6], M[7]);
	r4 = _mm256_unpacklo_epi16(M[8], M[9]);
	r5 = _mm256_unpacklo_epi16(M[10], M[11]);
	r6 = _mm256_unpacklo_epi16(M[12], M[13]);
	r7 = _mm256_unpacklo_epi16(M[14], M[15]);

	temp = _mm256_unpacklo_epi32(r0, r1);
	temp0 = _mm256_unpacklo_epi32(r2, r3);
	temp1 = _mm256_unpacklo_epi32(r4, r5);
	temp2 = _mm256_unpacklo_epi32(r6, r7);

	r8 = _mm256_unpackhi_epi32(r0, r1);
	r9 = _mm256_unpackhi_epi32(r2, r3);
	r10 = _mm256_unpackhi_epi32(r4, r5);
	r11 = _mm256_unpackhi_epi32(r6, r7);

	r0 = _mm256_unpacklo_epi64(temp, temp0);
	r2 = _mm256_unpackhi_epi64(temp, temp0);

	r1 = _mm256_unpacklo_epi64(temp1, temp2);
	r3 = _mm256_unpackhi_epi64(temp1, temp2);

	temp = _mm256_unpackhi_epi16(M[0], M[1]);
	temp0 = _mm256_unpackhi_epi16(M[2], M[3]);
	temp1 = _mm256_unpackhi_epi16(M[4], M[5]);
	temp2 = _mm256_unpackhi_epi16(M[6], M[7]);
	r4 = _mm256_unpackhi_epi16(M[8], M[9]);

	M[0] = _mm256_permute2f128_si256(r0, r1, 0x20);
	M[8] = _mm256_permute2f128_si256(r0, r1, 0x31);
	M[1] = _mm256_permute2f128_si256(r2, r3, 0x20);
	M[9] = _mm256_permute2f128_si256(r2, r3, 0x31);

	r5 = _mm256_unpackhi_epi16(M[10], M[11]);
	r6 = _mm256_unpackhi_epi16(M[12], M[13]);
	r7 = _mm256_unpackhi_epi16(M[14], M[15]);

	r0 = _mm256_unpacklo_epi64(r8, r9);
	r1 = _mm256_unpacklo_epi64(r10, r11);

	r2 = _mm256_unpackhi_epi64(r8, r9);
	r3 = _mm256_unpackhi_epi64(r10, r11);

	M[3] = _mm256_permute2f128_si256(r2, r3, 0x20);
	M[11] = _mm256_permute2f128_si256(r2, r3, 0x31);
	M[2] = _mm256_permute2f128_si256(r0, r1, 0x20);
	M[10] = _mm256_permute2f128_si256(r0, r1, 0x31);

	r0 = _mm256_unpacklo_epi32(temp, temp0);
	r1 = _mm256_unpacklo_epi32(temp1, temp2);
	r2 = _mm256_unpacklo_epi32(r4, r5);
	r3 = _mm256_unpacklo_epi32(r6, r7);

	r8 = _mm256_unpacklo_epi64(r0, r1);
	r10 = _mm256_unpackhi_epi64(r0, r1);

	r9 = _mm256_unpacklo_epi64(r2, r3);
	r11 = _mm256_unpackhi_epi64(r2, r3);

	M[4] = _mm256_permute2f128_si256(r8, r9, 0x20);
	M[12] = _mm256_permute2f128_si256(r8, r9, 0x31);
	M[5] = _mm256_permute2f128_si256(r10, r11, 0x20);
	M[13] = _mm256_permute2f128_si256(r10, r11, 0x31);

	r0 = _mm256_unpackhi_epi32(temp, temp0);
	r1 = _mm256_unpackhi_epi32(temp1, temp2);
	r2 = _mm256_unpackhi_epi32(r4, r5);
	r3 = _mm256_unpackhi_epi32(r6, r7);

	r4 = _mm256_unpacklo_epi64(r0, r1);
	r6 = _mm256_unpackhi_epi64(r0, r1);

	r5 = _mm256_unpacklo_epi64(r2, r3);
	r7 = _mm256_unpackhi_epi64(r2, r3);

	//-------------------------------------------------------

	M[6] = _mm256_permute2f128_si256(r4, r5, 0x20);
	M[14] = _mm256_permute2f128_si256(r4, r5, 0x31);
	M[7] = _mm256_permute2f128_si256(r6, r7, 0x20);
	M[15] = _mm256_permute2f128_si256(r6, r7, 0x31);
}