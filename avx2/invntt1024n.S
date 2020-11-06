#include "consts1024.h"
.include "shuffle.inc"
.include "fq.inc"

.macro update rln,rl0,rl1,rl2,rl3,rh0,rh1,rh2,rh3
vpaddw		%ymm\rh0,%ymm\rl0,%ymm\rln
vpsubw		%ymm\rh0,%ymm\rl0,%ymm\rh0

vpaddw		%ymm\rh1,%ymm\rl1,%ymm\rl0
vpsubw		%ymm\rh1,%ymm\rl1,%ymm\rh1

vpaddw		%ymm\rh2,%ymm\rl2,%ymm\rl1
vpsubw		%ymm\rh2,%ymm\rl2,%ymm\rh2

vpaddw		%ymm\rh3,%ymm\rl3,%ymm\rl2
vpsubw		%ymm\rh3,%ymm\rl3,%ymm\rh3
.endm

.macro levels0t4 off
vmovdqa		(128*\off+  0)*2(%rsi),%ymm6
vmovdqa		(128*\off+ 16)*2(%rsi),%ymm9
vmovdqa		(128*\off+ 32)*2(%rsi),%ymm8
vmovdqa		(128*\off+ 48)*2(%rsi),%ymm10
vmovdqa		(128*\off+ 64)*2(%rsi),%ymm7
vmovdqa		(128*\off+ 80)*2(%rsi),%ymm3
vmovdqa		(128*\off+ 96)*2(%rsi),%ymm4
vmovdqa		(128*\off+112)*2(%rsi),%ymm11

shuffle8	6,9,5,9
shuffle8	8,10,6,10
shuffle8	7,3,8,3
shuffle8	4,11,7,11

shuffle4	5,6,4,6
shuffle4	8,7,5,7
shuffle4	9,10,8,10
shuffle4	3,11,9,11

/* level0 */
update		3,4,6,8,10,5,7,9,11

vpbroadcastq	(_ZETAS_PINV+256-16*\off-4)*2(%rdx),%ymm1
vpbroadcastq	(_ZETAS+256-16*\off-4)*2(%rdx),%ymm2
vpshufb		%ymm15,%ymm1,%ymm1
vpshufb		%ymm15,%ymm2,%ymm2
fqmulprecomp	1,2,5,x=10,neg=1
vpbroadcastq	(_ZETAS_PINV+256-16*\off-8)*2(%rdx),%ymm1
vpbroadcastq	(_ZETAS+256-16*\off-8)*2(%rdx),%ymm2
vpshufb		%ymm15,%ymm1,%ymm1
vpshufb		%ymm15,%ymm2,%ymm2
fqmulprecomp	1,2,7,x=10,neg=1
vpbroadcastq	(_ZETAS_PINV+256-16*\off-12)*2(%rdx),%ymm1
vpbroadcastq	(_ZETAS+256-16*\off-12)*2(%rdx),%ymm2
vpshufb		%ymm15,%ymm1,%ymm1
vpshufb		%ymm15,%ymm2,%ymm2
fqmulprecomp	1,2,9,x=10,neg=1
vpbroadcastq	(_ZETAS_PINV+256-16*\off-16)*2(%rdx),%ymm1
vpbroadcastq	(_ZETAS+256-16*\off-16)*2(%rdx),%ymm2
vpshufb		%ymm15,%ymm1,%ymm1
vpshufb		%ymm15,%ymm2,%ymm2
fqmulprecomp	1,2,11,x=10,neg=1

pack		3,5,10,5
pack		4,7,3,7
pack		6,9,4,9
pack		8,11,6,11

/* level1 */
update		8,10,3,4,6,5,7,9,11

vpbroadcastd	(_ZETAS_PINV+128-8*\off-2)*2(%rdx),%ymm1
vpbroadcastd	(_ZETAS+128-8*\off-2)*2(%rdx),%ymm2
vpshufb		%ymm15,%ymm1,%ymm1
vpshufb		%ymm15,%ymm2,%ymm2
fqmulprecomp	1,2,5,x=6,neg=1
vpbroadcastd	(_ZETAS_PINV+128-8*\off-4)*2(%rdx),%ymm1
vpbroadcastd	(_ZETAS+128-8*\off-4)*2(%rdx),%ymm2
vpshufb		%ymm15,%ymm1,%ymm1
vpshufb		%ymm15,%ymm2,%ymm2
fqmulprecomp	1,2,7,x=6,neg=1
vpbroadcastd	(_ZETAS_PINV+128-8*\off-6)*2(%rdx),%ymm1
vpbroadcastd	(_ZETAS+128-8*\off-6)*2(%rdx),%ymm2
vpshufb		%ymm15,%ymm1,%ymm1
vpshufb		%ymm15,%ymm2,%ymm2
fqmulprecomp	1,2,9,x=6,neg=1
vpbroadcastd	(_ZETAS_PINV+128-8*\off-8)*2(%rdx),%ymm1
vpbroadcastd	(_ZETAS+128-8*\off-8)*2(%rdx),%ymm2
vpshufb		%ymm15,%ymm1,%ymm1
vpshufb		%ymm15,%ymm2,%ymm2
fqmulprecomp	1,2,11,x=6,neg=1

fqmulprecomp	12,13,8,x=6  // extra reduction
fqmulprecomp	12,13,10,x=6  // extra reduction
fqmulprecomp	12,13,3,x=6  // extra reduction
fqmulprecomp	12,13,4,x=6  // extra reduction

pack		8,5,6,5
pack		10,7,8,7
pack		3,9,10,9
pack		4,11,3,11

vmovdqa		cdecl(idxdata)+_UNPCKWIDX(%rip),%ymm1
vpermq		$0xD8,%ymm6,%ymm6
vpshufb		%ymm1,%ymm6,%ymm6
vpermq		$0xD8,%ymm5,%ymm5
vpshufb		%ymm1,%ymm5,%ymm5
vpermq		$0xD8,%ymm8,%ymm8
vpshufb		%ymm1,%ymm8,%ymm8
vpermq		$0xD8,%ymm7,%ymm7
vpshufb		%ymm1,%ymm7,%ymm7
vpermq		$0xD8,%ymm10,%ymm10
vpshufb		%ymm1,%ymm10,%ymm10
vpermq		$0xD8,%ymm9,%ymm9
vpshufb		%ymm1,%ymm9,%ymm9
vpermq		$0xD8,%ymm3,%ymm3
vpshufb		%ymm1,%ymm3,%ymm3
vpermq		$0xD8,%ymm11,%ymm11
vpshufb		%ymm1,%ymm11,%ymm11

/* level2 */
update		4,6,8,10,3,5,7,9,11

vpbroadcastw	(_ZETAS_PINV+63-4*\off)*2(%rdx),%ymm1
vpbroadcastw	(_ZETAS+63-4*\off)*2(%rdx),%ymm2
fqmulprecomp	1,2,5,x=3,neg=1
vpbroadcastw	(_ZETAS_PINV+63-4*\off-1)*2(%rdx),%ymm1
vpbroadcastw	(_ZETAS+63-4*\off-1)*2(%rdx),%ymm2
fqmulprecomp	1,2,7,x=3,neg=1
vpbroadcastw	(_ZETAS_PINV+63-4*\off-2)*2(%rdx),%ymm1
vpbroadcastw	(_ZETAS+63-4*\off-2)*2(%rdx),%ymm2
fqmulprecomp	1,2,9,x=3,neg=1
vpbroadcastw	(_ZETAS_PINV+63-4*\off-3)*2(%rdx),%ymm1
vpbroadcastw	(_ZETAS+63-4*\off-3)*2(%rdx),%ymm2
fqmulprecomp	1,2,11,x=3,neg=1

/* level3 */
update		3,4,5,8,9,6,7,10,11

vpbroadcastw	(_ZETAS_PINV+31-2*\off)*2(%rdx),%ymm1
vpbroadcastw	(_ZETAS+31-2*\off)*2(%rdx),%ymm2
fqmulprecomp	1,2,6,x=9,neg=1
fqmulprecomp	1,2,7,x=9,neg=1
vpbroadcastw	(_ZETAS_PINV+31-2*\off-1)*2(%rdx),%ymm1
vpbroadcastw	(_ZETAS+31-2*\off-1)*2(%rdx),%ymm2
fqmulprecomp	1,2,10,x=9,neg=1
fqmulprecomp	1,2,11,x=9,neg=1

fqmulprecomp	12,13,3,x=9  // extra reduction
fqmulprecomp	12,13,4,x=9  // extra reduction
fqmulprecomp	12,13,5,x=9  // extra reduction
fqmulprecomp	12,13,8,x=9  // extra reduction

/* level4 */
update		9,3,4,6,7,5,8,10,11

vpbroadcastw	(_ZETAS_PINV+15-\off)*2(%rdx),%ymm1
vpbroadcastw	(_ZETAS+15-\off)*2(%rdx),%ymm2
fqmulprecomp	1,2,5,x=7,neg=1
fqmulprecomp	1,2,8,x=7,neg=1
fqmulprecomp	1,2,10,x=7,neg=1
fqmulprecomp	1,2,11,x=7,neg=1

vmovdqa		%ymm9,(128*\off+  0)*2(%rdi)
vmovdqa		%ymm3,(128*\off+ 16)*2(%rdi)
vmovdqa		%ymm4,(128*\off+ 32)*2(%rdi)
vmovdqa		%ymm6,(128*\off+ 48)*2(%rdi)
vmovdqa		%ymm5,(128*\off+ 64)*2(%rdi)
vmovdqa		%ymm8,(128*\off+ 80)*2(%rdi)
vmovdqa		%ymm10,(128*\off+ 96)*2(%rdi)
vmovdqa		%ymm11,(128*\off+112)*2(%rdi)
.endm

.macro levels5t7 off
vmovdqa		(  0+16*\off)*2(%rdi),%ymm4
vmovdqa		(128+16*\off)*2(%rdi),%ymm5
vmovdqa		(256+16*\off)*2(%rdi),%ymm6
vmovdqa		(384+16*\off)*2(%rdi),%ymm7
vmovdqa		(512+16*\off)*2(%rdi),%ymm8
vmovdqa		(640+16*\off)*2(%rdi),%ymm9
vmovdqa		(768+16*\off)*2(%rdi),%ymm10
vmovdqa		(896+16*\off)*2(%rdi),%ymm11

/* level5 */
update		3,4,6,8,10,5,7,9,11

vpbroadcastw	(_ZETAS_PINV+7)*2(%rdx),%ymm1
vpbroadcastw	(_ZETAS+7)*2(%rdx),%ymm2
fqmulprecomp	1,2,5,x=10,neg=1
vpbroadcastw	(_ZETAS_PINV+6)*2(%rdx),%ymm1
vpbroadcastw	(_ZETAS+6)*2(%rdx),%ymm2
fqmulprecomp	1,2,7,x=10,neg=1
vpbroadcastw	(_ZETAS_PINV+5)*2(%rdx),%ymm1
vpbroadcastw	(_ZETAS+5)*2(%rdx),%ymm2
fqmulprecomp	1,2,9,x=10,neg=1
vpbroadcastw	(_ZETAS_PINV+4)*2(%rdx),%ymm1
vpbroadcastw	(_ZETAS+4)*2(%rdx),%ymm2
fqmulprecomp	1,2,11,x=10,neg=1

fqmulprecomp	12,13,3,x=10  // extra reduction
fqmulprecomp	12,13,4,x=10  // extra reduction

/* level6 */
update		10,3,5,6,9,4,7,8,11

vpbroadcastw	(_ZETAS_PINV+3)*2(%rdx),%ymm1
vpbroadcastw	(_ZETAS+3)*2(%rdx),%ymm2
fqmulprecomp	1,2,4,x=9,neg=1
fqmulprecomp	1,2,7,x=9,neg=1
vpbroadcastw	(_ZETAS_PINV+2)*2(%rdx),%ymm1
vpbroadcastw	(_ZETAS+2)*2(%rdx),%ymm2
fqmulprecomp	1,2,8,x=9,neg=1
fqmulprecomp	1,2,11,x=9,neg=1

/* level7 */
update		9,10,3,4,7,5,6,8,11

vpbroadcastw	(_ZETAS_PINV+1)*2(%rdx),%ymm1
vpbroadcastw	(_ZETAS+1)*2(%rdx),%ymm2
fqmulprecomp	1,2,5,x=7,neg=1
fqmulprecomp	1,2,6,x=7,neg=1
fqmulprecomp	1,2,8,x=7,neg=1
fqmulprecomp	1,2,11,x=7,neg=1

vmovdqa		%ymm9,(  0+16*\off)*2(%rdi)
vmovdqa		%ymm10,(128+16*\off)*2(%rdi)
vmovdqa		%ymm3,(256+16*\off)*2(%rdi)
vmovdqa		%ymm4,(384+16*\off)*2(%rdi)
vmovdqa		%ymm5,(512+16*\off)*2(%rdi)
vmovdqa		%ymm6,(640+16*\off)*2(%rdi)
vmovdqa		%ymm8,(768+16*\off)*2(%rdi)
vmovdqa		%ymm11,(896+16*\off)*2(%rdi)
.endm

.text
.global cdecl(poly_invntt_tomont)
cdecl(poly_invntt_tomont):
vmovdqa         _16XP*2(%rdx),%ymm0
vmovdqa		cdecl(idxdata)+_REVWIDX(%rip),%ymm15
vmovdqa		cdecl(idxdata)+_PACKWIDX(%rip),%ymm14
vmovdqa		(_16XMONT_PINV)*2(%rdx),%ymm12
vmovdqa		(_16XMONT)*2(%rdx),%ymm13

levels0t4	0
levels0t4	1
levels0t4	2
levels0t4	3
levels0t4	4
levels0t4	5
levels0t4	6
levels0t4	7

levels5t7	0
levels5t7	1
levels5t7	2
levels5t7	3
levels5t7	4
levels5t7	5
levels5t7	6
levels5t7	7

ret
