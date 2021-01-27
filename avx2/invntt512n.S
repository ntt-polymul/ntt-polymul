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

.macro levels0t5 off
vmovdqa		(128*\off+  0)*2(%rsi),%ymm8
vmovdqa		(128*\off+ 16)*2(%rsi),%ymm5
vmovdqa		(128*\off+ 32)*2(%rsi),%ymm10
vmovdqa		(128*\off+ 48)*2(%rsi),%ymm7
vmovdqa		(128*\off+ 64)*2(%rsi),%ymm3
vmovdqa		(128*\off+ 80)*2(%rsi),%ymm9
vmovdqa		(128*\off+ 96)*2(%rsi),%ymm4
vmovdqa		(128*\off+112)*2(%rsi),%ymm11

/* level 0 */
update		6,8,10,3,4,5,7,9,11

vbroadcasti128	(_ZETAS_PINV+256-32*\off-8)*2(%rdx),%ymm1
vbroadcasti128	(_ZETAS+256-32*\off-8)*2(%rdx),%ymm2
vpshufb		%ymm15,%ymm1,%ymm1
vpshufb		%ymm15,%ymm2,%ymm2
fqmulprecomp	1,2,5,x=4,neg=1
vbroadcasti128	(_ZETAS_PINV+256-32*\off-16)*2(%rdx),%ymm1
vbroadcasti128	(_ZETAS+256-32*\off-16)*2(%rdx),%ymm2
vpshufb		%ymm15,%ymm1,%ymm1
vpshufb		%ymm15,%ymm2,%ymm2
fqmulprecomp	1,2,7,x=4,neg=1
vbroadcasti128	(_ZETAS_PINV+256-32*\off-24)*2(%rdx),%ymm1
vbroadcasti128	(_ZETAS+256-32*\off-24)*2(%rdx),%ymm2
vpshufb		%ymm15,%ymm1,%ymm1
vpshufb		%ymm15,%ymm2,%ymm2
fqmulprecomp	1,2,9,x=4,neg=1
vbroadcasti128	(_ZETAS_PINV+256-32*\off-32)*2(%rdx),%ymm1
vbroadcasti128	(_ZETAS+256-32*\off-32)*2(%rdx),%ymm2
vpshufb		%ymm15,%ymm1,%ymm1
vpshufb		%ymm15,%ymm2,%ymm2
fqmulprecomp	1,2,11,x=4,neg=1

pack		6,5,4,5
pack		8,7,6,7
pack		10,9,8,9
pack		3,11,10,11

/* level 1 */
update		3,4,6,8,10,5,7,9,11

vpbroadcastq	(_ZETAS_PINV+128-16*\off-4)*2(%rdx),%ymm1
vpbroadcastq	(_ZETAS+128-16*\off-4)*2(%rdx),%ymm2
vpshufb		%ymm15,%ymm1,%ymm1
vpshufb		%ymm15,%ymm2,%ymm2
fqmulprecomp	1,2,5,x=10,neg=1
vpbroadcastq	(_ZETAS_PINV+128-16*\off-8)*2(%rdx),%ymm1
vpbroadcastq	(_ZETAS+128-16*\off-8)*2(%rdx),%ymm2
vpshufb		%ymm15,%ymm1,%ymm1
vpshufb		%ymm15,%ymm2,%ymm2
fqmulprecomp	1,2,7,x=10,neg=1
vpbroadcastq	(_ZETAS_PINV+128-16*\off-12)*2(%rdx),%ymm1
vpbroadcastq	(_ZETAS+128-16*\off-12)*2(%rdx),%ymm2
vpshufb		%ymm15,%ymm1,%ymm1
vpshufb		%ymm15,%ymm2,%ymm2
fqmulprecomp	1,2,9,x=10,neg=1
vpbroadcastq	(_ZETAS_PINV+128-16*\off-16)*2(%rdx),%ymm1
vpbroadcastq	(_ZETAS+128-16*\off-16)*2(%rdx),%ymm2
vpshufb		%ymm15,%ymm1,%ymm1
vpshufb		%ymm15,%ymm2,%ymm2
fqmulprecomp	1,2,11,x=10,neg=1

fqmulprecomp	12,13,3,x=10  // extra reduction
fqmulprecomp	12,13,4,x=10 // extra reduction
fqmulprecomp	12,13,6,x=10  // extra reduction
fqmulprecomp	12,13,8,x=10  // extra reduction

pack		3,5,10,5
pack		4,7,3,7
pack		6,9,4,9
pack		8,11,6,11

/* level 2 */
update		8,10,3,4,6,5,7,9,11

vpbroadcastd	(_ZETAS_PINV+64-8*\off-2)*2(%rdx),%ymm1
vpbroadcastd	(_ZETAS+64-8*\off-2)*2(%rdx),%ymm2
vpshufb		%ymm15,%ymm1,%ymm1
vpshufb		%ymm15,%ymm2,%ymm2
fqmulprecomp	1,2,5,x=6,neg=1
vpbroadcastd	(_ZETAS_PINV+64-8*\off-4)*2(%rdx),%ymm1
vpbroadcastd	(_ZETAS+64-8*\off-4)*2(%rdx),%ymm2
vpshufb		%ymm15,%ymm1,%ymm1
vpshufb		%ymm15,%ymm2,%ymm2
fqmulprecomp	1,2,7,x=6,neg=1
vpbroadcastd	(_ZETAS_PINV+64-8*\off-6)*2(%rdx),%ymm1
vpbroadcastd	(_ZETAS+64-8*\off-6)*2(%rdx),%ymm2
vpshufb		%ymm15,%ymm1,%ymm1
vpshufb		%ymm15,%ymm2,%ymm2
fqmulprecomp	1,2,9,x=6,neg=1
vpbroadcastd	(_ZETAS_PINV+64-8*\off-8)*2(%rdx),%ymm1
vpbroadcastd	(_ZETAS+64-8*\off-8)*2(%rdx),%ymm2
vpshufb		%ymm15,%ymm1,%ymm1
vpshufb		%ymm15,%ymm2,%ymm2
fqmulprecomp	1,2,11,x=6,neg=1

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

/* level 3 */
update		4,6,8,10,3,5,7,9,11

vpbroadcastw	(_ZETAS_PINV+32-4*\off-1)*2(%rdx),%ymm1
vpbroadcastw	(_ZETAS+32-4*\off-1)*2(%rdx),%ymm2
fqmulprecomp	1,2,5,x=3,neg=1
vpbroadcastw	(_ZETAS_PINV+32-4*\off-2)*2(%rdx),%ymm1
vpbroadcastw	(_ZETAS+32-4*\off-2)*2(%rdx),%ymm2
fqmulprecomp	1,2,7,x=3,neg=1
vpbroadcastw	(_ZETAS_PINV+32-4*\off-3)*2(%rdx),%ymm1
vpbroadcastw	(_ZETAS+32-4*\off-3)*2(%rdx),%ymm2
fqmulprecomp	1,2,9,x=3,neg=1
vpbroadcastw	(_ZETAS_PINV+32-4*\off-4)*2(%rdx),%ymm1
vpbroadcastw	(_ZETAS+32-4*\off-4)*2(%rdx),%ymm2
fqmulprecomp	1,2,11,x=3,neg=1

fqmulprecomp	12,13,4,x=3  // extra reduction
fqmulprecomp	12,13,6,x=3  // extra reduction
fqmulprecomp	12,13,8,x=3  // extra reduction
fqmulprecomp	12,13,10,x=3  // extra reduction

/* level 4 */
update		3,4,5,8,9,6,7,10,11

vpbroadcastw	(_ZETAS_PINV+16-2*\off-1)*2(%rdx),%ymm1
vpbroadcastw	(_ZETAS+16-2*\off-1)*2(%rdx),%ymm2
fqmulprecomp	1,2,6,x=9,neg=1
fqmulprecomp	1,2,7,x=9,neg=1
vpbroadcastw	(_ZETAS_PINV+16-2*\off-2)*2(%rdx),%ymm1
vpbroadcastw	(_ZETAS+16-2*\off-2)*2(%rdx),%ymm2
fqmulprecomp	1,2,10,x=9,neg=1
fqmulprecomp	1,2,11,x=9,neg=1

/* level 5 */
update		9,3,4,6,7,5,8,10,11

vpbroadcastw	(_ZETAS_PINV+8-\off-1)*2(%rdx),%ymm1
vpbroadcastw	(_ZETAS+8-\off-1)*2(%rdx),%ymm2
fqmulprecomp	1,2,5,x=7,neg=1
fqmulprecomp	1,2,8,x=7,neg=1
fqmulprecomp	1,2,10,x=7,neg=1
fqmulprecomp	1,2,11,x=7,neg=1

fqmulprecomp	12,13,9,x=7  // extra reduction
fqmulprecomp	12,13,3,x=7  // extra reduction
fqmulprecomp	12,13,4,x=7  // extra reduction
fqmulprecomp	12,13,6,x=7  // extra reduction

vmovdqa		%ymm9,(128*\off+  0)*2(%rdi)
vmovdqa		%ymm3,(128*\off+ 16)*2(%rdi)
vmovdqa		%ymm4,(128*\off+ 32)*2(%rdi)
vmovdqa		%ymm6,(128*\off+ 48)*2(%rdi)
vmovdqa		%ymm5,(128*\off+ 64)*2(%rdi)
vmovdqa		%ymm8,(128*\off+ 80)*2(%rdi)
vmovdqa		%ymm10,(128*\off+ 96)*2(%rdi)
vmovdqa		%ymm11,(128*\off+112)*2(%rdi)
.endm

.macro levels6t7 off
vmovdqa		(  0+16*\off)*2(%rdi),%ymm4
vmovdqa		( 64+16*\off)*2(%rdi),%ymm5
vmovdqa		(128+16*\off)*2(%rdi),%ymm6
vmovdqa		(192+16*\off)*2(%rdi),%ymm7
vmovdqa		(256+16*\off)*2(%rdi),%ymm8
vmovdqa		(320+16*\off)*2(%rdi),%ymm9
vmovdqa		(384+16*\off)*2(%rdi),%ymm10
vmovdqa		(448+16*\off)*2(%rdi),%ymm11

/* level6 */
update		3,4,5,8,9,6,7,10,11

vpbroadcastw	(_ZETAS_PINV+3)*2(%rdx),%ymm1
vpbroadcastw	(_ZETAS+3)*2(%rdx),%ymm2
fqmulprecomp	1,2,6,x=9,neg=1
fqmulprecomp	1,2,7,x=9,neg=1
vpbroadcastw	(_ZETAS_PINV+2)*2(%rdx),%ymm1
vpbroadcastw	(_ZETAS+2)*2(%rdx),%ymm2
fqmulprecomp	1,2,10,x=9,neg=1
fqmulprecomp	1,2,11,x=9,neg=1

/* level7 */
update		9,3,4,6,7,5,8,10,11

vpbroadcastw	(_ZETAS_PINV+1)*2(%rdx),%ymm1
vpbroadcastw	(_ZETAS+1)*2(%rdx),%ymm2
fqmulprecomp	1,2,5,x=7,neg=1
fqmulprecomp	1,2,8,x=7,neg=1
fqmulprecomp	1,2,10,x=7,neg=1
fqmulprecomp	1,2,11,x=7,neg=1

vmovdqa		%ymm9,(  0+16*\off)*2(%rdi)
vmovdqa		%ymm3,( 64+16*\off)*2(%rdi)
vmovdqa		%ymm4,(128+16*\off)*2(%rdi)
vmovdqa		%ymm6,(192+16*\off)*2(%rdi)
vmovdqa		%ymm5,(256+16*\off)*2(%rdi)
vmovdqa		%ymm8,(320+16*\off)*2(%rdi)
vmovdqa		%ymm10,(384+16*\off)*2(%rdi)
vmovdqa		%ymm11,(448+16*\off)*2(%rdi)
.endm

.text
.global cdecl(poly_invntt_tomont)
cdecl(poly_invntt_tomont):
vmovdqa         _16XP*2(%rdx),%ymm0
vmovdqa		cdecl(idxdata)+_REVWIDX(%rip),%ymm15
vmovdqa		cdecl(idxdata)+_PACKWIDX(%rip),%ymm14
vmovdqa		(_16XMONT_PINV)*2(%rdx),%ymm12
vmovdqa		(_16XMONT)*2(%rdx),%ymm13

levels0t5	0
levels0t5	1
levels0t5	2
levels0t5	3

levels6t7	0
levels6t7	1
levels6t7	2
levels6t7	3

ret
