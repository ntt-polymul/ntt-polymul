#include "consts1024.h"
.include "fq.inc"

.macro schoolbook off
vmovdqa		(64*\off+ 0)*2(%rsi),%ymm2	# a0 a1
vmovdqa		(64*\off+16)*2(%rsi),%ymm3	# b0 b1
vmovdqa		(64*\off+32)*2(%rsi),%ymm4	# c0 c1
vmovdqa		(64*\off+48)*2(%rsi),%ymm5	# d0 d1

vmovdqa		(_16XMONT_PINV)*2(%rcx),%ymm6
vmovdqa		(_16XMONT)*2(%rcx),%ymm7
fqmulprecomp	6,7,2,x=15
fqmulprecomp	6,7,3,x=15
fqmulprecomp	6,7,4,x=15
fqmulprecomp	6,7,5,x=15

vpmullw		%ymm1,%ymm2,%ymm6		# (a0 a1).pinv
vpmullw		%ymm1,%ymm3,%ymm7		# (b0 b1).pinv
vpmullw		%ymm1,%ymm4,%ymm8		# (c0 c1).pinv
vpmullw		%ymm1,%ymm5,%ymm9		# (d0 d1).pinv

vmovdqa		(64*\off+ 0)*2(%rdx),%ymm10	# e0 e1
vmovdqa		(64*\off+16)*2(%rdx),%ymm11	# f0 f1
vmovdqa		(64*\off+32)*2(%rdx),%ymm12	# g0 g1
vmovdqa		(64*\off+48)*2(%rdx),%ymm13	# h0 h1

vpermq		$0x4E,%ymm10,%ymm14		# e1 e0
fqmulprecomp	6,2,10,x=15			# a0*e0 a1*e1
fqmulprecomp	6,2,14,x=15			# a0*e1 a1*e0

vpermq		$0x4E,%ymm11,%ymm2		# f1 f0
vpermq		$0x4E,%ymm12,%ymm6		# g1 g0
fqmulprecomp	7,3,11,x=15			# b0*f0 b1*f1
fqmulprecomp	7,3,2,x=15			# b0*f1 b1*f0

vpermq		$0x4E,%ymm13,%ymm3		# h1 h0
fqmulprecomp	8,4,12,x=15			# c0*g0 c1*g1
fqmulprecomp	9,5,13,x=15			# d0*h0 d1*h1
fqmulprecomp	8,4,6,x=15			# c0*g1 c1*g0
fqmulprecomp	9,5,3,x=15			# d0*h1 d1*h0

vperm2i128	$0x31,%ymm12,%ymm10,%ymm4	# a1*e1 c1*g1
vperm2i128	$0x31,%ymm13,%ymm11,%ymm5	# b1*f1 d1*h1

vmovdqa		(_ZETAS_PINV+128+16*\off)*2(%rcx),%ymm7
vmovdqa		(_ZETAS+128+16*\off)*2(%rcx),%ymm8
fqmulprecomp	7,8,4,x=15			# a1*e1*r0 c1*g1*r1
fqmulprecomp	7,8,5,x=15,neg=1		# -b1*f1*r0 -d1*h1*r1

vperm2i128	$0x20,%ymm14,%ymm4,%ymm7	# a1*e1*r0 a0*e1
vpblendd	$0xF0,%ymm14,%ymm10,%ymm10	# a0*e0 a1*e0
vpaddw		%ymm10,%ymm7,%ymm7		# a0*e0+a1*e1*r0 a0*e1+a1*e0

vperm2i128	$0x21,%ymm6,%ymm4,%ymm8		# c1*g1*r1 c0*g1
vpblendd	$0xF0,%ymm6,%ymm12,%ymm12	# c0*g0 c1*g0
vpaddw		%ymm12,%ymm8,%ymm8		# c0*g0+c1*g1*r1 c0*g1+c1*g0

vperm2i128	$0x20,%ymm2,%ymm5,%ymm9		# -b1*f1*r0 b0*f1
vpblendd	$0xF0,%ymm2,%ymm11,%ymm11	# b0*f0 b1*f0
vpaddw		%ymm11,%ymm9,%ymm9		# b0*f0-b1*f1*r0 b0*f1+b1*f0

vperm2i128	$0x21,%ymm3,%ymm5,%ymm10	# -d1*h1*r1 d0*h1
vpblendd	$0xF0,%ymm3,%ymm13,%ymm13	# d0*h0 d1*h0
vpaddw		%ymm13,%ymm10,%ymm10		# d0*h0-d1*h1*r1 d0*h1+d1*h0

vmovdqa		(_16XF_PINV)*2(%rcx),%ymm2
vmovdqa		(_16XF)*2(%rcx),%ymm3
fqmulprecomp	2,3,7,x=15
fqmulprecomp	2,3,8,x=15
fqmulprecomp	2,3,9,x=15
fqmulprecomp	2,3,10,x=15

vmovdqa		%ymm7,(64*\off+ 0)*2(%rdi)
vmovdqa		%ymm9,(64*\off+16)*2(%rdi)
vmovdqa		%ymm8,(64*\off+32)*2(%rdi)
vmovdqa		%ymm10,(64*\off+48)*2(%rdi)
.endm

.text
.global cdecl(poly_basemul_montgomery)
cdecl(poly_basemul_montgomery):
vmovdqa		(_16XP)*2(%rcx),%ymm0
vmovdqa		(_16XPINV)*2(%rcx),%ymm1

schoolbook	0
schoolbook	1
schoolbook	2
schoolbook	3
schoolbook	4
schoolbook	5
schoolbook	6
schoolbook	7

ret
