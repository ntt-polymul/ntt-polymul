
/*
BSD 2-Clause License

Copyright (c) 2020, Bo-Yin Yang, Cheng-Jhih Shi, Chi-Ming Chung, Vincent Hwang
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this
   list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/
.macro montgomery_mul a, b, lower, upper, tmp, M_inv, M
    smull.w \lower, \upper, \a, \b
    mul.w \tmp, \lower, \M_inv
    smlal.w \lower, \upper, \tmp, \M
.endm

.macro add_sub a0, b0, a1, b1, a2, b2, a3, b3
    add \a0, \b0
    add \a1, \b1
    add \a2, \b2
    add \a3, \b3
    sub.w \b0, \a0, \b0, lsl #1        
    sub.w \b1, \a1, \b1, lsl #1        
    sub.w \b2, \a2, \b2, lsl #1        
    sub.w \b3, \a3, \b3, lsl #1        
.endm

.macro central_reduce target, Mhalf, M
    cmp \target, \Mhalf
    it gt
    subgt \target, \M
    cmn \target, \Mhalf
    it lt
    addlt \target, \M
.endm

.macro barret a, O_Mbar, O_M, tmp
    smmulr.w \tmp, \a, \O_Mbar
    mls.w \a, \tmp, \O_M, \a
.endm

.macro unsigned_char target, Mhalf, M
    cmn \target, \Mhalf
    it lt
    addlt \target, \M
.endm

.syntax unified
.cpu cortex-m4

.align 2
.global NTT_512
.type NTT_512, %function

NTT_512:
push.w {r4-r12, r14}
vpush.w {s16, s17}

tmp_root_table .req s8
outer_counter  .req s9
inner_counter  .req s10
ary_input1     .req s11
ary_input2     .req s12
tmp_counter    .req s13
tmp_inv_table  .req s14
tmp_ans        .req s15
q_bar          .req s16
tmp_save       .req s17

vldr.w ary_input2, [sp, #48]
vldr.w tmp_inv_table, [sp, #52]
vldr.w tmp_ans, [sp, #56]
vldr.w q_bar, [sp, #60]
vldm.w r1!, {s0-s6}
vmov.w tmp_root_table, r1
vmov.w ary_input1, r0

// allocate NTTd array using sp
sub.w sp, #4096
mov.w r0, sp
vmov.w tmp_save, r0

.align 2
_1_2_3:

add.w r12, r0, #256 // set counter
vmov.w inner_counter, r12

.align 2
loop_1_2_3:

    vmov.w r1, ary_input1
    ldrb.w r4, [r1, #0]
    ldrb.w r5, [r1, #64]
    ldrb.w r6, [r1, #128]
    ldrb.w r7, [r1, #192]
    ldrb.w r8, [r1, #256]
    ldrb.w r9, [r1, #320]
    ldrb.w r10, [r1, #384]
    ldrb.w r11, [r1, #448]
    add.w r1, #1
    vmov.w ary_input1, r1
    
    movw.w r12, 251
    movw.w r1, 125
    central_reduce r4, r1, r12
    central_reduce r5, r1, r12
    central_reduce r6, r1, r12
    central_reduce r7, r1, r12
    central_reduce r8, r1, r12
    central_reduce r9, r1, r12
    central_reduce r10, r1, r12
    central_reduce r11, r1, r12
    
    // level 1
    vmov.w r1, s0
    mul.w r8, r1
    mul.w r9, r1
    mul.w r10, r1
    mul.w r11, r1
    add_sub r4, r8, r5, r9, r6, r10, r7, r11
    
    // level 2
    vmov.w r1, s1
    montgomery_mul r6, r1, r12, r6, r14, r3, r2
    montgomery_mul r7, r1, r12, r7, r14, r3, r2
    vmov.w r1, s2
    montgomery_mul r10, r1, r12, r10, r14, r3, r2
    montgomery_mul r11, r1, r12, r11, r14, r3, r2
    
    add_sub r4, r6, r5, r7, r8, r10, r9, r11
    
    // level 3
    vmov.w r1, s3
    montgomery_mul r5, r1, r12, r5, r14, r3, r2 
    vmov.w r1, s4
    montgomery_mul r7, r1, r12, r7, r14, r3, r2 
    vmov.w r1, s5
    montgomery_mul r9, r1, r12, r9, r14, r3, r2
    vmov.w r1, s6
    montgomery_mul r11, r1, r12, r11, r14, r3, r2

    add_sub r4, r5, r6, r7, r8, r9, r10, r11
    
    str.w r4, [r0, #0]
    str.w r5, [r0, #256]
    str.w r6, [r0, #512]
    str.w r7, [r0, #768]
    str.w r8, [r0, #1024]
    str.w r9, [r0, #1280]
    str.w r10, [r0, #1536]
    str.w r11, [r0, #1792] 
    
    vmov.w r1, ary_input2
    ldrsb.w r4, [r1, #0]
    ldrsb.w r5, [r1, #64]
    ldrsb.w r6, [r1, #128]
    ldrsb.w r7, [r1, #192]
    ldrsb.w r8, [r1, #256]
    ldrsb.w r9, [r1, #320]
    ldrsb.w r10, [r1, #384]
    ldrsb.w r11, [r1, #448]
    add.w r1, #1
    vmov.w ary_input2, r1
    
    // level 1
    vmov.w r1, s0
    mul.w r8, r1
    mul.w r9, r1
    mul.w r10, r1
    mul.w r11, r1
    add_sub r4, r8, r5, r9, r6, r10, r7, r11
    
    // level 2
    vmov.w r1, s1
    montgomery_mul r6, r1, r12, r6, r14, r3, r2
    montgomery_mul r7, r1, r12, r7, r14, r3, r2
    vmov.w r1, s2
    montgomery_mul r10, r1, r12, r10, r14, r3, r2
    montgomery_mul r11, r1, r12, r11, r14, r3, r2
    
    add_sub r4, r6, r5, r7, r8, r10, r9, r11
    
    // level 3
    vmov.w r1, s3
    montgomery_mul r5, r1, r12, r5, r14, r3, r2 
    vmov.w r1, s4
    montgomery_mul r7, r1, r12, r7, r14, r3, r2 
    vmov.w r1, s5
    montgomery_mul r9, r1, r12, r9, r14, r3, r2
    vmov.w r1, s6
    montgomery_mul r11, r1, r12, r11, r14, r3, r2

    add_sub r4, r5, r6, r7, r8, r9, r10, r11
    
    str.w r4, [r0, #2048]
    str.w r5, [r0, #2304]
    str.w r6, [r0, #2560]
    str.w r7, [r0, #2816]
    str.w r8, [r0, #3072]
    str.w r9, [r0, #3328]
    str.w r10, [r0, #3584]
    str.w r11, [r0, #3840] 
    
    add.w r0, #4
    
    vmov.w r5, inner_counter
    cmp.w r5, r0
    bne.w loop_1_2_3      

.align 2
_4_5_6:

sub.w r0, #256 // reset pointer

add.w r1, r0, #2048   // outer iteration set counter
vmov.w outer_counter, r1


.align 2
normal_4_5_6_outer:
    vmov.w r1, tmp_root_table
    vldm.w r1!, {s0-s6}
    vmov.w tmp_root_table, r1

    add.w r4, r0, #32 // inner iteration set counter
    vmov.w inner_counter, r4
    normal_4_5_6_inner:

        ldr.w r4, [r0, #0]
        ldr.w r5, [r0, #32]
        ldr.w r6, [r0, #64]
        ldr.w r7, [r0, #96]
        ldr.w r8, [r0, #128]
        ldr.w r9, [r0, #160]
        ldr.w r10, [r0, #192]
        ldr.w r11, [r0, #224]
    
        // level 4
        vmov.w r1, s0
        montgomery_mul r8, r1, r12, r8, r14, r3, r2
        montgomery_mul r9, r1, r12, r9, r14, r3, r2
        montgomery_mul r10, r1, r12, r10, r14, r3, r2
        montgomery_mul r11, r1, r12, r11, r14, r3, r2
        
        add_sub r4, r8, r5, r9, r6, r10, r7, r11
        
        // level 5
        vmov.w r1, s1
        montgomery_mul r6, r1, r12, r6, r14, r3, r2   
        montgomery_mul r7, r1, r12, r7, r14, r3, r2
        vmov.w r1, s2
        montgomery_mul r10, r1, r12, r10, r14, r3, r2   
        montgomery_mul r11, r1, r12, r11, r14, r3, r2
        
        add_sub r4, r6, r5, r7, r8, r10, r9, r11
        
        // level 6
        vmov.w r1, s3
        montgomery_mul r5, r1, r12, r5, r14, r3, r2    
        vmov.w r1, s4
        montgomery_mul r7, r1, r12, r7, r14, r3, r2         
        vmov.w r1, s5
        montgomery_mul r9, r1, r12, r9, r14, r3, r2         
        vmov.w r1, s6
        montgomery_mul r11, r1, r12, r11, r14, r3, r2   
        
        add_sub r4, r5, r6, r7, r8, r9, r10, r11
    
        // save
        str.w r4, [r0, #0]
        str.w r5, [r0, #32]
        str.w r6, [r0, #64]
        str.w r7, [r0, #96]
        str.w r8, [r0, #128]
        str.w r9, [r0, #160]
        str.w r10, [r0, #192]
        str.w r11, [r0, #224]
    
        ldr.w r4, [r0, #2048]
        ldr.w r5, [r0, #2080]
        ldr.w r6, [r0, #2112]
        ldr.w r7, [r0, #2144]
        ldr.w r8, [r0, #2176]
        ldr.w r9, [r0, #2208]
        ldr.w r10, [r0, #2240]
        ldr.w r11, [r0, #2272]
    
        // level 4
        vmov.w r1, s0
        montgomery_mul r8, r1, r12, r8, r14, r3, r2
        montgomery_mul r9, r1, r12, r9, r14, r3, r2
        montgomery_mul r10, r1, r12, r10, r14, r3, r2
        montgomery_mul r11, r1, r12, r11, r14, r3, r2
        
        add_sub r4, r8, r5, r9, r6, r10, r7, r11
        
        // level 5
        vmov.w r1, s1
        montgomery_mul r6, r1, r12, r6, r14, r3, r2   
        montgomery_mul r7, r1, r12, r7, r14, r3, r2
        vmov.w r1, s2
        montgomery_mul r10, r1, r12, r10, r14, r3, r2   
        montgomery_mul r11, r1, r12, r11, r14, r3, r2
        
        add_sub r4, r6, r5, r7, r8, r10, r9, r11
        
        // level 6
        vmov.w r1, s3
        montgomery_mul r5, r1, r12, r5, r14, r3, r2    
        vmov.w r1, s4
        montgomery_mul r7, r1, r12, r7, r14, r3, r2         
        vmov.w r1, s5
        montgomery_mul r9, r1, r12, r9, r14, r3, r2         
        vmov.w r1, s6
        montgomery_mul r11, r1, r12, r11, r14, r3, r2   
        
        add_sub r4, r5, r6, r7, r8, r9, r10, r11
    
        // save
        str.w r4, [r0, #2048]
        str.w r5, [r0, #2080]
        str.w r6, [r0, #2112]
        str.w r7, [r0, #2144]
        str.w r8, [r0, #2176]
        str.w r9, [r0, #2208]
        str.w r10, [r0, #2240]
        str.w r11, [r0, #2272]
    
        add.w r0, #4
        vmov.w r4, inner_counter
        cmp.w r4, r0
        bne.w normal_4_5_6_inner

    add.w r0, #224
    vmov.w r4, outer_counter
    cmp.w r4, r0
    bne.w normal_4_5_6_outer

.align 2
_7_8_9:

sub.w r0, #2048 // reset pointer

add.w r12, r0, #2048   // set counter
vmov.w inner_counter, r12

.align 2
loop_7_8_9:
    vmov.w r1, tmp_root_table
    vldm.w r1!, {s0-s6}
    vmov.w tmp_root_table, r1

    ldr.w r4, [r0, #0]
    ldr.w r5, [r0, #4]
    ldr.w r6, [r0, #8]
    ldr.w r7, [r0, #12]
    ldr.w r8, [r0, #16]
    ldr.w r9, [r0, #20]
    ldr.w r10, [r0, #24]
    ldr.w r11, [r0, #28]
    
    // level 7
    vmov.w r1, s0
    montgomery_mul r8, r1, r12, r8, r14, r3, r2
    montgomery_mul r9, r1, r12, r9, r14, r3, r2
    montgomery_mul r10, r1, r12, r10, r14, r3, r2
    montgomery_mul r11, r1, r12, r11, r14, r3, r2
    
    add_sub r4, r8, r5, r9, r6, r10, r7, r11
    
    // level 8
    vmov.w r1, s1
    montgomery_mul r6, r1, r12, r6, r14, r3, r2   
    montgomery_mul r7, r1, r12, r7, r14, r3, r2
    vmov.w r1, s2
    montgomery_mul r10, r1, r12, r10, r14, r3, r2   
    montgomery_mul r11, r1, r12, r11, r14, r3, r2
    
    add_sub r4, r6, r5, r7, r8, r10, r9, r11
    
    // level 9
    vmov.w r1, s3
    montgomery_mul r5, r1, r12, r5, r14, r3, r2    
    vmov.w r1, s4
    montgomery_mul r7, r1, r12, r7, r14, r3, r2         
    vmov.w r1, s5
    montgomery_mul r9, r1, r12, r9, r14, r3, r2         
    vmov.w r1, s6
    montgomery_mul r11, r1, r12, r11, r14, r3, r2   
    
    add_sub r4, r5, r6, r7, r8, r9, r10, r11
    
    // save
    str.w r4, [r0, #0]
    str.w r5, [r0, #4]
    str.w r6, [r0, #8]
    str.w r7, [r0, #12]
    str.w r8, [r0, #16]
    str.w r9, [r0, #20]
    str.w r10, [r0, #24]
    str.w r11, [r0, #28]
    
    ldr.w r4, [r0, #2048]
    ldr.w r5, [r0, #2052]
    ldr.w r6, [r0, #2056]
    ldr.w r7, [r0, #2060]
    ldr.w r8, [r0, #2064]
    ldr.w r9, [r0, #2068]
    ldr.w r10, [r0, #2072]
    ldr.w r11, [r0, #2076]
    
    // level 7
    vmov.w r1, s0
    montgomery_mul r8, r1, r12, r8, r14, r3, r2
    montgomery_mul r9, r1, r12, r9, r14, r3, r2
    montgomery_mul r10, r1, r12, r10, r14, r3, r2
    montgomery_mul r11, r1, r12, r11, r14, r3, r2
    
    add_sub r4, r8, r5, r9, r6, r10, r7, r11
    
    // level 8
    vmov.w r1, s1
    montgomery_mul r6, r1, r12, r6, r14, r3, r2   
    montgomery_mul r7, r1, r12, r7, r14, r3, r2
    vmov.w r1, s2
    montgomery_mul r10, r1, r12, r10, r14, r3, r2   
    montgomery_mul r11, r1, r12, r11, r14, r3, r2
    
    add_sub r4, r6, r5, r7, r8, r10, r9, r11
    
    // level 9
    vmov.w r1, s3
    montgomery_mul r5, r1, r12, r5, r14, r3, r2    
    vmov.w r1, s4
    montgomery_mul r7, r1, r12, r7, r14, r3, r2         
    vmov.w r1, s5
    montgomery_mul r9, r1, r12, r9, r14, r3, r2         
    vmov.w r1, s6
    montgomery_mul r11, r1, r12, r11, r14, r3, r2   
    
    add_sub r4, r5, r6, r7, r8, r9, r10, r11
    
    // save
    str.w r4, [r0, #2048]
    str.w r5, [r0, #2052]
    str.w r6, [r0, #2056]
    str.w r7, [r0, #2060]
    str.w r8, [r0, #2064]
    str.w r9, [r0, #2068]
    str.w r10, [r0, #2072]
    str.w r11, [r0, #2076]
    
    add.w r0, #32
    vmov.w r4, inner_counter
    cmp.w r4, r0
    bne.w loop_7_8_9

sub.w r0, #2048

.align 2
my_mul:
    lower          .req r14
    tmp1           .req r12
    counter        .req r1

    M              .req r2
    M_inv          .req r3
    K0             .req r4
    K1             .req r5
    K2             .req r6
    K3             .req r7
    B0             .req r8
    B1             .req r9
    B2             .req r10
    B3             .req r11

    
    add.w counter, r0, #2048  // set counter
    
    my_multiply:

        ldr.w K0, [r0, #0]
        ldr.w K1, [r0, #4]
        ldr.w K2, [r0, #8]
        ldr.w K3, [r0, #12]
        ldr.w B0, [r0, #2048]
        ldr.w B1, [r0, #2052]
        ldr.w B2, [r0, #2056]
        ldr.w B3, [r0, #2060]
        
        smull.w lower, K0, K0, B0
        mul.w tmp1, lower, M_inv
        smlal.w lower, K0, tmp1, M 

        smull.w lower, K1, K1, B1
        mul.w tmp1, lower, M_inv
        smlal.w lower, K1, tmp1, M 

        smull.w lower, K2, K2, B2
        mul.w tmp1, lower, M_inv
        smlal.w lower, K2, tmp1, M

        smull.w lower, K3, K3, B3
        mul.w tmp1, lower, M_inv
        smlal.w lower, K3, tmp1, M  
        
        str.w K0, [r0, #0]
        str.w K1, [r0, #4]
        str.w K2, [r0, #8]
        str.w K3, [r0, #12]
    
        ldr.w K0, [r0, #16]
        ldr.w K1, [r0, #20]
        ldr.w K2, [r0, #24]
        ldr.w K3, [r0, #28]
        ldr.w B0, [r0, #2064]
        ldr.w B1, [r0, #2068]
        ldr.w B2, [r0, #2072]
        ldr.w B3, [r0, #2076]
        
        smull.w lower, K0, K0, B0
        mul.w tmp1, lower, M_inv
        smlal.w lower, K0, tmp1, M 

        smull.w lower, K1, K1, B1
        mul.w tmp1, lower, M_inv
        smlal.w lower, K1, tmp1, M 

        smull.w lower, K2, K2, B2
        mul.w tmp1, lower, M_inv
        smlal.w lower, K2, tmp1, M

        smull.w lower, K3, K3, B3
        mul.w tmp1, lower, M_inv
        smlal.w lower, K3, tmp1, M  
        
        str.w K0, [r0, #16]
        str.w K1, [r0, #20]
        str.w K2, [r0, #24]
        str.w K3, [r0, #28]
    
        ldr.w K0, [r0, #32]
        ldr.w K1, [r0, #36]
        ldr.w K2, [r0, #40]
        ldr.w K3, [r0, #44]
        ldr.w B0, [r0, #2080]
        ldr.w B1, [r0, #2084]
        ldr.w B2, [r0, #2088]
        ldr.w B3, [r0, #2092]
        
        smull.w lower, K0, K0, B0
        mul.w tmp1, lower, M_inv
        smlal.w lower, K0, tmp1, M 

        smull.w lower, K1, K1, B1
        mul.w tmp1, lower, M_inv
        smlal.w lower, K1, tmp1, M 

        smull.w lower, K2, K2, B2
        mul.w tmp1, lower, M_inv
        smlal.w lower, K2, tmp1, M

        smull.w lower, K3, K3, B3
        mul.w tmp1, lower, M_inv
        smlal.w lower, K3, tmp1, M  
        
        str.w K0, [r0, #32]
        str.w K1, [r0, #36]
        str.w K2, [r0, #40]
        str.w K3, [r0, #44]
    
        ldr.w K0, [r0, #48]
        ldr.w K1, [r0, #52]
        ldr.w K2, [r0, #56]
        ldr.w K3, [r0, #60]
        ldr.w B0, [r0, #2096]
        ldr.w B1, [r0, #2100]
        ldr.w B2, [r0, #2104]
        ldr.w B3, [r0, #2108]
        
        smull.w lower, K0, K0, B0
        mul.w tmp1, lower, M_inv
        smlal.w lower, K0, tmp1, M 

        smull.w lower, K1, K1, B1
        mul.w tmp1, lower, M_inv
        smlal.w lower, K1, tmp1, M 

        smull.w lower, K2, K2, B2
        mul.w tmp1, lower, M_inv
        smlal.w lower, K2, tmp1, M

        smull.w lower, K3, K3, B3
        mul.w tmp1, lower, M_inv
        smlal.w lower, K3, tmp1, M  
        
        str.w K0, [r0, #48]
        str.w K1, [r0, #52]
        str.w K2, [r0, #56]
        str.w K3, [r0, #60]
    
        ldr.w K0, [r0, #64]
        ldr.w K1, [r0, #68]
        ldr.w K2, [r0, #72]
        ldr.w K3, [r0, #76]
        ldr.w B0, [r0, #2112]
        ldr.w B1, [r0, #2116]
        ldr.w B2, [r0, #2120]
        ldr.w B3, [r0, #2124]
        
        smull.w lower, K0, K0, B0
        mul.w tmp1, lower, M_inv
        smlal.w lower, K0, tmp1, M 

        smull.w lower, K1, K1, B1
        mul.w tmp1, lower, M_inv
        smlal.w lower, K1, tmp1, M 

        smull.w lower, K2, K2, B2
        mul.w tmp1, lower, M_inv
        smlal.w lower, K2, tmp1, M

        smull.w lower, K3, K3, B3
        mul.w tmp1, lower, M_inv
        smlal.w lower, K3, tmp1, M  
        
        str.w K0, [r0, #64]
        str.w K1, [r0, #68]
        str.w K2, [r0, #72]
        str.w K3, [r0, #76]
    
        ldr.w K0, [r0, #80]
        ldr.w K1, [r0, #84]
        ldr.w K2, [r0, #88]
        ldr.w K3, [r0, #92]
        ldr.w B0, [r0, #2128]
        ldr.w B1, [r0, #2132]
        ldr.w B2, [r0, #2136]
        ldr.w B3, [r0, #2140]
        
        smull.w lower, K0, K0, B0
        mul.w tmp1, lower, M_inv
        smlal.w lower, K0, tmp1, M 

        smull.w lower, K1, K1, B1
        mul.w tmp1, lower, M_inv
        smlal.w lower, K1, tmp1, M 

        smull.w lower, K2, K2, B2
        mul.w tmp1, lower, M_inv
        smlal.w lower, K2, tmp1, M

        smull.w lower, K3, K3, B3
        mul.w tmp1, lower, M_inv
        smlal.w lower, K3, tmp1, M  
        
        str.w K0, [r0, #80]
        str.w K1, [r0, #84]
        str.w K2, [r0, #88]
        str.w K3, [r0, #92]
    
        ldr.w K0, [r0, #96]
        ldr.w K1, [r0, #100]
        ldr.w K2, [r0, #104]
        ldr.w K3, [r0, #108]
        ldr.w B0, [r0, #2144]
        ldr.w B1, [r0, #2148]
        ldr.w B2, [r0, #2152]
        ldr.w B3, [r0, #2156]
        
        smull.w lower, K0, K0, B0
        mul.w tmp1, lower, M_inv
        smlal.w lower, K0, tmp1, M 

        smull.w lower, K1, K1, B1
        mul.w tmp1, lower, M_inv
        smlal.w lower, K1, tmp1, M 

        smull.w lower, K2, K2, B2
        mul.w tmp1, lower, M_inv
        smlal.w lower, K2, tmp1, M

        smull.w lower, K3, K3, B3
        mul.w tmp1, lower, M_inv
        smlal.w lower, K3, tmp1, M  
        
        str.w K0, [r0, #96]
        str.w K1, [r0, #100]
        str.w K2, [r0, #104]
        str.w K3, [r0, #108]
    
        ldr.w K0, [r0, #112]
        ldr.w K1, [r0, #116]
        ldr.w K2, [r0, #120]
        ldr.w K3, [r0, #124]
        ldr.w B0, [r0, #2160]
        ldr.w B1, [r0, #2164]
        ldr.w B2, [r0, #2168]
        ldr.w B3, [r0, #2172]
        
        smull.w lower, K0, K0, B0
        mul.w tmp1, lower, M_inv
        smlal.w lower, K0, tmp1, M 

        smull.w lower, K1, K1, B1
        mul.w tmp1, lower, M_inv
        smlal.w lower, K1, tmp1, M 

        smull.w lower, K2, K2, B2
        mul.w tmp1, lower, M_inv
        smlal.w lower, K2, tmp1, M

        smull.w lower, K3, K3, B3
        mul.w tmp1, lower, M_inv
        smlal.w lower, K3, tmp1, M  
        
        str.w K0, [r0, #112]
        str.w K1, [r0, #116]
        str.w K2, [r0, #120]
        str.w K3, [r0, #124]
    
        add.w r0, #128
        cmp.w counter, r0
        bne.w my_multiply

sub.w r0, #2048 // reset pointer

_9_8_7:
add.w r12, r0, #2048   // set counter
vmov.w inner_counter, r12

.align 2
loop_9_8_7:

    vmov.w r1, tmp_inv_table
    vldm.w r1!, {s0-s6}
    vmov.w tmp_inv_table, r1

    ldr.w r4, [r0, #0]
    ldr.w r5, [r0, #4]
    ldr.w r6, [r0, #8]
    ldr.w r7, [r0, #12]
    ldr.w r8, [r0, #16]
    ldr.w r9, [r0, #20]
    ldr.w r10, [r0, #24]
    ldr.w r11, [r0, #28]
    
    // level 9
    add_sub r4, r5, r6, r7, r8, r9, r10, r11
    
    vmov.w r1, s3
    montgomery_mul r5, r1, r12, r5, r14, r3, r2    
    vmov.w r1, s4
    montgomery_mul r7, r1, r12, r7, r14, r3, r2         
    vmov.w r1, s5
    montgomery_mul r9, r1, r12, r9, r14, r3, r2         
    vmov.w r1, s6
    montgomery_mul r11, r1, r12, r11, r14, r3, r2   
    
    // level 8
    add_sub r4, r6, r5, r7, r8, r10, r9, r11
    
    vmov.w r1, s1
    montgomery_mul r6, r1, r12, r6, r14, r3, r2   
    montgomery_mul r7, r1, r12, r7, r14, r3, r2
    vmov.w r1, s2
    montgomery_mul r10, r1, r12, r10, r14, r3, r2   
    montgomery_mul r11, r1, r12, r11, r14, r3, r2
    
    // level 7
    add_sub r4, r8, r5, r9, r6, r10, r7, r11
    
    vmov.w r1, s0
    montgomery_mul r8, r1, r12, r8, r14, r3, r2
    montgomery_mul r9, r1, r12, r9, r14, r3, r2
    montgomery_mul r10, r1, r12, r10, r14, r3, r2
    montgomery_mul r11, r1, r12, r11, r14, r3, r2
    
    // save
    str.w r4, [r0, #0]
    str.w r5, [r0, #4]
    str.w r6, [r0, #8]
    str.w r7, [r0, #12]
    str.w r8, [r0, #16]
    str.w r9, [r0, #20]
    str.w r10, [r0, #24]
    str.w r11, [r0, #28]
    
    vmov.w r1, tmp_inv_table
    vldm.w r1!, {s0-s6}
    vmov.w tmp_inv_table, r1

    ldr.w r4, [r0, #32]
    ldr.w r5, [r0, #36]
    ldr.w r6, [r0, #40]
    ldr.w r7, [r0, #44]
    ldr.w r8, [r0, #48]
    ldr.w r9, [r0, #52]
    ldr.w r10, [r0, #56]
    ldr.w r11, [r0, #60]
    
    // level 9
    add_sub r4, r5, r6, r7, r8, r9, r10, r11
    
    vmov.w r1, s3
    montgomery_mul r5, r1, r12, r5, r14, r3, r2    
    vmov.w r1, s4
    montgomery_mul r7, r1, r12, r7, r14, r3, r2         
    vmov.w r1, s5
    montgomery_mul r9, r1, r12, r9, r14, r3, r2         
    vmov.w r1, s6
    montgomery_mul r11, r1, r12, r11, r14, r3, r2   
    
    // level 8
    add_sub r4, r6, r5, r7, r8, r10, r9, r11
    
    vmov.w r1, s1
    montgomery_mul r6, r1, r12, r6, r14, r3, r2   
    montgomery_mul r7, r1, r12, r7, r14, r3, r2
    vmov.w r1, s2
    montgomery_mul r10, r1, r12, r10, r14, r3, r2   
    montgomery_mul r11, r1, r12, r11, r14, r3, r2
    
    // level 7
    add_sub r4, r8, r5, r9, r6, r10, r7, r11
    
    vmov.w r1, s0
    montgomery_mul r8, r1, r12, r8, r14, r3, r2
    montgomery_mul r9, r1, r12, r9, r14, r3, r2
    montgomery_mul r10, r1, r12, r10, r14, r3, r2
    montgomery_mul r11, r1, r12, r11, r14, r3, r2
    
    // save
    str.w r4, [r0, #32]
    str.w r5, [r0, #36]
    str.w r6, [r0, #40]
    str.w r7, [r0, #44]
    str.w r8, [r0, #48]
    str.w r9, [r0, #52]
    str.w r10, [r0, #56]
    str.w r11, [r0, #60]
    
    add.w r0, #64
    vmov.w r4, inner_counter
    cmp.w r4, r0
    bne.w loop_9_8_7

sub.w r0, #2048 // reset pointer

.align 2
_6_5_4:

add.w r1, r0, #2048   // outer iteration set counter
vmov.w outer_counter, r1

.align 2
normal_6_5_4_outer:
    vmov.w r1, tmp_inv_table
    vldm.w r1!, {s0-s6}
    vmov.w tmp_inv_table, r1

    add.w r4, r0, #32 // inner iteration set counter
    vmov.w inner_counter, r4
 
    normal_6_5_4_inner:

        ldr.w r4, [r0, #0]
        ldr.w r5, [r0, #32]
        ldr.w r6, [r0, #64]
        ldr.w r7, [r0, #96]
        ldr.w r8, [r0, #128]
        ldr.w r9, [r0, #160]
        ldr.w r10, [r0, #192]
        ldr.w r11, [r0, #224]
    
        // level 6
        add_sub r4, r5, r6, r7, r8, r9, r10, r11
        
        vmov.w r1, s3
        montgomery_mul r5, r1, r12, r5, r14, r3, r2
        vmov.w r1, s4
        montgomery_mul r7, r1, r12, r7, r14, r3, r2
        vmov.w r1, s5
        montgomery_mul r9, r1, r12, r9, r14, r3, r2
        vmov.w r1, s6
        montgomery_mul r11, r1, r12, r11, r14, r3, r2
        
        // level 5
        add_sub r4, r6, r5, r7, r8, r10, r9, r11
        
        vmov.w r1, s1
        montgomery_mul r6, r1, r12, r6, r14, r3, r2
        montgomery_mul r7, r1, r12, r7, r14, r3, r2
        vmov.w r1, s2
        montgomery_mul r10, r1, r12, r10, r14, r3, r2
        montgomery_mul r11, r1, r12, r11, r14, r3, r2
        
        // level 4
        add_sub r4, r8, r5, r9, r6, r10, r7, r11
        
        vmov.w r1, s0
        montgomery_mul r8, r1, r12, r8, r14, r3, r2
        montgomery_mul r9, r1, r12, r9, r14, r3, r2
        montgomery_mul r10, r1, r12, r10, r14, r3, r2
        montgomery_mul r11, r1, r12, r11, r14, r3, r2
    
        // save
        str.w r4, [r0, #0]
        str.w r5, [r0, #32]
        str.w r6, [r0, #64]
        str.w r7, [r0, #96]
        str.w r8, [r0, #128]
        str.w r9, [r0, #160]
        str.w r10, [r0, #192]
        str.w r11, [r0, #224]
    
        ldr.w r4, [r0, #4]
        ldr.w r5, [r0, #36]
        ldr.w r6, [r0, #68]
        ldr.w r7, [r0, #100]
        ldr.w r8, [r0, #132]
        ldr.w r9, [r0, #164]
        ldr.w r10, [r0, #196]
        ldr.w r11, [r0, #228]
    
        // level 6
        add_sub r4, r5, r6, r7, r8, r9, r10, r11
        
        vmov.w r1, s3
        montgomery_mul r5, r1, r12, r5, r14, r3, r2
        vmov.w r1, s4
        montgomery_mul r7, r1, r12, r7, r14, r3, r2
        vmov.w r1, s5
        montgomery_mul r9, r1, r12, r9, r14, r3, r2
        vmov.w r1, s6
        montgomery_mul r11, r1, r12, r11, r14, r3, r2
        
        // level 5
        add_sub r4, r6, r5, r7, r8, r10, r9, r11
        
        vmov.w r1, s1
        montgomery_mul r6, r1, r12, r6, r14, r3, r2
        montgomery_mul r7, r1, r12, r7, r14, r3, r2
        vmov.w r1, s2
        montgomery_mul r10, r1, r12, r10, r14, r3, r2
        montgomery_mul r11, r1, r12, r11, r14, r3, r2
        
        // level 4
        add_sub r4, r8, r5, r9, r6, r10, r7, r11
        
        vmov.w r1, s0
        montgomery_mul r8, r1, r12, r8, r14, r3, r2
        montgomery_mul r9, r1, r12, r9, r14, r3, r2
        montgomery_mul r10, r1, r12, r10, r14, r3, r2
        montgomery_mul r11, r1, r12, r11, r14, r3, r2
    
        // save
        str.w r4, [r0, #4]
        str.w r5, [r0, #36]
        str.w r6, [r0, #68]
        str.w r7, [r0, #100]
        str.w r8, [r0, #132]
        str.w r9, [r0, #164]
        str.w r10, [r0, #196]
        str.w r11, [r0, #228]
    
        add.w r0, #8
        vmov.w r4, inner_counter
        cmp.w r4, r0
        bne.w normal_6_5_4_inner

    add.w r0, #224
    vmov.w r4, outer_counter
    cmp.w r4, r0
    bne.w normal_6_5_4_outer

sub.w r0, #2048 // reset pointer
vmov.w tmp_save, r0
add.w r12, r0, #256 // set counter
vmov.w inner_counter, r12

vmov.w r4, tmp_inv_table
vldm.w r4, {s0-s7}

.align 2
_3_2_1:
loop_3_2_1:

    ldr.w r4, [r0, #0]
    ldr.w r5, [r0, #256]
    ldr.w r6, [r0, #512]
    ldr.w r7, [r0, #768]
    ldr.w r8, [r0, #1024]
    ldr.w r9, [r0, #1280]
    ldr.w r10, [r0, #1536]
    ldr.w r11, [r0, #1792]
    add.w r0, #4 
    
    // level 3
    add_sub r4, r5, r6, r7, r8, r9, r10, r11
    
    vmov.w r1, s3
    montgomery_mul r5, r1, r12, r5, r14, r3, r2
    vmov.w r1, s4
    montgomery_mul r7, r1, r12, r7, r14, r3, r2
    vmov.w r1, s5
    montgomery_mul r9, r1, r12, r9, r14, r3, r2
    vmov.w r1, s6
    montgomery_mul r11, r1, r12, r11, r14, r3, r2 
    
    // level 2
    add_sub r4, r6, r5, r7, r8, r10, r9, r11
    
    vmov.w r1, s1
    montgomery_mul r6, r1, r12, r6, r14, r3, r2
    montgomery_mul r7, r1, r12, r7, r14, r3, r2
    vmov.w r1, s2
    montgomery_mul r10, r1, r12, r10, r14, r3, r2
    montgomery_mul r11, r1, r12, r11, r14, r3, r2
    
    // level 1
    add_sub r4, r8, r5, r9, r6, r10, r7, r11

    vmov.w r1, s7
    montgomery_mul r4, r1, r12, r4, r14, r3, r2
    montgomery_mul r5, r1, r12, r5, r14, r3, r2
    montgomery_mul r6, r1, r12, r6, r14, r3, r2
    montgomery_mul r7, r1, r12, r7, r14, r3, r2
    vmov.w r1, s0
    montgomery_mul r8, r1, r12, r8, r14, r3, r2
    montgomery_mul r9, r1, r12, r9, r14, r3, r2
    montgomery_mul r10, r1, r12, r10, r14, r3, r2
    montgomery_mul r11, r1, r12, r11, r14, r3, r2
    
    mov.w r1, r2, asr #1
    central_reduce r4, r1, r2
    central_reduce r5, r1, r2
    central_reduce r6, r1, r2
    central_reduce r7, r1, r2
    central_reduce r8, r1, r2
    central_reduce r9, r1, r2
    central_reduce r10, r1, r2
    central_reduce r11, r1, r2

    vmov.w r1, q_bar
    movw.w r14, 251
    barret r4, r1, r14, r12
    barret r5, r1, r14, r12
    barret r6, r1, r14, r12
    barret r7, r1, r14, r12
    barret r8, r1, r14, r12
    barret r9, r1, r14, r12
    barret r10, r1, r14, r12
    barret r11, r1, r14, r12

    mov.w r1, 0
    unsigned_char r4, r1, r14
    unsigned_char r5, r1, r14
    unsigned_char r6, r1, r14
    unsigned_char r7, r1, r14
    unsigned_char r8, r1, r14
    unsigned_char r9, r1, r14
    unsigned_char r10, r1, r14
    unsigned_char r11, r1, r14
    
    vmov.w r1, tmp_ans
    strb.w r4, [r1, #0]
    strb.w r5, [r1, #64]
    strb.w r6, [r1, #128]
    strb.w r7, [r1, #192]
    strb.w r8, [r1, #256]
    strb.w r9, [r1, #320]
    strb.w r10, [r1, #384]
    strb.w r11, [r1, #448]
    
    add.w r1, #1
    vmov.w tmp_ans, r1
    
    vmov.w r5, inner_counter
    cmp.w r5, r0
    bne.w loop_3_2_1               

add.w sp, #4096
vpop.w {s16, s17}
pop.w {r4-r12, pc}

