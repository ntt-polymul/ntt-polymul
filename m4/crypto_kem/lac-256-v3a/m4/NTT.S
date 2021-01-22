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
.global NTT_1024
.type NTT_1024, %function

NTT_1024:
push.w {r4-r12, r14}
vpush.w {s16, s17, s18}

tmp_root_table .req s8
outer_counter  .req s9
inner_counter  .req s10
ary_input1     .req s11
ary_input2     .req s12
tmp_counter    .req s13
tmp_inv_table  .req s14
tmp_ans        .req s15
q_bar          .req s16
tmp_mul_table  .req s17
tmp_lr         .req s18

vldr.w ary_input2, [sp, #52]
vldr.w tmp_inv_table, [sp, #56]
vldr.w tmp_ans, [sp, #60]
vldr.w q_bar, [sp, #64]
vldr.w tmp_mul_table, [sp, #68]
vldm.w r1!, {s0-s6}
vmov.w tmp_root_table, r1
vmov.w ary_input1, r0

// allocate NTTd array using sp
sub.w sp, #8192
mov.w r0, sp

bl.w _1_2_3_big
add.w r0, #3584 // adjust r0

bl.w _1_2_3_small
sub.w r0, #4608 // adjust r0

bl.w _4_5_6
vmov.w r1, tmp_root_table
sub.w r1, 224
vmov.w tmp_root_table, r1

bl.w _4_5_6
sub.w r0, #8192 // adjust r0

bl.w _7_8_9
vmov.w r1, tmp_root_table
sub.w r1, 1792
vmov.w tmp_root_table, r1

bl.w _7_8_9

sub.w r0, #8192
bl.w my_mul // do mul and inverse


.align 2
_1_2_3_big:
vmov.w tmp_lr, lr

add.w r12, r0, #512 // set counter
vmov.w inner_counter, r12

.align 2
loop_1_2_3_big:

    vmov.w r1, ary_input1
    ldrb.w r4, [r1, #0]
    ldrb.w r5, [r1, #128]
    ldrb.w r6, [r1, #256]
    ldrb.w r7, [r1, #384]
    ldrb.w r8, [r1, #512]
    ldrb.w r9, [r1, #640]
    ldrb.w r10, [r1, #768]
    ldrb.w r11, [r1, #896]
    
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
    str.w r5, [r0, #512]
    str.w r6, [r0, #1024]
    str.w r7, [r0, #1536]
    str.w r8, [r0, #2048]
    str.w r9, [r0, #2560]
    str.w r10, [r0, #3072]
    str.w r11, [r0, #3584] 
    
    vmov.w r1, ary_input1
    ldrb.w r4, [r1, #1]
    ldrb.w r5, [r1, #129]
    ldrb.w r6, [r1, #257]
    ldrb.w r7, [r1, #385]
    ldrb.w r8, [r1, #513]
    ldrb.w r9, [r1, #641]
    ldrb.w r10, [r1, #769]
    ldrb.w r11, [r1, #897]
    add.w r1, #2
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
    
    str.w r4, [r0, #4]
    str.w r5, [r0, #516]
    str.w r6, [r0, #1028]
    str.w r7, [r0, #1540]
    str.w r8, [r0, #2052]
    str.w r9, [r0, #2564]
    str.w r10, [r0, #3076]
    str.w r11, [r0, #3588] 
    
    add.w r0, #8
    
    vmov.w r5, inner_counter
    cmp.w r5, r0
    bne.w loop_1_2_3_big 

vmov.w lr, tmp_lr
bx lr     

.align 2
_1_2_3_small:
vmov.w tmp_lr, lr

add.w r12, r0, #512 // set counter
vmov.w inner_counter, r12

.align 2
loop_1_2_3_small:

    vmov.w r1, ary_input2
    ldrsb.w r4, [r1, #0]
    ldrsb.w r5, [r1, #128]
    ldrsb.w r6, [r1, #256]
    ldrsb.w r7, [r1, #384]
    ldrsb.w r8, [r1, #512]
    ldrsb.w r9, [r1, #640]
    ldrsb.w r10, [r1, #768]
    ldrsb.w r11, [r1, #896]
    
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
    str.w r5, [r0, #512]
    str.w r6, [r0, #1024]
    str.w r7, [r0, #1536]
    str.w r8, [r0, #2048]
    str.w r9, [r0, #2560]
    str.w r10, [r0, #3072]
    str.w r11, [r0, #3584] 
    
    vmov.w r1, ary_input2
    ldrsb.w r4, [r1, #1]
    ldrsb.w r5, [r1, #129]
    ldrsb.w r6, [r1, #257]
    ldrsb.w r7, [r1, #385]
    ldrsb.w r8, [r1, #513]
    ldrsb.w r9, [r1, #641]
    ldrsb.w r10, [r1, #769]
    ldrsb.w r11, [r1, #897]
    add.w r1, #2
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
    
    str.w r4, [r0, #4]
    str.w r5, [r0, #516]
    str.w r6, [r0, #1028]
    str.w r7, [r0, #1540]
    str.w r8, [r0, #2052]
    str.w r9, [r0, #2564]
    str.w r10, [r0, #3076]
    str.w r11, [r0, #3588] 
    
    add.w r0, #8
    
    vmov.w r5, inner_counter
    cmp.w r5, r0
    bne.w loop_1_2_3_small 

vmov.w lr, tmp_lr
bx lr     

.align 2
_4_5_6:
vmov.w tmp_lr, lr

add.w r1, r0, #4096   // outer iteration set counter
vmov.w outer_counter, r1


.align 2
normal_4_5_6_outer:
    vmov.w r1, tmp_root_table
    vldm.w r1!, {s0-s6}
    vmov.w tmp_root_table, r1

    add.w r4, r0, #64 // inner iteration set counter
    vmov.w inner_counter, r4
    normal_4_5_6_inner:

        ldr.w r4, [r0, #0]
        ldr.w r5, [r0, #64]
        ldr.w r6, [r0, #128]
        ldr.w r7, [r0, #192]
        ldr.w r8, [r0, #256]
        ldr.w r9, [r0, #320]
        ldr.w r10, [r0, #384]
        ldr.w r11, [r0, #448]
    
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
        str.w r5, [r0, #64]
        str.w r6, [r0, #128]
        str.w r7, [r0, #192]
        str.w r8, [r0, #256]
        str.w r9, [r0, #320]
        str.w r10, [r0, #384]
        str.w r11, [r0, #448]
    
        ldr.w r4, [r0, #4]
        ldr.w r5, [r0, #68]
        ldr.w r6, [r0, #132]
        ldr.w r7, [r0, #196]
        ldr.w r8, [r0, #260]
        ldr.w r9, [r0, #324]
        ldr.w r10, [r0, #388]
        ldr.w r11, [r0, #452]
    
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
        str.w r4, [r0, #4]
        str.w r5, [r0, #68]
        str.w r6, [r0, #132]
        str.w r7, [r0, #196]
        str.w r8, [r0, #260]
        str.w r9, [r0, #324]
        str.w r10, [r0, #388]
        str.w r11, [r0, #452]
    
        add.w r0, #8
        vmov.w r4, inner_counter
        cmp.w r4, r0
        bne.w normal_4_5_6_inner

    add.w r0, #448
    vmov.w r4, outer_counter
    cmp.w r4, r0
    bne.w normal_4_5_6_outer

vmov.w lr, tmp_lr
bx lr 

.align 2
_7_8_9:
vmov.w tmp_lr, lr

add.w r12, r0, #4096   // set counter
vmov.w inner_counter, r12

.align 2
loop_7_8_9:
    vmov.w r1, tmp_root_table
    vldm.w r1!, {s0-s6}
    vmov.w tmp_root_table, r1

    ldr r4, [r0, #0]
    ldr r5, [r0, #8]
    ldr r6, [r0, #16]
    ldr r7, [r0, #24]
    ldr.w r8, [r0, #32]
    ldr.w r9, [r0, #40]
    ldr.w r10, [r0, #48]
    ldr.w r11, [r0, #56]
    
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
    str r4, [r0, #0]
    str r5, [r0, #8]
    str r6, [r0, #16]
    str r7, [r0, #24]
    str.w r8, [r0, #32]
    str.w r9, [r0, #40]
    str.w r10, [r0, #48]
    str.w r11, [r0, #56]
    
    ldr r4, [r0, #4]
    ldr r5, [r0, #12]
    ldr r6, [r0, #20]
    ldr r7, [r0, #28]
    ldr.w r8, [r0, #36]
    ldr.w r9, [r0, #44]
    ldr.w r10, [r0, #52]
    ldr.w r11, [r0, #60]
    
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
    str r4, [r0, #4]
    str r5, [r0, #12]
    str r6, [r0, #20]
    str r7, [r0, #28]
    str.w r8, [r0, #36]
    str.w r9, [r0, #44]
    str.w r10, [r0, #52]
    str.w r11, [r0, #60]
    
    add.w r0, #64
    vmov.w r4, inner_counter
    cmp.w r4, r0
    bne.w loop_7_8_9

vmov.w lr, tmp_lr
bx lr

.align 2
my_mul:
    upper          .req r10
    lower          .req r14
    tmp1           .req r12
    root           .req r1

    M              .req r2
    M_inv          .req r3
    K0             .req r4
    K1             .req r5
    B0             .req r8
    B1             .req r9
    counter        .req r6
    mul_table      .req r7

    add.w r0, #128  // adjust r0 to access two polys
    add.w counter, r0, #4096  // set counter
    vmov.w mul_table, tmp_mul_table
    
    my_multiply:

        ldr.w root, [mul_table], #4
        ldr.w K0, [r0, #-128]
        ldr.w K1, [r0, #-124]
        ldr.w B0, [r0, #3968]
        ldr.w B1, [r0, #3972]
        
        smull.w lower, upper, K1, B0
        smlal.w lower, upper, K0, B1
        mul.w tmp1, lower, M_inv
        smlal.w lower, upper, tmp1, M  @ upper = K1B0 + K0B1

        smull.w lower, K1, K1, B1
        mul.w tmp1, lower, M_inv
        smlal.w lower, K1, tmp1, M
        smull.w lower, K1, K1, root
        smlal.w lower, K1, K0, B0
        mul.w tmp1, lower, M_inv
        smlal.w lower, K1, tmp1, M     @ K1 = K1B1*root + K0B0 
        
        str.w K1, [r0, #-128]
        str.w upper, [r0, #-124]
    
        ldr.w root, [mul_table], #4
        ldr.w K0, [r0, #-120]
        ldr.w K1, [r0, #-116]
        ldr.w B0, [r0, #3976]
        ldr.w B1, [r0, #3980]
        
        smull.w lower, upper, K1, B0
        smlal.w lower, upper, K0, B1
        mul.w tmp1, lower, M_inv
        smlal.w lower, upper, tmp1, M  @ upper = K1B0 + K0B1

        smull.w lower, K1, K1, B1
        mul.w tmp1, lower, M_inv
        smlal.w lower, K1, tmp1, M
        smull.w lower, K1, K1, root
        smlal.w lower, K1, K0, B0
        mul.w tmp1, lower, M_inv
        smlal.w lower, K1, tmp1, M     @ K1 = K1B1*root + K0B0 
        
        str.w K1, [r0, #-120]
        str.w upper, [r0, #-116]
    
        ldr.w root, [mul_table], #4
        ldr.w K0, [r0, #-112]
        ldr.w K1, [r0, #-108]
        ldr.w B0, [r0, #3984]
        ldr.w B1, [r0, #3988]
        
        smull.w lower, upper, K1, B0
        smlal.w lower, upper, K0, B1
        mul.w tmp1, lower, M_inv
        smlal.w lower, upper, tmp1, M  @ upper = K1B0 + K0B1

        smull.w lower, K1, K1, B1
        mul.w tmp1, lower, M_inv
        smlal.w lower, K1, tmp1, M
        smull.w lower, K1, K1, root
        smlal.w lower, K1, K0, B0
        mul.w tmp1, lower, M_inv
        smlal.w lower, K1, tmp1, M     @ K1 = K1B1*root + K0B0 
        
        str.w K1, [r0, #-112]
        str.w upper, [r0, #-108]
    
        ldr.w root, [mul_table], #4
        ldr.w K0, [r0, #-104]
        ldr.w K1, [r0, #-100]
        ldr.w B0, [r0, #3992]
        ldr.w B1, [r0, #3996]
        
        smull.w lower, upper, K1, B0
        smlal.w lower, upper, K0, B1
        mul.w tmp1, lower, M_inv
        smlal.w lower, upper, tmp1, M  @ upper = K1B0 + K0B1

        smull.w lower, K1, K1, B1
        mul.w tmp1, lower, M_inv
        smlal.w lower, K1, tmp1, M
        smull.w lower, K1, K1, root
        smlal.w lower, K1, K0, B0
        mul.w tmp1, lower, M_inv
        smlal.w lower, K1, tmp1, M     @ K1 = K1B1*root + K0B0 
        
        str.w K1, [r0, #-104]
        str.w upper, [r0, #-100]
    
        ldr.w root, [mul_table], #4
        ldr.w K0, [r0, #-96]
        ldr.w K1, [r0, #-92]
        ldr.w B0, [r0, #4000]
        ldr.w B1, [r0, #4004]
        
        smull.w lower, upper, K1, B0
        smlal.w lower, upper, K0, B1
        mul.w tmp1, lower, M_inv
        smlal.w lower, upper, tmp1, M  @ upper = K1B0 + K0B1

        smull.w lower, K1, K1, B1
        mul.w tmp1, lower, M_inv
        smlal.w lower, K1, tmp1, M
        smull.w lower, K1, K1, root
        smlal.w lower, K1, K0, B0
        mul.w tmp1, lower, M_inv
        smlal.w lower, K1, tmp1, M     @ K1 = K1B1*root + K0B0 
        
        str.w K1, [r0, #-96]
        str.w upper, [r0, #-92]
    
        ldr.w root, [mul_table], #4
        ldr.w K0, [r0, #-88]
        ldr.w K1, [r0, #-84]
        ldr.w B0, [r0, #4008]
        ldr.w B1, [r0, #4012]
        
        smull.w lower, upper, K1, B0
        smlal.w lower, upper, K0, B1
        mul.w tmp1, lower, M_inv
        smlal.w lower, upper, tmp1, M  @ upper = K1B0 + K0B1

        smull.w lower, K1, K1, B1
        mul.w tmp1, lower, M_inv
        smlal.w lower, K1, tmp1, M
        smull.w lower, K1, K1, root
        smlal.w lower, K1, K0, B0
        mul.w tmp1, lower, M_inv
        smlal.w lower, K1, tmp1, M     @ K1 = K1B1*root + K0B0 
        
        str.w K1, [r0, #-88]
        str.w upper, [r0, #-84]
    
        ldr.w root, [mul_table], #4
        ldr.w K0, [r0, #-80]
        ldr.w K1, [r0, #-76]
        ldr.w B0, [r0, #4016]
        ldr.w B1, [r0, #4020]
        
        smull.w lower, upper, K1, B0
        smlal.w lower, upper, K0, B1
        mul.w tmp1, lower, M_inv
        smlal.w lower, upper, tmp1, M  @ upper = K1B0 + K0B1

        smull.w lower, K1, K1, B1
        mul.w tmp1, lower, M_inv
        smlal.w lower, K1, tmp1, M
        smull.w lower, K1, K1, root
        smlal.w lower, K1, K0, B0
        mul.w tmp1, lower, M_inv
        smlal.w lower, K1, tmp1, M     @ K1 = K1B1*root + K0B0 
        
        str.w K1, [r0, #-80]
        str.w upper, [r0, #-76]
    
        ldr.w root, [mul_table], #4
        ldr.w K0, [r0, #-72]
        ldr.w K1, [r0, #-68]
        ldr.w B0, [r0, #4024]
        ldr.w B1, [r0, #4028]
        
        smull.w lower, upper, K1, B0
        smlal.w lower, upper, K0, B1
        mul.w tmp1, lower, M_inv
        smlal.w lower, upper, tmp1, M  @ upper = K1B0 + K0B1

        smull.w lower, K1, K1, B1
        mul.w tmp1, lower, M_inv
        smlal.w lower, K1, tmp1, M
        smull.w lower, K1, K1, root
        smlal.w lower, K1, K0, B0
        mul.w tmp1, lower, M_inv
        smlal.w lower, K1, tmp1, M     @ K1 = K1B1*root + K0B0 
        
        str.w K1, [r0, #-72]
        str.w upper, [r0, #-68]
    
        add.w r0, #64
        cmp.w counter, r0
        bne.w my_multiply

sub.w r0, #4224 // reset pointer

_9_8_7:
add.w r12, r0, #4096   // set counter
vmov.w inner_counter, r12

.align 2
loop_9_8_7:
    vmov.w r1, tmp_inv_table
    vldm.w r1!, {s0-s6}
    vmov.w tmp_inv_table, r1

    ldr r4, [r0, #0]
    ldr r5, [r0, #8]
    ldr r6, [r0, #16]
    ldr r7, [r0, #24]
    ldr.w r8, [r0, #32]
    ldr.w r9, [r0, #40]
    ldr.w r10, [r0, #48]
    ldr.w r11, [r0, #56]
    
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
    str r4, [r0, #0]
    str r5, [r0, #8]
    str r6, [r0, #16]
    str r7, [r0, #24]
    str.w r8, [r0, #32]
    str.w r9, [r0, #40]
    str.w r10, [r0, #48]
    str.w r11, [r0, #56]
    
    ldr r4, [r0, #4]
    ldr r5, [r0, #12]
    ldr r6, [r0, #20]
    ldr r7, [r0, #28]
    ldr.w r8, [r0, #36]
    ldr.w r9, [r0, #44]
    ldr.w r10, [r0, #52]
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
    str r4, [r0, #4]
    str r5, [r0, #12]
    str r6, [r0, #20]
    str r7, [r0, #28]
    str.w r8, [r0, #36]
    str.w r9, [r0, #44]
    str.w r10, [r0, #52]
    str.w r11, [r0, #60]
    
    add.w r0, #64
    vmov.w r4, inner_counter
    cmp.w r4, r0
    bne.w loop_9_8_7

sub.w r0, #4096 // reset pointer

.align 2
_6_5_4:

add.w r1, r0, #4096   // outer iteration set counter
vmov.w outer_counter, r1

.align 2
normal_6_5_4_outer:
    vmov.w r1, tmp_inv_table
    vldm.w r1!, {s0-s6}
    vmov.w tmp_inv_table, r1

    add.w r4, r0, #64 // inner iteration set counter
    vmov.w inner_counter, r4
 
    normal_6_5_4_inner:

        ldr.w r4, [r0, #0]
        ldr.w r5, [r0, #64]
        ldr.w r6, [r0, #128]
        ldr.w r7, [r0, #192]
        ldr.w r8, [r0, #256]
        ldr.w r9, [r0, #320]
        ldr.w r10, [r0, #384]
        ldr.w r11, [r0, #448]
    
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
        str.w r5, [r0, #64]
        str.w r6, [r0, #128]
        str.w r7, [r0, #192]
        str.w r8, [r0, #256]
        str.w r9, [r0, #320]
        str.w r10, [r0, #384]
        str.w r11, [r0, #448]
    
        ldr.w r4, [r0, #4]
        ldr.w r5, [r0, #68]
        ldr.w r6, [r0, #132]
        ldr.w r7, [r0, #196]
        ldr.w r8, [r0, #260]
        ldr.w r9, [r0, #324]
        ldr.w r10, [r0, #388]
        ldr.w r11, [r0, #452]
    
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
        str.w r5, [r0, #68]
        str.w r6, [r0, #132]
        str.w r7, [r0, #196]
        str.w r8, [r0, #260]
        str.w r9, [r0, #324]
        str.w r10, [r0, #388]
        str.w r11, [r0, #452]
    
        add.w r0, #8
        vmov.w r4, inner_counter
        cmp.w r4, r0
        bne.w normal_6_5_4_inner

    add.w r0, #448
    vmov.w r4, outer_counter
    cmp.w r4, r0
    bne.w normal_6_5_4_outer

sub.w r0, #4096 // reset pointer

add.w r12, r0, #512 // set counter
vmov.w inner_counter, r12

vmov.w r4, tmp_inv_table
vldm.w r4, {s0-s7}

.align 2
_3_2_1:
loop_3_2_1:

    ldr.w r4, [r0, #0]
    ldr.w r5, [r0, #512]
    ldr.w r6, [r0, #1024]
    ldr.w r7, [r0, #1536]
    ldr.w r8, [r0, #2048]
    ldr.w r9, [r0, #2560]
    ldr.w r10, [r0, #3072]
    ldr.w r11, [r0, #3584]
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
    strb.w r5, [r1, #128]
    strb.w r6, [r1, #256]
    strb.w r7, [r1, #384]
    strb.w r8, [r1, #512]
    strb.w r9, [r1, #640]
    strb.w r10, [r1, #768]
    strb.w r11, [r1, #896]
    
    add.w r1, #1
    vmov.w tmp_ans, r1
    
    vmov.w r5, inner_counter
    cmp.w r5, r0
    bne.w loop_3_2_1               

add.w sp, #8192
vpop.w {s16, s17, s18}
pop.w {r4-r12, pc}

