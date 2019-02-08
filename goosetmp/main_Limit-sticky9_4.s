	.file	"mainbody.bc"


	.section	.rodata.cst8,"aM",@progbits,8
	.align	8
.LCPI1_0:					
	.quad	13842939354630062080	# double value: -7.000000e+00
.LCPI1_1:					
	.quad	4618441417868443648	# double value: 6.000000e+00
.LCPI1_2:					
	.quad	4598175219545276416	# double value: 2.500000e-01
.LCPI1_3:					
	.quad	4602678819172646912	# double value: 5.000000e-01
	.text
	.align	16
	.globl	main_Limit_2D_sticky9_4
	.type	main_Limit_2D_sticky9_4,@function
main_Limit_2D_sticky9_4:
.Leh_func_begin1:
.Llabel1:
	subq	$120, %rsp
	movsd	%xmm7, 96(%rsp)
	movsd	%xmm6, 56(%rsp)
	movsd	%xmm5, 112(%rsp)
	movsd	%xmm4, 64(%rsp)
	movsd	%xmm3, 80(%rsp)
	movsd	%xmm2, 88(%rsp)
	movsd	%xmm1, 48(%rsp)
	movsd	%xmm0, 40(%rsp)
	movsd	x_0_1__Limit_2D_sticky9_4, %xmm0
	movsd	%xmm0, 104(%rsp)
	movapd	%xmm1, %xmm0
	subsd	104(%rsp), %xmm0
	mulsd	%xmm0, %xmm0
	movsd	x_0_0__Limit_2D_sticky9_4, %xmm1
	movsd	%xmm1, 32(%rsp)
	movsd	40(%rsp), %xmm1
	subsd	32(%rsp), %xmm1
	mulsd	%xmm1, %xmm1
	addsd	%xmm0, %xmm1
	movsd	x_0_2__Limit_2D_sticky9_4, %xmm0
	movsd	%xmm0, 24(%rsp)
	movapd	%xmm2, %xmm0
	subsd	24(%rsp), %xmm0
	mulsd	%xmm0, %xmm0
	addsd	%xmm1, %xmm0
	addsd	eps2_Limit_2D_sticky9_4, %xmm0
	call	rsqrt
	movsd	%xmm0, 16(%rsp)
	movsd	eps2_Limit_2D_sticky9_4, %xmm0
	movsd	%xmm0, 8(%rsp)
	movsd	128(%rsp), %xmm0
	subsd	104(%rsp), %xmm0
	mulsd	%xmm0, %xmm0
	movsd	96(%rsp), %xmm1
	subsd	32(%rsp), %xmm1
	mulsd	%xmm1, %xmm1
	addsd	%xmm0, %xmm1
	movsd	136(%rsp), %xmm0
	subsd	24(%rsp), %xmm0
	mulsd	%xmm0, %xmm0
	addsd	%xmm1, %xmm0
	addsd	8(%rsp), %xmm0
	call	rsqrt
	movsd	%xmm0, 32(%rsp)
	movsd	128(%rsp), %xmm0
	subsd	48(%rsp), %xmm0
	movsd	%xmm0, 72(%rsp)
	mulsd	%xmm0, %xmm0
	movsd	96(%rsp), %xmm1
	subsd	40(%rsp), %xmm1
	movsd	%xmm1, 96(%rsp)
	mulsd	%xmm1, %xmm1
	addsd	%xmm0, %xmm1
	movsd	136(%rsp), %xmm0
	subsd	88(%rsp), %xmm0
	movsd	%xmm0, 104(%rsp)
	mulsd	%xmm0, %xmm0
	addsd	%xmm1, %xmm0
	addsd	8(%rsp), %xmm0
	movsd	%xmm0, 88(%rsp)
	call	rsqrt
	movapd	%xmm0, %xmm1
	movsd	m_0__Limit_2D_sticky9_4, %xmm0
	movsd	160(%rsp), %xmm2
	movsd	72(%rsp), %xmm3
	mulsd	%xmm2, %xmm3
	movsd	152(%rsp), %xmm4
	movsd	96(%rsp), %xmm5
	mulsd	%xmm4, %xmm5
	addsd	%xmm3, %xmm5
	movsd	168(%rsp), %xmm3
	movsd	104(%rsp), %xmm6
	mulsd	%xmm3, %xmm6
	addsd	%xmm5, %xmm6
	movsd	72(%rsp), %xmm5
	mulsd	64(%rsp), %xmm5
	movsd	%xmm5, 72(%rsp)
	movsd	96(%rsp), %xmm5
	mulsd	80(%rsp), %xmm5
	addsd	72(%rsp), %xmm5
	movsd	%xmm5, 96(%rsp)
	movsd	104(%rsp), %xmm5
	mulsd	112(%rsp), %xmm5
	addsd	96(%rsp), %xmm5
	mulsd	%xmm6, %xmm5
	divsd	88(%rsp), %xmm5
	movsd	%xmm5, 104(%rsp)
	mulsd	64(%rsp), %xmm2
	mulsd	80(%rsp), %xmm4
	addsd	%xmm2, %xmm4
	mulsd	112(%rsp), %xmm3
	addsd	%xmm4, %xmm3
	mulsd	.LCPI1_0(%rip), %xmm3
	movsd	64(%rsp), %xmm2
	mulsd	%xmm2, %xmm2
	movsd	%xmm2, 64(%rsp)
	movsd	80(%rsp), %xmm2
	mulsd	%xmm2, %xmm2
	addsd	64(%rsp), %xmm2
	movsd	%xmm2, 80(%rsp)
	movsd	112(%rsp), %xmm2
	mulsd	%xmm2, %xmm2
	addsd	80(%rsp), %xmm2
	mulsd	.LCPI1_1(%rip), %xmm2
	addsd	%xmm3, %xmm2
	subsd	%xmm5, %xmm2
	movsd	%xmm2, 112(%rsp)
	movapd	%xmm1, %xmm2
	mulsd	.LCPI1_2(%rip), %xmm2
	mulsd	56(%rsp), %xmm2
	movsd	144(%rsp), %xmm3
	mulsd	%xmm3, %xmm2
	mulsd	112(%rsp), %xmm2
	movapd	%xmm0, %xmm4
	mulsd	56(%rsp), %xmm4
	mulsd	%xmm3, %xmm4
	mulsd	16(%rsp), %xmm4
	mulsd	%xmm1, %xmm4
	addsd	%xmm2, %xmm4
	mulsd	.LCPI1_3(%rip), %xmm0
	mulsd	56(%rsp), %xmm0
	mulsd	%xmm3, %xmm0
	mulsd	16(%rsp), %xmm0
	mulsd	32(%rsp), %xmm0
	addsd	%xmm4, %xmm0
	pxor	%xmm1, %xmm1
	addsd	%xmm1, %xmm0
	call	accum
	movsd	%xmm0, pot_pn2_i__Limit_2D_sticky9_4
	pxor	%xmm0, %xmm0
	addq	$120, %rsp
	ret
	.size	main_Limit_2D_sticky9_4, .-main_Limit_2D_sticky9_4
.Leh_func_end1:
	.type	pot_pn2_i__Limit_2D_sticky9_4,@object
	.bss
	.globl pot_pn2_i__Limit_2D_sticky9_4
	.align	8
pot_pn2_i__Limit_2D_sticky9_4:				# pot_pn2_i__Limit-sticky9_4
	.size	pot_pn2_i__Limit_2D_sticky9_4, 8
	.zero	8
	.section	.eh_frame,"aw",@progbits
.LEH_frame0:
.Lsection_eh_frame:
.Leh_frame_common:
	.long	.Leh_frame_common_end-.Leh_frame_common_begin
.Leh_frame_common_begin:
	.long	0x0
	.byte	0x1
	.asciz	"zR"
	.uleb128	1
	.sleb128	-8
	.byte	0x10
	.uleb128	1
	.byte	0x1B
	.byte	0xC
	.uleb128	7
	.uleb128	8
	.byte	0x90
	.uleb128	1
	.align	8
.Leh_frame_common_end:

.Lmain_Limit_2D_sticky9_4.eh:
	.long	.Leh_frame_end1-.Leh_frame_begin1
.Leh_frame_begin1:
	.long	.Leh_frame_begin1-.Leh_frame_common
	.long	.Leh_func_begin1-.
	.long	.Leh_func_end1-.Leh_func_begin1
	.uleb128	0
	.byte	0xE
	.uleb128	128
	.byte	0x4
	.long	.Llabel1-.Leh_func_begin1
	.byte	0xD
	.uleb128	7
	.align	8
.Leh_frame_end1:

	.section	.note.GNU-stack,"",@progbits
