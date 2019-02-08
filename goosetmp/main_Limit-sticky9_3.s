	.file	"mainbody.bc"


	.text
	.align	16
	.globl	main_Limit_2D_sticky9_3
	.type	main_Limit_2D_sticky9_3,@function
main_Limit_2D_sticky9_3:
.Leh_func_begin1:
.Llabel1:
	subq	$8, %rsp
	movsd	%xmm6, (%rsp)
	subsd	%xmm1, %xmm4
	mulsd	%xmm4, %xmm4
	subsd	%xmm0, %xmm3
	mulsd	%xmm3, %xmm3
	addsd	%xmm4, %xmm3
	subsd	%xmm2, %xmm5
	mulsd	%xmm5, %xmm5
	addsd	%xmm3, %xmm5
	movapd	%xmm5, %xmm0
	addsd	eps2_Limit_2D_sticky9_3, %xmm0
	call	rsqrt
	movapd	%xmm0, %xmm1
	mulsd	(%rsp), %xmm1
	pxor	%xmm0, %xmm0
	subsd	%xmm1, %xmm0
	call	accum
	movsd	%xmm0, pot_i__Limit_2D_sticky9_3
	pxor	%xmm0, %xmm0
	addq	$8, %rsp
	ret
	.size	main_Limit_2D_sticky9_3, .-main_Limit_2D_sticky9_3
.Leh_func_end1:
	.type	pot_i__Limit_2D_sticky9_3,@object
	.bss
	.globl pot_i__Limit_2D_sticky9_3
	.align	8
pot_i__Limit_2D_sticky9_3:				# pot_i__Limit-sticky9_3
	.size	pot_i__Limit_2D_sticky9_3, 8
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

.Lmain_Limit_2D_sticky9_3.eh:
	.long	.Leh_frame_end1-.Leh_frame_begin1
.Leh_frame_begin1:
	.long	.Leh_frame_begin1-.Leh_frame_common
	.long	.Leh_func_begin1-.
	.long	.Leh_func_end1-.Leh_func_begin1
	.uleb128	0
	.byte	0xE
	.uleb128	16
	.byte	0x4
	.long	.Llabel1-.Leh_func_begin1
	.byte	0xD
	.uleb128	7
	.align	8
.Leh_frame_end1:

	.section	.note.GNU-stack,"",@progbits
