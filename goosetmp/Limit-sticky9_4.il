il_ps_2_0
dcl_input_interp(linear) v0.xy
;
;################################################################################
; generated on 2017-04-20 20:51:49
; generated from Limit-sticky9_4.q
; VARI  x_i_0_ x_i_1_ x_i_2_ v_i_0_ v_i_1_ v_i_2_ m_i_
; VARJ  x_j_0_ x_j_1_ x_j_2_ m_j_ v_j_0_ v_j_1_ v_j_2_
; VARF  pot_pn2_i_
; CONST x_0_0_ x_0_1_ x_0_2_ eps2 m_0_
;################################################################################
;x_i_0_ -> res0
;x_i_1_ -> res1
;x_i_2_ -> res2
;v_i_0_ -> res3
;v_i_1_ -> res4
;v_i_2_ -> res5
;m_i_ -> res6
;x_j_0_ -> res7
;x_j_1_ -> res8
;x_j_2_ -> res9
;m_j_ -> res10
;v_j_0_ -> res11
;v_j_1_ -> res12
;v_j_2_ -> res13
;pot_pn2_i_ -> o0
dcl_resource_id(0)_type(2d,unnorm)_fmtx(float)_fmty(float)_fmtz(float)_fmtw(float)
dcl_resource_id(1)_type(2d,unnorm)_fmtx(float)_fmty(float)_fmtz(float)_fmtw(float)
dcl_resource_id(2)_type(2d,unnorm)_fmtx(float)_fmty(float)_fmtz(float)_fmtw(float)
dcl_resource_id(3)_type(2d,unnorm)_fmtx(float)_fmty(float)_fmtz(float)_fmtw(float)
dcl_resource_id(4)_type(2d,unnorm)_fmtx(float)_fmty(float)_fmtz(float)_fmtw(float)
dcl_resource_id(5)_type(2d,unnorm)_fmtx(float)_fmty(float)_fmtz(float)_fmtw(float)
dcl_resource_id(6)_type(2d,unnorm)_fmtx(float)_fmty(float)_fmtz(float)_fmtw(float)
dcl_resource_id(7)_type(2d,unnorm)_fmtx(float)_fmty(float)_fmtz(float)_fmtw(float)
dcl_resource_id(8)_type(2d,unnorm)_fmtx(float)_fmty(float)_fmtz(float)_fmtw(float)
dcl_resource_id(9)_type(2d,unnorm)_fmtx(float)_fmty(float)_fmtz(float)_fmtw(float)
dcl_resource_id(10)_type(2d,unnorm)_fmtx(float)_fmty(float)_fmtz(float)_fmtw(float)
dcl_resource_id(11)_type(2d,unnorm)_fmtx(float)_fmty(float)_fmtz(float)_fmtw(float)
dcl_resource_id(12)_type(2d,unnorm)_fmtx(float)_fmty(float)_fmtz(float)_fmtw(float)
dcl_resource_id(13)_type(2d,unnorm)_fmtx(float)_fmty(float)_fmtz(float)_fmtw(float)
dcl_output_generic o0
; const literal 0.000000e+00
dcl_literal l20, 0x0, 0x0, 0x0, 0x0
; const literal 2.500000e-01
dcl_literal l21, 0x0, 0x3fd00000, 0x0, 0x3fd00000
; const literal 5.000000e-01
dcl_literal l22, 0x0, 0x3fe00000, 0x0, 0x3fe00000
; const literal 6.000000e+00
dcl_literal l23, 0x0, 0x40180000, 0x0, 0x40180000
; const literal 7.000000e+00
dcl_literal l24, 0x0, 0x401c0000, 0x0, 0x401c0000
dcl_cb cb0[5]
; const x_0_0_ -> r505 cb0[0]
; const x_0_1_ -> r506 cb0[1]
; const x_0_2_ -> r507 cb0[2]
; const eps2 -> r508 cb0[3]
; const m_0_ -> r509 cb0[4]

; nj
dcl_cb cb10[1]

; float zero
dcl_literal l10, 0.0, 0.0, 0.0, 0.0
; for loop counter
dcl_literal l11, 0.0, 0.0, 1.0, 0x1
; for rsq
dcl_literal l4, 0x0, 0x3fe00000, 0x0, 0x3ff80000
; for negation
dcl_literal l8, 0x0, 0x80000000, 0x0, 0x80000000
;;
dcl_literal l50, 0x0, 0x0, 0x0, 0x0
;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;load VARI
sample_resource(0)_sampler(0) r100, v0.xy
sample_resource(1)_sampler(1) r101, v0.xy
sample_resource(2)_sampler(2) r102, v0.xy
sample_resource(3)_sampler(3) r103, v0.xy
sample_resource(4)_sampler(4) r104, v0.xy
sample_resource(5)_sampler(5) r105, v0.xy
sample_resource(6)_sampler(6) r106, v0.xy
;load CB
mov r505, cb0[0]
mov r506, cb0[1]
mov r507, cb0[2]
mov r508, cb0[3]
mov r509, cb0[4]
mov r500, l20
mov r501, l21
mov r502, l22
mov r503, l23
mov r504, l24
; zero accum
mov r400, l10



mov r200.xy, l11.xy
mov r201.w, cb10[0].w
mov r201.z, cb10[0].z
ixor r201.x, r201.x, r201.x
itof r203.x, r201.z

whileloop
  sample_resource(7)_sampler(7) r300, r200.xy
  sample_resource(8)_sampler(8) r301, r200.xy
  sample_resource(9)_sampler(9) r302, r200.xy
  sample_resource(10)_sampler(10) r303, r200.xy
  sample_resource(11)_sampler(11) r304, r200.xy
  sample_resource(12)_sampler(12) r305, r200.xy
  sample_resource(13)_sampler(13) r306, r200.xy

  call 0
  call 1
  call 0


  add r200.x, r200.x, l11.z
  iadd r201.x, r201.x, l11.w
  ige r201.y, r201.x, r201.w
  break_logicalnz r201.y

;  umod r202.x, r201.x, r201.z
;  if_logicalz r202.x
;    add r200.0y, r200.0y, l11.z
;  endif

  eq r202.x, r200.x, r203.x
  if_logicalnz r202.x
    add r200.0y, r200.0y, l11.z
  endif
endloop

mov o0, r400


ret_dyn
endmain
; end of the main function

;; force
func 0
  dadd r768.xy, r100.xy, r505_neg(yw).xy
  dadd r768.zw, r100.zw, r505_neg(yw).zw
  dadd r772.xy, r300.xy, r505_neg(yw).xy
  dadd r772.zw, r300.zw, r505_neg(yw).zw
  dadd r773.xy, r300.xy, r100_neg(yw).xy
  dadd r773.zw, r300.zw, r100_neg(yw).zw
  dadd r774.xy, r101.xy, r506_neg(yw).xy
  dadd r774.zw, r101.zw, r506_neg(yw).zw
  dadd r777.xy, r301.xy, r506_neg(yw).xy
  dadd r777.zw, r301.zw, r506_neg(yw).zw
  dadd r778.xy, r301.xy, r101_neg(yw).xy
  dadd r778.zw, r301.zw, r101_neg(yw).zw
  dadd r769.xy, r102.xy, r507_neg(yw).xy
  dadd r769.zw, r102.zw, r507_neg(yw).zw
  dadd r770.xy, r302.xy, r507_neg(yw).xy
  dadd r770.zw, r302.zw, r507_neg(yw).zw
  dadd r771.xy, r302.xy, r102_neg(yw).xy
  dadd r771.zw, r302.zw, r102_neg(yw).zw
  dmul r730.xy, r768.xy, r768.xy
  dmul r730.zw, r768.zw, r768.zw
  dmul r731.xy, r774.xy, r774.xy
  dmul r731.zw, r774.zw, r774.zw
  dadd r701.xy, r730.xy, r731.xy
  dadd r701.zw, r730.zw, r731.zw
  dmul r732.xy, r769.xy, r769.xy
  dmul r732.zw, r769.zw, r769.zw
  dadd r702.xy, r701.xy, r732.xy
  dadd r702.zw, r701.zw, r732.zw
  dadd r703.xy, r702.xy, r508.xy
  dadd r703.zw, r702.zw, r508.zw
  mov r10.xy, r703.xy
  call 15
  mov r721.xy, r11.xy
  mov r10.xy, r703.zw
  call 15
  mov r721.__zw, r11.xyxy
  dmul r733.xy, r772.xy, r772.xy
  dmul r733.zw, r772.zw, r772.zw
  dmul r734.xy, r777.xy, r777.xy
  dmul r734.zw, r777.zw, r777.zw
  dadd r704.xy, r733.xy, r734.xy
  dadd r704.zw, r733.zw, r734.zw
  dmul r735.xy, r770.xy, r770.xy
  dmul r735.zw, r770.zw, r770.zw
  dadd r705.xy, r704.xy, r735.xy
  dadd r705.zw, r704.zw, r735.zw
  dadd r706.xy, r705.xy, r508.xy
  dadd r706.zw, r705.zw, r508.zw
  dmul r736.xy, r773.xy, r773.xy
  dmul r736.zw, r773.zw, r773.zw
  dmul r737.xy, r778.xy, r778.xy
  dmul r737.zw, r778.zw, r778.zw
  dadd r707.xy, r736.xy, r737.xy
  dadd r707.zw, r736.zw, r737.zw
  dmul r738.xy, r771.xy, r771.xy
  dmul r738.zw, r771.zw, r771.zw
  dadd r708.xy, r707.xy, r738.xy
  dadd r708.zw, r707.zw, r738.zw
  dadd r709.xy, r708.xy, r508.xy
  dadd r709.zw, r708.zw, r508.zw
  mov r10.xy, r706.xy
  call 15
  mov r722.xy, r11.xy
  mov r10.xy, r706.zw
  call 15
  mov r722.__zw, r11.xyxy
  mov r10.xy, r709.xy
  call 15
  mov r723.xy, r11.xy
  mov r10.xy, r709.zw
  call 15
  mov r723.__zw, r11.xyxy
  dmul r739.xy, r103.xy, r103.xy
  dmul r739.zw, r103.zw, r103.zw
  dmul r740.xy, r104.xy, r104.xy
  dmul r740.zw, r104.zw, r104.zw
  dadd r710.xy, r739.xy, r740.xy
  dadd r710.zw, r739.zw, r740.zw
  dmul r741.xy, r105.xy, r105.xy
  dmul r741.zw, r105.zw, r105.zw
  dadd r711.xy, r710.xy, r741.xy
  dadd r711.zw, r710.zw, r741.zw
  dmul r742.xy, r723.xy, r501.xy
  dmul r742.zw, r723.zw, r501.zw
  dmul r743.xy, r742.xy, r106.xy
  dmul r743.zw, r742.zw, r106.zw
  dmul r744.xy, r743.xy, r303.xy
  dmul r744.zw, r743.zw, r303.zw
  dmul r745.xy, r711.xy, r503.xy
  dmul r745.zw, r711.zw, r503.zw
  dmul r746.xy, r103.xy, r304.xy
  dmul r746.zw, r103.zw, r304.zw
  dmul r747.xy, r104.xy, r305.xy
  dmul r747.zw, r104.zw, r305.zw
  dadd r712.xy, r746.xy, r747.xy
  dadd r712.zw, r746.zw, r747.zw
  dmul r748.xy, r105.xy, r306.xy
  dmul r748.zw, r105.zw, r306.zw
  dadd r713.xy, r712.xy, r748.xy
  dadd r713.zw, r712.zw, r748.zw
  dmul r749.xy, r713.xy, r504.xy
  dmul r749.zw, r713.zw, r504.zw
  dadd r775.xy, r745.xy, r749_neg(yw).xy
  dadd r775.zw, r745.zw, r749_neg(yw).zw
  dmul r750.xy, r773.xy, r103.xy
  dmul r750.zw, r773.zw, r103.zw
  dmul r751.xy, r778.xy, r104.xy
  dmul r751.zw, r778.zw, r104.zw
  dadd r714.xy, r750.xy, r751.xy
  dadd r714.zw, r750.zw, r751.zw
  dmul r752.xy, r771.xy, r105.xy
  dmul r752.zw, r771.zw, r105.zw
  dadd r715.xy, r714.xy, r752.xy
  dadd r715.zw, r714.zw, r752.zw
  dmul r753.xy, r773.xy, r304.xy
  dmul r753.zw, r773.zw, r304.zw
  dmul r754.xy, r778.xy, r305.xy
  dmul r754.zw, r778.zw, r305.zw
  dadd r716.xy, r753.xy, r754.xy
  dadd r716.zw, r753.zw, r754.zw
  dmul r755.xy, r771.xy, r306.xy
  dmul r755.zw, r771.zw, r306.zw
  dadd r717.xy, r716.xy, r755.xy
  dadd r717.zw, r716.zw, r755.zw
  dmul r756.xy, r715.xy, r717.xy
  dmul r756.zw, r715.zw, r717.zw
  ddiv r729.xy, r756.xy, r709.xy
  ddiv r729.zw, r756.zw, r709.zw
  dadd r776.xy, r775.xy, r729_neg(yw).xy
  dadd r776.zw, r775.zw, r729_neg(yw).zw
  dmul r757.xy, r744.xy, r776.xy
  dmul r757.zw, r744.zw, r776.zw
  dmul r758.xy, r509.xy, r106.xy
  dmul r758.zw, r509.zw, r106.zw
  dmul r759.xy, r758.xy, r303.xy
  dmul r759.zw, r758.zw, r303.zw
  dmul r760.xy, r759.xy, r721.xy
  dmul r760.zw, r759.zw, r721.zw
  dmul r761.xy, r760.xy, r723.xy
  dmul r761.zw, r760.zw, r723.zw
  dadd r718.xy, r757.xy, r761.xy
  dadd r718.zw, r757.zw, r761.zw
  dmul r762.xy, r509.xy, r502.xy
  dmul r762.zw, r509.zw, r502.zw
  dmul r763.xy, r762.xy, r106.xy
  dmul r763.zw, r762.zw, r106.zw
  dmul r764.xy, r763.xy, r303.xy
  dmul r764.zw, r763.zw, r303.zw
  dmul r765.xy, r764.xy, r721.xy
  dmul r765.zw, r764.zw, r721.zw
  dmul r766.xy, r765.xy, r722.xy
  dmul r766.zw, r765.zw, r722.zw
  dadd r719.xy, r718.xy, r766.xy
  dadd r719.zw, r718.zw, r766.zw
  dadd r720.xy, r719.xy, r500.xy
  dadd r720.zw, r719.zw, r500.zw
  dadd r400.xy, r400.xy, r720.xy
  dadd r400.zw, r400.zw, r720.zw

  ret
endfunc

;; permutation
func 1
  mov r300, r300.zwxy
  mov r301, r301.zwxy
  mov r302, r302.zwxy
  mov r303, r303.zwxy
  mov r304, r304.zwxy
  mov r305, r305.zwxy
  mov r306, r306.zwxy

  ret
endfunc

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; rsq double
; in r10.xy
; out r11.xy
func 15
  ; x0 = rsq(x)
  d2f r11, r10.xy
  rsq r11, r11
  f2d r12, r11

  ; x*0.5
  dmul r13.xy, r10.xy, l4.xy

  ; b = 1.5*x0
  dmul r14.xy, r12.xy, l4.zw

  ; a = x0*x0*x0*x*0.5
  dmul r16.xy, r12.xy, r12.xy
  dmul r16.xy, r16.xy, r12.xy
  dmul r16.xy, r16.xy, r13.xy

  ; x1 = b - a
  ixor r999, r16, l8
  dadd r11.xy, r14.xy, r999.xy
  ret
endfunc

end

