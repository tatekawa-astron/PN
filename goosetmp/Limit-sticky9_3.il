il_ps_2_0
dcl_input_interp(linear) v0.xy
;
;################################################################################
; generated on 2017-04-20 20:51:49
; generated from Limit-sticky9_3.q
; VARI  x_i_0_ x_i_1_ x_i_2_
; VARJ  x_j_0_ x_j_1_ x_j_2_ m_j_
; VARF  pot_i_
; CONST eps2
;################################################################################
;x_i_0_ -> res0
;x_i_1_ -> res1
;x_i_2_ -> res2
;x_j_0_ -> res3
;x_j_1_ -> res4
;x_j_2_ -> res5
;m_j_ -> res6
;pot_i_ -> o0
dcl_resource_id(0)_type(2d,unnorm)_fmtx(float)_fmty(float)_fmtz(float)_fmtw(float)
dcl_resource_id(1)_type(2d,unnorm)_fmtx(float)_fmty(float)_fmtz(float)_fmtw(float)
dcl_resource_id(2)_type(2d,unnorm)_fmtx(float)_fmty(float)_fmtz(float)_fmtw(float)
dcl_resource_id(3)_type(2d,unnorm)_fmtx(float)_fmty(float)_fmtz(float)_fmtw(float)
dcl_resource_id(4)_type(2d,unnorm)_fmtx(float)_fmty(float)_fmtz(float)_fmtw(float)
dcl_resource_id(5)_type(2d,unnorm)_fmtx(float)_fmty(float)_fmtz(float)_fmtw(float)
dcl_resource_id(6)_type(2d,unnorm)_fmtx(float)_fmty(float)_fmtz(float)_fmtw(float)
dcl_output_generic o0
dcl_cb cb0[1]
; const eps2 -> r500 cb0[0]

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
;load CB
mov r500, cb0[0]
; zero accum
mov r400, l10



mov r200.xy, l11.xy
mov r201.w, cb10[0].w
mov r201.z, cb10[0].z
ixor r201.x, r201.x, r201.x
itof r203.x, r201.z

whileloop
  sample_resource(3)_sampler(3) r300, r200.xy
  sample_resource(4)_sampler(4) r301, r200.xy
  sample_resource(5)_sampler(5) r302, r200.xy
  sample_resource(6)_sampler(6) r303, r200.xy

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
  dadd r710.xy, r300.xy, r100_neg(yw).xy
  dadd r710.zw, r300.zw, r100_neg(yw).zw
  dadd r711.xy, r301.xy, r101_neg(yw).xy
  dadd r711.zw, r301.zw, r101_neg(yw).zw
  dadd r712.xy, r302.xy, r102_neg(yw).xy
  dadd r712.zw, r302.zw, r102_neg(yw).zw
  dmul r705.xy, r710.xy, r710.xy
  dmul r705.zw, r710.zw, r710.zw
  dmul r706.xy, r711.xy, r711.xy
  dmul r706.zw, r711.zw, r711.zw
  dadd r701.xy, r705.xy, r706.xy
  dadd r701.zw, r705.zw, r706.zw
  dmul r707.xy, r712.xy, r712.xy
  dmul r707.zw, r712.zw, r712.zw
  dadd r702.xy, r701.xy, r707.xy
  dadd r702.zw, r701.zw, r707.zw
  dadd r703.xy, r702.xy, r500.xy
  dadd r703.zw, r702.zw, r500.zw
  mov r10.xy, r703.xy
  call 15
  mov r704.xy, r11.xy
  mov r10.xy, r703.zw
  call 15
  mov r704.__zw, r11.xyxy
  dmul r708.xy, r704.xy, r303.xy
  dmul r708.zw, r704.zw, r303.zw
  ixor r713, r708, l8
  dadd r400.xy, r400.xy, r713.xy
  dadd r400.zw, r400.zw, r713.zw

  ret
endfunc

;; permutation
func 1
  mov r300, r300.zwxy
  mov r301, r301.zwxy
  mov r302, r302.zwxy
  mov r303, r303.zwxy

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

