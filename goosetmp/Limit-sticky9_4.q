CONST x_0_0_, x_0_1_, x_0_2_, eps2, m_0_;
VARI x_i_0_, x_i_1_, x_i_2_, v_i_0_, v_i_1_, v_i_2_, m_i_;
VARJ x_j_0_, x_j_1_, x_j_2_, m_j_, v_j_0_, v_j_1_, v_j_2_;
VARF pot_pn2_i_;
pot_pn2_i__ = 0.0;

dxa_0_ = 0.0;  // the value 0.0 is meaningless.
dxb_0_ = 0.0;  // the value 0.0 is meaningless.
dxab_0_ = 0.0;  // the value 0.0 is meaningless.
dxa_1_ = 0.0;  // the value 0.0 is meaningless.
dxb_1_ = 0.0;  // the value 0.0 is meaningless.
dxab_1_ = 0.0;  // the value 0.0 is meaningless.
dxa_2_ = 0.0;  // the value 0.0 is meaningless.
dxb_2_ = 0.0;  // the value 0.0 is meaningless.
dxab_2_ = 0.0;  // the value 0.0 is meaningless.
r1a2e = 0.0;  // the value 0.0 is meaningless.
r1ainv = 0.0;  // the value 0.0 is meaningless.
r1b2e = 0.0;  // the value 0.0 is meaningless.
rab2e = 0.0;  // the value 0.0 is meaningless.
r1binv = 0.0;  // the value 0.0 is meaningless.
rabinv = 0.0;  // the value 0.0 is meaningless.
vi2 = 0.0;  // the value 0.0 is meaningless.
dxa_0_ = x_i_0_ - x_0_0_;
dxb_0_ = x_j_0_ - x_0_0_;
dxab_0_ = x_j_0_ - x_i_0_;
dxa_1_ = x_i_1_ - x_0_1_;
dxb_1_ = x_j_1_ - x_0_1_;
dxab_1_ = x_j_1_ - x_i_1_;
dxa_2_ = x_i_2_ - x_0_2_;
dxb_2_ = x_j_2_ - x_0_2_;
dxab_2_ = x_j_2_ - x_i_2_;
r1a2e = dxa_0_ * dxa_0_ + dxa_1_ * dxa_1_ + dxa_2_ * dxa_2_ + eps2;
r1ainv = rsqrt(r1a2e);
r1b2e = dxb_0_ * dxb_0_ + dxb_1_ * dxb_1_ + dxb_2_ * dxb_2_ + eps2;
rab2e = dxab_0_ * dxab_0_ + dxab_1_ * dxab_1_ + dxab_2_ * dxab_2_ + eps2;
r1binv = rsqrt(r1b2e);
rabinv = rsqrt(rab2e);
vi2 = v_i_0_ * v_i_0_ + v_i_1_ * v_i_1_ + v_i_2_ * v_i_2_;
pot_pn2_i__ = pot_pn2_i__ + (0.25 * rabinv * m_i_ * m_j_ * (6.0 * vi2 - 7.0 * (v_i_0_ * v_j_0_ + v_i_1_ * v_j_1_ + v_i_2_ * v_j_2_) - (dxab_0_ * v_i_0_ + dxab_1_ * v_i_1_ + dxab_2_ * v_i_2_) * (dxab_0_ * v_j_0_ + dxab_1_ * v_j_1_ + dxab_2_ * v_j_2_) / rab2e) + m_0_ * m_i_ * m_j_ * r1ainv * rabinv + 0.5 * m_0_ * m_i_ * m_j_ * r1ainv * r1binv);
pot_pn2_i_ += pot_pn2_i__;
