CONST eps2;
VARI x_i_0_, x_i_1_, x_i_2_;
VARJ x_j_0_, x_j_1_, x_j_2_, m_j_;
VARF a_i_0_, a_i_1_, a_i_2_;
a_i_2__ = 0.0;
a_i_1__ = 0.0;
a_i_0__ = 0.0;

dx_0_ = 0.0;  // the value 0.0 is meaningless.
dx_1_ = 0.0;  // the value 0.0 is meaningless.
dx_2_ = 0.0;  // the value 0.0 is meaningless.
r2 = 0.0;  // the value 0.0 is meaningless.
rinv = 0.0;  // the value 0.0 is meaningless.
mrinv = 0.0;  // the value 0.0 is meaningless.
mr3inv = 0.0;  // the value 0.0 is meaningless.
dx_0_ = x_j_0_ - x_i_0_;
dx_1_ = x_j_1_ - x_i_1_;
dx_2_ = x_j_2_ - x_i_2_;
r2 = dx_0_ * dx_0_ + dx_1_ * dx_1_ + dx_2_ * dx_2_ + eps2;
rinv = rsqrt(r2);
mrinv = m_j_ * rinv;
mr3inv = mrinv * rinv * rinv;
a_i_0__ = a_i_0__ + mr3inv * dx_0_;
a_i_1__ = a_i_1__ + mr3inv * dx_1_;
a_i_2__ = a_i_2__ + mr3inv * dx_2_;
a_i_0_ += a_i_0__;
a_i_1_ += a_i_1__;
a_i_2_ += a_i_2__;
