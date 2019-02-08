CONST eps2;
VARI x_i_0_, x_i_1_, x_i_2_;
VARJ x_j_0_, x_j_1_, x_j_2_, m_j_;
VARF pot_i_;
pot_i__ = 0.0;

dx_0_ = 0.0;  // the value 0.0 is meaningless.
dx_1_ = 0.0;  // the value 0.0 is meaningless.
dx_2_ = 0.0;  // the value 0.0 is meaningless.
r2 = 0.0;  // the value 0.0 is meaningless.
rinv = 0.0;  // the value 0.0 is meaningless.
mrinv = 0.0;  // the value 0.0 is meaningless.
dx_0_ = x_j_0_ - x_i_0_;
dx_1_ = x_j_1_ - x_i_1_;
dx_2_ = x_j_2_ - x_i_2_;
r2 = dx_0_ * dx_0_ + dx_1_ * dx_1_ + dx_2_ * dx_2_ + eps2;
rinv = rsqrt(r2);
mrinv = rinv * m_j_;
pot_i__ = pot_i__ - mrinv;
pot_i_ += pot_i__;
