                      struct f01_jp_t {
double x_j_0_;
double x_j_1_;
double x_j_2_;
double m_j_;
double v_j_0_;
double v_j_1_;
double v_j_2_;
double padding0__;
};

struct f01_ip_t {
double x_i_0_;
double x_i_1_;
double x_i_2_;
double m_i_;
double v_i_0_;
double v_i_1_;
double v_i_2_;
double padding0__;
};

struct f01_result_t {
double a1c_0_;
double a1c_1_;
double a1c_2_;
double padding0__;
};


                      __global__ void f01_calculator(int ioff_, int ni_, int nj_, double x_0_0_, double x_0_1_, double x_0_2_, double eps2, double m_0_, f01_jp_t *f01_jp_, f01_ip_t *f01_ip_, f01_result_t *f01_result_);
                      __global__ void f01_reducer(int njdiv_, int njdiv_ru_, f01_result_t *f01_result_, f01_result_t *f01_result_sub_);
