                      struct f02_jp_t {
double x_j_0_;
double v_j_0_;
double x_j_1_;
double v_j_1_;
double x_j_2_;
double v_j_2_;
double m_j_;
double padding0__;
};

struct f02_ip_t {
double x_i_0_;
double v_i_0_;
double x_i_1_;
double v_i_1_;
double x_i_2_;
double v_i_2_;
};

struct f02_result_t {
double ac_i_0_;
double ac_i_1_;
double ac_i_2_;
double padding0__;
};


                      __global__ void f02_calculator(int ioff_, int ni_, int nj_, double x_0_0_, double x_0_1_, double x_0_2_, double eps2, double m_0_, f02_jp_t *f02_jp_, f02_ip_t *f02_ip_, f02_result_t *f02_result_);
                      __global__ void f02_reducer(int njdiv_, int njdiv_ru_, f02_result_t *f02_result_, f02_result_t *f02_result_sub_);
