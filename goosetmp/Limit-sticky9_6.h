                      struct f06_jp_t {
double x_j_0_;
double x_j_1_;
double x_j_2_;
double m_j_;
double v_j_0_;
double v_j_1_;
double v_j_2_;
double padding0__;
};

struct f06_ip_t {
double x_i_0_;
double x_i_1_;
double x_i_2_;
double v_i_0_;
double v_i_1_;
double v_i_2_;
double m_i_;
double padding0__;
};

struct f06_result_t {
double pot_pn2_i_;
double padding0__;
};


                      __global__ void f06_calculator(int ioff_, int ni_, int nj_, double x_0_0_, double x_0_1_, double x_0_2_, double eps2, double m_0_, f06_jp_t *f06_jp_, f06_ip_t *f06_ip_, f06_result_t *f06_result_);
                      __global__ void f06_reducer(int njdiv_, int njdiv_ru_, f06_result_t *f06_result_, f06_result_t *f06_result_sub_);
