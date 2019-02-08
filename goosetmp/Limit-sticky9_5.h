                      struct f05_jp_t {
double x_j_0_;
double x_j_1_;
double x_j_2_;
double m_j_;
};

struct f05_ip_t {
double x_i_0_;
double x_i_1_;
double x_i_2_;
double padding0__;
};

struct f05_result_t {
double pot_i_;
double padding0__;
};


                      __global__ void f05_calculator(int ioff_, int ni_, int nj_, double eps2, f05_jp_t *f05_jp_, f05_ip_t *f05_ip_, f05_result_t *f05_result_);
                      __global__ void f05_reducer(int njdiv_, int njdiv_ru_, f05_result_t *f05_result_, f05_result_t *f05_result_sub_);
