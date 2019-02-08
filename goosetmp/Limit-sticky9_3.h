                      struct f03_jp_t {
double x_j_0_;
double x_j_1_;
double x_j_2_;
double m_j_;
};

struct f03_ip_t {
double x_i_0_;
double x_i_1_;
double x_i_2_;
double padding0__;
};

struct f03_result_t {
double pot_i_;
double padding0__;
};


                      __global__ void f03_calculator(int ioff_, int ni_, int nj_, double eps2, f03_jp_t *f03_jp_, f03_ip_t *f03_ip_, f03_result_t *f03_result_);
                      __global__ void f03_reducer(int njdiv_, int njdiv_ru_, f03_result_t *f03_result_, f03_result_t *f03_result_sub_);
