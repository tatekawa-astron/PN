                      struct f00_jp_t {
double x_j_0_;
double x_j_1_;
double x_j_2_;
double m_j_;
};

struct f00_ip_t {
double x_i_0_;
double x_i_1_;
double x_i_2_;
double padding0__;
};

struct f00_result_t {
double a_i_0_;
double a_i_1_;
double a_i_2_;
double padding0__;
};


                      __global__ void f00_calculator(int ioff_, int ni_, int nj_, double eps2, f00_jp_t *f00_jp_, f00_ip_t *f00_ip_, f00_result_t *f00_result_);
                      __global__ void f00_reducer(int njdiv_, int njdiv_ru_, f00_result_t *f00_result_, f00_result_t *f00_result_sub_);
