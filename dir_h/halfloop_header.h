#include <fcntl.h>
#include <math.h>
#include <omp.h>
#include <inttypes.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/mman.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <sys/time.h>
#include <sys/resource.h>
#include <time.h>
#include <unistd.h>
#include <errno.h>
#include <assert.h>
#include <emmintrin.h>
#include <immintrin.h>
#ifdef _CBLAS
#include <cblas.h>
#endif /* _CBLAS */
#ifdef _COMPLEX
#include <complex.h>
#endif /* _COMPLEX */

#define PI_LF 3.141592653589793
#define FNAMESIZE (4096+256)
#define PNAMESIZE (8192+512)
#define BIT8 8
#define bget_off(bmr,nr) (((bmr[(nr)/BIT8] >> (7-((nr)%BIT8))) & 1) ^ 1)
#define bget_0on(bmr,nr) (((bmr[(nr)/BIT8] >> (7-((nr)%BIT8))) & 1) ^ 0)
#define bget____(bmr,nr) ((int)(((bmr[(nr)/BIT8] >> (7-((nr)%BIT8))) & 1) << 1)-(int)1)
#define bset_off(bmr,nr) (bmr[(nr)/BIT8] &= ~(1 << (7 - ((nr)%BIT8))))
#define bset_0on(bmr,nr) (bmr[(nr)/BIT8] |=  (1 << (7 - ((nr)%BIT8))))
#define uchar_mask_size(A) (((A) + ((BIT8 - ((A) % BIT8)) % BIT8))/BIT8)

/* ---------------------------------------------------------------- */

/* global variables used for timing */
#define GLOBAL_NTICKS 8
extern clock_t GLOBAL_t_start[GLOBAL_NTICKS];
extern clock_t GLOBAL_t_final[GLOBAL_NTICKS];
extern struct timeval GLOBAL_d_start[GLOBAL_NTICKS];
extern struct timeval GLOBAL_d_final[GLOBAL_NTICKS];
extern long GLOBAL_l_msec[GLOBAL_NTICKS],GLOBAL_l_ssec[GLOBAL_NTICKS],GLOBAL_l_usec[GLOBAL_NTICKS];
extern double GLOBAL_elct[GLOBAL_NTICKS],GLOBAL_elrt[GLOBAL_NTICKS];

/* global variables used for some routines */
extern int GLOBAL_verbose;
extern char GLOBAL_mode[FNAMESIZE];
extern double GLOBAL_tolerance;
extern unsigned int GLOBAL_quicksort_recursion_limit;
extern unsigned int GLOBAL_halfloop_recursion_limit;
extern int addressable_1;
extern int addressable_0;
extern int addressable_int_length;
extern int addressable_int[128];
extern int GLOBAL_flag_orth_brute;
extern char GLOBAL_E_base_mda_r4[PNAMESIZE];
extern int GLOBAL_flag_r0drop_vs_rcdrop;
extern double GLOBAL_gamma;
extern int GLOBAL_n_shuffle;
extern double GLOBAL_p_set;
extern int GLOBAL_n_member_lob;
extern char GLOBAL_dir_trunk[PNAMESIZE];
extern char GLOBAL_prefix_base[FNAMESIZE];
extern int GLOBAL_flag_force_create;
extern int GLOBAL_flag_omp_use;

#define rup(A,B) ((A) + !!((A)%(B))*((B) - ((A)%(B))))
#define maximum(A,B) ((A) > (B) ? (A) : (B))
#define minimum(A,B) ((A) < (B) ? (A) : (B))
#define periodize(A,B,C) ((A) < (B) ? (A) + (C) - (B) : ((A) >= (C) ? (A) - (C) + (B) : (A)))
#define crop(A,B,C) ((A) < (B) ? (B) : ((A) > (C) ? (C) : (A)))
#define rand01 ((double)rand()/(double)RAND_MAX)
#define bsize(A) ((rup((A) + ((BITJ - ((A) % BITJ)) % BITJ),POPLENGTH))/BIT8)

/* ---------------------------------------------------------------- */
extern int GLOBAL_malloc1_notupdate;
extern unsigned long long int GLOBAL_n_malloc1;
extern int GLOBAL_n_malloc1_[GLOBAL_NTICKS];

/* ---------------------------------------------------------------- */

/* RAND functions */
extern unsigned long int POW2RPOWPLUSRADD;
extern unsigned long int POW22RPOWMINUSONE;

/* ---------------------------------------------------------------- */

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

void array_extract_i_from_i(int n_r,int n_c,int *A_rc__,int n_r_rtn,int *r_rtn_,int n_c_rtn,int *c_rtn_,int **B_rc_p_,int **B_cr_p_);
void array_extract_i_from_i_test();
void array_extract_f_from_d(int n_r,int n_c,double *A_rc__,int n_r_rtn,int *r_rtn_,int n_c_rtn,int *c_rtn_,float **B_rc_p_,float **B_cr_p_);
void array_extract_f_from_d_test();
void array_extract_d_from_d(int n_r,int n_c,double *A_rc__,int n_r_rtn,int *r_rtn_,int n_c_rtn,int *c_rtn_,double **B_rc_p_,double **B_cr_p_);
void array_extract_d_from_d_test();
void array_extract_f_from_f(int n_r,int n_c,float *A_rc__,int n_r_rtn,int *r_rtn_,int n_c_rtn,int *c_rtn_,float **B_rc_p_,float **B_cr_p_);
void array_extract_f_from_f_test();
void array_extract_d_from_f(int n_r,int n_c,float *A_rc__,int n_r_rtn,int *r_rtn_,int n_c_rtn,int *c_rtn_,double **B_rc_p_,double **B_cr_p_);
void array_extract_d_from_f_test();
void array_extract(int n_r,int n_c,void *A_rc__,const char *type_A,int n_r_rtn,int *r_rtn_,int n_c_rtn,int *c_rtn_,void *B_rc_p_,void *B_cr_p_,const char *type_B);
void array_extract_test();
void array_mean_center_row_f(int n_r,int n_c,float *f_rc_0in__,float *f_cr_0in__,float **f_rc_out_p_,float **f_cr_out_p_);
void array_mean_center_row_d(int n_r,int n_c,double *d_rc_0in__,double *d_cr_0in__,double **d_rc_out_p_,double **d_cr_out_p_);
void array_mean_center_row(int n_r,int n_c,void *v_rc_0in__,void *v_cr_0in__,const char *type,void *v_rc_out_p_,void *v_cr_out_p_);
void array_mean_center_row_test();
void array_normalize_row_f(int n_r,int n_c,float *f_rc_0in__,float *f_cr_0in__,float **f_rc_out_p_,float **f_cr_out_p_);
void array_normalize_row_d(int n_r,int n_c,double *d_rc_0in__,double *d_cr_0in__,double **d_rc_out_p_,double **d_cr_out_p_);
void array_normalize_row(int n_r,int n_c,void *v_rc_0in__,void *v_cr_0in__,const char *type,void *v_rc_out_p_,void *v_cr_out_p_);
void array_normalize_row_test();
void array_gram_schmidt_inplace_f(int n_r,int n_c,float *Q_rc__);
void array_orth_f(int n_r,int n_c,float **Q_rc_p_,unsigned long int *rseed_p);
void array_orth_f_test_error();
void array_orth_f_test_speed();
void array_gram_schmidt_inplace_d(int n_r,int n_c,double *Q_rc__);
void array_orth_d(int n_r,int n_c,double **Q_rc_p_,unsigned long int *rseed_p);
void array_orth_d_test_error();
void array_orth_d_test_speed();
void iarray_printf_margin(int *i_,int n_r,int n_c,const char *prefix);
void farray_printf_margin(float *f_,int n_r,int n_c,const char *prefix);
void darray_printf_margin(double *d_,int n_r,int n_c,const char *prefix);
void array_printf_margin(void *v_,const char *type,int n_row,int n_col,const char *prefix);
void array_printf(void *v_,const char *type,int n_row,int n_col,const char *prefix);
void array_fprintf(const char *fname,void *v_,const char *type,int n_row,int n_col,const char *prefix);
void bitstring_from_uchar(unsigned char *w_, char *str_, int k);
void bitstring_from_uchar_printf(unsigned char *w_,int nrows,int ncols,const char *prefix);
void icumsum(unsigned long long int ulli_length,int *i_0in_,int **i_out_p_);
double ifnormn(unsigned long long int ulli_length,int *i_0_,int *i_1_);
double ullifnormn(unsigned long long int ulli_length,unsigned long long int *ulli_0_,unsigned long long int *ulli_1_);
float ffnorm(unsigned long long int ulli_length,float *f_0_,float *f_1_);
float ffnormn(unsigned long long int ulli_length,float *f_0_,float *f_1_);
double dfnormn(unsigned long long int ulli_length,double *d_0_,double *d_1_);
double dfnorm(unsigned long long int ulli_length,double *d_0_,double *d_1_);
#ifdef _COMPLEX
double cfnorm(unsigned long long int ulli_length,float complex *c_0_,float complex *c_1_);
double zfnorm(unsigned long long int ulli_length,double complex *z_0_,double complex *z_1_);
#endif /* _COMPLEX */
void array_maximum_minimum(void *v_,const char *type,unsigned long long int ulli_length,void *max_p,int *index_max_p,void *min_p,int *index_min_p);
void array_stats(void *v_,const char *type,unsigned long long int ulli_length,void *max_p,void *min_p,double *mean_p,double *stdev_p);
void dtranspose_bruteforce(int n_row_A,int n_col_A,double *d_0in__,double *d_out__);
void dtranspose_block_AtoB(const int n_row_A,const int n_col_A,const double* A_,double* B_,const int block_size);
void dtranspose_block_BtoA(const int n_row_A,const int n_col_A,const double* A_,double* B_,const int block_size);
void dtranspose(const int n_row_A,const int n_col_A,const double* A_,const double* B_);
void dtranspose_test();
#ifdef _COMPLEX
void cntranspose_bruteforce(int n_row_A,int n_col_A,float complex *c_0in__,float complex *c_out__);
void cntranspose_block_AtoB(const int n_row_A,const int n_col_A,const float complex* A_,float complex* B_,const int block_size);
void cntranspose_block_BtoA(const int n_row_A,const int n_col_A,const float complex* A_,float complex* B_,const int block_size);
void cntranspose(const int n_row_A,const int n_col_A,const float complex* A_,float complex* B_);
void cntranspose_test();
void cctranspose_bruteforce(int n_row_A,int n_col_A,float complex *c_0in__,float complex *c_out__);
void cctranspose_block_AtoB(const int n_row_A,const int n_col_A,const float complex* A_,float complex* B_,const int block_size);
void cctranspose_block_BtoA(const int n_row_A,const int n_col_A,const float complex* A_,float complex* B_,const int block_size);
void cctranspose(const int n_row_A,const int n_col_A,const float complex* A_,float complex* B_);
void cctranspose_test();
#endif /* _COMPLEX */
unsigned long int lrand();
double randn();
unsigned long int RGET(unsigned long int *rseed_p);
double R01GET(unsigned long int *rseed_p);
void RSEED_adv8(unsigned long int *rseed_p);
double RNGET(unsigned long int *rseed_p);
double RISIGET(unsigned long int *rseed_p,double rate);
void R01GET_test();
void dp_pd_bruteforce(int n_col_X,double *d_A_,double *d_B_,double *d_C_);
void dp_pd_mult_bruteforce(int n_row_A,int n_col_X,double *d_A_trn__,int n_row_B,double *d_B_trn__,double **d_C_p_);
void dp_pd_immintrin_loadu_wrap(int n_col_X,double *d_A_,double *d_B_,double *d_C_);
void dp_pd_mult_immintrin_loadu_wrap(int n_row_A,int n_col_X,double *d_A_trn__,int n_row_B,double *d_B_trn__,double **d_C_p_);
void dp_pd_immintrin_loadu_fma(int n_col_X,double *d_A_,double *d_B_,double *d_C_);
void dp_pd_mult_immintrin_loadu_fma(int n_row_A,int n_col_X,double *d_A_trn__,int n_row_B,double *d_B_trn__,double **d_C_p_);
void dp_pd_mult_immintrin_test();
void dp_ps_mult_cblas_sgemm(int n_row_A,int n_col_X,float *f_A_trn__,int n_row_B,float *f_B_trn__,float **f_C_p_);
void dp_ps_mult_cblas_sdot(int n_row_A,int n_col_X,float *f_A_trn__,int n_row_B,float *f_B_trn__,float **f_C_p_);
void dp_ps_bruteforce(int n_col_X,float *f_A_,float *f_B_,float *f_C_);
void dp_ps_mult_bruteforce(int n_row_A,int n_col_X,float *f_A_trn__,int n_row_B,float *f_B_trn__,float **f_C_p_);
void dp_ps_immintrin_loadu_wrap(int n_col_X,float *f_A_,float *f_B_,float *f_C_);
void dp_ps_mult_immintrin_loadu_wrap(int n_row_A,int n_col_X,float *f_A_trn__,int n_row_B,float *f_B_trn__,float **f_C_p_);
void dp_ps_immintrin_loadu_avx(int n_col_X,float *f_A_,float *f_B_,float *f_C_);
void dp_ps_mult_immintrin_loadu_avx(int n_row_A,int n_col_X,float *f_A_trn__,int n_row_B,float *f_B_trn__,float **f_C_p_);
void dp_ps_immintrin_loadu_fma(int n_col_X,float *f_A_,float *f_B_,float *f_C_);
void dp_ps_mult_immintrin_loadu_fma(int n_row_A,int n_col_X,float *f_A_trn__,int n_row_B,float *f_B_trn__,float **f_C_p_);
void dp_ps_mult_immintrin_load1_fma(int n_row_A,int n_col_X,float *f_A_trn__,int n_row_B,float *f_B_trn__,float **f_C_p_);
void dp_ps_mult_immintrin_fma2(int n_row_A,int n_col_X,__m256 *ps_A_trn__,int n_row_B,__m256 *ps_B_trn__,float **f_C_p_);
void dp_ps_mult_immintrin_fma(int n_row_A,int n_col_X,__m256 *ps_A_trn__,int n_row_B,__m256 *ps_B_trn__,float **f_C_p_);
void dp_ps_mult_immintrin_avx(int n_row_A,int n_col_X,__m256 *ps_A_trn__,int n_row_B,__m256 *ps_B_trn__,float **f_C_p_);
void dp_ps_mult_immintrin_test();
void dp_ps_single_test();
void hp_interleave_mult_bruteforce(int n_row_A,int n_col_X,float complex *c_A_trn__,int n_row_B,float complex *c_B_trn__,float complex **c_C_p_);
void hp_interleave_mult_cblas_cgemm(int n_row_A,int n_col_X,float complex *c_A_trn__,int n_row_B,float complex *c_B_trn__,float complex **c_C_p_);
void hp_segregated_mult_bruteforce(int n_row_A,int n_col_X,float *f_AR_trn__,float *f_AI_trn__,int n_row_B,float *f_BR_trn__,float *f_BI_trn__,float complex **c_C_p_);
void hp_ps_mult_immintrin_test();
float erfcln_single_f(float f_0in);
void erfcln_f(int n_r,float *f_0in_,float **f_out_p_);
void erfcln_f_test();
double erfcln_single_d(double d_0in);
void erfcln_d(int n_r,double *d_0in_,double **d_out_p_);
void erfcln_d_test();
double z_to_lp_single_d(double d_0in);
void z_to_lp_d(int n_r,double *d_0in_,double **d_out_p_);
void z_to_lp_d_test();
int is_internal_maximum(int n_Z,double *Z_,int index);
void find_local_maxima(int n_Z,double *Z_,int *n_index_p,int **index_p_);
void find_internal_maximum(int verbose,int n_Z,double *Z_,double Z_min,double *zone_max_p,int *zone_max_index_p);
void find_internal_maximum_test();
void get_xdrop_logscale(double n_row,double n_col,double gamma,int *rdrop_p,int *cdrop_p);
int get_xdrop_logscale_length(double n_row,double n_col,double gamma);
void get_xdrop_logscale_array(double n_row,double n_col,double gamma,int *length_p,int **rdrop_p_,int **cdrop_p_,int **rkeep_p_,int **ckeep_p_);
void get_xdrop_logscale_array_test();
void GLOBAL_tic(int nx);
void GLOBAL_toc(int nx,int verbose,const char *prefix);
void gumbel_nll(int n_x,double *x_,double *g_,double tol,double **nll_p_,double *nll_sum_p);
void gumbel_nll_wrap(int n_g,double *g_,double *nll_sum_p,void *args);
void gumbel_nll_test();
void gumbel_pdf(int n_x,double *x_,double *g_,double tol,double **pdf_p_);
void gumbel_cdf(int n_x,double *x_,double *g_,double tol,double **cdf_p_);
double gumbel_cdf_single(double x,double *g_,double tol);
void gumbel_fit(int n_x,double *x_,double x,double **g_opt_p_,double *nlp_opt_p,double *nlp_emp_p,double *p_opt_p,double *p_emp_p);
void gumbel_fit_test();
void halfloop_nonbinary_f(int n_r,int n_c,float *A1_rc__,float *A1_cr__,double gamma,int **out_xdrop_p_,int *n_iteration_p,double **trace_p_);
void halfloop_nonbinary_d(int n_r,int n_c,double *B1_rc__,double *B1_cr__,double gamma,int **out_xdrop_p_,int *n_iteration_p,double **trace_p_);
void halfloop_nonbinary_test_error();
void halfloop_nonbinary_test_speed();
void halfloop_nonbinary_rdrop_f_bkp(int n_r,int n_c,float *A1_rc__,float *A1_cr__,double gamma,int **out_xdrop_p_,int *n_iteration_p,double **trace_p_);
void halfloop_nonbinary_rdrop_f(int n_r,int n_c,float *A1_rc__,float *A1_cr__,double gamma,int **out_xdrop_p_,int *n_iteration_p,double **trace_p_);
void halfloop_nonbinary_rdrop_d_bkp(int n_r,int n_c,double *A1_rc__,double *A1_cr__,double gamma,int **out_xdrop_p_,int *n_iteration_p,double **trace_p_);
void halfloop_nonbinary_rdrop_d(int n_r,int n_c,double *A1_rc__,double *A1_cr__,double gamma,int **out_xdrop_p_,int *n_iteration_p,double **trace_p_);
void halfloop_nonbinary_rdrop_test_error();
void halfloop_nonbinary_rdrop_test_speed();
void halfloop_nonbinary_f_recursive_helper_QR_helper_orth
(
  int verbose
 ,int n_r_index
 ,int n_c_index
 ,float *E_rc__
 ,float *QE_rc__
 ,unsigned long long int *rseed_p
);
void halfloop_nonbinary_f_recursive_helper_QR_helper_orth_test_speed();
void halfloop_nonbinary_f_recursive_helper_QR__
(
  int verbose
 ,int n_r
 ,int n_c
 ,float * E_base_rc__
 ,int n_r_index
 ,int *r_index_
 ,int n_c_index
 ,int *c_index_
 ,int flag_r0drop_vs_rcdrop
 ,double gamma
 ,int n_shuffle
 ,int flag_force_create
 ,char *fname_trace__
 ,char *fname_xdrop__
 ,char *fname_QR__
);
void halfloop_nonbinary_recursive_helper_ZR__
(
 int verbose
 ,int flag_r0drop_vs_rcdrop
 ,char *fname_trace__
 ,char *fname_xdrop__
 ,char *fname_QR__
 ,double p_use
 ,int n_member_lob
 ,double *nlp_ZR_max_p
 ,int *nlp_ZR_index_p
 ,int *n_r_rtn_index_p
 ,int **r_rtn_index_p_
 ,int *n_r_rmv_index_p
 ,int **r_rmv_index_p_
 ,int *n_c_rtn_index_p
 ,int **c_rtn_index_p_
 ,int *n_c_rmv_index_p
 ,int **c_rmv_index_p_
 ,double *nlp_gumb_opt_p
 ,double *nlp_gumb_emp_p
);
void halfloop_nonbinary_f_recursive_omp_helper_QR_
(
 int verbose
 ,int n_r_index
 ,int n_c_index
 ,double *E_rc__
 ,double *E_cr__
 ,int flag_r0drop_vs_rcdrop
 ,double gamma
 ,int nshuffle
 ,int n_iteration
 ,double *trace__
 ,double *QR_
 ,int n_xdrop
 ,double *xdrop__
);
void halfloop_nonbinary_f_recursive_omp_helper_QR__
(
  int verbose
 ,int n_r
 ,int n_c
 ,float * E_base_rc__
 ,int n_r_index
 ,int *r_index_
 ,int n_c_index
 ,int *c_index_
 ,int flag_r0drop_vs_rcdrop
 ,double gamma
 ,int n_shuffle
 ,int flag_force_create
 ,char *fname_trace__
 ,char *fname_xdrop__
 ,char *fname_QR__
);
void halfloop_nonbinary_f_recursive
(
 int verbose
 ,int flag_omp
 ,int ndepth
 ,int recursion_limit
 ,int n_r
 ,int n_c
 ,float *E_base_rc__
 ,int n_r_index_0in
 ,int *r_index_0in_
 ,int n_c_index_0in
 ,int *c_index_0in_
 ,int flag_r0drop_vs_rcdrop
 ,double gamma_0in
 ,int n_shuffle_0in
 ,double p_set_0in
 ,int n_member_lob_0in
 ,double p_prev_0in
 ,char *dir_trunk_0in
 ,char *dir_out_0in
 ,char *prefix_base_0in
 ,int flag_force_create_0in
 ,unsigned long long int **binary_label_p_
 ,char ***output_label_p_
 ,char ***nlpbra_label_p_
 ,char ***nlpnex_label_p_
);
void halfloop_nonbinary_f_recursive_test();
void halfloop_nonbinary_f_recursive_test_speed();
void halfloop_nonbinary_f_gateway_shell();
void * malloc1(size_t size);
void free1(void **vp);
void malloc1_char_FNAMESIZE__(int n_l,char ***str_p_);
void malloc1_char__(int n_l,int n_s,char ***str_p_);
void free1_char__(int n_l,char ***str_p_);
void MDA_write_i4(int n_dim,int *dim_,int *i4_,const char *fname);
void MDA_read_i4(int *n_dim_p,int **dim_p_,int **i4_p_,const char *fname);
void MDA_printf_r4_margin(const char *fname);
void MDA_write_r4(int n_dim,int *dim_,float *r4_,const char *fname);
void MDA_read_r4(int *n_dim_p,int **dim_p_,float **r4_p_,const char *fname);
void MDA_printf_i4_margin(const char *fname);
void MDA_write_r8(int n_dim,int *dim_,double *r8_,const char *fname);
void MDA_read_r8(int *n_dim_p,int **dim_p_,double **r8_p_,const char *fname);
void MDA_printf_r8_margin(const char *fname);
void MDA_write_ulli(int n_dim,int *dim_,unsigned long long int *ulli_,const char *fname);
void MDA_read_ulli(int *n_dim_p,int **dim_p_,unsigned long long int **ulli_p_,const char *fname);
void MDA_printf_ulli_margin(const char *fname);
void MDA_io_test();
void nelder_mead_terminate(int n_d,double *simplex_point__,int *index_simplex_point_,double *simplex_cost_,double option_tolx,double option_tolf,int *flag_continue_p);
void nelder_mead_simplex_centroid(int n_d,double *simplex_point__,int *index_simplex_point_,double *simplex_centroid_);
void nelder_mead_update_point(int n_d,double *simplex_point_,double *simplex_centroid_,double lambda,double *point_);
void nelder_mead_optimization(int n_d,double *point_start_,double *point_final_, void cost_function(int,double *,double *,void *),void *cost_function_args,double option_tolx,double option_tolf,int option_maxiter,int option_maxfeval);
void nelder_mead_test();
int iquicksort_partition_index(int *i_,int stride,int *index_,int l,int r);
unsigned int iquicksort_index(unsigned int recursion_level,int *i_,int stride,int *index_,int l,int r);
unsigned int iquicksort_index_driver(int n_i,int *i_,int stride,int *i_workspace_,int *index_);
void iquicksort_index_driver_test();
unsigned int iquicksort_index_index_driver(int n_i,int *i_,int stride,int *i_workspace_,int *index_orig_from_sort_,int *index_sort_from_orig_,int *index_workspace_);
void iquicksort_index_index_driver_test();
void irandperm(int n_i,int **index_p_,unsigned long long int *rseed);
void irandperm_test();
int fquicksort_partition_index(float *f_,int stride,int *index_,int l,int r);
unsigned int fquicksort_index(unsigned int recursion_level,float *f_,int stride,int *index_,int l,int r);
unsigned int fquicksort_index_driver(int n_f,float *f_,int stride,float *f_workspace_,int *index_);
void fquicksort_index_driver_test();
unsigned int fquicksort_index_index_driver(int n_f,float *f_,int stride,float *f_workspace_,int *index_orig_from_sort_,int *index_sort_from_orig_,int *index_workspace_);
void fquicksort_index_index_driver_test();
int dquicksort_partition_index(double *d_,int stride,int *index_,int l,int r);
unsigned int dquicksort_index(unsigned int recursion_level,double *d_,int stride,int *index_,int l,int r);
unsigned int dquicksort_index_driver(int n_d,double *d_,int stride,double *d_workspace_,int *index_);
void dquicksort_index_driver_test();
unsigned int dquicksort_index_index_driver(int n_d,double *d_,int stride,double *d_workspace_,int *index_orig_from_sort_,int *index_sort_from_orig_,int *index_workspace_);
void dquicksort_index_index_driver_test();
void update_global(char *vname);
void read_input();
void ping();
void pong();

#ifdef __cplusplus
}
#endif /* __cplusplus */




