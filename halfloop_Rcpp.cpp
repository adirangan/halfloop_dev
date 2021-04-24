/*
 Sys.setenv("PKG_CXXFLAGS"=" -w -O2 -mfma -I\"/home/rangan/dir_bcc/dir_halfloop_dev/dir_h/\" " , "PKG_LIBS"=" -fopenmp -L\"/home/rangan/dir_bcc/dir_halfloop_dev/\" -lm -lhalfloop")
 Rcpp::sourceCpp('dir_bcc/dir_halfloop_dev/halfloop_Rcpp.cpp' , showOutput = TRUE , rebuild = TRUE )
 */

#include "halfloop_header.h"
#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
IntegerVector halfloop_test() {
  dp_ps_mult_immintrin_test();
  halfloop_nonbinary_f_recursive_omp_test();
  return 0;
}

// [[Rcpp::export]]
IntegerVector halfloop_nonbinary_f_gateway_Rcpp(NumericMatrix x__,IntegerVector flag_r0drop_vs_rcdrop_0in,NumericVector gamma_0in) {
  int verbose=0;
  int nr=0,n_r=0;
  int nc=0,n_c=0;
  float *E_base_rc__=NULL;
  int flag_r0drop_vs_rcdrop = flag_r0drop_vs_rcdrop_0in[0];
  double gamma = gamma_0in[0];
  int n_shuffle = 0;
  double p_set = 0;
  int n_member_lob = 3;
  char *dir_trunk = NULL;
  char *prefix_base = NULL;
  int flag_force_create = 1;
  int flag_omp_use = 1;
  unsigned long long int *binary_label_out_ = NULL;
  n_r = x__.nrow();
  IntegerVector integer_label_out_(n_r);
  n_c = x__.ncol();
  E_base_rc__ = (float *) malloc((size_t)n_r*(size_t)n_c*sizeof(float));
  for (nc=0;nc<n_c;nc++){ for (nr=0;nr<n_r;nr++){ E_base_rc__[nr+nc*n_r] = x__(nr,nc);}}
  if (verbose){ farray_printf_margin(E_base_rc__,n_r,n_c," %% E_base_rc__: ");}
  binary_label_out_ = (unsigned long long int *) malloc((size_t)n_r*sizeof(unsigned long long int));
  memset(binary_label_out_,0,(size_t)n_r*sizeof(unsigned long long int));
  halfloop_nonbinary_f_gateway_matlab
  (
   verbose
  ,n_r
  ,n_c
  ,E_base_rc__
  ,flag_r0drop_vs_rcdrop
  ,gamma
  ,n_shuffle
  ,p_set
  ,n_member_lob
  ,dir_trunk
  ,prefix_base
  ,flag_force_create
  ,flag_omp_use
  ,binary_label_out_
  );
  for (nr=0;nr<n_r;nr++){ 
    if (binary_label_out_[nr]> (unsigned long long int)2147483647){ printf(" %% Warning, label overflow in halfloop_nonbinary_f_gateway_Rcpp\n");}
    integer_label_out_[nr]=binary_label_out_[nr];
    /* for (nr=0;nr<n_r;nr++){ } */}
  free(binary_label_out_);
  free(E_base_rc__);
  return integer_label_out_;
}


/*** R
n_r <- 128
n_c <- 512
A__ <- matrix(c(1:n_r*n_c), nrow = n_r, ncol = n_c)
for (nc in 0:(n_c-1)){
  for (nr in 0:(n_r-1)){
    x <- 2.0*nr/(n_r-1) - 1.0
    y <- 2.0*nc/(n_c-1) - 1.0
    z <- (x+y)/4.0
    A__[1+nr,1+nc] = sin(2*pi*x) + cos(2*pi*2*y) + x*x + y*y*y + cos(2*pi*4*z)      
}
}
binary_label_ <- halfloop_nonbinary_f_gateway_Rcpp(A__,1,0.02)
binary_label_ans_ <- c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,7,7,7,7,7,7,7,7,7,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,3,3,3,3,3,3,3,3,3,3,3,3,3,3,7,7,7,7,7,7,7,0,0,0,0)
binary_label_dif_ <- binary_label_ans_ - binary_label_
binary_label_err <- sqrt(sum(binary_label_dif_^2))
sprintf(" %% r0drop: binary_label_err: %f ",binary_label_err)
binary_label_ <- halfloop_nonbinary_f_gateway_Rcpp(A__,0,0.02)
binary_label_ans_ <- c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,3,3,3,3,3,3,3,3,3,3,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,3,3,3,3,3,3,3,3,0,0,0,0,0)
binary_label_dif_ <- binary_label_ans_ - binary_label_
binary_label_err <- sqrt(sum(binary_label_dif_^2))
sprintf(" %% rcdrop: binary_label_err: %f ",binary_label_err)
#halfloop_test()
*/
