/*
 Sys.setenv("PKG_CXXFLAGS"=" -w -O2 -mfma -I\"/home/rangan/dir_bcc/dir_halfloop_dev/dir_h/\" " , "PKG_LIBS"=" -fopenmp -L\"/home/rangan/dir_bcc/dir_halfloop_dev/\" -lm -lhalfloop")
 Rcpp::sourceCpp('dir_bcc/dir_halfloop_dev/halfloop_Rcpp.cpp' , showOutput = TRUE , rebuild = TRUE )
 */

#include "halfloop_header.h"
#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
IntegerVector halfloop_test() {
  //dp_ps_mult_immintrin_test();
  halfloop_nonbinary_f_recursive_omp_test();
  return 0;
}


// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically 
// run after the compilation.
//

/*** R
halfloop_test()
*/
