// Written by Aaditya Rangan, with thanks to Chris Chang (see https://www.cog-genomics.org/plink2/). ;
// This program is free software: ;
// you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
// This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  
// See the GNU General Public License for more details.
// You should have received a copy of the GNU General Public License along with this program.  
// If not, see <http://www.gnu.org/licenses/>.

#ifndef _MONOLITH
#include "halfloop_header.h"
#endif /* _MONOLITH */
#ifdef _MONOLITH
#include "halfloop_header.h"
#endif /* _MONOLITH */

/* global variables used for timing */
clock_t GLOBAL_t_start[GLOBAL_NTICKS],GLOBAL_t_final[GLOBAL_NTICKS];
struct timeval GLOBAL_d_start[GLOBAL_NTICKS],GLOBAL_d_final[GLOBAL_NTICKS];
long GLOBAL_l_msec[GLOBAL_NTICKS],GLOBAL_l_ssec[GLOBAL_NTICKS],GLOBAL_l_usec[GLOBAL_NTICKS]; 
double GLOBAL_elct[GLOBAL_NTICKS],GLOBAL_elrt[GLOBAL_NTICKS];

int GLOBAL_verbose=0;
double GLOBAL_tolerance=0.000000000001; //%<-- 1e-12;
unsigned int GLOBAL_recursion_limit=1024*32; //%<-- 2^15;
int addressable_1=1;
int addressable_0=0;
int addressable_int_length=128;
int addressable_int[128] = {0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52,53,54,55,56,57,58,59,60,61,62,63,64,65,66,67,68,69,70,71,72,73,74,75,76,77,78,79,80,81,82,83,84,85,86,87,88,89,90,91,92,93,94,95,96,97,98,99,100,101,102,103,104,105,106,107,108,109,110,111,112,113,114,115,116,117,118,119,120,121,122,123,124,125,126,127};

/* ---------------------------------------------------------------- */
int GLOBAL_malloc1_notupdate=0;
unsigned long long int GLOBAL_n_malloc1=0;
int GLOBAL_n_malloc1_[GLOBAL_NTICKS];

/* ---------------------------------------------------------------- */

/* RAND functions */
unsigned long int POW2RPOWPLUSRADD=35L;
unsigned long int POW22RPOWMINUSONE=2147483647LL;

/*---------------------------------------------------------------- */

#ifdef _MONOLITH
#include "halfloop_function.c"
#endif /* _MONOLITH */

inline void ping(){ printf(" %% ping\n");}
inline void pong(){ printf(" %% pong\n");}

int main(int argc, char** argv) {
  int cur_arg = 1,ii=0;
  char* argptr;
  int retval=0;
  int cmdline_param=0;
  int flag_test=0;
  while (cur_arg < argc) {
    argptr = argv[cur_arg];
    if (0){ /* do nothing */}
    else if (!strcmp(argptr, "--test")) {
      if (cur_arg == argc - 1) { printf("Error: Missing --test parameter."); exit(EXIT_FAILURE);}
      ii = atoi(argv[cur_arg + 1]);
      if (ii < 1) { printf("Error: Invalid --test parameter.\n"); exit(EXIT_FAILURE);}
      flag_test = ii; printf(" %% flag_test %d\n",flag_test);
      cur_arg += 2; /* else if flag_test */}
    else { printf("Error: Invalid argument (%s).", argv[cur_arg]); exit(EXIT_FAILURE);}
    /* while (cur_arg < argc) { } */}
  //MDA_io_test();
  //R01GET_test();
  //dtranspose_test();
  //dp_ps_single_test();
  //dp_ps_mult_immintrin_test();
  //dp_pd_mult_immintrin_test();
  //get_xdrop_logscale_array_test();
  //iquicksort_index_driver_test();
  //fquicksort_index_driver_test();
  //dquicksort_index_driver_test();
  //iquicksort_index_index_driver_test();
  //fquicksort_index_index_driver_test();
  //dquicksort_index_index_driver_test();
  //irandperm_test();
  //halfloop_nonbinary_test_error();
  //halfloop_nonbinary_rdrop_test_error();
  //gumbel_nll_test();
  //nelder_mead_test();
  //gumbel_fit_test();
  //array_extract_i_from_i_test();
  //array_extract_f_from_d_test();
  //array_extract_d_from_d_test();
  //array_extract_f_from_f_test();
  //array_extract_d_from_f_test();
  //array_extract_test();
  //array_mean_center_row_test();
  //array_normalize_row_test();
  //array_orth_f_test_error();
  //array_orth_f_test_speed();
  //array_orth_d_test_error();
  //array_orth_d_test_speed();
  //erfcln_f_test();
  //erfcln_d_test();
  //z_to_lp_d_test();
  //find_internal_maximum_test();
  halfloop_nonbinary_f_recursive_test();
  halfloop_nonbinary_f_recursive_omp_test();
  if (GLOBAL_verbose>-1){ printf("exiting successfully\n");}
  return 0;
}
