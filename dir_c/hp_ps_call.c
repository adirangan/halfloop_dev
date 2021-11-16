#ifndef _MONOLITH
#include "halfloop_header.h"
#endif /* _MONOLITH */
/* Trying to build a faster complex multiply -- not there yet! */
/* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */

void hp_interleave_mult_bruteforce(int n_row_A,int n_col_X,float complex *c_A_trn__,int n_row_B,float complex *c_B_trn__,float complex **c_C_p_)
{
  /* unaligned */
  int nrow_A=0,nrow_B=0,ncol_X=0;
  float complex *c_C__=NULL;
  c_C__=NULL;
  if (c_C_p_!=NULL){
    if ( (*c_C_p_)==NULL ){ (*c_C_p_) = (float complex *) malloc1((unsigned long long int)n_row_A*(unsigned long long int)n_row_B*sizeof(float complex));}
    c_C__ = *c_C_p_;
    /* if (c_C_p_!=NULL){ } */}
  if (c_C__!=NULL){
    for (nrow_B=0;nrow_B<n_row_B;nrow_B++){
      for (nrow_A=0;nrow_A<n_row_A;nrow_A++){
	for (ncol_X=0;ncol_X<n_col_X;ncol_X++){
	  c_C__[nrow_A + nrow_B*n_row_A] += conjf(c_A_trn__[ncol_X + nrow_A*n_col_X]) * c_B_trn__[ncol_X + nrow_B*n_col_X];
	  /* for (ncol_X=0;ncol_X<n_col_X;ncol_X++){ } */}
	/* for (nrow_A=0;nrow_A<n_row_A;nrow_A++){ } */}
      /* for (nrow_B=0;nrow_B<n_row_B;nrow_B++){ } */}
    /* if (c_C__!=NULL){ } */}
}

void hp_interleave_mult_cblas_cgemm(int n_row_A,int n_col_X,float complex *c_A_trn__,int n_row_B,float complex *c_B_trn__,float complex **c_C_p_)
{
  /* unaligned */
  float complex *c_C__=NULL;
  float complex calpha = (float complex)1.0;
  float complex cbeta = (float complex)0.0;
  c_C__=NULL;
  if (c_C_p_!=NULL){
    if ( (*c_C_p_)==NULL ){ (*c_C_p_) = (float complex *) malloc1((unsigned long long int)n_row_A*(unsigned long long int)n_row_B*sizeof(float complex));}
    c_C__ = *c_C_p_;
    /* if (c_C_p_!=NULL){ } */}
  if (c_C__!=NULL){
#ifdef _CBLAS
    cblas_cgemm(CblasColMajor,CblasConjTrans,CblasNoTrans,n_row_A,n_row_B,n_col_X,&calpha,c_A_trn__,n_col_X,c_B_trn__,n_col_X,&cbeta,c_C__,n_row_A);
#endif /* _CBLAS */
#ifndef _CBLAS
    printf(" %% Warning, _CBLAS not defined in hp_interleave_mult_cblas_cgemm, using bruteforce instead.\n");
    hp_interleave_mult_bruteforce(n_row_A,n_col_X,c_A_trn__,n_row_B,c_B_trn__,c_C_p_);
#endif /* _CBLAS */
  /* if (c_C__!=NULL){ } */}
}

void hp_segregated_mult_bruteforce(int n_row_A,int n_col_X,float *f_AR_trn__,float *f_AI_trn__,int n_row_B,float *f_BR_trn__,float *f_BI_trn__,float complex **c_C_p_)
{
  /* unaligned */
  int nrow_A=0,nrow_B=0,ncol_X=0;
  float complex *c_C__=NULL;
  c_C__=NULL;
  if (c_C_p_!=NULL){
    if ( (*c_C_p_)==NULL ){ (*c_C_p_) = (float complex *) malloc1((unsigned long long int)n_row_A*(unsigned long long int)n_row_B*sizeof(float complex));}
    c_C__ = *c_C_p_;
    /* if (c_C_p_!=NULL){ } */}
  if (c_C__!=NULL){
    for (nrow_B=0;nrow_B<n_row_B;nrow_B++){
      for (nrow_A=0;nrow_A<n_row_A;nrow_A++){
	c_C__[nrow_A + nrow_B*n_row_A] = (float complex) 0.0;
	for (ncol_X=0;ncol_X<n_col_X;ncol_X++){
	  c_C__[nrow_A + nrow_B*n_row_A] +=
	    ((float complex)f_AR_trn__[ncol_X + nrow_A*n_col_X] - _Complex_I * (float complex)f_AI_trn__[ncol_X + nrow_A*n_col_X])
	    *
	    ((float complex)f_BR_trn__[ncol_X + nrow_B*n_col_X] + _Complex_I * (float complex)f_BI_trn__[ncol_X + nrow_B*n_col_X])
	    ;
	  /* for (ncol_X=0;ncol_X<n_col_X;ncol_X++){ } */}
	/* for (nrow_A=0;nrow_A<n_row_A;nrow_A++){ } */}
      /* for (nrow_B=0;nrow_B<n_row_B;nrow_B++){ } */}
    /* if (c_C__!=NULL){ } */}
}

void hp_segregated_mult_immintrin_loadu_fma(int n_row_A,int n_col_X,float *f_AR_trn__,float *f_AI_trn__,int n_row_B,float *f_BR_trn__,float *f_BI_trn__,float complex **c_C_p_)
{
  /* does not assume alignment */
  int nrow_A=0,nrow_B=0;
  float *tmp_f_AR_point0_=NULL;
  float *tmp_f_AI_point0_=NULL;
  float *tmp_f_BR_point0_=NULL;
  float *tmp_f_BI_point0_=NULL;
  float complex *c_C__=NULL;
  float f_ARBR=0.0,f_ARBI=0.0,f_AIBR=0.0,f_AIBI=0.0;
  int na=0;
  c_C__=NULL;
  if (c_C_p_!=NULL){
    if ( (*c_C_p_)==NULL ){ (*c_C_p_) = (float complex *) malloc1((unsigned long long int)n_row_A*(unsigned long long int)n_row_B*sizeof(float complex));}
    c_C__ = *c_C_p_;
    /* if (c_C_p_!=NULL){ } */}
  if (c_C__!=NULL){
    na=0;
    for (nrow_B=0;nrow_B<n_row_B;nrow_B++){
      for (nrow_A=0;nrow_A<n_row_A;nrow_A++){
	tmp_f_AR_point0_ = &(f_AR_trn__[nrow_A*n_col_X]);
	tmp_f_AI_point0_ = &(f_AI_trn__[nrow_A*n_col_X]);
	tmp_f_BR_point0_ = &(f_BR_trn__[nrow_B*n_col_X]);
	tmp_f_BI_point0_ = &(f_BI_trn__[nrow_B*n_col_X]);
	dp_ps_immintrin_loadu_fma(n_col_X,tmp_f_AR_point0_,tmp_f_BR_point0_,&f_ARBR);
	dp_ps_immintrin_loadu_fma(n_col_X,tmp_f_AR_point0_,tmp_f_BI_point0_,&f_ARBI);
	dp_ps_immintrin_loadu_fma(n_col_X,tmp_f_AI_point0_,tmp_f_BI_point0_,&f_AIBI);
	dp_ps_immintrin_loadu_fma(n_col_X,tmp_f_AI_point0_,tmp_f_BR_point0_,&f_AIBR);
	c_C__[na] = (float complex)(f_ARBR + f_AIBI) + _Complex_I * (float complex)(f_ARBI - f_AIBR);
	na += 1;
	/* for (nrow_A=0;nrow_A<n_row_A;nrow_A++){ } */}
      /* for (nrow_B=0;nrow_B<n_row_B;nrow_B++){ } */}
    /* if (c_C__!=NULL){ } */}
}

void hp_ps_mult_immintrin_test()
{
  int n_row_A = 153;
  int n_col_X = 1152;
  int n_col_X_rup = rup(n_col_X,8);
  int n_col_X_256 = n_col_X_rup/8; //%<-- 8 floats per __m256. ;
  int n_row_B = 147;
  int n_row_A_sub = 5;
  int n_col_X_sub = 15;
  int n_row_B_sub = 4;
  int nrow_A=0,ncol_X=0,nrow_B=0;
  unsigned long long int ulli_A_total = (unsigned long long int) n_row_A * (unsigned long long int) n_col_X_rup ;
  unsigned long long int ulli_B_total = (unsigned long long int) n_row_B * (unsigned long long int) n_col_X_rup ;
  unsigned long long int ulli_C_total = (unsigned long long int) n_row_A * (unsigned long long int) n_row_B ;
  unsigned long long int ulli=0;
  __m256 *ps_AR_trn__ = NULL;
  __m256 *ps_AI_trn__ = NULL;
  __m256 *ps_BR_trn__ = NULL;
  __m256 *ps_BI_trn__ = NULL;
  float *f_AR_trn__ = NULL;
  float *f_AI_trn__ = NULL;
  float *f_AR_u_trn__ = NULL;
  float *f_AI_u_trn__ = NULL;
  float *f_BR_trn__ = NULL;
  float *f_BI_trn__ = NULL;
  float *f_BR_u_trn__ = NULL;
  float *f_BI_u_trn__ = NULL;
  float *f_CR_bf__ = NULL;
  float *f_CI_bf__ = NULL;
  float *f_CR_ps__ = NULL;
  float *f_CI_ps__ = NULL;
  float *f_AR_sub__ = NULL;
  float *f_AI_sub__ = NULL;
  float *f_BR_sub__ = NULL;
  float *f_BI_sub__ = NULL;
  float *f_CR_sub__ = NULL;
  float *f_CI_sub__ = NULL;
  float complex *c_A_u_trn__ = NULL;
  float complex *c_B_u_trn__ = NULL;
  float complex *c_C_bf__ = NULL;
  float complex *c_C_al__ = NULL;
  float complex *c_A_sub__ = NULL;
  float complex *c_B_sub__ = NULL;
  float complex *c_C_sub__ = NULL;
  float ferror=0;
  ps_AR_trn__ = (__m256 *) _mm_malloc(n_col_X_256*n_row_A*sizeof(__m256),32);
  ps_AI_trn__ = (__m256 *) _mm_malloc(n_col_X_256*n_row_A*sizeof(__m256),32);
  ps_BR_trn__ = (__m256 *) _mm_malloc(n_col_X_256*n_row_A*sizeof(__m256),32);
  ps_BI_trn__ = (__m256 *) _mm_malloc(n_col_X_256*n_row_A*sizeof(__m256),32);
  f_AR_trn__ = (float *) ps_AR_trn__;
  f_AI_trn__ = (float *) ps_AI_trn__;
  f_BR_trn__ = (float *) ps_BR_trn__;
  f_BI_trn__ = (float *) ps_BI_trn__;
  f_AR_u_trn__ = (float *) malloc1(n_col_X*n_row_A*sizeof(float));
  f_AI_u_trn__ = (float *) malloc1(n_col_X*n_row_A*sizeof(float));
  f_BR_u_trn__ = (float *) malloc1(n_col_X*n_row_A*sizeof(float));
  f_BI_u_trn__ = (float *) malloc1(n_col_X*n_row_A*sizeof(float));
  f_CR_bf__ = (float *) malloc1(n_row_A*n_row_B*sizeof(float));
  f_CI_bf__ = (float *) malloc1(n_row_A*n_row_B*sizeof(float));
  f_CR_ps__ = (float *) malloc1(n_row_A*n_row_B*sizeof(float));
  f_CI_ps__ = (float *) malloc1(n_row_A*n_row_B*sizeof(float));
  f_AR_sub__ = (float *) malloc1(n_row_A_sub*n_col_X_sub*sizeof(float));
  f_AI_sub__ = (float *) malloc1(n_row_A_sub*n_col_X_sub*sizeof(float));
  f_BR_sub__ = (float *) malloc1(n_row_B_sub*n_col_X_sub*sizeof(float));
  f_BI_sub__ = (float *) malloc1(n_row_B_sub*n_col_X_sub*sizeof(float));
  f_CR_sub__ = (float *) malloc1(n_row_A_sub*n_row_B_sub*sizeof(float));
  f_CI_sub__ = (float *) malloc1(n_row_A_sub*n_row_B_sub*sizeof(float));
  c_A_u_trn__ = (float complex *) malloc1(n_col_X*n_row_A*sizeof(float complex));
  c_B_u_trn__ = (float complex *) malloc1(n_col_X*n_row_A*sizeof(float complex));
  c_C_bf__ = (float complex *) malloc1(n_row_A*n_row_B*sizeof(float complex));
  c_C_al__ = (float complex *) malloc1(n_row_A*n_row_B*sizeof(float complex));
  c_A_sub__ = (float complex *) malloc1(n_row_A_sub*n_col_X_sub*sizeof(float complex));
  c_B_sub__ = (float complex *) malloc1(n_row_B_sub*n_col_X_sub*sizeof(float complex));
  c_C_sub__ = (float complex *) malloc1(n_row_A_sub*n_row_B_sub*sizeof(float complex));
  GLOBAL_tic(0);
  ulli=0;
  for (nrow_A=0;nrow_A<n_row_A;nrow_A++){ for (ncol_X=0;ncol_X<n_col_X;ncol_X++){
      f_AR_trn__[ncol_X + nrow_A*n_col_X_rup] = (float)(((int)ulli%7)-3);
      f_AI_trn__[ncol_X + nrow_A*n_col_X_rup] = (float)(((int)ulli%7)-2);
      f_AR_u_trn__[ncol_X + nrow_A*n_col_X] = (float)(((int)ulli%7)-3);
      f_AI_u_trn__[ncol_X + nrow_A*n_col_X] = (float)(((int)ulli%7)-2);
      c_A_u_trn__[ncol_X + nrow_A*n_col_X] = (float complex) f_AR_u_trn__[ncol_X + nrow_A*n_col_X] + _Complex_I * (float complex) f_AI_u_trn__[ncol_X + nrow_A*n_col_X];
      ulli++;
      /* for (nrow_A=0;nrow_A<n_row_A;nrow_A++){ for (ncol_X=0;ncol_X<n_col_X;ncol_X++){ }} */}}  
  for (nrow_A=0;nrow_A<n_row_A_sub;nrow_A++){ for (ncol_X=0;ncol_X<n_col_X_sub;ncol_X++){ 
      f_AR_sub__[nrow_A+ncol_X*n_row_A_sub] = f_AR_trn__[ncol_X + nrow_A*n_col_X_rup];
      f_AI_sub__[nrow_A+ncol_X*n_row_A_sub] = f_AI_trn__[ncol_X + nrow_A*n_col_X_rup];
      c_A_sub__[nrow_A+ncol_X*n_row_A_sub] = (float complex) f_AR_sub__[nrow_A+ncol_X*n_row_A_sub] + _Complex_I * (float complex) f_AI_sub__[nrow_A+ncol_X*n_row_A_sub];
      /* for (nrow_A=0;nrow_A<n_row_A_sub;nrow_A++){ for (ncol_X=0;ncol_X<n_col_X_sub;ncol_X++){ }} */}}
  printf(" %% upper corner of f_AR_trn__: \n");
  array_printf(f_AR_sub__,"float",n_row_A_sub,n_col_X_sub," % f_AR_sub__: ");
  printf(" %% upper corner of f_AI_trn__: \n");
  array_printf(f_AI_sub__,"float",n_row_A_sub,n_col_X_sub," % f_AI_sub__: ");
  printf(" %% upper corner of c_A_trn__: \n");
  array_printf(c_A_sub__,"float complex",n_row_A_sub,n_col_X_sub," % c_A_sub__: ");
  ulli=0;
  for (nrow_B=0;nrow_B<n_row_B;nrow_B++){ for (ncol_X=0;ncol_X<n_col_X;ncol_X++){
      f_BR_trn__[ncol_X + nrow_B*n_col_X_rup] = (float)(((int)ulli%11)-2);
      f_BI_trn__[ncol_X + nrow_B*n_col_X_rup] = (float)(((int)ulli%11)-1);
      f_BR_u_trn__[ncol_X + nrow_B*n_col_X] = (float)(((int)ulli%11)-2);
      f_BI_u_trn__[ncol_X + nrow_B*n_col_X] = (float)(((int)ulli%11)-1);
      c_B_u_trn__[ncol_X + nrow_B*n_col_X] = (float complex) f_BR_u_trn__[ncol_X + nrow_B*n_col_X] + _Complex_I * (float complex) f_BI_u_trn__[ncol_X + nrow_B*n_col_X];
      ulli++;
      /* for (nrow_B=0;nrow_B<n_row_B;nrow_B++){ for (ncol_X=0;ncol_X<n_col_X;ncol_X++){ }} */}}
  for (nrow_B=0;nrow_B<n_row_B_sub;nrow_B++){ for (ncol_X=0;ncol_X<n_col_X_sub;ncol_X++){
      f_BR_sub__[nrow_B+ncol_X*n_row_B_sub] = f_BR_trn__[ncol_X + nrow_B*n_col_X_rup];
      f_BI_sub__[nrow_B+ncol_X*n_row_B_sub] = f_BI_trn__[ncol_X + nrow_B*n_col_X_rup];
      c_B_sub__[nrow_B+ncol_X*n_row_B_sub] = (float complex) f_BR_sub__[nrow_B+ncol_X*n_row_B_sub] + _Complex_I * (float complex) f_BI_sub__[nrow_B+ncol_X*n_row_B_sub];
      /* for (nrow_B=0;nrow_B<n_row_B_sub;nrow_B++){ for (ncol_X=0;ncol_X<n_col_X_sub;ncol_X++){ }} */}}
  printf(" %% upper corner of f_BR_trn__: \n");
  array_printf(f_BR_sub__,"float",n_row_B_sub,n_col_X_sub," % f_BR_sub__: ");
  printf(" %% upper corner of f_BI_trn__: \n");
  array_printf(f_BI_sub__,"float",n_row_B_sub,n_col_X_sub," % f_BI_sub__: ");
  printf(" %% upper corner of c_B_trn__: \n");
  array_printf(c_B_sub__,"float complex",n_row_B_sub,n_col_X_sub," % c_B_sub__: ");
  GLOBAL_toc(0,1," initialize: ");
  /* %%%%%%%% */
  GLOBAL_tic(0);
  memset(c_C_bf__,0,ulli_C_total*sizeof(float complex));
  hp_interleave_mult_bruteforce(n_row_A,n_col_X,c_A_u_trn__,n_row_B,c_B_u_trn__,&c_C_bf__);
  GLOBAL_toc(0,1," hp_ps_mult_bruteforce: ");
  printf("Gops %0.6f\n",(double)n_row_A*(double)n_row_B*(double)n_col_X_rup/GLOBAL_elrt[0]/1e9);
  for (nrow_A=0;nrow_A<n_row_A_sub;nrow_A++){ for (nrow_B=0;nrow_B<n_row_B_sub;nrow_B++){
      c_C_sub__[nrow_A+nrow_B*n_row_A_sub] = c_C_bf__[nrow_A+nrow_B*n_row_A];
      /* for (nrow_A=0;nrow_A<n_row_A_sub;nrow_A++){ for (nrow_B=0;nrow_B<n_row_B_sub;nrow_B++){ }} */}}
  printf(" %% upper corner of c_C_bf__: \n");
  array_printf(c_C_sub__,"float complex",n_row_A_sub,n_row_B_sub," % c_C_bf__: ");
  ferror = cfnorm(ulli_C_total,c_C_bf__,c_C_bf__);
  printf(" %% ferror %0.16f\n",ferror);
  /* %%%%%%%% */
  GLOBAL_tic(0);
  memset(c_C_al__,0,ulli_C_total*sizeof(float complex));
  hp_segregated_mult_bruteforce(n_row_A,n_col_X,f_AR_u_trn__,f_AI_u_trn__,n_row_B,f_BR_u_trn__,f_BI_u_trn__,&c_C_al__);
  GLOBAL_toc(0,1," hp_segregated_mult_bruteforce: ");
  printf("Gops %0.6f\n",(double)n_row_A*(double)n_row_B*(double)n_col_X_rup/GLOBAL_elrt[0]/1e9);
  for (nrow_A=0;nrow_A<n_row_A_sub;nrow_A++){ for (nrow_B=0;nrow_B<n_row_B_sub;nrow_B++){
      c_C_sub__[nrow_A+nrow_B*n_row_A_sub] = c_C_al__[nrow_A+nrow_B*n_row_A];
      /* for (nrow_A=0;nrow_A<n_row_A_sub;nrow_A++){ for (nrow_B=0;nrow_B<n_row_B_sub;nrow_B++){ }} */}}
  printf(" %% upper corner of c_C_al__: \n");
  array_printf(c_C_sub__,"float complex",n_row_A_sub,n_row_B_sub," % c_C_al__: ");
  ferror = cfnorm(ulli_C_total,c_C_bf__,c_C_al__);
  printf(" %% ferror %0.16f\n",ferror);
  /* %%%%%%%% */
  GLOBAL_tic(0);
  memset(c_C_al__,0,ulli_C_total*sizeof(float complex));
  hp_segregated_mult_immintrin_loadu_fma(n_row_A,n_col_X,f_AR_u_trn__,f_AI_u_trn__,n_row_B,f_BR_u_trn__,f_BI_u_trn__,&c_C_al__);
  GLOBAL_toc(0,1," hp_segregated_mult_immintrin_loadu_fma: ");
  printf("Gops %0.6f\n",(double)n_row_A*(double)n_row_B*(double)n_col_X_rup/GLOBAL_elrt[0]/1e9);
  for (nrow_A=0;nrow_A<n_row_A_sub;nrow_A++){ for (nrow_B=0;nrow_B<n_row_B_sub;nrow_B++){
      c_C_sub__[nrow_A+nrow_B*n_row_A_sub] = c_C_al__[nrow_A+nrow_B*n_row_A];
      /* for (nrow_A=0;nrow_A<n_row_A_sub;nrow_A++){ for (nrow_B=0;nrow_B<n_row_B_sub;nrow_B++){ }} */}}
  printf(" %% upper corner of c_C_al__: \n");
  array_printf(c_C_sub__,"float complex",n_row_A_sub,n_row_B_sub," % c_C_al__: ");
  ferror = cfnorm(ulli_C_total,c_C_bf__,c_C_al__);
  printf(" %% ferror %0.16f\n",ferror);
  /* %%%%%%%% */
#ifdef _CBLAS
  GLOBAL_tic(0);
  memset(c_C_al__,0,ulli_C_total*sizeof(float complex));
  hp_interleave_mult_cblas_cgemm(n_row_A,n_col_X,c_A_u_trn__,n_row_B,c_B_u_trn__,&c_C_al__);
  GLOBAL_toc(0,1," hp_interleave_mult_cblas_cgemm: ");
  printf("Gops %0.6f\n",(double)n_row_A*(double)n_row_B*(double)n_col_X_rup/GLOBAL_elrt[0]/1e9);
  for (nrow_A=0;nrow_A<n_row_A_sub;nrow_A++){ for (nrow_B=0;nrow_B<n_row_B_sub;nrow_B++){
      c_C_sub__[nrow_A+nrow_B*n_row_A_sub] = c_C_al__[nrow_A+nrow_B*n_row_A];
      /* for (nrow_A=0;nrow_A<n_row_A_sub;nrow_A++){ for (nrow_B=0;nrow_B<n_row_B_sub;nrow_B++){ }} */}}
  printf(" %% upper corner of c_C_al__: \n");
  array_printf(c_C_sub__,"float complex",n_row_A_sub,n_row_B_sub," % c_C_al__: ");
  ferror = cfnorm(ulli_C_total,c_C_bf__,c_C_al__);
  printf(" %% ferror %0.16f\n",ferror);
#endif /* _CBLAS */
  /* %%%%%%%% */
  _mm_free(ps_AR_trn__); ps_AR_trn__ = NULL;
  _mm_free(ps_AI_trn__); ps_AI_trn__ = NULL;
  _mm_free(ps_BR_trn__); ps_BR_trn__ = NULL;
  _mm_free(ps_BI_trn__); ps_BI_trn__ = NULL;
  free1(&f_AR_u_trn__);
  free1(&f_AI_u_trn__);
  free1(&f_BR_u_trn__);
  free1(&f_BI_u_trn__);
  free1(&f_CR_bf__);
  free1(&f_CI_bf__);
  free1(&f_CR_ps__);
  free1(&f_CI_ps__);
  free1(&f_AR_sub__);
  free1(&f_AI_sub__);
  free1(&f_BR_sub__);
  free1(&f_BI_sub__);
  free1(&f_CR_sub__);
  free1(&f_CI_sub__);
  free1(&c_A_u_trn__);
  free1(&c_B_u_trn__);
  free1(&c_C_bf__);
  free1(&c_C_al__);
  free1(&c_A_sub__);
  free1(&c_B_sub__);
  free1(&c_C_sub__);
  //wkspace_printf();
}


