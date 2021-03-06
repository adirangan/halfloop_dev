#ifndef _MONOLITH
#include "halfloop_header.h"
#endif /* _MONOLITH */

void halfloop_nonbinary_f_recursive_helper_QR_helper_orth
(
  int verbose
 ,int n_r_index
 ,int n_c_index
 ,float *E_rc__
 ,float *QE_rc__
 ,unsigned long long int *rseed_p
)
{
  int niteration=0,n_iteration=0;
  if (verbose){ printf(" %% [entering halfloop_nonbinary_f_recursive_helper_QR_helper_orth]\n");}
  float *Q_rc__=NULL;
  float *QE_0_rc__=NULL;
  float *QE_1_rc__=NULL;
  int *r_perm_index_=NULL;
  int *r_perm_index_thrice_=NULL;
  int *r_perm_index_sub_=NULL;
  int n_m_index = 24; //%<-- 8 floats per __m256. ;
  int n_m_index_2 = n_m_index/2;
  int n_loop=0,nloop=0;
  if ( (GLOBAL_flag_orth_brute==1) || (n_r_index<=n_m_index) ){
    Q_rc__ = (float *) malloc1((unsigned long long int)n_r_index*(unsigned long long int)n_r_index*sizeof(float));
    array_orth_f(n_r_index,n_r_index,&Q_rc__,rseed_p);
    if (verbose>9){ array_printf_margin(Q_rc__,"float",n_r_index,n_r_index," % Q_rc__: "); printf(" %% %% %% %%\n");}
    dp_ps_mult_immintrin_loadu_wrap(n_r_index,n_r_index,Q_rc__,n_c_index,E_rc__,&QE_rc__);
    free1(&Q_rc__);
    /* if (n_r_index<=n_m_index){ } */}
  if ( (GLOBAL_flag_orth_brute==0) && (n_r_index> n_m_index) ){
    Q_rc__ = (float *) malloc1((unsigned long long int)n_m_index*(unsigned long long int)n_m_index*sizeof(float));
    QE_0_rc__ = (float *) malloc1((unsigned long long int)n_m_index*(unsigned long long int)n_c_index*sizeof(float));
    QE_1_rc__ = (float *) malloc1((unsigned long long int)n_m_index*(unsigned long long int)n_c_index*sizeof(float));
    memcpy(QE_rc__,E_rc__,(unsigned long long int)n_r_index*(unsigned long long int)n_c_index*sizeof(float));
    n_loop = maximum(1,ceil(log(n_r_index)/maximum(1,log(n_m_index))));
    for (nloop=0;nloop<n_loop;nloop++){
      irandperm(n_r_index,&r_perm_index_,rseed_p);
      r_perm_index_thrice_ = (int *) malloc1((unsigned long long int)n_r_index*(unsigned long long int)3*sizeof(int));
      memcpy(r_perm_index_thrice_+0*n_r_index,r_perm_index_,(unsigned long long int)n_r_index*sizeof(int));
      memcpy(r_perm_index_thrice_+1*n_r_index,r_perm_index_,(unsigned long long int)n_r_index*sizeof(int));
      memcpy(r_perm_index_thrice_+2*n_r_index,r_perm_index_,(unsigned long long int)n_r_index*sizeof(int));
      n_iteration = 2*ceil((double)n_r_index/(double)n_m_index);
      for (niteration=0;niteration<n_iteration;niteration++){
	if (verbose>9){ printf(" %% %d/%d\n",niteration,n_iteration);}
	r_perm_index_sub_ = r_perm_index_thrice_ + (int)floor(niteration*n_m_index_2);
	if (verbose>9){ array_printf_margin(r_perm_index_sub_,"int",1,n_m_index," %% r_perm_index_sub_: ");}
	array_extract_f_from_f(n_r_index,n_c_index,QE_rc__,n_m_index,r_perm_index_sub_,0,NULL,&QE_0_rc__,NULL);
	if (verbose>9){ array_printf_margin(QE_0_rc__,"float",n_m_index,n_c_index," %% QE_0_rc__: ");}
	array_orth_f(n_m_index,n_m_index,&Q_rc__,rseed_p);
	if (verbose>9){ array_printf_margin(Q_rc__,"float",n_m_index,n_m_index," %% Q_rc__: ");}
	dp_ps_mult_immintrin_loadu_wrap(n_m_index,n_m_index,Q_rc__,n_c_index,QE_0_rc__,&QE_1_rc__);
	array_implant_f_from_f(n_r_index,n_c_index,QE_rc__,n_m_index,r_perm_index_sub_,0,NULL,QE_1_rc__,NULL);
	/* for (niteration=0;niteration<n_iteration;niteration++){ } */}
      free1(&r_perm_index_thrice_);
      free1(&r_perm_index_);
      /* for (nloop=0;nloop<n_loop;nloop++){ } */}
    free1(&Q_rc__);
    free1(&QE_0_rc__);
    free1(&QE_1_rc__);
    /* if (n_r_index> n_m_index){ } */}
  if (verbose){ printf(" %% [finished halfloop_nonbinary_f_recursive_helper_QR_helper_orth]\n");}
}

void halfloop_nonbinary_f_recursive_helper_QR_helper_orth_test_speed()
{
  int verbose=0;
  int nr=0,n_r = 6000;
  int nc=0,n_c = 8000;
  float *E_rc__=NULL;
  float *QE_rc__=NULL;
  unsigned long long int ulli=0;
  float x=0,y=0,z=0;
  unsigned long long int rseed=0;
  E_rc__ = (float *) malloc1((unsigned long long int)n_r*(unsigned long long int)n_c*sizeof(float));
  rseed=1; RSEED_adv8(&rseed);
  ulli=0;
  for (nc=0;nc<n_c;nc++){ for (nr=0;nr<n_r;nr++){ 
      //E_rc__[ulli] = RNGET(&rseed);
      //if ( (nc<floor(sqrt((double)n_c))) && (nr<floor(sqrt((double)n_r))) ){ E_rc__[ulli] += 0.5; }
      x = 2.0*(float)nr/(float)(n_r-1) - 1.0;
      y = 2.0*(float)nc/(float)(n_c-1) - 1.0;
      z = (x+y)/4.0;
      E_rc__[ulli] = sin(2*PI_LF*x) + cos(2*PI_LF*2*y) + x*x + y*y*y + cos(2*PI_LF*4*z);
      ulli++;
      /* for (nc=0;nc<n_c;nc++){ for (nr=0;nr<n_r;nr++){ }} */}}
  array_printf_margin(E_rc__,"float",n_r,n_c," %  E_r__: ");
  QE_rc__ = (float *) malloc1((unsigned long long int)n_r*(unsigned long long int)n_c*sizeof(float));
  /* %%%%%%%% */
  GLOBAL_tic(0);
  halfloop_nonbinary_f_recursive_helper_QR_helper_orth(verbose,n_r,n_c,E_rc__,QE_rc__,&rseed);
  GLOBAL_toc(0,1," % halfloop_nonbinary_f_recursive: ");
  printf("Effective Multiplication Gops %0.6f\n",(double)n_r*(double)n_r*(double)n_c/GLOBAL_elrt[0]/1e9);
  printf("Effective Output Gops %0.6f\n",(double)n_r*(double)n_c/GLOBAL_elrt[0]/1e9);
  /* %%%%%%%% */
  array_printf_margin(QE_rc__,"float",n_r,n_c," % QE_r__: ");
  free1(&E_rc__);
  free1(&QE_rc__);
}

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
)
{
  struct stat stat_file = {0};
  int flag_not_exist=0;
  int flag_rcdrop = 0;
  int flag_r0drop = 1;
  int trace_QR_index = 3;
  int n_iteration=0,tmp_n_iteration=0;
  int nshuffle=0;
  unsigned long int rseed=0;
  double *trace__=NULL;
  double *trace_shuffle__=NULL;
  double *QR_=NULL;
  double *QR__=NULL;
  int n_xdrop=0;
  int *xdrop__=NULL;
  int *xdrop_shuffle__=NULL;
  float *E_cr__=NULL;
  float *E_rc__=NULL;
  float *QE_rc__=NULL;
  int MDA_n_dim=0;
  int *MDA_dim_=NULL;
  char MDA_fname[PNAMESIZE];
  if (verbose){ printf(" %% [entering halfloop_nonbinary_f_recursive_helper_QR__]\n");}
  flag_not_exist = (stat(fname_trace__,&stat_file)==-1) || (stat(fname_QR__,&stat_file)==-1) || (stat(fname_xdrop__,&stat_file)==-1);
  if ( flag_force_create || flag_not_exist ){
    if (verbose){ printf(" %% %s not found, creating\n",fname_trace__);}
    MDA_dim_ = (int *) malloc1(2*sizeof(int));
    n_iteration = get_xdrop_logscale_length(n_r_index,flag_r0drop_vs_rcdrop==flag_rcdrop ? n_c_index : maximum(n_r_index,n_c_index),gamma);
    trace__ = (double *) malloc1((unsigned long long int)6*(unsigned long long int)n_iteration*sizeof(double));
    trace_shuffle__ = (double *) malloc1((unsigned long long int)6*(unsigned long long int)n_iteration*sizeof(double));
    QR__ = (double *) malloc1((unsigned long long int)n_iteration*(unsigned long long int)(1+n_shuffle)*sizeof(double));
    if (flag_r0drop_vs_rcdrop==flag_rcdrop){ n_xdrop = n_r_index + n_c_index;}
    if (flag_r0drop_vs_rcdrop==flag_r0drop){ n_xdrop = n_r_index;}
    xdrop__ = (int *) malloc1((unsigned long long int)2*(unsigned long long int)n_xdrop*sizeof(double));
    xdrop_shuffle__ = (int *) malloc1((unsigned long long int)2*(unsigned long long int)n_xdrop*sizeof(double));
    E_rc__ = (float *) malloc1((unsigned long long int)n_r_index*(unsigned long long int)n_c_index*sizeof(float));
    E_cr__ = (float *) malloc1((unsigned long long int)n_c_index*(unsigned long long int)n_r_index*sizeof(float));
    array_extract_f_from_f(n_r,n_c,E_base_rc__,n_r_index,r_index_,n_c_index,c_index_,&E_rc__,NULL);
    if (verbose>9){ array_printf_margin(E_rc__,"float",n_r_index,n_c_index," % E_rc__: "); printf(" %% %% %% %%\n");}
    array_mean_center_row(n_r_index,n_c_index,E_rc__,NULL,"float",&E_rc__,&E_cr__);
    if (verbose>9){ array_printf_margin(E_rc__,"float",n_r_index,n_c_index," % E_rc__: "); printf(" %% %% %% %%\n");}
    if (verbose>9){ array_printf_margin(E_cr__,"float",n_c_index,n_r_index," % E_cr__: "); printf(" %% %% %% %%\n");}
    if (flag_r0drop_vs_rcdrop==flag_rcdrop){ halfloop_nonbinary_f(n_r_index,n_c_index,E_rc__,E_cr__,gamma,&xdrop__,&tmp_n_iteration,&trace__);}
    if (flag_r0drop_vs_rcdrop==flag_r0drop){ halfloop_nonbinary_rdrop_f(n_r_index,n_c_index,E_rc__,E_cr__,gamma,&xdrop__,&tmp_n_iteration,&trace__);}
    if (verbose>9){ array_printf_margin(xdrop__,"int",2,n_xdrop," % xdrop__: "); printf(" %% %% %% %%\n");}
    if (verbose>9){ array_printf_margin(trace__,"double",6,n_iteration," % trace__: "); printf(" %% %% %% %%\n");}
    if (tmp_n_iteration!=n_iteration){ printf(" %% Warning, tmp_n_iteration %d not equal n_iteration %d\n",tmp_n_iteration,n_iteration);}
    sprintf(MDA_fname,"%s",fname_trace__); MDA_n_dim = 2; MDA_dim_[0] = 6; MDA_dim_[1] = n_iteration;
    MDA_write_r8(MDA_n_dim,MDA_dim_,trace__,MDA_fname);
    if (verbose>9){ MDA_printf_r8_margin(MDA_fname);}
    sprintf(MDA_fname,"%s",fname_xdrop__); MDA_n_dim = 2; MDA_dim_[0] = 2; MDA_dim_[1] = n_xdrop;
    MDA_write_i4(MDA_n_dim,MDA_dim_,xdrop__,MDA_fname);
    if (verbose>9){ MDA_printf_i4_margin(MDA_fname);}
    QR_ = QR__ + (unsigned long long int)0*(unsigned long long int)n_iteration;
    array_extract_d_from_d(6,n_iteration,trace__,1,&trace_QR_index,0,NULL,&QR_,NULL);
    if (verbose>9){ array_printf_margin(QR_,"double",1,n_iteration," % QR_: ");}
    QE_rc__ = (float *) malloc1((unsigned long long int)n_r_index*(unsigned long long int)n_c_index*sizeof(float));
    for (nshuffle=0;nshuffle<n_shuffle;nshuffle++){
      if (verbose>9){ printf(" %% nshuffle %d/%d\n",nshuffle,n_shuffle);}
      rseed = (1+nshuffle); RSEED_adv8(&rseed);
      array_extract_f_from_f(n_r,n_c,E_base_rc__,n_r_index,r_index_,n_c_index,c_index_,&E_rc__,NULL);
      array_mean_center_row(n_r_index,n_c_index,E_rc__,NULL,"float",&E_rc__,NULL);
      if (verbose>9){ array_printf_margin(E_rc__,"float",n_r_index,n_c_index," % E_rc__: "); printf(" %% %% %% %%\n");}
      halfloop_nonbinary_f_recursive_helper_QR_helper_orth(verbose,n_r_index,n_c_index,E_rc__,QE_rc__,&rseed);
      if (verbose>9){ array_printf_margin(QE_rc__,"float",n_r_index,n_c_index," % QE_rc__: "); printf(" %% %% %% %%\n");}
      array_mean_center_row(n_r_index,n_c_index,QE_rc__,NULL,"float",&E_rc__,&E_cr__);
      if (verbose>9){ array_printf_margin(E_rc__,"float",n_r_index,n_c_index," % E_rc__: "); printf(" %% %% %% %%\n");}
      if (flag_r0drop_vs_rcdrop==flag_rcdrop){ halfloop_nonbinary_f(n_r_index,n_c_index,E_rc__,E_cr__,gamma,&xdrop_shuffle__,&tmp_n_iteration,&trace_shuffle__);}
      if (flag_r0drop_vs_rcdrop==flag_r0drop){ halfloop_nonbinary_rdrop_f(n_r_index,n_c_index,E_rc__,E_cr__,gamma,&xdrop_shuffle__,&tmp_n_iteration,&trace_shuffle__);}
      if (verbose>9){ array_printf_margin(xdrop_shuffle__,"int",2,n_xdrop," % xdrop_shuffle__: "); printf(" %% %% %% %%\n");}
      if (verbose>9){ array_printf_margin(trace_shuffle__,"double",6,n_iteration," % trace_shuffle__: "); printf(" %% %% %% %%\n");}
      if (tmp_n_iteration!=n_iteration){ printf(" %% Warning, tmp_n_iteration %d not equal n_iteration %d\n",tmp_n_iteration,n_iteration);}
      QR_ = QR__ + (unsigned long long int)(1+nshuffle)*(unsigned long long int)n_iteration;
      array_extract_d_from_d(6,n_iteration,trace_shuffle__,1,&trace_QR_index,0,NULL,&QR_,NULL);
      if (verbose>9){ array_printf_margin(QR_,"double",1,n_iteration," % QR_: ");}
      /* for (nshuffle=0;nshuffle<n_shuffle;nshuffle++){ } */}
    if (verbose>9){ array_printf_margin(QR__,"double",n_iteration,(1+n_shuffle)," % QR__: ");}
    sprintf(MDA_fname,"%s",fname_QR__); MDA_n_dim = 2; MDA_dim_[0] = n_iteration; MDA_dim_[1] = (1+n_shuffle);
    MDA_write_r8(MDA_n_dim,MDA_dim_,QR__,MDA_fname);
    if (verbose>9){ MDA_printf_r8_margin(MDA_fname);}
    free1(&QE_rc__);
    free1(&E_cr__);
    free1(&E_rc__);
    free1(&trace__);
    free1(&trace_shuffle__);
    free1(&QR__);
    free1(&xdrop__);
    free1(&xdrop_shuffle__);
    free1(&MDA_dim_);
    /* not found */}
  else{ if (verbose){ printf(" %% %s found, not creating\n",fname_trace__);}}
  if (verbose){ printf(" %% [finished halfloop_nonbinary_f_recursive_helper_QR__]\n");}
}

/* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */

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
)
{
  int flag_rcdrop = 0;
  int flag_r0drop = 1;
  int trace_QR_index = 3;
  int tmp_n_iteration=0;
  double *trace_shuffle__=NULL;
  int *xdrop_shuffle__=NULL;
  unsigned long int rseed=0;
  float *tmp_E_cr__=NULL;
  float *tmp_E_rc__=NULL;
  float *tmp_QE_rc__=NULL;
  if (verbose){ printf(" %% [entering halfloop_nonbinary_f_recursive_omp_helper_QR_]\n");}
  if (nshuffle==0){
    if (flag_r0drop_vs_rcdrop==flag_rcdrop){ halfloop_nonbinary_f(n_r_index,n_c_index,E_rc__,E_cr__,gamma,&xdrop__,&tmp_n_iteration,&trace__);}
    if (flag_r0drop_vs_rcdrop==flag_r0drop){ halfloop_nonbinary_rdrop_f(n_r_index,n_c_index,E_rc__,E_cr__,gamma,&xdrop__,&tmp_n_iteration,&trace__);}
    if (tmp_n_iteration!=n_iteration){ printf(" %% Warning, tmp_n_iteration %d not equal n_iteration %d\n",tmp_n_iteration,n_iteration);}
    array_extract_d_from_d(6,n_iteration,trace__,1,&trace_QR_index,0,NULL,&QR_,NULL);
    /* if (nshuffle==0){ } */}
  if (nshuffle> 0){
    trace_shuffle__ = (double *) malloc1((unsigned long long int)6*(unsigned long long int)n_iteration*sizeof(double));
    xdrop_shuffle__ = (int *) malloc1((unsigned long long int)2*(unsigned long long int)n_xdrop*sizeof(double));
    tmp_E_rc__ = (float *) malloc1((unsigned long long int)n_r_index*(unsigned long long int)n_c_index*sizeof(float));
    tmp_E_cr__ = (float *) malloc1((unsigned long long int)n_c_index*(unsigned long long int)n_r_index*sizeof(float));
    tmp_QE_rc__ = (float *) malloc1((unsigned long long int)n_r_index*(unsigned long long int)n_c_index*sizeof(float));
    rseed = (0+nshuffle); RSEED_adv8(&rseed);
    halfloop_nonbinary_f_recursive_helper_QR_helper_orth(verbose,n_r_index,n_c_index,E_rc__,tmp_QE_rc__,&rseed);
    array_mean_center_row(n_r_index,n_c_index,tmp_QE_rc__,NULL,"float",&tmp_E_rc__,&tmp_E_cr__);
    if (flag_r0drop_vs_rcdrop==flag_rcdrop){ halfloop_nonbinary_f(n_r_index,n_c_index,tmp_E_rc__,tmp_E_cr__,gamma,&xdrop_shuffle__,&tmp_n_iteration,&trace_shuffle__);}
    if (flag_r0drop_vs_rcdrop==flag_r0drop){ halfloop_nonbinary_rdrop_f(n_r_index,n_c_index,tmp_E_rc__,tmp_E_cr__,gamma,&xdrop_shuffle__,&tmp_n_iteration,&trace_shuffle__);}
    array_extract_d_from_d(6,n_iteration,trace_shuffle__,1,&trace_QR_index,0,NULL,&QR_,NULL);
    free1(&trace_shuffle__);
    free1(&xdrop_shuffle__);
    free1(&tmp_E_rc__);
    free1(&tmp_E_cr__);
    free1(&tmp_QE_rc__);
    /* if (nshuffle> 0){ } */}
  if (verbose){ printf(" %% [finished halfloop_nonbinary_f_recursive_omp_helper_QR_]\n");}
}
 
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
)
{
  struct stat stat_file = {0};
  int flag_not_exist=0;
  int flag_rcdrop = 0;
  int flag_r0drop = 1;
  int trace_QR_index = 3;
  int n_iteration=0,tmp_n_iteration=0;
  int nshuffle=0;
  unsigned long int rseed=0;
  double *trace__=NULL;
  double *QR_=NULL;
  double *QR__=NULL;
  int n_xdrop=0;
  int *xdrop__=NULL;
  float *E_cr__=NULL;
  float *E_rc__=NULL;
  int MDA_n_dim=0;
  int *MDA_dim_=NULL;
  char MDA_fname[PNAMESIZE];
  if (verbose){ printf(" %% [entering halfloop_nonbinary_f_recursive_omp_helper_QR__]\n");}
  flag_not_exist = (stat(fname_trace__,&stat_file)==-1) || (stat(fname_QR__,&stat_file)==-1) || (stat(fname_xdrop__,&stat_file)==-1);
  if ( flag_force_create || flag_not_exist ){
    if (verbose){ printf(" %% %s not found, creating\n",fname_trace__);}
    MDA_dim_ = (int *) malloc1(2*sizeof(int));
    n_iteration = get_xdrop_logscale_length(n_r_index,flag_r0drop_vs_rcdrop==flag_rcdrop ? n_c_index : maximum(n_r_index,n_c_index),gamma);
    trace__ = (double *) malloc1((unsigned long long int)6*(unsigned long long int)n_iteration*sizeof(double));
    QR__ = (double *) malloc1((unsigned long long int)n_iteration*(unsigned long long int)(1+n_shuffle)*sizeof(double));
    if (flag_r0drop_vs_rcdrop==flag_rcdrop){ n_xdrop = n_r_index + n_c_index;}
    if (flag_r0drop_vs_rcdrop==flag_r0drop){ n_xdrop = n_r_index;}
    xdrop__ = (int *) malloc1((unsigned long long int)2*(unsigned long long int)n_xdrop*sizeof(double));
    E_rc__ = (float *) malloc1((unsigned long long int)n_r_index*(unsigned long long int)n_c_index*sizeof(float));
    E_cr__ = (float *) malloc1((unsigned long long int)n_c_index*(unsigned long long int)n_r_index*sizeof(float));
    array_extract_f_from_f(n_r,n_c,E_base_rc__,n_r_index,r_index_,n_c_index,c_index_,&E_rc__,NULL);
    if (verbose>9){ array_printf_margin(E_rc__,"float",n_r_index,n_c_index," % E_rc__: "); printf(" %% %% %% %%\n");}
    array_mean_center_row(n_r_index,n_c_index,E_rc__,NULL,"float",&E_rc__,&E_cr__);
    if (verbose>9){ array_printf_margin(E_rc__,"float",n_r_index,n_c_index," % E_rc__: "); printf(" %% %% %% %%\n");}
    if (verbose>9){ array_printf_margin(E_cr__,"float",n_c_index,n_r_index," % E_cr__: "); printf(" %% %% %% %%\n");}
    GLOBAL_malloc1_notupdate=1;
#pragma omp parallel private(nshuffle,QR_)
    { /* begin omp parallel */
      nshuffle=0; QR_ = NULL;
#pragma omp for schedule(dynamic)
      for (nshuffle=0;nshuffle<=n_shuffle;nshuffle++){
	QR_ = QR__ + (unsigned long long int)nshuffle*(unsigned long long int)n_iteration;
	halfloop_nonbinary_f_recursive_omp_helper_QR_
	  (
	   verbose
	   ,n_r_index
	   ,n_c_index
	   ,E_rc__
	   ,E_cr__
	   ,flag_r0drop_vs_rcdrop
	   ,gamma
	   ,nshuffle
	   ,n_iteration
	   ,trace__
	   ,QR_
	   ,n_xdrop
	   ,xdrop__
	   );
	/* for (nshuffle=0;nshuffle<=n_shuffle;nshuffle++){ } */}
      /* end omp parallel */}
    GLOBAL_malloc1_notupdate=0;
    sprintf(MDA_fname,"%s",fname_trace__); MDA_n_dim = 2; MDA_dim_[0] = 6; MDA_dim_[1] = n_iteration;
    MDA_write_r8(MDA_n_dim,MDA_dim_,trace__,MDA_fname);
    if (verbose>9){ MDA_printf_r8_margin(MDA_fname);}
    sprintf(MDA_fname,"%s",fname_xdrop__); MDA_n_dim = 2; MDA_dim_[0] = 2; MDA_dim_[1] = n_xdrop;
    MDA_write_i4(MDA_n_dim,MDA_dim_,xdrop__,MDA_fname);
    if (verbose>9){ MDA_printf_i4_margin(MDA_fname);}
    sprintf(MDA_fname,"%s",fname_QR__); MDA_n_dim = 2; MDA_dim_[0] = n_iteration; MDA_dim_[1] = (1+n_shuffle);
    MDA_write_r8(MDA_n_dim,MDA_dim_,QR__,MDA_fname);
    if (verbose>9){ MDA_printf_r8_margin(MDA_fname);}
    free1(&E_cr__);
    free1(&E_rc__);
    free1(&trace__);
    free1(&QR__);
    free1(&xdrop__);
    free1(&MDA_dim_);
    /* not found */}
  else{ if (verbose){ printf(" %% %s found, not creating\n",fname_trace__);}}
  if (verbose){ printf(" %% [finished halfloop_nonbinary_f_recursive_omp_helper_QR__]\n");}
}

/* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */

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
)
{
  int flag_rcdrop = 0;
  int flag_r0drop = 1;
  int trace_r_rtn_index = 1;
  int trace_c_rtn_index = 2;
  int MDA_n_dim=0;
  int *MDA_dim_=NULL;
  int niteration=0,n_iteration=0;
  double *trace__=NULL;
  int nshuffle=0,n_shuffle=0;
  double QR_avg=0,QR_std=0,*QR_=NULL,*QR__=NULL;
  int n_xdrop=0,*xdrop_=NULL,*xdrop__=NULL;
  double ZR=0,nlp_ZR=0,*nlp_ZR_=NULL;
  int nsub=0,n_sub=0,*index_sub_=NULL;
  double *nlp_ZR_sub_=NULL;
  double nlp_ZR_max=0;
  int nlp_ZR_index;
  int r_rtn_max=0,*r_rtn_=NULL;
  int c_rtn_max=0,*c_rtn_=NULL;
  int n_r_rtn_index=0,*r_rtn_index_=NULL;
  int n_r_rmv_index=0,*r_rmv_index_=NULL;
  int n_c_rtn_index=0,*c_rtn_index_=NULL;
  int n_c_rmv_index=0,*c_rmv_index_=NULL;
  double *ZR__=NULL;
  double *ZR_sub_max_s_=NULL;
  double nlp_gumb_opt=0,nlp_gumb_emp=0,p_gumb_opt=0,p_gumb_emp=0;
  int nl=0,n_l=0;
  MDA_dim_ = (int *) malloc1(2*sizeof(int));
  if (verbose){ printf(" %% [entering halfloop_nonbinary_recursive_helper_ZR__]\n");}
  /* %%%%%%%% */
  MDA_read_r8(&MDA_n_dim,&MDA_dim_,&trace__,fname_trace__);
  if (MDA_n_dim!=2){ printf(" %% Warning, improper n_dim %d in fanme_trace__\n",MDA_n_dim);}
  if (MDA_dim_[0]!=6){ printf(" %% Warning, improper dim_[0] %d in fname_trace__\n",MDA_dim_[0]);}
  n_iteration = MDA_dim_[1];
  if (verbose>2){ array_printf_margin(trace__,"double",6,n_iteration," % trace__: ");}
  r_rtn_ = (int *) malloc1((unsigned long long int)n_iteration*sizeof(int));
  for (niteration=0;niteration<n_iteration;niteration++){ r_rtn_[niteration] = (int)round(trace__[trace_r_rtn_index + niteration*6]);}
  r_rtn_max = r_rtn_[0];
  c_rtn_ = (int *) malloc1((unsigned long long int)n_iteration*sizeof(int));
  for (niteration=0;niteration<n_iteration;niteration++){ c_rtn_[niteration] = (int)round(trace__[trace_c_rtn_index + niteration*6]);}
  c_rtn_max = c_rtn_[0];
  if (flag_r0drop_vs_rcdrop==flag_rcdrop){ n_xdrop = r_rtn_max + c_rtn_max;}
  if (flag_r0drop_vs_rcdrop==flag_r0drop){ n_xdrop = r_rtn_max;}
  MDA_read_r8(&MDA_n_dim,&MDA_dim_,&QR__,fname_QR__);
  if (MDA_n_dim!=2){ printf(" %% Warning, improper n_dim %d in fname_QR__\n",MDA_n_dim);}
  if (MDA_dim_[0]!=n_iteration){ printf(" %% Warning, improper dim_[0] %d in fname_QR__\n",MDA_dim_[0]);}
  n_shuffle = MDA_dim_[1]-1;
  if (verbose>2){ array_printf_margin(QR__,"double",n_iteration,1+n_shuffle," % QR__: ");}
  MDA_read_i4(&MDA_n_dim,&MDA_dim_,&xdrop__,fname_xdrop__);
  if (MDA_n_dim!=2){ printf(" %% Warning, improper n_dim %d in fname_xdrop__\n",MDA_n_dim);}
  if (MDA_dim_[0]!=2){ printf(" %% Warning, improper dim_[0] %d in fname_xdrop__\n",MDA_dim_[0]);}
  if (MDA_dim_[1]!=n_xdrop){ printf(" %% Warning, improper dim_[1] %d in fname_xdrop__\n",MDA_dim_[1]);}
  /* %%%%%%%% */
  QR_ = QR__;
  ZR__ = (double *) malloc1((unsigned long long int)n_iteration*(unsigned long long int)(1+n_shuffle)*sizeof(double));
  ZR_sub_max_s_ = (double *) malloc1((unsigned long long int)(1+n_shuffle)*sizeof(double));
  if (verbose>2){ array_printf_margin(QR_,"double",1,n_iteration," % QR_: ");}
  nlp_ZR_ = (double *) malloc1((unsigned long long int)n_iteration*sizeof(double));
  for (niteration=0;niteration<n_iteration;niteration++){
    QR_avg = 0; for (nshuffle=1;nshuffle<1+n_shuffle;nshuffle++){ QR_avg += QR__[niteration+nshuffle*n_iteration];}
    QR_avg /= maximum(1,n_shuffle);
    QR_std = 0; for (nshuffle=1;nshuffle<1+n_shuffle;nshuffle++){ QR_std += pow(QR__[niteration+nshuffle*n_iteration] - QR_avg,2);}
    QR_std /= maximum(1,n_shuffle); QR_std = sqrt(QR_std);
    ZR = (QR_[niteration] - QR_avg)/maximum(1e-12,QR_std);
    nlp_ZR = -z_to_lp_single_d(ZR);
    nlp_ZR_[niteration] = nlp_ZR;
    for (nshuffle=0;nshuffle<1+n_shuffle;nshuffle++){
      ZR = (QR__[niteration + nshuffle*n_iteration] - QR_avg)/maximum(1e-12,QR_std);
      ZR__[niteration + nshuffle*n_iteration] = ZR;
      /* for (nshuffle=0;nshuffle<1+n_shuffle;nshuffle++){ } */}
    /* for (niteration=0;niteration<n_iteration;niteration++){ } */}
  if (verbose>2){ array_printf_margin(nlp_ZR_,"double",1,n_iteration," % nlp_ZR_: ");}
  index_sub_ = (int *) malloc1((unsigned long long int)n_iteration*sizeof(int));
  nlp_ZR_sub_ = (double *) malloc1((unsigned long long int)n_iteration*sizeof(double));
  n_sub=0;
  for (niteration=0;niteration<n_iteration;niteration++){
    if ( (r_rtn_[niteration]>=n_member_lob) && (r_rtn_max - r_rtn_[niteration]>=n_member_lob) ){
      index_sub_[n_sub] = niteration; nlp_ZR_sub_[n_sub] = nlp_ZR_[niteration]; n_sub++;
      /* if sufficiently many members */}
    /* for (niteration=0;niteration<n_iteration;niteration++){ } */}
  if (n_sub<1){
    if (verbose>2){ printf(" %% n_sub %d in halfloop_nonbinary_recursive_helper_ZR__, using all of nlp_ZR_\n",n_sub);}
    n_sub = n_iteration;
    for (niteration=0;niteration<n_iteration;niteration++){
      index_sub_[niteration] = niteration;
      nlp_ZR_sub_[niteration] = nlp_ZR_[niteration];
      /* for (niteration=0;niteration<n_iteration;niteration++){ } */}
    /* if (n_sub<1){ } */}
  for (nshuffle=0;nshuffle<1+n_shuffle;nshuffle++){
    ZR_sub_max_s_[nshuffle] = ZR__[index_sub_[0] + nshuffle*n_iteration];
    for (nsub=0;nsub<n_sub;nsub++){
      ZR = ZR__[index_sub_[nsub] + nshuffle*n_iteration];
      ZR_sub_max_s_[nshuffle] = maximum(ZR,ZR_sub_max_s_[nshuffle]);
      /* for (nsub=0;nsub<n_sub;nsub++){ } */}
    /* for (nshuffle=0;nshuffle<1+n_shuffle;nshuffle++){ } */}
  gumbel_fit(n_shuffle,ZR_sub_max_s_+1,ZR_sub_max_s_[0],NULL,&nlp_gumb_opt,&nlp_gumb_emp,&p_gumb_opt,&p_gumb_emp);
  if (nlp_gumb_opt_p!=NULL){ *nlp_gumb_opt_p = nlp_gumb_opt;}
  if (nlp_gumb_emp_p!=NULL){ *nlp_gumb_emp_p = nlp_gumb_emp;}
  if (verbose>2){ printf(" %% nlp_gumb_opt %f p_gumb_opt %f nlp_gumb_emp %f p_gumb_emp %f\n",nlp_gumb_opt,p_gumb_opt,nlp_gumb_emp,p_gumb_emp);}
  find_internal_maximum(0*verbose,n_sub,nlp_ZR_sub_,0,&nlp_ZR_max,&nlp_ZR_index);
  if (nlp_ZR_max>=-log(p_use)){ find_internal_maximum(0*verbose,n_sub,nlp_ZR_sub_,-log(p_use),&nlp_ZR_max,&nlp_ZR_index);}
  if (nlp_ZR_index<=-1){
    if ( (verbose>0) && (n_sub>=5) ){
      printf(" %% Warning, no maximum found, n_sub %d\n",n_sub);
      array_printf_margin(nlp_ZR_sub_,"double",1,n_sub," % nlp_ZR_sub: ");
      /* if ( (verbose>0) && (n_sub>=5) ){ } */}
    nlp_ZR_index = 0; nlp_ZR_max = nlp_ZR_sub_[nlp_ZR_index];
    /* if (nlp_ZR_index<=-1){ } */}
  if (verbose>2){ printf(" %% index_sub_[%d] = %d\n",nlp_ZR_index,index_sub_[nlp_ZR_index]);}
  nlp_ZR_index = index_sub_[nlp_ZR_index];
  n_r_rtn_index = trace__[trace_r_rtn_index + nlp_ZR_index*6]; n_r_rmv_index = r_rtn_max - n_r_rtn_index;
  n_c_rtn_index = trace__[trace_c_rtn_index + nlp_ZR_index*6]; n_c_rmv_index = c_rtn_max - n_c_rtn_index;
  if (verbose){ printf(" %% nlp_ZR_max %f nlp_ZR_index %d n_r_rtn_index %d/%d n_c_rtn_index %d/%d\n",nlp_ZR_max,nlp_ZR_index,n_r_rtn_index,n_r_rmv_index,n_c_rtn_index,n_c_rmv_index);}
  if (nlp_ZR_max_p!=NULL){ *nlp_ZR_max_p = nlp_ZR_max;}
  if (nlp_ZR_index_p!=NULL){ *nlp_ZR_index_p = nlp_ZR_index;}
  if (n_r_rtn_index_p!=NULL){ *n_r_rtn_index_p = n_r_rtn_index;}
  if (n_r_rmv_index_p!=NULL){ *n_r_rmv_index_p = n_r_rmv_index;}
  if (n_c_rtn_index_p!=NULL){ *n_c_rtn_index_p = n_c_rtn_index;}
  if (n_c_rmv_index_p!=NULL){ *n_c_rmv_index_p = n_c_rmv_index;}
  r_rtn_index_=NULL;
  if (r_rtn_index_p_!=NULL){
    if ( (*r_rtn_index_p_)==NULL ){ (*r_rtn_index_p_) = (int *) malloc1((unsigned long long int)n_r_rtn_index*sizeof(int));}
    r_rtn_index_ = *r_rtn_index_p_;
    /* if (r_rtn_index_p_!=NULL){ } */}
  if (r_rtn_index_!=NULL){
    xdrop_ = &(xdrop__[0 + (n_xdrop-1)*2]);
    nl = 0; while (nl<n_r_rtn_index){ if ((*xdrop_)>-1){ r_rtn_index_[nl++] = *xdrop_;} xdrop_-=2;}
    if (verbose>2){ array_printf_margin(r_rtn_index_,"int",1,n_r_rtn_index," % r_rtn_index_: ");}
    /* if (r_rtn_index_!=NULL){ } */}
  r_rmv_index_=NULL;
  if (r_rmv_index_p_!=NULL){
    if ( (*r_rmv_index_p_)==NULL ){ (*r_rmv_index_p_) = (int *) malloc1((unsigned long long int)n_r_rmv_index*sizeof(int));}
    r_rmv_index_ = *r_rmv_index_p_;
    /* if (r_rmv_index_p_!=NULL){ } */}
  if (r_rmv_index_!=NULL){
    xdrop_ = &(xdrop__[0 + (0)*2]);
    nl = 0; while (nl<n_r_rmv_index){ if ((*xdrop_)>-1){ r_rmv_index_[nl++] = *xdrop_;} xdrop_+=2;}
    if (verbose>2){ array_printf_margin(r_rmv_index_,"int",1,n_r_rmv_index," % r_rmv_index_: ");}
    /* if (r_rtn_index_!=NULL){ } */}
  if (flag_r0drop_vs_rcdrop==flag_rcdrop){
    c_rtn_index_=NULL;
    if (c_rtn_index_p_!=NULL){
      if ( (*c_rtn_index_p_)==NULL ){ (*c_rtn_index_p_) = (int *) malloc1((unsigned long long int)n_c_rtn_index*sizeof(int));}
      c_rtn_index_ = *c_rtn_index_p_;
      /* if (c_rtn_index_p_!=NULL){ } */}
    if (c_rtn_index_!=NULL){
      xdrop_ = &(xdrop__[1 + (n_xdrop-1)*2]);
      nl = 0; while (nl<n_c_rtn_index){ if ((*xdrop_)>-1){ c_rtn_index_[nl++] = *xdrop_;} xdrop_-=2;}
      if (verbose>2){ array_printf_margin(c_rtn_index_,"int",1,n_c_rtn_index," % c_rtn_index_: ");}
      /* if (c_rtn_index_!=NULL){ } */}
    c_rmv_index_=NULL;
    if (c_rmv_index_p_!=NULL){
      if ( (*c_rmv_index_p_)==NULL ){ (*c_rmv_index_p_) = (int *) malloc1((unsigned long long int)n_c_rmv_index*sizeof(int));}
      c_rmv_index_ = *c_rmv_index_p_;
      /* if (c_rmv_index_p_!=NULL){ } */}
    if (c_rmv_index_!=NULL){
      xdrop_ = &(xdrop__[1 + (0)*2]);
      nl = 0; while (nl<n_c_rmv_index){ if ((*xdrop_)>-1){ c_rmv_index_[nl++] = *xdrop_;} xdrop_+=2;}
      if (verbose>2){ array_printf_margin(c_rmv_index_,"int",1,n_c_rmv_index," % c_rmv_index_: ");}
      /* if (c_rmv_index_!=NULL){ } */}    
    /* if (flag_r0drop_vs_rcdrop==flag_rcdrop){ } */}
  free1(&index_sub_);
  free1(&nlp_ZR_sub_);
  free1(&nlp_ZR_);
  free1(&ZR__);
  free1(&ZR_sub_max_s_);
  free1(&r_rtn_);
  free1(&c_rtn_);
  free1(&trace__);
  free1(&QR__);
  free1(&xdrop__);
  free1(&MDA_dim_);
  /* %%%%%%%% */
  if (verbose){ printf(" %% [finished halfloop_nonbinary_recursive_helper_ZR__]\n");}
}

/* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */

void halfloop_nonbinary_f_recursive_prefix_base1(char *prefix_base0,double gamma,int n_member_lob,char *prefix_base1)
{
  char prefix_n_member_lob[FNAMESIZE];
  char prefix_gamma[FNAMESIZE];
  if (gamma> 0){ sprintf(prefix_gamma,"_g%.3d",(int)floor(1000*gamma));} else{ sprintf(prefix_gamma,"");}
  if (n_member_lob> 2){ sprintf(prefix_n_member_lob,"_n%.2d",n_member_lob);} else{ sprintf(prefix_n_member_lob,"");}
  sprintf(prefix_base1,"%s%s%s",prefix_base0,prefix_gamma,prefix_n_member_lob);
}

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
)
{
  int verbose_t=0;
  int flag_rcdrop = 0;
  int flag_r0drop = 1;
  char *cwd=NULL;
  int n_str_path=0;
  char dir_trunk[PNAMESIZE];
  char dir_out[PNAMESIZE],dir_A_out[PNAMESIZE],dir_B_out[PNAMESIZE];
  char prefix_base0[FNAMESIZE];
  char prefix_base1[FNAMESIZE];
  char dir_0in[PNAMESIZE];
  struct stat stat_dir = {0};
  struct stat stat_file = {0};
  int na=0;
  int nr=0,flag_free_r_index=0,nr_index=0,n_r_index=0,*r_index_=NULL;
  int nc=0,flag_free_c_index=0,nc_index=0,n_c_index=0,*c_index_=NULL;
  double gamma=0;
  int nshuffle=0,n_shuffle=0;
  int n_member_lob=0;
  double p_set=0,p_use=0,p_prev=0;
  int flag_force_create=0;
  int MDA_n_dim=0;
  int *MDA_dim_=NULL;
  char MDA_fname[PNAMESIZE];
  unsigned long long int *binary_label_=NULL;
  int length_output_label = 0;
  int length_nlpbra_label = 0;
  int length_nlpnex_label = 0;
  char **output_label__=NULL;
  char **nlpbra_label__=NULL;
  char **nlpnex_label__=NULL;
  unsigned long long int *binary_label_A_=NULL;
  char **output_label_A__=NULL;
  char **nlpbra_label_A__=NULL;
  char **nlpnex_label_A__=NULL;
  unsigned long long int *binary_label_B_=NULL;
  char **output_label_B__=NULL;
  char **nlpbra_label_B__=NULL;
  char **nlpnex_label_B__=NULL;
  char fname_trace_E__[PNAMESIZE];
  char fname_QR_E__[PNAMESIZE];
  char fname_xdrop_E__[PNAMESIZE];
  char fname_trace_F__[PNAMESIZE];
  char fname_QR_F__[PNAMESIZE];
  char fname_xdrop_F__[PNAMESIZE];
  double nlp_ZR_E_max=0; int nlp_ZR_E_index=0;
  double nlp_ZR_F_max=0; int nlp_ZR_F_index=0;
  int n_c_rtn_index_E=0,*c_rtn_index_E_=NULL,*c_index_rtn_sub_=NULL;
  double nlp_gumb_opt=0,nlp_gumb_emp=0;
  int nr_rtn_index_F=0,n_r_rtn_index_F=0,*r_rtn_index_F_=NULL,*r_index_rtn_sub_=NULL;
  int nr_rmv_index_F=0,n_r_rmv_index_F=0,*r_rmv_index_F_=NULL,*r_index_rmv_sub_=NULL;
  double *nlp_p_split_=NULL;
  double p_branch=0,p_next=0;
  int flag_split=0;
  char fname_output_label[PNAMESIZE];
  char fname_nlpbra_label[PNAMESIZE];
  char fname_nlpnex_label[PNAMESIZE];
  FILE *fp=NULL;
  if (verbose>0){ printf(" %% [entering halfloop_nonbinary_f_recursive], flag_omp %d\n",flag_omp);}
  /* %%%%%%%%%%%%%%%% */
  n_str_path = pathconf(".",_PC_PATH_MAX);
  MDA_dim_ = (int *) malloc1(2*sizeof(int)); //%<-- should only ever need 2. ;
  nlp_p_split_ = (double *) malloc1((1+14)*sizeof(double)); //%<-- gumb_emp, gumb_opt, F, E, p_branch, p_next, p_set, p_use, n_r_rtn_index_F, n_r_rmv_index_F, n_member_lob, ndepth, recursion_limit, flag_split. ;
  cwd = (char *) malloc1((size_t)n_str_path); getcwd(cwd,n_str_path);
  if (n_str_path> FNAMESIZE){ printf(" %% Warning, n_str_path %d in halfloop_nonbinary_f_recursive\n",n_str_path);}
  if ( (dir_trunk_0in!=NULL) && (strlen(dir_trunk_0in)> FNAMESIZE) ){ printf(" %% Warning, dir_trunk_0in %s too long in halfloop_nonbinary_f_recursive\n",dir_trunk_0in);}
  if ( (dir_out_0in!=NULL) && (strlen(dir_out_0in)> FNAMESIZE) ){ printf(" %% Warning, dir_out_0in %s too long in halfloop_nonbinary_f_recursive\n",dir_out_0in);}
  if ( (prefix_base_0in!=NULL) && (strlen(prefix_base_0in)> FNAMESIZE) ){ printf(" %% Warning, prefix_base_0in %s too long in halfloop_nonbinary_f_recursive\n",prefix_base_0in);}
  if ( dir_trunk_0in==NULL ){ if (verbose>0){ printf(" %% cwd: [%s]\n",cwd);} sprintf(dir_trunk,"%s",cwd);} else{ sprintf(dir_trunk,"%s",dir_trunk_0in);}
  if ( prefix_base_0in==NULL ){ sprintf(prefix_base0,"%s","test");} else{ sprintf(prefix_base0,"%s",prefix_base_0in);}
  if (verbose>1){
    printf(" %% dir_trunk: %s\n",dir_trunk);
    printf(" %% prefix_base0: %s\n",prefix_base0);
    /* if (verbose>1){ } */}
  if ( (n_r_index_0in<=0) || (r_index_0in_==NULL) ){
    n_r_index = n_r;
    r_index_ = (int *) malloc1((unsigned long long int)n_r*sizeof(int));
    for (nr=0;nr<n_r;nr++){ r_index_[nr] = nr;}
    flag_free_r_index = 1;
    /* if ( (n_r_index_0in<=0) || (r_index_0in_==NULL) ){ } */}
  else /* if ( (n_r_index_0in>0) && (r_index_0in_!=NULL) ) */{
    n_r_index = n_r_index_0in;
    r_index_ = r_index_0in_;
    flag_free_r_index = 0;
    /* exists */}
  if ( (n_c_index_0in<=0) || (c_index_0in_==NULL) ){
    n_c_index = n_c;
    c_index_ = (int *) malloc1((unsigned long long int)n_c*sizeof(int));
    for (nc=0;nc<n_c;nc++){ c_index_[nc] = nc;}
    flag_free_c_index = 1;
    /* if ( (n_r_index_0in<=0) || (r_index_0in_==NULL) ){ } */}
  else /* if ( (n_c_index_0in>0) && (c_index_0in_!=NULL) ) */{
    n_c_index = n_c_index_0in;
    c_index_ = c_index_0in_;
    flag_free_c_index = 0;
    /* exists */}
  if (dir_out_0in!=NULL){
    sprintf(dir_out,"%s",dir_out_0in);
    /* if (dir_out_0in!=NULL){ } */}
  if ( gamma_0in<=0 ){ gamma = 0.0; } else{ gamma = minimum(1,gamma_0in);}
  if (n_member_lob_0in<=0){ n_member_lob = 2;} else{ n_member_lob = n_member_lob_0in;}
  if (dir_out_0in==NULL){
    sprintf(dir_0in,"%s/dir_%s",dir_trunk,prefix_base0);
    if (stat(dir_0in,&stat_dir)==-1){ printf(" %% %s not found, creating\n",dir_0in); mkdir(dir_0in,0755);} else{ if (verbose){ printf(" %% %s found, not creating\n",dir_0in);}}
    halfloop_nonbinary_f_recursive_prefix_base1(prefix_base0,gamma,n_member_lob,prefix_base1);
    sprintf(dir_out,"%s/dir_%s",dir_0in,prefix_base1);
    if (verbose>1){
      printf(" %% dir_0in: %s\n",dir_0in);
      printf(" %% prefix_base1: %s\n",prefix_base1);
      printf(" %% dir_out: %s\n",dir_out);
      /* if (verbose>1){ } */}
    /* if (dir_out_0in==NULL){ } */}
  if (stat(dir_out,&stat_dir)==-1){ printf(" %% %s not found, creating\n",dir_out); mkdir(dir_out,0755);} else{ if (verbose){ printf(" %% %s found, not creating\n",dir_out);}}
  if (n_shuffle_0in<=0){ n_shuffle = 64;} else{ n_shuffle = n_shuffle_0in;}
  if (p_set_0in<=0){ p_set = 0.05;} else{ p_set = p_set_0in;}
  p_use = (double)p_set / (double)(1+2*p_set);
  if (p_prev_0in<=0){ p_prev = 0.00;} else{ p_prev = p_prev_0in;}
  if (flag_force_create_0in<=0){ flag_force_create = 0;} else{ flag_force_create = flag_force_create_0in;}
  if (verbose>0){
    printf(" %% n_r %d n_c %d --> n_r_index %d n_c_index %d\n",n_r,n_c,n_r_index,n_c_index);
    printf(" %% gamma %0.3f n_shuffle %d p_set %0.3f (p_use %0.3f) n_member_lob %d p_prev %0.3f flag_force_create %d  --> dir_out %s\n",gamma,n_shuffle,p_set,p_use,n_member_lob,p_prev,flag_force_create,dir_out);
    /* if (verbose>0){ } */}
  sprintf(MDA_fname,"%s/r_index_.mda",dir_out); MDA_n_dim = 2; MDA_dim_[0] = n_r_index; MDA_dim_[1] = 1;
  MDA_write_i4(MDA_n_dim,MDA_dim_,r_index_,MDA_fname);
  if (verbose>2){ MDA_printf_i4_margin(MDA_fname);}
  //sprintf(MDA_fname,"%s/c_index_.mda",dir_out); MDA_n_dim = 2; MDA_dim_[0] = n_c_index; MDA_dim_[1] = 1;
  //MDA_write_i4(MDA_n_dim,MDA_dim_,c_index_,MDA_fname);
  //if (verbose>2){ MDA_printf_i4_margin(MDA_fname);}
  /* %%%%%%%%%%%%%%%% */
  /* initialize labels */
  /* %%%%%%%%%%%%%%%% */
  binary_label_ = NULL;
  *binary_label_p_ = (unsigned  long long int *) malloc1((unsigned long long int)n_r_index*sizeof(unsigned long long int));
  binary_label_ = *binary_label_p_; for (nr_index=0;nr_index<n_r_index;nr_index++){ binary_label_[nr_index]=(unsigned long long int)0;}
  length_output_label = 8 + 1 + 1*(recursion_limit-ndepth) ; //%<-- allows for 1 character per level, as well as root label. ;
  output_label__ = NULL; malloc1_char__(n_r_index,length_output_label,output_label_p_); output_label__ = *output_label_p_;
  if (output_label__!=NULL){ for (nr_index=0;nr_index<n_r_index;nr_index++){ sprintf(output_label__[nr_index],"0");}}
  length_nlpbra_label = 8 + 1 + 8*(recursion_limit-ndepth+1) ; //%<-- allows for 4+2 digits and space, as well as leaf label. ;
  nlpbra_label__ = NULL; malloc1_char__(n_r_index,length_nlpbra_label,nlpbra_label_p_); nlpbra_label__ = *nlpbra_label_p_;
  if (nlpbra_label__!=NULL){ for (nr_index=0;nr_index<n_r_index;nr_index++){ sprintf(nlpbra_label__[nr_index],"");}}
  length_nlpnex_label = 8 + 1 + 8*(recursion_limit-ndepth+1) ; //%<-- allows for 4+2 digits and space, as well as leaf label. ;
  nlpnex_label__ = NULL; malloc1_char__(n_r_index,length_nlpnex_label,nlpnex_label_p_); nlpnex_label__ = *nlpnex_label_p_;
  if (nlpnex_label__!=NULL){ for (nr_index=0;nr_index<n_r_index;nr_index++){ sprintf(nlpnex_label__[nr_index],"");}}
  /* %%%%%%%%%%%%%%%% */
  /* perform calculation */
  /* %%%%%%%%%%%%%%%% */
  if (flag_r0drop_vs_rcdrop==flag_rcdrop){
    if (verbose>1){ printf(" %% n_iteration_E %d\n",get_xdrop_logscale_length(n_r_index,flag_r0drop_vs_rcdrop==flag_rcdrop ? n_c_index : maximum(n_r_index,n_c_index),gamma));}
    sprintf(fname_trace_E__,"%s/trace_E__.mda",dir_out);
    sprintf(fname_QR_E__,"%s/QR_E__.mda",dir_out);
    sprintf(fname_xdrop_E__,"%s/xdrop_E__.mda",dir_out);
    if (flag_omp==0){
      GLOBAL_tic(1);
      halfloop_nonbinary_f_recursive_helper_QR__
	(
	 maximum(0,verbose-1)
	 ,n_r
	 ,n_c
	 ,E_base_rc__
	 ,n_r_index
	 ,r_index_
	 ,n_c_index
	 ,c_index_
	 ,flag_rcdrop
	 ,gamma
	 ,n_shuffle
	 ,flag_force_create
	 ,fname_trace_E__
	 ,fname_xdrop_E__
	 ,fname_QR_E__
	 );
      GLOBAL_toc(1,verbose_t," % halfloop_nonbinary_f_recursive_helper_QR__: ");
    /* if (flag_omp==0){ } */}
    if (flag_omp==1){
      GLOBAL_tic(1);
      halfloop_nonbinary_f_recursive_omp_helper_QR__
	(
	 maximum(0,verbose-1)
	 ,n_r
	 ,n_c
	 ,E_base_rc__
	 ,n_r_index
	 ,r_index_
	 ,n_c_index
	 ,c_index_
	 ,flag_rcdrop
	 ,gamma
	 ,n_shuffle
	 ,flag_force_create
	 ,fname_trace_E__
	 ,fname_xdrop_E__
	 ,fname_QR_E__
	 );
      GLOBAL_toc(1,verbose_t," % halfloop_nonbinary_f_recursive_omp_helper_QR__: ");
    /* if (flag_omp==1){ } */}
    GLOBAL_tic(2);
    halfloop_nonbinary_recursive_helper_ZR__
      (
       maximum(0,verbose-1)
       ,flag_rcdrop
       ,fname_trace_E__
       ,fname_xdrop_E__
       ,fname_QR_E__
       ,p_use
       ,0*n_member_lob
       ,&nlp_ZR_E_max
       ,&nlp_ZR_E_index
       ,NULL
       ,NULL
       ,NULL
       ,NULL
       ,&n_c_rtn_index_E
       ,&c_rtn_index_E_
       ,NULL
       ,NULL
       ,&nlp_gumb_opt
       ,&nlp_gumb_emp
       );
    GLOBAL_toc(2,verbose_t," % halfloop_nonbinary_f_recursive_helper_ZR__: ");
    /* if (flag_r0drop_vs_rcdrop==flag_rcdrop){ } */}
  if (flag_r0drop_vs_rcdrop==flag_r0drop){
    if (verbose>1){ printf(" %% using all columns, skipping variable-selection.\n");}
    n_c_rtn_index_E = n_c_index;
    c_rtn_index_E_ = (int *) malloc1((unsigned long long int)n_c_rtn_index_E*sizeof(int));
    memcpy(c_rtn_index_E_,c_index_,(unsigned long long int)n_c_index*sizeof(int));
    /* if (flag_r0drop_vs_rcdrop==flag_r0drop){ } */}
  array_extract_i_from_i(1,n_c_index,c_index_,0,NULL,n_c_rtn_index_E,c_rtn_index_E_,&c_index_rtn_sub_,NULL);
  iquicksort_index(0,c_index_rtn_sub_,1,NULL,0,n_c_rtn_index_E-1);
  if (verbose>2){ array_printf_margin(c_index_rtn_sub_,"int",1,n_c_rtn_index_E," % c_index_rtn_sub_: ");}
  /* %%%% */
  sprintf(fname_trace_F__,"%s/trace_F__.mda",dir_out);
  sprintf(fname_QR_F__,"%s/QR_F__.mda",dir_out);
  sprintf(fname_xdrop_F__,"%s/xdrop_F__.mda",dir_out);
  if (flag_r0drop_vs_rcdrop==flag_rcdrop){
    if (flag_omp==0){
      GLOBAL_tic(1);
      halfloop_nonbinary_f_recursive_helper_QR__
	(
	 maximum(0,verbose-1)
	 ,n_r
	 ,n_c
	 ,E_base_rc__
	 ,n_r_index
	 ,r_index_
	 ,n_c_rtn_index_E
	 ,c_index_rtn_sub_
	 ,flag_r0drop
	 ,0*gamma
	 ,n_shuffle
	 ,flag_force_create
	 ,fname_trace_F__
	 ,fname_xdrop_F__
	 ,fname_QR_F__
	 );
      GLOBAL_toc(1,verbose_t," % halfloop_nonbinary_f_recursive_helper_QR__: ");
      /* if (flag_omp==0){ } */}
    if (flag_omp==1){
      GLOBAL_tic(1);
      halfloop_nonbinary_f_recursive_omp_helper_QR__
	(
	 maximum(0,verbose-1)
	 ,n_r
	 ,n_c
	 ,E_base_rc__
	 ,n_r_index
	 ,r_index_
	 ,n_c_rtn_index_E
	 ,c_index_rtn_sub_
	 ,flag_r0drop
	 ,0*gamma
	 ,n_shuffle
	 ,flag_force_create
	 ,fname_trace_F__
	 ,fname_xdrop_F__
	 ,fname_QR_F__
	 );
      GLOBAL_toc(1,verbose_t," % halfloop_nonbinary_f_recursive_omp_helper_QR__: ");
      /* if (flag_omp==1){ } */}
    GLOBAL_tic(2);
    halfloop_nonbinary_recursive_helper_ZR__
      (
       maximum(0,verbose-1)
       ,flag_r0drop
       ,fname_trace_F__
       ,fname_xdrop_F__
       ,fname_QR_F__
       ,p_use
       ,1*n_member_lob
       ,&nlp_ZR_F_max
       ,&nlp_ZR_F_index
       ,&n_r_rtn_index_F
       ,&r_rtn_index_F_
       ,&n_r_rmv_index_F
       ,&r_rmv_index_F_
       ,NULL
       ,NULL
       ,NULL
       ,NULL
       ,NULL
       ,NULL
       );
    GLOBAL_toc(2,verbose_t," % halfloop_nonbinary_f_recursive_helper_ZR__: ");
    /* if (flag_r0drop_vs_rcdrop==flag_rcdrop){ } */}
  if (flag_r0drop_vs_rcdrop==flag_r0drop){
    if (flag_omp==0){
      GLOBAL_tic(1);
      halfloop_nonbinary_f_recursive_helper_QR__
	(
	 maximum(0,verbose-1)
	 ,n_r
	 ,n_c
	 ,E_base_rc__
	 ,n_r_index
	 ,r_index_
	 ,n_c_rtn_index_E
	 ,c_index_rtn_sub_
	 ,flag_r0drop
	 ,1*gamma
	 ,n_shuffle
	 ,flag_force_create
	 ,fname_trace_F__
	 ,fname_xdrop_F__
	 ,fname_QR_F__
	 );
      GLOBAL_toc(1,verbose_t," % halfloop_nonbinary_f_recursive_helper_QR__: ");
      /* if (flag_omp==0){ } */}
    if (flag_omp==1){
      GLOBAL_tic(1);
      halfloop_nonbinary_f_recursive_omp_helper_QR__
	(
	 maximum(0,verbose-1)
	 ,n_r
	 ,n_c
	 ,E_base_rc__
	 ,n_r_index
	 ,r_index_
	 ,n_c_rtn_index_E
	 ,c_index_rtn_sub_
	 ,flag_r0drop
	 ,1*gamma
	 ,n_shuffle
	 ,flag_force_create
	 ,fname_trace_F__
	 ,fname_xdrop_F__
	 ,fname_QR_F__
	 );
      GLOBAL_toc(1,verbose_t," % halfloop_nonbinary_f_recursive_omp_helper_QR__: ");
      /* if (flag_omp==1){ } */}
    GLOBAL_tic(2);
    halfloop_nonbinary_recursive_helper_ZR__
      (
       maximum(0,verbose-1)
       ,flag_r0drop
       ,fname_trace_F__
       ,fname_xdrop_F__
       ,fname_QR_F__
       ,p_use
       ,1*n_member_lob
       ,&nlp_ZR_F_max
       ,&nlp_ZR_F_index
       ,&n_r_rtn_index_F
       ,&r_rtn_index_F_
       ,&n_r_rmv_index_F
       ,&r_rmv_index_F_
       ,NULL
       ,NULL
       ,NULL
       ,NULL
       ,&nlp_gumb_opt
       ,&nlp_gumb_emp
       );
    GLOBAL_toc(2,verbose_t," % halfloop_nonbinary_f_recursive_helper_ZR__: ");
    nlp_ZR_E_max = nlp_ZR_F_max;
    /* if (flag_r0drop_vs_rcdrop==flag_r0drop){ } */}
  iquicksort_index(0,r_rtn_index_F_,1,NULL,0,n_r_rtn_index_F-1);
  iquicksort_index(0,r_rmv_index_F_,1,NULL,0,n_r_rmv_index_F-1);
  array_extract_i_from_i(1,n_r_index,r_index_,0,NULL,n_r_rtn_index_F,r_rtn_index_F_,&r_index_rtn_sub_,NULL);
  iquicksort_index(0,r_index_rtn_sub_,1,NULL,0,n_r_rtn_index_F-1);
  if (verbose>2){ array_printf_margin(r_index_rtn_sub_,"int",1,n_r_rtn_index_F," % r_index_rtn_sub_: ");}
  array_extract_i_from_i(1,n_r_index,r_index_,0,NULL,n_r_rmv_index_F,r_rmv_index_F_,&r_index_rmv_sub_,NULL);
  iquicksort_index(0,r_index_rmv_sub_,1,NULL,0,n_r_rmv_index_F-1);
  if (verbose>2){ array_printf_margin(r_index_rmv_sub_,"int",1,n_r_rmv_index_F," % r_index_rmv_sub_: ");}
  /* %%%% */
  p_branch = exp(-nlp_gumb_emp); if (p_branch<(double)1.0/(double)n_shuffle){ p_branch = exp(-nlp_gumb_opt);}
  p_next = 1.0 - (1.0-p_prev)*(1-p_branch); //%<-- q_next = q_prev*q_branch;
  flag_split = (ndepth< recursion_limit) && (n_r_rtn_index_F>=n_member_lob) && (n_r_rmv_index_F>=n_member_lob) && (p_next<=p_use);
  if (verbose>0){ printf(" %% p_branch %f p_next %f n_r_rtn_index_F %d n_r_rmv_index_F %d flag_split %d\n",p_branch,p_next,n_r_rtn_index_F,n_r_rmv_index_F,flag_split);} 
  na=0;
  nlp_p_split_[na++] = nlp_gumb_emp;
  nlp_p_split_[na++] = nlp_gumb_opt;
  nlp_p_split_[na++] = nlp_ZR_E_max;
  nlp_p_split_[na++] = nlp_ZR_F_max;
  nlp_p_split_[na++] = p_branch;
  nlp_p_split_[na++] = p_next;
  nlp_p_split_[na++] = p_set;
  nlp_p_split_[na++] = p_use;
  nlp_p_split_[na++] = (double)n_r_rtn_index_F;
  nlp_p_split_[na++] = (double)n_r_rmv_index_F;
  nlp_p_split_[na++] = (double)n_member_lob;
  nlp_p_split_[na++] = (double)ndepth;
  nlp_p_split_[na++] = (double)recursion_limit;
  nlp_p_split_[na++] = (double)flag_split;
  sprintf(MDA_fname,"%s/nlp_p_split_.mda",dir_out); MDA_n_dim = 2; MDA_dim_[0] = 9; MDA_dim_[1] = 1;
  MDA_write_r8(MDA_n_dim,MDA_dim_,nlp_p_split_,MDA_fname);
  if (verbose>2){ MDA_printf_r8_margin(MDA_fname);}
  if (flag_split==1){
    /* %%%% */
    sprintf(dir_A_out,"%s/A",dir_out);
    if (stat(dir_A_out,&stat_dir)==-1){ printf(" %% %s not found, creating\n",dir_A_out); mkdir(dir_A_out,0755);} else{ if (verbose){ printf(" %% %s found, not creating\n",dir_A_out);}}
    halfloop_nonbinary_f_recursive
      (
       verbose
       ,flag_omp
       ,1+ndepth
       ,recursion_limit
       ,n_r
       ,n_c
       ,E_base_rc__
       ,n_r_rtn_index_F
       ,r_index_rtn_sub_
       ,n_c_index
       ,c_index_
       ,flag_r0drop_vs_rcdrop
       ,gamma
       ,n_shuffle
       ,p_set
       ,n_member_lob
       ,p_next
       ,dir_trunk
       ,dir_A_out
       ,prefix_base0
       ,flag_force_create
       ,&binary_label_A_
       ,&output_label_A__
       ,&nlpbra_label_A__
       ,&nlpnex_label_A__
       );      
    /* %%%% */
    sprintf(dir_B_out,"%s/B",dir_out);
    if (stat(dir_B_out,&stat_dir)==-1){ printf(" %% %s not found, creating\n",dir_B_out); mkdir(dir_B_out,0755);} else{ if (verbose){ printf(" %% %s found, not creating\n",dir_B_out);}}
    halfloop_nonbinary_f_recursive
      (
       verbose
       ,flag_omp
       ,1+ndepth
       ,recursion_limit
       ,n_r
       ,n_c
       ,E_base_rc__
       ,n_r_rmv_index_F
       ,r_index_rmv_sub_
       ,n_c_index
       ,c_index_
       ,flag_r0drop_vs_rcdrop
       ,gamma
       ,n_shuffle
       ,p_set
       ,n_member_lob
       ,p_next
       ,dir_trunk
       ,dir_B_out
       ,prefix_base0
       ,flag_force_create
       ,&binary_label_B_
       ,&output_label_B__
       ,&nlpbra_label_B__
       ,&nlpnex_label_B__
       );      
    /* %%%% */
    /* update labels */
    /* %%%% */
    for (nr_rtn_index_F=0;nr_rtn_index_F<n_r_rtn_index_F;nr_rtn_index_F++){
      nr_index = r_rtn_index_F_[nr_rtn_index_F];
      if (binary_label_!=NULL){ if (binary_label_A_[nr_rtn_index_F]>(unsigned long long int)(4.295e9)){ printf(" %% Warning, too many clusters for binary_label_, binary_label overflow\n");} binary_label_[nr_index] += (unsigned long long int)0+(unsigned long long int)2*binary_label_A_[nr_rtn_index_F];}
      if (output_label__!=NULL){ if (strlen(output_label_A__[nr_rtn_index_F])>FNAMESIZE-16){ printf(" %% Warning, too many clusters for output_label__, increase label length\n");} sprintf(output_label__[nr_index],"%s%s","A",output_label_A__[nr_rtn_index_F]);}
      if (nlpbra_label__!=NULL){ if (strlen(nlpbra_label_A__[nr_rtn_index_F])>FNAMESIZE-16){ printf(" %% Warning, too many clusters for nlpbra_label__, increase label length\n");} sprintf(nlpbra_label__[nr_index],"%0.2f %s",minimum(9999,-log(p_branch)),nlpbra_label_A__[nr_rtn_index_F]);} //%<-- restrict to 4+2 digits. ;
      if (nlpnex_label__!=NULL){ if (strlen(nlpnex_label_A__[nr_rtn_index_F])>FNAMESIZE-16){ printf(" %% Warning, too many clusters for nlpnex_label__, increase label length\n");} sprintf(nlpnex_label__[nr_index],"%0.2f %s",minimum(9999,-log(p_next)),nlpnex_label_A__[nr_rtn_index_F]);} //%<-- restrict to 4+2 digits. ;
      /* for (nr_rtn_index_F=0;nr_rtn_index_F<n_r_rtn_index_F;nr_rtn_index_F++){ } */}
    for (nr_rmv_index_F=0;nr_rmv_index_F<n_r_rmv_index_F;nr_rmv_index_F++){
      nr_index = r_rmv_index_F_[nr_rmv_index_F];
      if (binary_label_!=NULL){ if (binary_label_B_[nr_rmv_index_F]>(unsigned long long int)(4.295e9)){ printf(" %% Warning, too many clusters for binary_label_, binary_label overflow\n");} binary_label_[nr_index] += (unsigned long long int)1+(unsigned long long int)2*binary_label_B_[nr_rmv_index_F];}
      if (output_label__!=NULL){ if (strlen(output_label_B__[nr_rmv_index_F])>FNAMESIZE-16){ printf(" %% Warning, too many clusters for output_label__, increase label length\n");} sprintf(output_label__[nr_index],"%s%s","B",output_label_B__[nr_rmv_index_F]);}
      if (nlpbra_label__!=NULL){ if (strlen(nlpbra_label_B__[nr_rmv_index_F])>FNAMESIZE-16){ printf(" %% Warning, too many clusters for nlpbra_label__, increase label length\n");} sprintf(nlpbra_label__[nr_index],"%0.2f %s",minimum(9999,-log(p_branch)),nlpbra_label_B__[nr_rmv_index_F]);} //%<-- restrict to 4+2 digits. ;
      if (nlpnex_label__!=NULL){ if (strlen(nlpnex_label_B__[nr_rmv_index_F])>FNAMESIZE-16){ printf(" %% Warning, too many clusters for nlpnex_label__, increase label length\n");} sprintf(nlpnex_label__[nr_index],"%0.2f %s",minimum(9999,-log(p_next)),nlpnex_label_B__[nr_rmv_index_F]);} //%<-- restrict to 4+2 digits. ;
      /* for (nr_rmv_index_F=0;nr_rmv_index_F<n_r_rmv_index_F;nr_rmv_index_F++){ } */}
    /* %%%% */
    /* free labels */
    /* %%%% */
    free1(&binary_label_A_);
    free1_char__(n_r_rtn_index_F,&output_label_A__);
    free1_char__(n_r_rtn_index_F,&nlpbra_label_A__);
    free1_char__(n_r_rtn_index_F,&nlpnex_label_A__);
    free1(&binary_label_B_);
    free1_char__(n_r_rmv_index_F,&output_label_B__);
    free1_char__(n_r_rmv_index_F,&nlpbra_label_B__);
    free1_char__(n_r_rmv_index_F,&nlpnex_label_B__);
    /* %%%% */
    /* if (flag_split==1){ } */}
  if (flag_split==0){
    /* %%%% */
    /* update labels with nlp */
    /* %%%% */
    for (nr_index=0;nr_index<n_r_index;nr_index++){
      sprintf(nlpbra_label__[nr_index],"%0.2f",minimum(9999,-log(p_branch))); //%<-- restrict to 4+2 digits. ;
      sprintf(nlpnex_label__[nr_index],"%0.2f",minimum(9999,-log(p_next))); //%<-- restrict to 4+2 digits. ;
      /* for (nr_index=0;nr_index<n_r_index;nr_index++){ } */}
    /* if (flag_split==0){ } */}
  /* %%%% */
  /* if ndepth==0 then print labels to file */
  /* %%%% */
  if (ndepth==0){
    sprintf(fname_output_label,"%s/output_label__.txt",dir_out);
    if ((fp=fopen(fname_output_label,"w"))==NULL){ printf(" %% Warning! could not open %s.\n",fname_output_label); exit(EXIT_FAILURE);}
    for (nr_index=0;nr_index<n_r_index;nr_index++){
      fprintf(fp,"%s\n",output_label__[nr_index]);
      /* for (nr_index=0;nr_index<n_r_index;nr_index++){ } */}
    fclose(fp);fp=NULL;
    sprintf(fname_nlpbra_label,"%s/nlpbra_label__.txt",dir_out);
    if ((fp=fopen(fname_nlpbra_label,"w"))==NULL){ printf(" %% Warning! could not open %s.\n",fname_nlpbra_label); exit(EXIT_FAILURE);}
    for (nr_index=0;nr_index<n_r_index;nr_index++){
      fprintf(fp,"%s\n",nlpbra_label__[nr_index]);
      /* for (nr_index=0;nr_index<n_r_index;nr_index++){ } */}
    fclose(fp);fp=NULL;
    sprintf(fname_nlpnex_label,"%s/nlpnex_label__.txt",dir_out);
    if ((fp=fopen(fname_nlpnex_label,"w"))==NULL){ printf(" %% Warning! could not open %s.\n",fname_nlpnex_label); exit(EXIT_FAILURE);}
    for (nr_index=0;nr_index<n_r_index;nr_index++){
      fprintf(fp,"%s\n",nlpnex_label__[nr_index]);
      /* for (nr_index=0;nr_index<n_r_index;nr_index++){ } */}
    fclose(fp);fp=NULL;
    /* if (ndepth==0){ } */}
  /* %%%% */
  free1(&r_index_rmv_sub_);
  free1(&r_index_rtn_sub_);
  free1(&c_index_rtn_sub_);
  free1(&c_rtn_index_E_);
  free1(&r_rmv_index_F_);
  free1(&r_rtn_index_F_);
  /* %%%%%%%%%%%%%%%% */
  if (flag_free_r_index){ free1(&r_index_);}
  if (flag_free_c_index){ free1(&c_index_);}
  free1(&nlp_p_split_);
  free1(&cwd);
  free1(&MDA_dim_);
  /* %%%%%%%%%%%%%%%% */
  if (verbose>0){ printf(" %% [finished halfloop_nonbinary_f_recursive], flag_omp %d\n",flag_omp);}
}

void halfloop_nonbinary_f_recursive_test()
{
  int verbose=0;
  int flag_omp=GLOBAL_flag_omp_use;
  int recursion_limit=GLOBAL_halfloop_recursion_limit;
  int flag_r0drop_vs_rcdrop=0;
  int flag_rcdrop = 0;
  int flag_r0drop = 1;
  double gamma = 0.02;
  int n_shuffle = 64;
  double p_set = 0.05;
  int n_member_lob = 3;
  int nr=0,n_r = 128;
  int nc=0,n_c = 512;
  int nr_index=0,n_r_index=0,*r_index_=NULL;
  int nc_index=0,n_c_index=0,*c_index_=NULL;
  float x=0,y=0,z=0;
  unsigned long int rseed=0;
  unsigned long long int ulli=0;
  float *E_base_rc__=NULL;
  int flag_force_create=1;
  unsigned long long int *binary_label_=NULL;
  char binary_label_str_[FNAMESIZE];
  char **output_label__=NULL;
  char **nlpbra_label__=NULL;
  char **nlpnex_label__=NULL;
  int flag_error=0;
  char output_label_rc_ans__[128][32] = { "A0" , "A0" , "A0" , "A0" , "A0" , "A0" , "A0" , "A0" , "A0" , "A0" , "A0" , "A0" , "A0" , "A0" , "A0" , "A0" , "A0" , "A0" , "A0" , "A0" , "A0" , "A0" , "A0" , "A0" , "A0" , "A0" , "A0" , "A0" , "A0" , "A0" , "A0" , "A0" , "A0" , "BA0" , "BA0" , "BA0" , "BA0" , "BA0" , "BA0" , "BA0" , "BA0" , "BA0" , "BA0" , "BA0" , "BA0" , "BA0" , "BA0" , "BA0" , "BA0" , "BA0" , "BA0" , "BA0" , "BA0" , "BA0" , "BA0" , "BB0" , "BB0" , "BB0" , "BB0" , "BB0" , "BB0" , "BB0" , "BB0" , "BB0" , "BB0" , "A0" , "A0" , "A0" , "A0" , "A0" , "A0" , "A0" , "A0" , "A0" , "A0" , "A0" , "A0" , "A0" , "A0" , "A0" , "A0" , "A0" , "A0" , "A0" , "A0" , "A0" , "A0" , "A0" , "A0" , "A0" , "A0" , "A0" , "A0" , "A0" , "A0" , "A0" , "BA0" , "BA0" , "BA0" , "BA0" , "BA0" , "BA0" , "BA0" , "BA0" , "BA0" , "BA0" , "BA0" , "BA0" , "BA0" , "BA0" , "BA0" , "BA0" , "BA0" , "BA0" , "BA0" , "BB0" , "BB0" , "BB0" , "BB0" , "BB0" , "BB0" , "BB0" , "BB0" , "A0" , "A0" , "A0" , "A0" , "A0" };
  char nlpbra_label_rc_ans__[128][32] = { "5.97 " , "5.97 " , "5.97 " , "5.97 " , "5.97 " , "5.97 " , "5.97 " , "5.97 " , "5.97 " , "5.97 " , "5.97 " , "5.97 " , "5.97 " , "5.97 " , "5.97 " , "5.97 " , "5.97 " , "5.97 " , "5.97 " , "5.97 " , "5.97 " , "5.97 " , "5.97 " , "5.97 " , "5.97 " , "5.97 " , "5.97 " , "5.97 " , "5.97 " , "5.97 " , "5.97 " , "5.97 " , "5.97 " , "5.97 3.47 " , "5.97 3.47 " , "5.97 3.47 " , "5.97 3.47 " , "5.97 3.47 " , "5.97 3.47 " , "5.97 3.47 " , "5.97 3.47 " , "5.97 3.47 " , "5.97 3.47 " , "5.97 3.47 " , "5.97 3.47 " , "5.97 3.47 " , "5.97 3.47 " , "5.97 3.47 " , "5.97 3.47 " , "5.97 3.47 " , "5.97 3.47 " , "5.97 3.47 " , "5.97 3.47 " , "5.97 3.47 " , "5.97 3.47 " , "5.97 3.47 " , "5.97 3.47 " , "5.97 3.47 " , "5.97 3.47 " , "5.97 3.47 " , "5.97 3.47 " , "5.97 3.47 " , "5.97 3.47 " , "5.97 3.47 " , "5.97 3.47 " , "5.97 " , "5.97 " , "5.97 " , "5.97 " , "5.97 " , "5.97 " , "5.97 " , "5.97 " , "5.97 " , "5.97 " , "5.97 " , "5.97 " , "5.97 " , "5.97 " , "5.97 " , "5.97 " , "5.97 " , "5.97 " , "5.97 " , "5.97 " , "5.97 " , "5.97 " , "5.97 " , "5.97 " , "5.97 " , "5.97 " , "5.97 " , "5.97 " , "5.97 " , "5.97 " , "5.97 " , "5.97 3.47 " , "5.97 3.47 " , "5.97 3.47 " , "5.97 3.47 " , "5.97 3.47 " , "5.97 3.47 " , "5.97 3.47 " , "5.97 3.47 " , "5.97 3.47 " , "5.97 3.47 " , "5.97 3.47 " , "5.97 3.47 " , "5.97 3.47 " , "5.97 3.47 " , "5.97 3.47 " , "5.97 3.47 " , "5.97 3.47 " , "5.97 3.47 " , "5.97 3.47 " , "5.97 3.47 " , "5.97 3.47 " , "5.97 3.47 " , "5.97 3.47 " , "5.97 3.47 " , "5.97 3.47 " , "5.97 3.47 " , "5.97 3.47 " , "5.97 " , "5.97 " , "5.97 " , "5.97 " , "5.97 " };
  char nlpnex_label_rc_ans__[128][32] = { "5.97 " , "5.97 " , "5.97 " , "5.97 " , "5.97 " , "5.97 " , "5.97 " , "5.97 " , "5.97 " , "5.97 " , "5.97 " , "5.97 " , "5.97 " , "5.97 " , "5.97 " , "5.97 " , "5.97 " , "5.97 " , "5.97 " , "5.97 " , "5.97 " , "5.97 " , "5.97 " , "5.97 " , "5.97 " , "5.97 " , "5.97 " , "5.97 " , "5.97 " , "5.97 " , "5.97 " , "5.97 " , "5.97 " , "5.97 3.39 " , "5.97 3.39 " , "5.97 3.39 " , "5.97 3.39 " , "5.97 3.39 " , "5.97 3.39 " , "5.97 3.39 " , "5.97 3.39 " , "5.97 3.39 " , "5.97 3.39 " , "5.97 3.39 " , "5.97 3.39 " , "5.97 3.39 " , "5.97 3.39 " , "5.97 3.39 " , "5.97 3.39 " , "5.97 3.39 " , "5.97 3.39 " , "5.97 3.39 " , "5.97 3.39 " , "5.97 3.39 " , "5.97 3.39 " , "5.97 3.39 " , "5.97 3.39 " , "5.97 3.39 " , "5.97 3.39 " , "5.97 3.39 " , "5.97 3.39 " , "5.97 3.39 " , "5.97 3.39 " , "5.97 3.39 " , "5.97 3.39 " , "5.97 " , "5.97 " , "5.97 " , "5.97 " , "5.97 " , "5.97 " , "5.97 " , "5.97 " , "5.97 " , "5.97 " , "5.97 " , "5.97 " , "5.97 " , "5.97 " , "5.97 " , "5.97 " , "5.97 " , "5.97 " , "5.97 " , "5.97 " , "5.97 " , "5.97 " , "5.97 " , "5.97 " , "5.97 " , "5.97 " , "5.97 " , "5.97 " , "5.97 " , "5.97 " , "5.97 " , "5.97 3.39 " , "5.97 3.39 " , "5.97 3.39 " , "5.97 3.39 " , "5.97 3.39 " , "5.97 3.39 " , "5.97 3.39 " , "5.97 3.39 " , "5.97 3.39 " , "5.97 3.39 " , "5.97 3.39 " , "5.97 3.39 " , "5.97 3.39 " , "5.97 3.39 " , "5.97 3.39 " , "5.97 3.39 " , "5.97 3.39 " , "5.97 3.39 " , "5.97 3.39 " , "5.97 3.39 " , "5.97 3.39 " , "5.97 3.39 " , "5.97 3.39 " , "5.97 3.39 " , "5.97 3.39 " , "5.97 3.39 " , "5.97 3.39 " , "5.97 " , "5.97 " , "5.97 " , "5.97 " , "5.97 "};
  char output_label_r0_ans__[128][32] = { "A0" , "A0" , "A0" , "A0" , "A0" , "A0" , "A0" , "A0" , "A0" , "A0" , "A0" , "A0" , "A0" , "A0" , "A0" , "A0" , "A0" , "A0" , "A0" , "A0" , "A0" , "A0" , "A0" , "A0" , "A0" , "A0" , "A0" , "A0" , "A0" , "A0" , "A0" , "BA0" , "BA0" , "BA0" , "BA0" , "BA0" , "BA0" , "BA0" , "BA0" , "BBA0" , "BBA0" , "BBA0" , "BBA0" , "BBA0" , "BBA0" , "BBA0" , "BBA0" , "BBA0" , "BBA0" , "BBA0" , "BBA0" , "BBA0" , "BBA0" , "BBA0" , "BBA0" , "BBA0" , "BBA0" , "BBA0" , "BBB0" , "BBB0" , "BBB0" , "BBB0" , "BBB0" , "BBB0" , "BBB0" , "BBB0" , "BBB0" , "A0" , "A0" , "A0" , "A0" , "A0" , "A0" , "A0" , "A0" , "A0" , "A0" , "A0" , "A0" , "A0" , "A0" , "A0" , "A0" , "A0" , "A0" , "A0" , "A0" , "A0" , "A0" , "A0" , "A0" , "A0" , "A0" , "A0" , "BA0" , "BA0" , "BA0" , "BA0" , "BA0" , "BA0" , "BA0" , "BA0" , "BA0" , "BBA0" , "BBA0" , "BBA0" , "BBA0" , "BBA0" , "BBA0" , "BBA0" , "BBA0" , "BBA0" , "BBA0" , "BBA0" , "BBA0" , "BBA0" , "BBA0" , "BBB0" , "BBB0" , "BBB0" , "BBB0" , "BBB0" , "BBB0" , "BBB0" , "A0" , "A0" , "A0" , "A0" };
  char nlpbra_label_r0_ans__[128][32] = { "7.21 " , "7.21 " , "7.21 " , "7.21 " , "7.21 " , "7.21 " , "7.21 " , "7.21 " , "7.21 " , "7.21 " , "7.21 " , "7.21 " , "7.21 " , "7.21 " , "7.21 " , "7.21 " , "7.21 " , "7.21 " , "7.21 " , "7.21 " , "7.21 " , "7.21 " , "7.21 " , "7.21 " , "7.21 " , "7.21 " , "7.21 " , "7.21 " , "7.21 " , "7.21 " , "7.21 " , "7.21 4.23 " , "7.21 4.23 " , "7.21 4.23 " , "7.21 4.23 " , "7.21 4.23 " , "7.21 4.23 " , "7.21 4.23 " , "7.21 4.23 " , "7.21 4.23 4.78 " , "7.21 4.23 4.78 " , "7.21 4.23 4.78 " , "7.21 4.23 4.78 " , "7.21 4.23 4.78 " , "7.21 4.23 4.78 " , "7.21 4.23 4.78 " , "7.21 4.23 4.78 " , "7.21 4.23 4.78 " , "7.21 4.23 4.78 " , "7.21 4.23 4.78 " , "7.21 4.23 4.78 " , "7.21 4.23 4.78 " , "7.21 4.23 4.78 " , "7.21 4.23 4.78 " , "7.21 4.23 4.78 " , "7.21 4.23 4.78 " , "7.21 4.23 4.78 " , "7.21 4.23 4.78 " , "7.21 4.23 4.78 " , "7.21 4.23 4.78 " , "7.21 4.23 4.78 " , "7.21 4.23 4.78 " , "7.21 4.23 4.78 " , "7.21 4.23 4.78 " , "7.21 4.23 4.78 " , "7.21 4.23 4.78 " , "7.21 4.23 4.78 " , "7.21 " , "7.21 " , "7.21 " , "7.21 " , "7.21 " , "7.21 " , "7.21 " , "7.21 " , "7.21 " , "7.21 " , "7.21 " , "7.21 " , "7.21 " , "7.21 " , "7.21 " , "7.21 " , "7.21 " , "7.21 " , "7.21 " , "7.21 " , "7.21 " , "7.21 " , "7.21 " , "7.21 " , "7.21 " , "7.21 " , "7.21 " , "7.21 4.23 " , "7.21 4.23 " , "7.21 4.23 " , "7.21 4.23 " , "7.21 4.23 " , "7.21 4.23 " , "7.21 4.23 " , "7.21 4.23 " , "7.21 4.23 " , "7.21 4.23 4.78 " , "7.21 4.23 4.78 " , "7.21 4.23 4.78 " , "7.21 4.23 4.78 " , "7.21 4.23 4.78 " , "7.21 4.23 4.78 " , "7.21 4.23 4.78 " , "7.21 4.23 4.78 " , "7.21 4.23 4.78 " , "7.21 4.23 4.78 " , "7.21 4.23 4.78 " , "7.21 4.23 4.78 " , "7.21 4.23 4.78 " , "7.21 4.23 4.78 " , "7.21 4.23 4.78 " , "7.21 4.23 4.78 " , "7.21 4.23 4.78 " , "7.21 4.23 4.78 " , "7.21 4.23 4.78 " , "7.21 4.23 4.78 " , "7.21 4.23 4.78 " , "7.21 " , "7.21 " , "7.21 " , "7.21 " };
  char nlpnex_label_r0_ans__[128][32] = { "7.21 " , "7.21 " , "7.21 " , "7.21 " , "7.21 " , "7.21 " , "7.21 " , "7.21 " , "7.21 " , "7.21 " , "7.21 " , "7.21 " , "7.21 " , "7.21 " , "7.21 " , "7.21 " , "7.21 " , "7.21 " , "7.21 " , "7.21 " , "7.21 " , "7.21 " , "7.21 " , "7.21 " , "7.21 " , "7.21 " , "7.21 " , "7.21 " , "7.21 " , "7.21 " , "7.21 " , "7.21 4.18 " , "7.21 4.18 " , "7.21 4.18 " , "7.21 4.18 " , "7.21 4.18 " , "7.21 4.18 " , "7.21 4.18 " , "7.21 4.18 " , "7.21 4.18 3.75 " , "7.21 4.18 3.75 " , "7.21 4.18 3.75 " , "7.21 4.18 3.75 " , "7.21 4.18 3.75 " , "7.21 4.18 3.75 " , "7.21 4.18 3.75 " , "7.21 4.18 3.75 " , "7.21 4.18 3.75 " , "7.21 4.18 3.75 " , "7.21 4.18 3.75 " , "7.21 4.18 3.75 " , "7.21 4.18 3.75 " , "7.21 4.18 3.75 " , "7.21 4.18 3.75 " , "7.21 4.18 3.75 " , "7.21 4.18 3.75 " , "7.21 4.18 3.75 " , "7.21 4.18 3.75 " , "7.21 4.18 3.75 " , "7.21 4.18 3.75 " , "7.21 4.18 3.75 " , "7.21 4.18 3.75 " , "7.21 4.18 3.75 " , "7.21 4.18 3.75 " , "7.21 4.18 3.75 " , "7.21 4.18 3.75 " , "7.21 4.18 3.75 " , "7.21 " , "7.21 " , "7.21 " , "7.21 " , "7.21 " , "7.21 " , "7.21 " , "7.21 " , "7.21 " , "7.21 " , "7.21 " , "7.21 " , "7.21 " , "7.21 " , "7.21 " , "7.21 " , "7.21 " , "7.21 " , "7.21 " , "7.21 " , "7.21 " , "7.21 " , "7.21 " , "7.21 " , "7.21 " , "7.21 " , "7.21 " , "7.21 4.18 " , "7.21 4.18 " , "7.21 4.18 " , "7.21 4.18 " , "7.21 4.18 " , "7.21 4.18 " , "7.21 4.18 " , "7.21 4.18 " , "7.21 4.18 " , "7.21 4.18 3.75 " , "7.21 4.18 3.75 " , "7.21 4.18 3.75 " , "7.21 4.18 3.75 " , "7.21 4.18 3.75 " , "7.21 4.18 3.75 " , "7.21 4.18 3.75 " , "7.21 4.18 3.75 " , "7.21 4.18 3.75 " , "7.21 4.18 3.75 " , "7.21 4.18 3.75 " , "7.21 4.18 3.75 " , "7.21 4.18 3.75 " , "7.21 4.18 3.75 " , "7.21 4.18 3.75 " , "7.21 4.18 3.75 " , "7.21 4.18 3.75 " , "7.21 4.18 3.75 " , "7.21 4.18 3.75 " , "7.21 4.18 3.75 " , "7.21 4.18 3.75 " , "7.21 " , "7.21 " , "7.21 " , "7.21 " };
  GLOBAL_flag_orth_brute = 1;
  E_base_rc__ = (float *) malloc1((unsigned long long int)n_r*(unsigned long long int)n_c*sizeof(float));
  rseed=1; RSEED_adv8(&rseed);
  ulli=0;
  for (nc=0;nc<n_c;nc++){ for (nr=0;nr<n_r;nr++){ 
      //E_base_rc__[ulli] = RNGET(&rseed);
      //if ( (nc<floor(sqrt((double)n_c))) && (nr<floor(sqrt((double)n_r))) ){ E_base_rc__[ulli] += 0.5; }
      x = 2.0*(float)nr/(float)(n_r-1) - 1.0;
      y = 2.0*(float)nc/(float)(n_c-1) - 1.0;
      z = (x+y)/4.0;
      E_base_rc__[ulli] = sin(2*PI_LF*x) + cos(2*PI_LF*2*y) + x*x + y*y*y + cos(2*PI_LF*4*z);
      ulli++;
      /* for (nc=0;nc<n_c;nc++){ for (nr=0;nr<n_r;nr++){ }} */}}
  array_printf_margin(E_base_rc__,"float",n_r,n_c," % E_base_r__: ");
  r_index_ = (int *) malloc1((unsigned long long int)n_r*sizeof(int));
  n_r_index = 0; for (nr=0;nr<n_r;nr++){ if (nr%1==0){ r_index_[n_r_index++]=nr;}}
  c_index_ = (int *) malloc1((unsigned long long int)n_c*sizeof(int));
  n_c_index = 0; for (nc=0;nc<n_c;nc++){ if (nc%1==0){ c_index_[n_c_index++]=nc;}}
  /* %%%%%%%% */
  GLOBAL_tic(0);
  halfloop_nonbinary_f_recursive(
  verbose
 ,flag_omp
 ,0
 ,recursion_limit
 ,n_r
 ,n_c
 ,E_base_rc__
 ,n_r_index
 ,r_index_
 ,n_c_index
 ,c_index_
 ,flag_rcdrop
 ,gamma
 ,n_shuffle
 ,p_set
 ,n_member_lob
 ,-1
 ,NULL
 ,NULL
 ,"test_rc"
 ,flag_force_create
 ,&binary_label_
 ,&output_label__
 ,&nlpbra_label__
 ,&nlpnex_label__
 );
  if (verbose){
    for (nr_index=0;nr_index<n_r_index;nr_index++){
      bitstring_from_ulli(binary_label_ + nr_index,binary_label_str_,64);
      printf(" %% %lld --> %s --> %s --> %s --> %s\n",binary_label_[nr_index],binary_label_str_,output_label__[nr_index],nlpbra_label__[nr_index],nlpnex_label__[nr_index]);
      /* for (nr_index=0;nr_index<n_r_index;nr_index++){ } */}
  /* if (verbose){ } */}
  flag_error = 0;
  for (nr_index=0;nr_index<n_r_index;nr_index++){
    flag_error += strcmp(output_label__[nr_index],output_label_rc_ans__[nr_index]);
    flag_error += (strstr(nlpbra_label__[nr_index],nlpbra_label_rc_ans__[nr_index]) == nlpbra_label__);
    flag_error += (strstr(nlpnex_label__[nr_index],nlpnex_label_rc_ans__[nr_index]) == nlpnex_label__);
    if (verbose){ printf(" %% nr_index %d/%d: %d %s %s %s\n",nr_index,n_r_index,binary_label_[nr_index],output_label__[nr_index],nlpbra_label__[nr_index],nlpnex_label__[nr_index]);}
    /* for (nr_index=0;nr_index<n_r_index;nr_index++){ } */}
  printf(" %% output_label_rc_ans__, nlpbra_label_rc_ans__ and nlpnex_label_rc_ans__ vs output_label__, nlpbra_label__ and nlpex_label__:  flag_error %d\n",flag_error);
  free1(&binary_label_);
  free1_char__(n_r_index,&output_label__);
  free1_char__(n_r_index,&nlpbra_label__);
  free1_char__(n_r_index,&nlpnex_label__);
  GLOBAL_toc(0,1," % halfloop_nonbinary_f_recursive: ");
  /* %%%%%%%% */
  GLOBAL_tic(0);
  halfloop_nonbinary_f_recursive(
  verbose
 ,flag_omp
 ,0
 ,recursion_limit
 ,n_r
 ,n_c
 ,E_base_rc__
 ,n_r_index
 ,r_index_
 ,n_c_index
 ,c_index_
 ,flag_r0drop
 ,gamma
 ,n_shuffle
 ,p_set
 ,n_member_lob
 ,-1
 ,NULL
 ,NULL
 ,"test_r0"
 ,flag_force_create
 ,&binary_label_
 ,&output_label__
 ,&nlpbra_label__
 ,&nlpnex_label__
 );
  if (verbose){
    for (nr_index=0;nr_index<n_r_index;nr_index++){
      bitstring_from_ulli(binary_label_ + nr_index,binary_label_str_,64);
      printf(" %% %lld --> %s --> %s --> %s --> %s\n",binary_label_[nr_index],binary_label_str_,output_label__[nr_index],nlpbra_label__[nr_index],nlpnex_label__[nr_index]);
      /* for (nr_index=0;nr_index<n_r_index;nr_index++){ } */}
  /* if (verbose){ } */}
  flag_error = 0;
  for (nr_index=0;nr_index<n_r_index;nr_index++){
    flag_error += strcmp(output_label__[nr_index],output_label_r0_ans__[nr_index]);
    flag_error += (strstr(nlpbra_label__[nr_index],nlpbra_label_r0_ans__[nr_index]) == nlpbra_label__);
    flag_error += (strstr(nlpnex_label__[nr_index],nlpnex_label_r0_ans__[nr_index]) == nlpnex_label__);
    if (verbose){ printf(" %% nr_index %d/%d: %d %s %s %s\n",nr_index,n_r_index,binary_label_[nr_index],output_label__[nr_index],nlpbra_label__[nr_index],nlpnex_label__[nr_index]);}
    /* for (nr_index=0;nr_index<n_r_index;nr_index++){ } */}
  printf(" %% output_label_r0_ans__, nlpbra_label_r0_ans__ and nlpnex_label_r0_ans__ vs output_label__, nlpbra_label__ and nlpex_label__:  flag_error %d\n",flag_error);
  free1(&binary_label_);
  free1_char__(n_r_index,&output_label__);
  free1_char__(n_r_index,&nlpbra_label__);
  free1_char__(n_r_index,&nlpnex_label__);
  GLOBAL_toc(0,1," % halfloop_nonbinary_f_recursive: ");
  /* %%%%%%%% */
  free1(&E_base_rc__);
  free1(&r_index_);
  free1(&c_index_);
}

void halfloop_nonbinary_f_recursive_test_speed()
{
  int verbose=1;
  int flag_omp=GLOBAL_flag_omp_use;
  int recursion_limit=GLOBAL_halfloop_recursion_limit;
  int flag_r0drop_vs_rcdrop=0;
  int flag_rcdrop = 0;
  int flag_r0drop = 1;
  double gamma = 0.01;
  int n_shuffle = 64;
  double p_set = 0.05;
  int n_member_lob = 3;
  int nr=0,n_r = 178*2;
  int nc=0,n_c = 2e4;
  int nr_index=0,n_r_index=0,*r_index_=NULL;
  int nc_index=0,n_c_index=0,*c_index_=NULL;
  float x=0,y=0,z=0;
  unsigned long int rseed=0;
  unsigned long long int ulli=0;
  float *E_base_rc__=NULL;
  int flag_force_create=1;
  unsigned long long int *binary_label_=NULL;
  char **output_label__=NULL;
  char **nlpbra_label__=NULL;
  char **nlpnex_label__=NULL;
  E_base_rc__ = (float *) malloc1((unsigned long long int)n_r*(unsigned long long int)n_c*sizeof(float));
  rseed=1; RSEED_adv8(&rseed);
  ulli=0;
  for (nc=0;nc<n_c;nc++){ for (nr=0;nr<n_r;nr++){ 
      //E_base_rc__[ulli] = RNGET(&rseed);
      //if ( (nc<floor(sqrt((double)n_c))) && (nr<floor(sqrt((double)n_r))) ){ E_base_rc__[ulli] += 0.5; }
      x = 2.0*(float)nr/(float)(n_r-1) - 1.0;
      y = 2.0*(float)nc/(float)(n_c-1) - 1.0;
      z = (x+y)/4.0;
      E_base_rc__[ulli] = sin(2*PI_LF*x) + cos(2*PI_LF*2*y) + x*x + y*y*y + cos(2*PI_LF*4*z);
      ulli++;
      /* for (nc=0;nc<n_c;nc++){ for (nr=0;nr<n_r;nr++){ }} */}}
  array_printf_margin(E_base_rc__,"float",n_r,n_c," % E_base_r__: ");
  r_index_ = (int *) malloc1((unsigned long long int)n_r*sizeof(int));
  n_r_index = 0; for (nr=0;nr<n_r;nr++){ if (nr%1==0){ r_index_[n_r_index++]=nr;}}
  c_index_ = (int *) malloc1((unsigned long long int)n_c*sizeof(int));
  n_c_index = 0; for (nc=0;nc<n_c;nc++){ if (nc%1==0){ c_index_[n_c_index++]=nc;}}
  /* %%%%%%%% */
  GLOBAL_tic(0);
  halfloop_nonbinary_f_recursive(
  verbose
 ,flag_omp
 ,0
 ,recursion_limit
 ,n_r
 ,n_c
 ,E_base_rc__
 ,n_r_index
 ,r_index_
 ,n_c_index
 ,c_index_
 ,flag_rcdrop
 ,gamma
 ,n_shuffle
 ,p_set
 ,n_member_lob
 ,-1
 ,NULL
 ,NULL
 ,"test_speed_rc"
 ,flag_force_create
 ,&binary_label_
 ,&output_label__
 ,&nlpbra_label__
 ,&nlpnex_label__
 );
  free1(&binary_label_);
  free1_char__(n_r_index,&output_label__);
  free1_char__(n_r_index,&nlpbra_label__);
  free1_char__(n_r_index,&nlpnex_label__);
  GLOBAL_toc(0,1," % halfloop_nonbinary_f_recursive: ");
  /* %%%%%%%% */
  free1(&E_base_rc__);
  free1(&r_index_);
  free1(&c_index_);
}

/* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */

void halfloop_nonbinary_f_gateway_shell()
{
  int verbose=GLOBAL_verbose;
  int n_d=0;
  int *d_=NULL;
  float *E_base_rc__=NULL;
  int n_r=0,n_c=0,nr=0,nc=0;
  unsigned long long int *binary_label_=NULL;
  char **output_label__=NULL;
  char **nlpbra_label__=NULL;
  char **nlpnex_label__=NULL;
  if (verbose){ printf(" %% [entering halfloop_nonbinary_f_gateway_shell]\n");}
  MDA_read_r4(&n_d,&d_,&E_base_rc__,GLOBAL_E_base_mda_r4);
  if (verbose){ array_printf(d_,"int",1,n_d," % d_: ");}
  n_r = d_[0]; n_c = d_[1];
  if (verbose>1){
    for (nr=0;nr<16;nr++){
      printf(" %%");
      for (nc=0;nc<16;nc++){
	printf(" %0.4f",E_base_rc__[nr+nc*n_r]);
	/* for (nc=0;nc<16;nc++){ } */}
      printf(" \n");
      /* for (nr=0;nr<16;nr++){} */}
    /* if (verbose>1){ } */}
  /* %%%%%%%% */
  if (verbose){
    printf(" %% calling halfloop_nonbinary_f_recursive with parameters: \n");
    printf(" %% verbose %d\n",GLOBAL_verbose);
    printf(" %% flag_omp_use %d\n",GLOBAL_flag_omp_use);
    printf(" %% flag_split %d\n",GLOBAL_flag_split);
    printf(" %% recursion_limit %d\n",GLOBAL_halfloop_recursion_limit);
    printf(" %% E_base_mda_r4 %s\n",GLOBAL_E_base_mda_r4);
    printf(" %% flag_r0drop_vs_rcdrop %d\n",GLOBAL_flag_r0drop_vs_rcdrop);
    printf(" %% gamma %0.3f\n",GLOBAL_gamma);
    printf(" %% n_shuffle %d\n",GLOBAL_n_shuffle);
    printf(" %% p_set %0.3f\n",GLOBAL_p_set);
    printf(" %% n_member_lob %d\n",GLOBAL_n_member_lob);
    printf(" %% dir_trunk %s\n",GLOBAL_dir_trunk);
    printf(" %% prefix_base %s\n",GLOBAL_prefix_base);
    printf(" %% flag_force_create %d\n",GLOBAL_flag_force_create);
    /* if (verbose){ } */}
  /* %%%%%%%% */
  if (GLOBAL_flag_split==0){
  halfloop_nonbinary_f_recursive(
  GLOBAL_verbose
 ,GLOBAL_flag_omp_use
 ,0
 ,GLOBAL_halfloop_recursion_limit
 ,n_r
 ,n_c
 ,E_base_rc__
 ,0
 ,NULL
 ,0
 ,NULL
 ,GLOBAL_flag_r0drop_vs_rcdrop
 ,GLOBAL_gamma
 ,GLOBAL_n_shuffle
 ,GLOBAL_p_set
 ,GLOBAL_n_member_lob
 ,-1
 ,GLOBAL_dir_trunk
 ,NULL
 ,GLOBAL_prefix_base
 ,GLOBAL_flag_force_create
 ,&binary_label_
 ,&output_label__
 ,&nlpbra_label__
 ,&nlpnex_label__
 );
  /* if (GLOBAL_flag_split==0){ } */}
  /* %%%%%%%% */
  if (GLOBAL_flag_split==1){
  halfloop_split_nonbinary_f_recursive(
  GLOBAL_verbose
 ,GLOBAL_flag_omp_use
 ,0
 ,GLOBAL_halfloop_recursion_limit
 ,n_r
 ,n_c
 ,E_base_rc__
 ,0
 ,NULL
 ,0
 ,NULL
 ,GLOBAL_flag_r0drop_vs_rcdrop
 ,GLOBAL_gamma
 ,GLOBAL_n_shuffle
 ,GLOBAL_p_set
 ,GLOBAL_n_member_lob
 ,-1
 ,GLOBAL_dir_trunk
 ,NULL
 ,GLOBAL_prefix_base
 ,GLOBAL_flag_force_create
 ,&binary_label_
 ,&output_label__
 ,&nlpbra_label__
 ,&nlpnex_label__
 );
  /* if (GLOBAL_flag_split==1){ } */}
  /* %%%%%%%% */
  free1(&d_);
  free1(&E_base_rc__);
  free1(&binary_label_);
  free1_char__(n_r,&output_label__);
  free1_char__(n_r,&nlpbra_label__);
  free1_char__(n_r,&nlpnex_label__);
  if (verbose){ printf(" %% [finished halfloop_nonbinary_f_gateway_shell]\n");}
}
