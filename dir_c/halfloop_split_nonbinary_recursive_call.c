#ifndef _MONOLITH
#include "halfloop_header.h"
#endif /* _MONOLITH */

void trace_trn_from_data(int n_iteration,int n_r,int n_c,int *rdrop_,int *cdrop_,double *QR_,double *QC_,double **trace_p_)
{
  int niteration=0;
  int n_r_rem=0,n_c_rem=0;
  int rdrop=0,cdrop=0;
  double *trace_trn__=NULL;
  double QR=0,QC=0;
  if (trace_p_!=NULL){
    if ( (*trace_p_)==NULL ){ (*trace_p_) = (double *) malloc1((unsigned long long int)6*(unsigned long long int)n_iteration*sizeof(double)); }
    trace_trn__ = *trace_p_;
    /* if (trace_p_!=NULL){ } */}
  if (trace_trn__!=NULL){
    n_r_rem = n_r; n_c_rem = n_c;
    for (niteration=0;niteration<n_iteration;niteration++){
      rdrop = rdrop_[niteration];
      cdrop = cdrop_[niteration];
      QR = QR_[niteration];
      QC = QC_[niteration];
      trace_trn__[0+6*niteration] = 1+niteration;
      trace_trn__[1+6*niteration] = n_r_rem;
      trace_trn__[2+6*niteration] = n_c_rem;
      trace_trn__[3+6*niteration] = QR;
      trace_trn__[4+6*niteration] = QC;
      trace_trn__[5+6*niteration] = 1.0;
      n_r_rem -= rdrop;
      n_c_rem -= cdrop;
      /* for (niteration=0;niteration<n_iteration;niteration++){ } */}
    /* if (trace_trn__!=NULL){ } */}
}

void halfloop_split_nonbinary_f_recursive_helper_QR__
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
  int trace_QC_index = 4;
  int n_iteration=0,tmp_n_iteration=0;
  int nshuffle=0;
  unsigned long int rseed=0;
  double *trace__=NULL;
  double *trace_shuffle__=NULL;
  double *QR_=NULL;
  double *QR__=NULL;
  double *QC_=NULL;
  double *QC__=NULL;
  int n_xdrop=0;
  int *xdrop__=NULL;
  int *xdrop_shuffle__=NULL;
  float *E_cr__=NULL;
  float *E_rc__=NULL;
  float *QE_rc__=NULL;
  int MDA_n_dim=0;
  int *MDA_dim_=NULL;
  char MDA_fname[PNAMESIZE];
  float *E_split_cr__=NULL;
  float *E_split_rc__=NULL;
  int n_c_split_index=0;
  double *trace_split__=NULL;
  double *trace_split_shuffle__=NULL;
  int n_xdrop_split=0;
  int *xdrop_split__=NULL;
  int *xdrop_split_shuffle__=NULL;
  int *rdrop_=NULL;
  int *cdrop_=NULL;
  int *rkeep_=NULL;
  int *ckeep_=NULL;
  if (verbose){ printf(" %% [entering halfloop_split_nonbinary_f_recursive_helper_QR__]\n");}
  flag_not_exist = (stat(fname_trace__,&stat_file)==-1) || (stat(fname_QR__,&stat_file)==-1) || (stat(fname_xdrop__,&stat_file)==-1);
  if ( flag_force_create || flag_not_exist ){
    if (verbose){ printf(" %% %s not found, creating\n",fname_trace__);}
    MDA_dim_ = (int *) malloc1(2*sizeof(int));
    /* %%%%%%%% */
    n_c_split_index = 2*n_c_index;
    n_iteration = get_xdrop_logscale_length(n_r_index,flag_r0drop_vs_rcdrop==flag_rcdrop ? n_c_split_index : maximum(n_r_index,n_c_split_index),gamma);
    trace_split__ = (double *) malloc1((unsigned long long int)6*(unsigned long long int)n_iteration*sizeof(double));
    trace_split_shuffle__ = (double *) malloc1((unsigned long long int)6*(unsigned long long int)n_iteration*sizeof(double));
    if (flag_r0drop_vs_rcdrop==flag_rcdrop){ n_xdrop_split = n_r_index + n_c_split_index;}
    if (flag_r0drop_vs_rcdrop==flag_r0drop){ n_xdrop_split = n_r_index;}
    xdrop_split__ = (int *) malloc1((unsigned long long int)2*(unsigned long long int)n_xdrop_split*sizeof(double));
    xdrop_split_shuffle__ = (int *) malloc1((unsigned long long int)2*(unsigned long long int)n_xdrop_split*sizeof(double));
    /* %%%%%%%% */    
    trace__ = (double *) malloc1((unsigned long long int)6*(unsigned long long int)n_iteration*sizeof(double));
    trace_shuffle__ = (double *) malloc1((unsigned long long int)6*(unsigned long long int)n_iteration*sizeof(double));
    QR__ = (double *) malloc1((unsigned long long int)n_iteration*(unsigned long long int)(1+n_shuffle)*sizeof(double));
    QC__ = (double *) malloc1((unsigned long long int)n_iteration*(unsigned long long int)(1+n_shuffle)*sizeof(double));
    if (flag_r0drop_vs_rcdrop==flag_rcdrop){ n_xdrop = n_r_index + n_c_index;}
    if (flag_r0drop_vs_rcdrop==flag_r0drop){ n_xdrop = n_r_index;}
    xdrop__ = (int *) malloc1((unsigned long long int)2*(unsigned long long int)n_xdrop*sizeof(double));
    xdrop_shuffle__ = (int *) malloc1((unsigned long long int)2*(unsigned long long int)n_xdrop*sizeof(double));
    /* %%%%%%%% */
    E_rc__ = (float *) malloc1((unsigned long long int)n_r_index*(unsigned long long int)n_c_index*sizeof(float));
    E_cr__ = (float *) malloc1((unsigned long long int)n_c_index*(unsigned long long int)n_r_index*sizeof(float));
    E_split_rc__ = (float *) malloc1((unsigned long long int)n_r_index*(unsigned long long int)n_c_split_index*sizeof(float));
    E_split_cr__ = (float *) malloc1((unsigned long long int)n_c_split_index*(unsigned long long int)n_r_index*sizeof(float));
    array_extract_f_from_f(n_r,n_c,E_base_rc__,n_r_index,r_index_,n_c_index,c_index_,&E_rc__,NULL);
    if (verbose>9){ array_printf_margin(E_rc__,"float",n_r_index,n_c_index," % E_rc__: "); printf(" %% %% %% %%\n");}
    array_mean_center_row(n_r_index,n_c_index,E_rc__,NULL,"float",&E_rc__,&E_cr__);
    if (verbose>9){ array_printf_margin(E_rc__,"float",n_r_index,n_c_index," % E_rc__: "); printf(" %% %% %% %%\n");}
    if (verbose>9){ array_printf_margin(E_cr__,"float",n_c_index,n_r_index," % E_cr__: "); printf(" %% %% %% %%\n");}
    array_split_f_from_f(n_r_index,n_c_index,E_rc__,&E_split_rc__,&E_split_cr__);
    if (verbose>9){ array_printf_margin(E_split_rc__,"float",n_r_index,n_c_split_index," % E_split_rc__: "); printf(" %% %% %% %%\n");}
    if (verbose>9){ array_printf_margin(E_split_cr__,"float",n_c_split_index,n_r_index," % E_split_cr__: "); printf(" %% %% %% %%\n");}
    if (flag_r0drop_vs_rcdrop==flag_rcdrop){ halfloop_nonbinary_f(n_r_index,n_c_split_index,E_split_rc__,E_split_cr__,gamma,&xdrop_split__,&tmp_n_iteration,&trace_split__);}
    if (flag_r0drop_vs_rcdrop==flag_r0drop){ halfloop_nonbinary_rdrop_f(n_r_index,n_c_split_index,E_split_rc__,E_split_cr__,gamma,&xdrop_split__,&tmp_n_iteration,&trace_split__);}
    if (verbose>9){ array_printf_margin(xdrop_split__,"int",2,n_xdrop_split," % xdrop_split__: "); printf(" %% %% %% %%\n");}
    if (verbose>9){ array_printf_margin(trace_split__,"double",6,n_iteration," % trace_split__: "); printf(" %% %% %% %%\n");}
    if (tmp_n_iteration!=n_iteration){ printf(" %% Warning, tmp_n_iteration %d not equal n_iteration %d\n",tmp_n_iteration,n_iteration);}
    QR_ = QR__ + (unsigned long long int)0*(unsigned long long int)n_iteration;
    array_extract_d_from_d(6,n_iteration,trace_split__,1,&trace_QR_index,0,NULL,&QR_,NULL);
    QC_ = QC__ + (unsigned long long int)0*(unsigned long long int)n_iteration;
    array_extract_d_from_d(6,n_iteration,trace_split__,1,&trace_QC_index,0,NULL,&QC_,NULL);
    rdrop_=NULL;cdrop_=NULL;rkeep_=NULL;ckeep_=NULL;
    xdrop_from_xdrop_split
      (
        n_r_index
       ,n_c_split_index
       ,gamma
       ,xdrop_split__
       ,&xdrop__
       ,&tmp_n_iteration
       ,&rdrop_
       ,&cdrop_
       ,&rkeep_
       ,&ckeep_
       );
    trace_trn_from_data(n_iteration,n_r_index,n_c_index,rdrop_,cdrop_,QR_,QC_,&trace__);
    if (verbose>9){ array_printf_margin(xdrop__,"int",2,n_xdrop," % xdrop__: "); printf(" %% %% %% %%\n");}
    if (verbose>9){ array_printf_margin(trace__,"double",6,n_iteration," % trace__: "); printf(" %% %% %% %%\n");}
    if (tmp_n_iteration!=n_iteration){ printf(" %% Warning, tmp_n_iteration %d not equal n_iteration %d\n",tmp_n_iteration,n_iteration);}
    free1(&rdrop_); free1(&cdrop_); free1(&rkeep_); free1(&ckeep_);
    sprintf(MDA_fname,"%s",fname_trace__); MDA_n_dim = 2; MDA_dim_[0] = 6; MDA_dim_[1] = n_iteration;
    MDA_write_r8(MDA_n_dim,MDA_dim_,trace__,MDA_fname);
    if (verbose>9){ MDA_printf_r8_margin(MDA_fname);}
    sprintf(MDA_fname,"%s",fname_xdrop__); MDA_n_dim = 2; MDA_dim_[0] = 2; MDA_dim_[1] = n_xdrop;
    MDA_write_i4(MDA_n_dim,MDA_dim_,xdrop__,MDA_fname);
    if (verbose>9){ MDA_printf_i4_margin(MDA_fname);}
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
      array_split_f_from_f(n_r_index,n_c_index,E_rc__,&E_split_rc__,&E_split_cr__);
      if (verbose>9){ array_printf_margin(E_split_rc__,"float",n_r_index,n_c_split_index," % E_split_rc__: "); printf(" %% %% %% %%\n");}
      if (flag_r0drop_vs_rcdrop==flag_rcdrop){ halfloop_nonbinary_f(n_r_index,n_c_split_index,E_split_rc__,E_split_cr__,gamma,&xdrop_split_shuffle__,&tmp_n_iteration,&trace_split_shuffle__);}
      if (flag_r0drop_vs_rcdrop==flag_r0drop){ halfloop_nonbinary_rdrop_f(n_r_index,n_c_split_index,E_split_rc__,E_split_cr__,gamma,&xdrop_split_shuffle__,&tmp_n_iteration,&trace_split_shuffle__);}
      if (verbose>9){ array_printf_margin(xdrop_split_shuffle__,"int",2,n_xdrop_split," % xdrop_split_shuffle__: "); printf(" %% %% %% %%\n");}
      if (verbose>9){ array_printf_margin(trace_split_shuffle__,"double",6,n_iteration," % trace_split_shuffle__: "); printf(" %% %% %% %%\n");}
      if (tmp_n_iteration!=n_iteration){ printf(" %% Warning, tmp_n_iteration %d not equal n_iteration %d\n",tmp_n_iteration,n_iteration);}
      QR_ = QR__ + (unsigned long long int)(1+nshuffle)*(unsigned long long int)n_iteration;
      array_extract_d_from_d(6,n_iteration,trace_split_shuffle__,1,&trace_QR_index,0,NULL,&QR_,NULL);
      if (verbose>9){ array_printf_margin(QR_,"double",1,n_iteration," % QR_: ");}
      QC_ = QC__ + (unsigned long long int)(1+nshuffle)*(unsigned long long int)n_iteration;
      array_extract_d_from_d(6,n_iteration,trace_split_shuffle__,1,&trace_QC_index,0,NULL,&QC_,NULL);
      if (verbose>9){ array_printf_margin(QC_,"double",1,n_iteration," % QC_: ");}
      /* for (nshuffle=0;nshuffle<n_shuffle;nshuffle++){ } */}
    if (verbose>9){ array_printf_margin(QR__,"double",n_iteration,(1+n_shuffle)," % QR__: ");}
    sprintf(MDA_fname,"%s",fname_QR__); MDA_n_dim = 2; MDA_dim_[0] = n_iteration; MDA_dim_[1] = (1+n_shuffle);
    MDA_write_r8(MDA_n_dim,MDA_dim_,QR__,MDA_fname);
    if (verbose>9){ MDA_printf_r8_margin(MDA_fname);}
    free1(&QE_rc__);
    free1(&E_cr__);
    free1(&E_rc__);
    free1(&E_split_cr__);
    free1(&E_split_rc__);
    free1(&trace_split__);
    free1(&trace_split_shuffle__);
    free1(&xdrop_split__);
    free1(&xdrop_split_shuffle__);
    free1(&trace__);
    free1(&trace_shuffle__);
    free1(&QR__);
    free1(&QC__);
    free1(&xdrop__);
    free1(&xdrop_shuffle__);
    free1(&MDA_dim_);
    /* not found */}
  else{ if (verbose){ printf(" %% %s found, not creating\n",fname_trace__);}}
  if (verbose){ printf(" %% [finished halfloop_split_nonbinary_f_recursive_helper_QR__]\n");}
}

/* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */

void halfloop_split_nonbinary_f_recursive_omp_helper_QR_
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
 ,double *xdrop__
)
{
  int flag_rcdrop = 0;
  int flag_r0drop = 1;
  int trace_QR_index = 3;
  int trace_QC_index = 4;
  int n_c_split_index = 2*n_c_index;
  int n_xdrop_split = 0;
  int n_xdrop = 0;
  int *xdrop_split__=NULL;
  double *trace_split__=NULL;
  int tmp_n_iteration=0;
  double *trace_split_shuffle__=NULL;
  int *xdrop_split_shuffle__=NULL;
  unsigned long int rseed=0;
  float *tmp_E_split_cr__=NULL;
  float *tmp_E_split_rc__=NULL;
  float *tmp_E_cr__=NULL;
  float *tmp_E_rc__=NULL;
  float *tmp_QE_rc__=NULL;
  double *QC_ = NULL;
  int *rdrop_=NULL,*cdrop_=NULL;
  int *rkeep_=NULL,*ckeep_=NULL;
  if (verbose){ printf(" %% [entering halfloop_split_nonbinary_f_recursive_omp_helper_QR_]\n");}
  tmp_E_split_rc__ = (float *) malloc1((unsigned long long int)n_r_index*(unsigned long long int)n_c_split_index*sizeof(float));
  tmp_E_split_cr__ = (float *) malloc1((unsigned long long int)n_c_split_index*(unsigned long long int)n_r_index*sizeof(float));
  if (flag_r0drop_vs_rcdrop==flag_rcdrop){ n_xdrop_split = n_r_index + n_c_split_index;}
  if (flag_r0drop_vs_rcdrop==flag_r0drop){ n_xdrop_split = n_r_index;}
  if (nshuffle==0){
    trace_split__ = (double *) malloc1((unsigned long long int)6*(unsigned long long int)n_iteration*sizeof(double));
    xdrop_split__ = (int *) malloc1((unsigned long long int)2*(unsigned long long int)n_xdrop_split*sizeof(double));
    QC_ = (double *) malloc1((unsigned long long int)n_iteration*sizeof(double));
    array_split_f_from_f(n_r_index,n_c_index,E_rc__,&tmp_E_split_rc__,&tmp_E_split_cr__);
    if (flag_r0drop_vs_rcdrop==flag_rcdrop){ halfloop_nonbinary_f(n_r_index,n_c_split_index,tmp_E_split_rc__,tmp_E_split_cr__,gamma,&xdrop_split__,&tmp_n_iteration,&trace_split__);}
    if (flag_r0drop_vs_rcdrop==flag_r0drop){ halfloop_nonbinary_rdrop_f(n_r_index,n_c_split_index,tmp_E_split_rc__,tmp_E_split_cr__,gamma,&xdrop_split__,&tmp_n_iteration,&trace_split__);}
    if (tmp_n_iteration!=n_iteration){ printf(" %% Warning, tmp_n_iteration %d not equal n_iteration %d\n",tmp_n_iteration,n_iteration);}
    array_extract_d_from_d(6,n_iteration,trace_split__,1,&trace_QR_index,0,NULL,&QR_,NULL);
    array_extract_d_from_d(6,n_iteration,trace_split__,1,&trace_QR_index,0,NULL,&QR_,NULL);
    rdrop_=NULL;cdrop_=NULL;rkeep_=NULL;ckeep_=NULL;
    xdrop_from_xdrop_split
      (
        n_r_index
       ,n_c_split_index
       ,gamma
       ,xdrop_split__
       ,&xdrop__
       ,&tmp_n_iteration
       ,&rdrop_
       ,&cdrop_
       ,&rkeep_
       ,&ckeep_
       );
    trace_trn_from_data(n_iteration,n_r_index,n_c_index,rdrop_,cdrop_,QR_,QC_,&trace__);
    free1(&rdrop_);
    free1(&cdrop_);
    free1(&rkeep_);
    free1(&ckeep_);
    free1(&QC_);
    free1(&trace_split__);
    free1(&xdrop_split__);
    /* if (nshuffle==0){ } */}
  if (nshuffle> 0){
    trace_split_shuffle__ = (double *) malloc1((unsigned long long int)6*(unsigned long long int)n_iteration*sizeof(double));
    xdrop_split_shuffle__ = (int *) malloc1((unsigned long long int)2*(unsigned long long int)n_xdrop_split*sizeof(double));
    tmp_E_rc__ = (float *) malloc1((unsigned long long int)n_r_index*(unsigned long long int)n_c_index*sizeof(float));
    tmp_E_cr__ = (float *) malloc1((unsigned long long int)n_c_index*(unsigned long long int)n_r_index*sizeof(float));
    tmp_QE_rc__ = (float *) malloc1((unsigned long long int)n_r_index*(unsigned long long int)n_c_index*sizeof(float));
    rseed = (0+nshuffle); RSEED_adv8(&rseed);
    halfloop_nonbinary_f_recursive_helper_QR_helper_orth(verbose,n_r_index,n_c_index,E_rc__,tmp_QE_rc__,&rseed);
    array_mean_center_row(n_r_index,n_c_index,tmp_QE_rc__,NULL,"float",&tmp_E_rc__,&tmp_E_cr__);
    array_split_f_from_f(n_r_index,n_c_index,tmp_E_rc__,&tmp_E_split_rc__,&tmp_E_split_cr__);
    if (flag_r0drop_vs_rcdrop==flag_rcdrop){ halfloop_nonbinary_f(n_r_index,n_c_split_index,tmp_E_split_rc__,tmp_E_split_cr__,gamma,&xdrop_split_shuffle__,&tmp_n_iteration,&trace_split_shuffle__);}
    if (flag_r0drop_vs_rcdrop==flag_r0drop){ halfloop_nonbinary_rdrop_f(n_r_index,n_c_split_index,tmp_E_split_rc__,tmp_E_split_cr__,gamma,&xdrop_split_shuffle__,&tmp_n_iteration,&trace_split_shuffle__);}
    array_extract_d_from_d(6,n_iteration,trace_split_shuffle__,1,&trace_QR_index,0,NULL,&QR_,NULL);
    free1(&trace_split_shuffle__);
    free1(&xdrop_split_shuffle__);
    free1(&tmp_E_rc__);
    free1(&tmp_E_cr__);
    free1(&tmp_QE_rc__);
    /* if (nshuffle> 0){ } */}
  free1(&tmp_E_split_rc__);
  free1(&tmp_E_split_cr__);
  if (verbose){ printf(" %% [finished halfloop_split_nonbinary_f_recursive_omp_helper_QR_]\n");}
}
 
void halfloop_split_nonbinary_f_recursive_omp_helper_QR__
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
  int n_c_split_index=0;
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
  if (verbose){ printf(" %% [entering halfloop_split_nonbinary_f_recursive_omp_helper_QR__]\n");}
  flag_not_exist = (stat(fname_trace__,&stat_file)==-1) || (stat(fname_QR__,&stat_file)==-1) || (stat(fname_xdrop__,&stat_file)==-1);
  if ( flag_force_create || flag_not_exist ){
    if (verbose){ printf(" %% %s not found, creating\n",fname_trace__);}
    MDA_dim_ = (int *) malloc1(2*sizeof(int));
    n_c_split_index = 2*n_c_index;
    n_iteration = get_xdrop_logscale_length(n_r_index,flag_r0drop_vs_rcdrop==flag_rcdrop ? n_c_split_index : maximum(n_r_index,n_c_split_index),gamma);
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
	halfloop_split_nonbinary_f_recursive_omp_helper_QR_
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
  if (verbose){ printf(" %% [finished halfloop_split_nonbinary_f_recursive_omp_helper_QR__]\n");}
}

/* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */

void halfloop_split_nonbinary_f_recursive
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
  if (verbose>0){ printf(" %% [entering halfloop_split_nonbinary_f_recursive], flag_omp %d\n",flag_omp);}
  /* %%%%%%%%%%%%%%%% */
  n_str_path = pathconf(".",_PC_PATH_MAX);
  MDA_dim_ = (int *) malloc1(2*sizeof(int)); //%<-- should only ever need 2. ;
  nlp_p_split_ = (double *) malloc1((1+14)*sizeof(double)); //%<-- gumb_emp, gumb_opt, F, E, p_branch, p_next, p_set, p_use, n_r_rtn_index_F, n_r_rmv_index_F, n_member_lob, ndepth, recursion_limit, flag_split. ;
  cwd = (char *) malloc1((size_t)n_str_path); getcwd(cwd,n_str_path);
  if (n_str_path> FNAMESIZE){ printf(" %% Warning, n_str_path %d in halfloop_split_nonbinary_f_recursive\n",n_str_path);}
  if ( (dir_trunk_0in!=NULL) && (strlen(dir_trunk_0in)> FNAMESIZE) ){ printf(" %% Warning, dir_trunk_0in %s too long in halfloop_split_nonbinary_f_recursive\n",dir_trunk_0in);}
  if ( (dir_out_0in!=NULL) && (strlen(dir_out_0in)> FNAMESIZE) ){ printf(" %% Warning, dir_out_0in %s too long in halfloop_split_nonbinary_f_recursive\n",dir_out_0in);}
  if ( (prefix_base_0in!=NULL) && (strlen(prefix_base_0in)> FNAMESIZE) ){ printf(" %% Warning, prefix_base_0in %s too long in halfloop_split_nonbinary_f_recursive\n",prefix_base_0in);}
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
      halfloop_split_nonbinary_f_recursive_helper_QR__
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
      GLOBAL_toc(1,verbose_t," % halfloop_split_nonbinary_f_recursive_helper_QR__: ");
    /* if (flag_omp==0){ } */}
    if (flag_omp==1){
      GLOBAL_tic(1);
      halfloop_split_nonbinary_f_recursive_omp_helper_QR__
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
      GLOBAL_toc(1,verbose_t," % halfloop_split_nonbinary_f_recursive_omp_helper_QR__: ");
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
      halfloop_split_nonbinary_f_recursive_helper_QR__
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
      GLOBAL_toc(1,verbose_t," % halfloop_split_nonbinary_f_recursive_helper_QR__: ");
      /* if (flag_omp==0){ } */}
    if (flag_omp==1){
      GLOBAL_tic(1);
      halfloop_split_nonbinary_f_recursive_omp_helper_QR__
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
      GLOBAL_toc(1,verbose_t," % halfloop_split_nonbinary_f_recursive_omp_helper_QR__: ");
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
      halfloop_split_nonbinary_f_recursive_helper_QR__
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
      GLOBAL_toc(1,verbose_t," % halfloop_split_nonbinary_f_recursive_helper_QR__: ");
      /* if (flag_omp==0){ } */}
    if (flag_omp==1){
      GLOBAL_tic(1);
      halfloop_split_nonbinary_f_recursive_omp_helper_QR__
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
      GLOBAL_toc(1,verbose_t," % halfloop_split_nonbinary_f_recursive_omp_helper_QR__: ");
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
    halfloop_split_nonbinary_f_recursive
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
    halfloop_split_nonbinary_f_recursive
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
  if (verbose>0){ printf(" %% [finished halfloop_split_nonbinary_f_recursive], flag_omp %d\n",flag_omp);}
}

void halfloop_split_nonbinary_f_recursive_test_speed()
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
  halfloop_split_nonbinary_f_recursive(
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
 ,"test_split_speed_rc"
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
  GLOBAL_toc(0,1," % halfloop_split_nonbinary_f_recursive: ");
  /* %%%%%%%% */
  free1(&E_base_rc__);
  free1(&r_index_);
  free1(&c_index_);
}

/* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
