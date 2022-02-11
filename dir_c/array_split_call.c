#ifndef _MONOLITH
#include "halfloop_header.h"
#endif /* _MONOLITH */

void array_f_from_split_f(int n_r,int n_c_split,float *A_rc__,float **B_rc_p_,float **B_cr_p_)
{
  int n_c = (n_c_split/2) + (n_c_split%2);
  int nr=0,nc=0,nc_split=0;
  float *A_r0_=NULL;
  float *A_r1_=NULL;
  float *B_rc__=NULL,*B_cr__=NULL;
  float tmp_f=0;
  B_rc__=NULL;
  if (B_rc_p_!=NULL){
    if ( (*B_rc_p_)==NULL ){ (*B_rc_p_) = (float *) malloc1((unsigned long long int)n_r*(unsigned long long int)n_c*sizeof(float)); }
    B_rc__ = *B_rc_p_;
    /* if (B_rc_p_!=NULL){ } */}
  B_cr__=NULL;
  if (B_cr_p_!=NULL){
    if ( (*B_cr_p_)==NULL ){ (*B_cr_p_) = (float *) malloc1((unsigned long long int)n_c*(unsigned long long int)n_r*sizeof(float)); }
    B_cr__ = *B_cr_p_;
    /* if (B_cr_p_!=NULL){ } */}
  nc=0;
  for (nc=0;nc<n_c;nc++){
    A_r0_ = A_rc__ + (0 + 2*(unsigned long long int)nc) * (unsigned long long int)n_r;
    A_r1_ = A_rc__ + (1 + 2*(unsigned long long int)nc) * (unsigned long long int)n_r;
    for (nr=0;nr<n_r;nr++){
      tmp_f = A_r0_[nr] - A_r1_[nr];
      if (B_rc__!=NULL){ B_rc__[nr + nc*n_r] = tmp_f;}
      if (B_cr__!=NULL){ B_cr__[nc + nr*n_c] = tmp_f;}
      /* for (nr=0;nr<n_r;nr++){ } */}
    /* for (nc=0;nc<n_c;nc++){ } */}
}

void array_split_f_from_f(int n_r,int n_c,float *A_rc__,float **B_rc_p_,float **B_cr_p_)
{
  int n_c_split = 2*n_c;
  int nr=0,nc=0,nc_split=0;
  float *A_r_=NULL;
  float *B_rc__=NULL,*B_cr__=NULL;
  float tmp_f=0;
  B_rc__=NULL;
  if (B_rc_p_!=NULL){
    if ( (*B_rc_p_)==NULL ){ (*B_rc_p_) = (float *) malloc1((unsigned long long int)n_r*(unsigned long long int)n_c_split*sizeof(float)); }
    B_rc__ = *B_rc_p_;
    /* if (B_rc_p_!=NULL){ } */}
  B_cr__=NULL;
  if (B_cr_p_!=NULL){
    if ( (*B_cr_p_)==NULL ){ (*B_cr_p_) = (float *) malloc1((unsigned long long int)n_c_split*(unsigned long long int)n_r*sizeof(float)); }
    B_cr__ = *B_cr_p_;
    /* if (B_cr_p_!=NULL){ } */}
  nc_split=0;
  for (nc=0;nc<n_c;nc++){
    A_r_ = A_rc__ + (unsigned long long int)nc * (unsigned long long int)n_r;
    for (nr=0;nr<n_r;nr++){
      tmp_f = maximum(0,+A_r_[nr]);
      if (B_rc__!=NULL){ B_rc__[nr       + nc_split*n_r] = tmp_f;}
      if (B_cr__!=NULL){ B_cr__[nc_split + nr*n_c_split] = tmp_f;}
      /* for (nr=0;nr<n_r;nr++){ } */}
    nc_split += 1;
    for (nr=0;nr<n_r;nr++){
      tmp_f = maximum(0,-A_r_[nr]);
      if (B_rc__!=NULL){ B_rc__[nr       + nc_split*n_r] = tmp_f;}
      if (B_cr__!=NULL){ B_cr__[nc_split + nr*n_c_split] = tmp_f;}
      /* for (nr=0;nr<n_r;nr++){ } */}
    nc_split += 1;
    /* for (nc=0;nc<n_c;nc++){ } */}
}

void array_split_f_from_f_test()
{
  int n_r = 5;
  int n_c = 8;
  int n_c_split = 2*n_c;
  float *A_rc__=NULL;
  float *A_cr__=NULL;
  float *B_rc__=NULL;
  float *B_cr__=NULL;
  float *C_rc__=NULL;
  float *C_cr__=NULL;
  int nr=0,nc=0;
  A_rc__ = (float *) malloc1(n_r*n_c*sizeof(float));
  A_cr__ = (float *) malloc1(n_c*n_r*sizeof(float));
  for (nr=0;nr<n_r;nr++){ for (nc=0;nc<n_c;nc++){ 
      A_rc__[nr+nc*n_r] = ((nr+1) + 0.1*(nc+1)) * ( ((nr+nc) % 2) ? -1 : +1);
      A_cr__[nc+nr*n_c] = A_rc__[nr+nc*n_r];
      /* for (nr=0;nr<n_r;nr++){ for (nc=0;nc<n_c;nc++){  }} */}};
  array_printf(A_rc__,"float",n_r,n_c," % A_rc__: ");
  /* %%%%%%%% */
  GLOBAL_tic(0);
  array_split_f_from_f(n_r,n_c,A_rc__,&B_rc__,&B_cr__);
  array_printf(B_rc__,"float",n_r,n_c_split," % B_rc__: ");
  array_printf(B_cr__,"float",n_c_split,n_r," % B_cr__: ");
  GLOBAL_toc(0,1," array_split_f_from_f: ");
  /* %%%%%%%% */
  GLOBAL_tic(0);
  array_f_from_split_f(n_r,n_c_split,B_rc__,&C_rc__,&C_cr__);
  array_printf(C_rc__,"float",n_r,n_c," % C_rc__: ");
  array_printf(C_cr__,"float",n_c,n_r," % C_cr__: ");
  GLOBAL_toc(0,1," array_split_f_from_f: ");
  /* %%%%%%%% */
  printf(" %% A_rc__ vs C_rc__: %0.16f\n",ffnormn(n_r*n_c,A_rc__,C_rc__));
  printf(" %% A_cr__ vs C_cr__: %0.16f\n",ffnormn(n_r*n_c,A_cr__,C_cr__));
  /* %%%%%%%% */
  free1(&C_rc__);
  free1(&C_cr__);
  free1(&B_rc__);
  free1(&B_cr__);
  free1(&A_rc__);
  free1(&A_cr__);
}

/* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */

void xdrop_from_xdrop_split
(
  int n_r
 ,int n_c_split
 ,double gamma
 ,int *xdrop_split_trn__
 ,int **xdrop_trn_p_
 ,int *n_iteration_p_
 ,int **rdrop_p_
 ,int **cdrop_p_
 ,int **rkeep_p_
 ,int **ckeep_p_
 )
{
  int n_c=0;
  int n_xdrop_split=0,n_xdrop=0;
  int nr=0,nc=0,nc_split=0;
  int nx=0,nx_split=0;
  int n_iteration=0,n_iteration_split=0;
  int niteration=0;
  int *rdrop_split_=NULL;
  int *cdrop_split_=NULL;
  int *rkeep_split_=NULL;
  int *ckeep_split_=NULL;
  int *flag_c_=NULL;
  int rdrop_split=0,cdrop_split=0;
  int rdrop=0,cdrop=0;
  int rkeep=0,ckeep=0;
  int n_l=0,nl=0;
  int i_r=0,i_c_split=0,i_c=0;
  n_xdrop_split = n_r + n_c_split;
  n_c = (n_c_split/2) + (n_c_split%2);
  n_xdrop = n_r + n_c;
  get_xdrop_logscale_array(n_r,n_c_split,gamma,&n_iteration_split,&rdrop_split_,&cdrop_split_,&rkeep_split_,&ckeep_split_);
  n_iteration = n_iteration_split;
  if (n_iteration_p_!=NULL){ *n_iteration_p_ = n_iteration;}
  if (xdrop_trn_p_!=NULL){ if (*xdrop_trn_p_==NULL){ (*xdrop_trn_p_) = (int *) malloc1((unsigned long long int)2*(unsigned long long int)n_xdrop*sizeof(int));}}
  if (rdrop_p_!=NULL){ if (*rdrop_p_==NULL){ (*rdrop_p_) = (int *) malloc1(n_iteration*sizeof(int));}}
  if (cdrop_p_!=NULL){ if (*cdrop_p_==NULL){ (*cdrop_p_) = (int *) malloc1(n_iteration*sizeof(int));}}
  if (rkeep_p_!=NULL){ if (*rkeep_p_==NULL){ (*rkeep_p_) = (int *) malloc1(n_iteration*sizeof(int));}}
  if (ckeep_p_!=NULL){ if (*ckeep_p_==NULL){ (*ckeep_p_) = (int *) malloc1(n_iteration*sizeof(int));}}
  flag_c_ = (int *) malloc1((unsigned long long int)n_c*sizeof(int));
  nx_split=0; nx=0; rkeep = n_r; ckeep = n_c;
  for (niteration=0;niteration<n_iteration;niteration++){
    rdrop_split = rdrop_split_[niteration]; cdrop_split = cdrop_split_[niteration]; n_l = rdrop_split + cdrop_split;
    rdrop=0; cdrop=0;
    for (nl=0;nl<n_l;nl++){
      i_r = xdrop_split_trn__[ 0 + nx_split*2 ];
      i_c_split = xdrop_split_trn__[ 1 + nx_split*2 ];
      if ( (i_r>=0) && (i_c_split< 0) ){
	rdrop += 1;
	if (xdrop_trn_p_!=NULL){ (*xdrop_trn_p_)[ 0 + nx*2 ] = i_r; (*xdrop_trn_p_)[ 1 + nx*2 ] = -1;}
	nx += 1;
	nx_split += 1;
	/* if ( (i_r>=0) && (i_c_split< 0) ){ } */}
      if ( (i_r< 0) && (i_c_split>=0) ){
	i_c = i_c_split/2;
	flag_c_[i_c] += 1;
	if (flag_c_[i_c]< 2){ /* do nothing */}
	if (flag_c_[i_c]==2){ if (xdrop_trn_p_!=NULL){ (*xdrop_trn_p_)[ 0 + nx*2 ] = -1; (*xdrop_trn_p_)[ 1 + nx*2 ] = i_c;} cdrop += 1; nx += 1;}
	if (flag_c_[i_c]> 2){ printf(" %% Warning, flag_c> 2 in xdrop_from_xdrop_split\n");}
	nx_split += 1;
	/* if ( (i_r< 0) && (i_c_split>=0) ){ } */}
      /* for (nl=0;nl<n_l;nl++){ } */}
      if (rdrop_p_!=NULL){ (*rdrop_p_)[niteration] = rdrop;}
      if (cdrop_p_!=NULL){ (*cdrop_p_)[niteration] = cdrop;}
      rkeep-=rdrop; ckeep-=cdrop;
      if (rkeep_p_!=NULL){ (*rkeep_p_)[niteration] = rkeep;}
      if (ckeep_p_!=NULL){ (*ckeep_p_)[niteration] = ckeep;}
      /* for (niteration=0;niteration<n_iteration;niteration++){ } */}
  free1(&flag_c_);
  free1(&rdrop_split_);
  free1(&cdrop_split_);
  free1(&rkeep_split_);
  free1(&ckeep_split_);
}

void xdrop_from_xdrop_split_test()
{
  int n_r = 8;
  int n_c = 7;
  int n_c_split = 2*n_c;
  int nr=0,nc=0,nc_split=0,nx_split=0,nl=0;
  int *index_r_=NULL,*index_c_split_=NULL;
  unsigned long long int rseed=67892;
  double gamma = 0.25;
  int n_xdrop = n_r + n_c;
  int n_xdrop_split = n_r + n_c_split;
  int niteration_split=0,n_iteration_split=0;
  int *xdrop_0_split_trn__=NULL;
  int *xdrop_1_split_trn__=NULL;
  int *rdrop_split_=NULL;  
  int *cdrop_split_=NULL;  
  int *rkeep_split_=NULL;  
  int *ckeep_split_=NULL;  
  int rdrop_split=0,cdrop_split=0;
  int n_iteration=0;
  int *xdrop_0_trn__=NULL;
  int *rdrop_0_=NULL;  
  int *cdrop_0_=NULL;  
  int *rkeep_0_=NULL;  
  int *ckeep_0_=NULL;
  int *xdrop_1_trn__=NULL;
  int *rdrop_1_=NULL;  
  int *cdrop_1_=NULL;  
  int *rkeep_1_=NULL;  
  int *ckeep_1_=NULL;
  irandperm(n_r,&index_r_,&rseed);
  irandperm(n_c_split,&index_c_split_,&rseed);
  /* %%%%%%%% */
  get_xdrop_logscale_array(n_r,n_c_split,gamma,&n_iteration_split,&rdrop_split_,&cdrop_split_,&rkeep_split_,&ckeep_split_);
  xdrop_0_split_trn__ = (int *) malloc1((unsigned long long int)2*(unsigned long long int)n_xdrop_split*sizeof(int));
  xdrop_1_split_trn__ = (int *) malloc1((unsigned long long int)2*(unsigned long long int)n_xdrop_split*sizeof(int));
  nr=0;nc_split=0;nx_split=0;
  for (niteration_split=0;niteration_split<n_iteration_split;niteration_split++){
    rdrop_split = rdrop_split_[niteration_split];
    cdrop_split = cdrop_split_[niteration_split];
    printf(" %% rdrop_split %d cdrop_split %d\n",rdrop_split,cdrop_split);
    for (nl=0;nl<rdrop_split;nl++){ 
      xdrop_0_split_trn__[ 0 + nx_split*2 ] = index_r_[nr];
      xdrop_0_split_trn__[ 1 + nx_split*2 ] = -1;
      xdrop_1_split_trn__[ 0 + nx_split*2 ] = index_r_[nr];
      xdrop_1_split_trn__[ 1 + nx_split*2 ] = -1;
      nx_split += 1; nr += 1;
      /* for (nl=0;nl<rdrop_split;nl++){ } */}
    for (nl=0;nl<cdrop_split;nl++){ 
      xdrop_0_split_trn__[ 0 + nx_split*2 ] = -1;
      xdrop_0_split_trn__[ 1 + nx_split*2 ] = index_c_split_[nc_split];
      xdrop_1_split_trn__[ 0 + nx_split*2 ] = -1;
      xdrop_1_split_trn__[ 1 + nx_split*2 ] = 2*(index_c_split_[nc_split]/2) + (1-(index_c_split_[nc_split]%2));
      nx_split += 1; nc_split += 1;
      /* for (nl=0;nl<cdrop_split;nl++){ } */}
    /* for (niteration_split=0;niteration_split<n_iteration_split;niteration_split++){ } */}
  array_printf(xdrop_0_split_trn__,"int",2,n_xdrop_split," %% xdrop_0_split_trn__: ");
  array_printf(xdrop_1_split_trn__,"int",2,n_xdrop_split," %% xdrop_1_split_trn__: ");
  array_printf(rdrop_split_,"int",1,n_iteration_split," %% rdrop_split_: ");
  array_printf(cdrop_split_,"int",1,n_iteration_split," %% cdrop_split_: ");
  array_printf(rkeep_split_,"int",1,n_iteration_split," %% rkeep_split_: ");
  array_printf(ckeep_split_,"int",1,n_iteration_split," %% ckeep_split_: ");
  /* %%%%%%%% */
  xdrop_from_xdrop_split
    (
      n_r
     ,n_c_split
     ,gamma
     ,xdrop_0_split_trn__
     ,&xdrop_0_trn__
     ,&n_iteration
     ,&rdrop_0_
     ,&cdrop_0_
     ,&rkeep_0_
     ,&ckeep_0_
     );
  array_printf(xdrop_0_trn__,"int",2,n_xdrop," %% xdrop_0_trn__: ");
  array_printf(rdrop_0_,"int",1,n_iteration," %% rdrop_0_: ");
  array_printf(cdrop_0_,"int",1,n_iteration," %% cdrop_0_: ");
  array_printf(rkeep_0_,"int",1,n_iteration," %% rkeep_0_: ");
  array_printf(ckeep_0_,"int",1,n_iteration," %% ckeep_0_: ");
  /* %%%%%%%% */
  xdrop_from_xdrop_split
    (
      n_r
     ,n_c_split
     ,gamma
     ,xdrop_1_split_trn__
     ,&xdrop_1_trn__
     ,&n_iteration
     ,&rdrop_1_
     ,&cdrop_1_
     ,&rkeep_1_
     ,&ckeep_1_
     );
  array_printf(xdrop_1_trn__,"int",2,n_xdrop," %% xdrop_1_trn__: ");
  array_printf(rdrop_1_,"int",1,n_iteration," %% rdrop_1_: ");
  array_printf(cdrop_1_,"int",1,n_iteration," %% cdrop_1_: ");
  array_printf(rkeep_1_,"int",1,n_iteration," %% rkeep_1_: ");
  array_printf(ckeep_1_,"int",1,n_iteration," %% ckeep_1_: ");
  /* %%%%%%%% */
  printf(" %% xdrop_0_trn__ vs xdrop_1_trn__: %0.16f\n",ifnormn(n_xdrop*2,xdrop_0_trn__,xdrop_1_trn__));
  printf(" %% rdrop_0_ vs rdrop_1_: %0.16f\n",ifnormn(n_iteration,rdrop_0_,rdrop_1_));
  printf(" %% cdrop_0_ vs cdrop_1_: %0.16f\n",ifnormn(n_iteration,cdrop_0_,cdrop_1_));
  printf(" %% rkeep_0_ vs rkeep_1_: %0.16f\n",ifnormn(n_iteration,rkeep_0_,rkeep_1_));
  printf(" %% ckeep_0_ vs ckeep_1_: %0.16f\n",ifnormn(n_iteration,ckeep_0_,ckeep_1_));
  /* %%%%%%%% */
  free1(&rdrop_0_);
  free1(&cdrop_0_);
  free1(&rkeep_0_);
  free1(&ckeep_0_);
  free1(&xdrop_0_trn__);
  free1(&rdrop_1_);
  free1(&cdrop_1_);
  free1(&rkeep_1_);
  free1(&ckeep_1_);
  free1(&xdrop_1_trn__);
  free1(&rdrop_split_);
  free1(&cdrop_split_);
  free1(&rkeep_split_);
  free1(&ckeep_split_);
  free1(&xdrop_0_split_trn__);
  free1(&xdrop_1_split_trn__);
  free1(&index_r_);
  free1(&index_c_split_);
}

/* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */

