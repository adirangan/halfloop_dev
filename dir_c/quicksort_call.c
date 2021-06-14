#ifndef _MONOLITH
#include "halfloop_header.h"
#endif /* _MONOLITH */

/* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */

int iquicksort_partition_index_bkp(int *i_,int stride,int *index_,int l,int r) 
{
  int pivot=0,tmpl=0;
  int i=0,j=0,tmpi=0;
  pivot = i_[stride*l];
  i = l; j = r+1;
  do{
    do{ i++;} while( i<=r && i_[stride*i] <= pivot );
    do{ j--;} while( i_[stride*j] > pivot );
    if( i >= j ) break;
    tmpl = i_[stride*i]; i_[stride*i] = i_[stride*j]; i_[stride*j] = tmpl;
    if (index_!=NULL){ tmpi = index_[i]; index_[i] = index_[j]; index_[j] = tmpi;}
  }while(1);
  tmpl = i_[stride*l]; i_[stride*l] = i_[stride*j]; i_[stride*j] = tmpl;
  if (index_!=NULL){ tmpi = index_[l]; index_[l] = index_[j]; index_[j] = tmpi;}
  return j;
}

int iquicksort_partition_index(int *i_full_,int stride,int *index_full_,int l,int r)
{
  int verbose=0;
  int *i_ = i_full_ + stride*l;
  int *index_ = index_full_ + l;
  int *g_wkspace_=NULL;
  int *index_wkspace_ = NULL;
  int n_i = 1+r-l;
  int m = n_i/2;
  int i_pivot=0;
  int index_pivot=0;
  int g_index_lower=0,g_index_upper=0;
  int ni=0;
  int flag_alt=0;
  int i=0;
  int index=0;
  int m_upd=0;
  if (verbose){ printf(" %% [entering iquicksort_partition_index]\n");}
  if (verbose){ array_printf(i_,"int",1,stride*n_i," %% pre i_: ");}
  if (verbose){ array_printf(index_,"int",1,n_i," %% pre index_: ");}
  g_wkspace_ = (int *) malloc1((unsigned long long int)(2*n_i)*sizeof(int));
  index_wkspace_ = (int *) malloc1((unsigned long long int)(2*n_i)*sizeof(int));
  i_pivot = i_[stride*m];
  index_pivot = index_[m];
  if (verbose){ printf(" %% m %d, i_pivot %d, index_pivot %d\n",m,i_pivot,index_pivot);}
  g_wkspace_[n_i] = i_pivot;
  index_wkspace_[n_i] = index_pivot;
  g_index_lower = n_i-1;
  g_index_upper = n_i+1;
  flag_alt = 0;
  for (ni=0;ni<n_i;ni++){
    if (ni!=m){
      i = i_[stride*ni];
      index = index_[ni];
      if (i< i_pivot){ g_wkspace_[g_index_lower] = i; index_wkspace_[g_index_lower] = index; g_index_lower--;}
      if (i> i_pivot){ g_wkspace_[g_index_upper] = i; index_wkspace_[g_index_upper] = index; g_index_upper++;}
      if (i==i_pivot){
	if (flag_alt==0){ g_wkspace_[g_index_lower] = i; index_wkspace_[g_index_lower] = index; g_index_lower--;}
	if (flag_alt==1){ g_wkspace_[g_index_upper] = i; index_wkspace_[g_index_upper] = index; g_index_upper++;}
	flag_alt = 1-flag_alt;
	/* if (i==i_pivot){ } */}
      /* if (ni!=m){ } */}
    /* for (ni=0;ni<n_i;ni++){ } */}
  for (ni=0;ni<n_i;ni++){ i_[stride*ni] = g_wkspace_[1+g_index_lower+ni]; index_[ni] = index_wkspace_[1+g_index_lower+ni];}
  if (verbose){ array_printf(g_wkspace_,"int",1,2*n_i," %% g_wkspace_: ");}
  if (verbose){ array_printf(index_wkspace_,"int",1,2*n_i," %% index_wkspace_: ");}
  if (verbose){ array_printf(i_,"int",1,stride*n_i," %% pos i_: ");}
  if (verbose){ array_printf(index_,"int",1,n_i," %% pos index_: ");}
  m_upd = n_i - 1 - g_index_lower;
  if (verbose){ printf(" %% local pivot position %d, full pivot position %d\n",m_upd,l+m_upd);}
  free1(&g_wkspace_);
  free1(&index_wkspace_);
  if (verbose){ printf(" %% [finished iquicksort_partition_index]\n");}
  return l+m_upd;
}

unsigned int iquicksort_index(unsigned int recursion_level,int *i_,int stride,int *index_,int l,int r)
{
  /* modified from http://www.comp.dit.ie/rlawlor/Alg_DS/sorting/quicksort.c */
  int j=0; unsigned int n1=recursion_level,n2=recursion_level;
  if( l < r ) { 
    j = iquicksort_partition_index(i_,stride,index_,l,r); 
    if (recursion_level<GLOBAL_quicksort_recursion_limit){ 
      n1 = iquicksort_index(recursion_level+1,i_,stride,index_,l,j-1); 
      n2 = iquicksort_index(recursion_level+1,i_,stride,index_,j+1,r); 
      /* if (recursion_level<GLOBAL_quicksort_recursion_limit){ } */}
    else /* (recursion_level>=GLOBAL_quicksort_recursion_limit) */{ 
      printf(" %% Warning! recursion limit %d reached in iquicksort_index.\n",recursion_level);
      /* recursion_limit breached */}
    /* if( l < r ) { } */}
  return maximum(n1,n2);
}

unsigned int iquicksort_index_driver(int n_i,int *i_,int stride,int *i_workspace_,int *index_)
{
  /* finds index listing index_ so that i_[stride*index_[.]] is sorted. */
  int ni=0;
  unsigned int ndepth=0;
  for (ni=0;ni<n_i;ni++){ i_workspace_[ni] = i_[stride*ni]; index_[ni]=ni;}
  ndepth = iquicksort_index(0,i_workspace_,1,index_,0,n_i-1);
  return ndepth;
}

void iquicksort_index_driver_test()
{
  int n_i = 5;
  int index_[n_i];
  int i_[2*n_i];
  int i_workspace_[1*n_i];
  int index_ans_[5] = { 2 , 0 , 3 , 4 , 1 };
  int n_g = 1024*8;
  int g_stride = 2;
  int flag_g_error=0;
  int *g_index_=NULL;
  int *g_val_=NULL;
  int *g_workspace_=NULL;
  int ng=0,ndepth=0;
  GLOBAL_tic(0);
  printf(" %% testing length 5\n");
  memset(i_,0,2*n_i*sizeof(int));
  i_[0] = 15; i_[2] = 25; i_[4] = 13; i_[6] = 18; i_[8] = 20;
  array_printf(i_,"int",1,2*n_i," %% i_: ");
  iquicksort_index_driver(n_i,i_,2,i_workspace_,index_);
  array_printf(i_workspace_,"int",1,n_i," %% i_workspace_: ");
  array_printf(index_,"int",1,n_i," %% index_: ");
  printf(" %% index_ans_ vs index_: relative error %0.16f\n",ifnormn(5,index_ans_,index_));
  g_index_ = (int *) malloc1(n_g*sizeof(int));
  g_val_ = (int *) malloc1(g_stride*n_g*sizeof(int));
  g_workspace_ = (int *) malloc1(n_g*sizeof(int));
  printf(" %% testing length %d: constant\n",n_g);
  for (ng=0;ng<n_g;ng++){ g_val_[g_stride*ng] = 0; };
  ndepth = iquicksort_index_driver(n_g,g_val_,g_stride,g_workspace_,g_index_);
  flag_g_error=0; for (ng=0;ng<n_g-1;ng++){ flag_g_error += (g_workspace_[ng]> g_workspace_[ng+1]);} for (ng=0;ng<n_g;ng++){ flag_g_error += g_workspace_[ng]!=g_val_[g_stride*g_index_[ng]];}
  printf(" %% ndepth %d, flag_g_error %d\n",ndepth,flag_g_error);
  printf(" %% testing length %d: mod3\n",n_g);
  for (ng=0;ng<n_g;ng++){ g_val_[g_stride*ng] = ng%3; };
  ndepth = iquicksort_index_driver(n_g,g_val_,g_stride,g_workspace_,g_index_);
  flag_g_error=0; for (ng=0;ng<n_g-1;ng++){ flag_g_error += (g_workspace_[ng]> g_workspace_[ng+1]);} for (ng=0;ng<n_g;ng++){ flag_g_error += g_workspace_[ng]!=g_val_[g_stride*g_index_[ng]];}
  printf(" %% ndepth %d, flag_g_error %d\n",ndepth,flag_g_error);
  printf(" %% testing length %d: sorted\n",n_g);
  for (ng=0;ng<n_g;ng++){ g_val_[g_stride*ng] = ng; };
  ndepth = iquicksort_index_driver(n_g,g_val_,g_stride,g_workspace_,g_index_);
  flag_g_error=0; for (ng=0;ng<n_g-1;ng++){ flag_g_error += (g_workspace_[ng]> g_workspace_[ng+1]);} for (ng=0;ng<n_g;ng++){ flag_g_error += g_workspace_[ng]!=g_val_[g_stride*g_index_[ng]];}
  printf(" %% ndepth %d, flag_g_error %d\n",ndepth,flag_g_error);
  printf(" %% testing length %d: reverse sorted\n",n_g);
  for (ng=0;ng<n_g;ng++){ g_val_[g_stride*ng] = n_g-ng; };
  ndepth = iquicksort_index_driver(n_g,g_val_,g_stride,g_workspace_,g_index_);
  flag_g_error=0; for (ng=0;ng<n_g-1;ng++){ flag_g_error += (g_workspace_[ng]> g_workspace_[ng+1]);} for (ng=0;ng<n_g;ng++){ flag_g_error += g_workspace_[ng]!=g_val_[g_stride*g_index_[ng]];}
  printf(" %% ndepth %d, flag_g_error %d\n",ndepth,flag_g_error);
  free1(&g_index_);
  free1(&g_workspace_);
  free1(&g_val_);
  GLOBAL_toc(0,1," iquicksort_index_driver_test: ");
}

unsigned int iquicksort_index_index_driver(int n_i,int *i_,int stride,int *i_workspace_,int *index_orig_from_sort_,int *index_sort_from_orig_,int *index_workspace_)
{
  /* finds index listing index_orig_from_sort_ so that i_[stride*index_orig_from_sort_[.]] is sorted. */
  int ni=0;
  unsigned int ndepth0=0,ndepth1=0;
  for (ni=0;ni<n_i;ni++){ i_workspace_[ni] = i_[stride*ni]; if (index_orig_from_sort_!=NULL){ index_orig_from_sort_[ni]=ni;}}
  ndepth0 = iquicksort_index(0,i_workspace_,1,index_orig_from_sort_,0,n_i-1);
  if (index_sort_from_orig_!=NULL){
    for (ni=0;ni<n_i;ni++){ index_workspace_[ni] = index_orig_from_sort_[ni]; index_sort_from_orig_[ni]=ni;}
    ndepth1 = iquicksort_index(0,index_workspace_,1,index_sort_from_orig_,0,n_i-1);
    /* if (index_sort_from_orig_!=NULL){ } */}
  return maximum(ndepth0,ndepth1);
}

void iquicksort_index_index_driver_test()
{
  int n_i = 5;
  int index_orig_from_sort_[n_i];
  int index_sort_from_orig_[n_i];
  int i_[2*n_i];
  int i_workspace_[1*n_i];
  int index_workspace_[1*n_i];
  int index_orig_from_sort_ans_[5] = { 2 , 0 , 3 , 4 , 1 };
  int index_sort_from_orig_ans_[5] = { 1 , 4 , 0 , 2 , 3 };
  memset(i_,0,2*n_i*sizeof(int));
  i_[0] = 15; i_[2] = 25; i_[4] = 13; i_[6] = 18; i_[8] = 20;
  iquicksort_index_index_driver(n_i,i_,2,i_workspace_,index_orig_from_sort_,index_sort_from_orig_,index_workspace_);
  array_printf(i_,"int",1,2*n_i," %% i_: ");
  array_printf(i_workspace_,"int",1,n_i," %% i_workspace_: ");
  array_printf(index_orig_from_sort_,"int",1,n_i," %% index_orig_from_sort_: ");
  array_printf(index_sort_from_orig_,"int",1,n_i," %% index_sort_from_orig_: ");
  printf(" %% index_orig_from_sort_ans_ vs index_orig_from_sort_: relative error %0.16f\n",ifnormn(5,index_orig_from_sort_ans_,index_orig_from_sort_));
  printf(" %% index_sort_from_orig_ans_ vs index_sort_from_orig_: relative error %0.16f\n",ifnormn(5,index_sort_from_orig_ans_,index_sort_from_orig_));
}

void irandperm(int n_i,int **index_p_,unsigned long long int *rseed_p)
{
  double *d_=NULL;
  int nd=0;
  double *d_workspace_=NULL;
  int *index_;
  d_ = (double *) malloc1(n_i*sizeof(double));
  d_workspace_ = (double *) malloc1(n_i*sizeof(double));
  for (nd=0;nd<n_i;nd++){ d_[nd] = R01GET(rseed_p);}
  if (*index_p_==NULL){ (*index_p_) = (int *) malloc1(n_i*sizeof(int));}
  index_ = *index_p_;
  dquicksort_index_driver(n_i,d_,1,d_workspace_,index_);
  free1(&d_workspace_);
  free1(&d_);
}

void irandperm_test()
{
  int n_i = 10;
  int *index_ = NULL;
  unsigned long long int rseed = 67892;
  GLOBAL_tic(0);
  irandperm(n_i,&index_,&rseed); array_printf(index_,"int",1,n_i," irandperm index_: ");
  irandperm(n_i,&index_,&rseed); array_printf(index_,"int",1,n_i," irandperm index_: ");
  irandperm(n_i,&index_,&rseed); array_printf(index_,"int",1,n_i," irandperm index_: ");
  irandperm(n_i,&index_,&rseed); array_printf(index_,"int",1,n_i," irandperm index_: ");
  irandperm(n_i,&index_,&rseed); array_printf(index_,"int",1,n_i," irandperm index_: ");
  irandperm(n_i,&index_,&rseed); array_printf(index_,"int",1,n_i," irandperm index_: ");
  irandperm(n_i,&index_,&rseed); array_printf(index_,"int",1,n_i," irandperm index_: ");
  free1(&index_);
  GLOBAL_toc(0,1," irandperm_test: ");
}


/* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */

int fquicksort_partition_index_bkp(float *f_,int stride,int *index_,int l,int r) 
{
  float pivot=0,tmpd=0;
  int i=0,j=0,tmpi=0;
  pivot = f_[stride*l];
  i = l; j = r+1;
  do{
    do{ i++;} while( i<=r && f_[stride*i] <= pivot );
    do{ j--;} while( f_[stride*j] > pivot );
    if( i >= j ) break;
    tmpd = f_[stride*i]; f_[stride*i] = f_[stride*j]; f_[stride*j] = tmpd;
    if (index_!=NULL){ tmpi = index_[i]; index_[i] = index_[j]; index_[j] = tmpi;}
  }while(1);
  tmpd = f_[stride*l]; f_[stride*l] = f_[stride*j]; f_[stride*j] = tmpd;
  if (index_!=NULL){ tmpi = index_[l]; index_[l] = index_[j]; index_[j] = tmpi;}
  return j;
}

int fquicksort_partition_index(float *f_full_,int stride,int *index_full_,int l,int r)
{
  int verbose=0;
  float *f_ = f_full_ + stride*l;
  int *index_ = index_full_ + l;
  float *g_wkspace_=NULL;
  int *index_wkspace_ = NULL;
  int n_f = 1+r-l;
  int m = n_f/2;
  float f_pivot=0;
  int index_pivot=0;
  int g_index_lower=0,g_index_upper=0;
  int nf=0;
  int flag_alt=0;
  float f=0;
  int index=0;
  int m_upd=0;
  if (verbose){ printf(" %% [entering fquicksort_partition_index]\n");}
  if (verbose){ array_printf(f_,"float",1,stride*n_f," %% pre f_: ");}
  if (verbose){ array_printf(index_,"int",1,n_f," %% pre index_: ");}
  g_wkspace_ = (float *) malloc1((unsigned long long int)(2*n_f)*sizeof(float));
  index_wkspace_ = (int *) malloc1((unsigned long long int)(2*n_f)*sizeof(int));
  f_pivot = f_[stride*m];
  index_pivot = index_[m];
  if (verbose){ printf(" %% m %d, f_pivot %f, index_pivot %d\n",m,f_pivot,index_pivot);}
  g_wkspace_[n_f] = f_pivot;
  index_wkspace_[n_f] = index_pivot;
  g_index_lower = n_f-1;
  g_index_upper = n_f+1;
  flag_alt = 0;
  for (nf=0;nf<n_f;nf++){
    if (nf!=m){
      f = f_[stride*nf];
      index = index_[nf];
      if (f< f_pivot){ g_wkspace_[g_index_lower] = f; index_wkspace_[g_index_lower] = index; g_index_lower--;}
      if (f> f_pivot){ g_wkspace_[g_index_upper] = f; index_wkspace_[g_index_upper] = index; g_index_upper++;}
      if (f==f_pivot){
	if (flag_alt==0){ g_wkspace_[g_index_lower] = f; index_wkspace_[g_index_lower] = index; g_index_lower--;}
	if (flag_alt==1){ g_wkspace_[g_index_upper] = f; index_wkspace_[g_index_upper] = index; g_index_upper++;}
	flag_alt = 1-flag_alt;
	/* if (f==f_pivot){ } */}
      /* if (nf!=m){ } */}
    /* for (nf=0;nf<n_f;nf++){ } */}
  for (nf=0;nf<n_f;nf++){ f_[stride*nf] = g_wkspace_[1+g_index_lower+nf]; index_[nf] = index_wkspace_[1+g_index_lower+nf];}
  if (verbose){ array_printf(g_wkspace_,"float",1,2*n_f," %% g_wkspace_: ");}
  if (verbose){ array_printf(index_wkspace_,"int",1,2*n_f," %% index_wkspace_: ");}
  if (verbose){ array_printf(f_,"float",1,stride*n_f," %% pos f_: ");}
  if (verbose){ array_printf(index_,"int",1,n_f," %% pos index_: ");}
  m_upd = n_f - 1 - g_index_lower;
  if (verbose){ printf(" %% local pivot position %d, full pivot position %d\n",m_upd,l+m_upd);}
  free1(&g_wkspace_);
  free1(&index_wkspace_);
  if (verbose){ printf(" %% [finished fquicksort_partition_index]\n");}
  return l+m_upd;
}

unsigned int fquicksort_index(unsigned int recursion_level,float *f_,int stride,int *index_,int l,int r)
{
  /* modified from http://www.comp.dit.ie/rlawlor/Alg_DS/sorting/quicksort.c */
  int j=0; unsigned int n1=recursion_level,n2=recursion_level;
  if( l < r ) { 
    j = fquicksort_partition_index(f_,stride,index_,l,r); 
    if (recursion_level<GLOBAL_quicksort_recursion_limit){ 
      n1 = fquicksort_index(recursion_level+1,f_,stride,index_,l,j-1); 
      n2 = fquicksort_index(recursion_level+1,f_,stride,index_,j+1,r); 
      /* if (recursion_level<GLOBAL_quicksort_recursion_limit){ } */}
    else /* (recursion_level>=GLOBAL_quicksort_recursion_limit) */{ 
      printf(" %% Warning! recursion limit %d reached in fquicksort_index.\n %% This means that many row- and col-scores are identical.\n %% It is likely that your data is mostly zero or constant.\n",recursion_level);
      /* recursion_limit breached */}
    /* if( l < r ) { } */}
  return maximum(n1,n2);
}

unsigned int fquicksort_index_driver(int n_f,float *f_,int stride,float *f_workspace_,int *index_)
{
  /* finds index listing index_ so that f_[stride*index_[.]] is sorted. */
  int nf=0;
  unsigned int ndepth=0;
  for (nf=0;nf<n_f;nf++){ f_workspace_[nf] = f_[stride*nf]; index_[nf]=nf;}
  ndepth = fquicksort_index(0,f_workspace_,1,index_,0,n_f-1);
  return ndepth;
}

void fquicksort_index_driver_test()
{
  int n_f = 5;
  int index_[n_f];
  float f_[2*n_f];
  float f_workspace_[1*n_f];
  int index_ans_[5] = { 2 , 0 , 3 , 4 , 1 };
  int n_g = 1024*8;
  int g_stride = 2;
  int flag_g_error=0;
  int *g_index_=NULL;
  float *g_val_=NULL;
  float *g_workspace_=NULL;
  int ng=0,ndepth=0;
  GLOBAL_tic(0);
  printf(" %% testing length 5\n");
  memset(f_,0,2*n_f*sizeof(float));
  f_[0] = 15; f_[2] = 25; f_[4] = 13; f_[6] = 18; f_[8] = 20;
  array_printf(f_,"float",1,2*n_f," %% f_: ");
  fquicksort_index_driver(n_f,f_,2,f_workspace_,index_);
  array_printf(f_workspace_,"float",1,n_f," %% f_workspace_: ");
  array_printf(index_,"int",1,n_f," %% index_: ");
  printf(" %% index_ans_ vs index_: relative error %0.16f\n",ifnormn(5,index_ans_,index_));
  g_index_ = (int *) malloc1(n_g*sizeof(int));
  g_val_ = (float *) malloc1(g_stride*n_g*sizeof(float));
  g_workspace_ = (float *) malloc1(n_g*sizeof(float));
  printf(" %% testing length %d: constant\n",n_g);
  for (ng=0;ng<n_g;ng++){ g_val_[g_stride*ng] = 0; };
  ndepth = fquicksort_index_driver(n_g,g_val_,g_stride,g_workspace_,g_index_);
  flag_g_error=0; for (ng=0;ng<n_g-1;ng++){ flag_g_error += (g_workspace_[ng]> g_workspace_[ng+1]);} for (ng=0;ng<n_g;ng++){ flag_g_error += g_workspace_[ng]!=g_val_[g_stride*g_index_[ng]];}
  printf(" %% ndepth %d, flag_g_error %d\n",ndepth,flag_g_error);
  printf(" %% testing length %d: mod3\n",n_g);
  for (ng=0;ng<n_g;ng++){ g_val_[g_stride*ng] = ng%3; };
  ndepth = fquicksort_index_driver(n_g,g_val_,g_stride,g_workspace_,g_index_);
  flag_g_error=0; for (ng=0;ng<n_g-1;ng++){ flag_g_error += (g_workspace_[ng]> g_workspace_[ng+1]);} for (ng=0;ng<n_g;ng++){ flag_g_error += g_workspace_[ng]!=g_val_[g_stride*g_index_[ng]];}
  printf(" %% ndepth %d, flag_g_error %d\n",ndepth,flag_g_error);
  printf(" %% testing length %d: sorted\n",n_g);
  for (ng=0;ng<n_g;ng++){ g_val_[g_stride*ng] = ng; };
  ndepth = fquicksort_index_driver(n_g,g_val_,g_stride,g_workspace_,g_index_);
  flag_g_error=0; for (ng=0;ng<n_g-1;ng++){ flag_g_error += (g_workspace_[ng]> g_workspace_[ng+1]);} for (ng=0;ng<n_g;ng++){ flag_g_error += g_workspace_[ng]!=g_val_[g_stride*g_index_[ng]];}
  printf(" %% ndepth %d, flag_g_error %d\n",ndepth,flag_g_error);
  printf(" %% testing length %d: reverse sorted\n",n_g);
  for (ng=0;ng<n_g;ng++){ g_val_[g_stride*ng] = n_g-ng; };
  ndepth = fquicksort_index_driver(n_g,g_val_,g_stride,g_workspace_,g_index_);
  flag_g_error=0; for (ng=0;ng<n_g-1;ng++){ flag_g_error += (g_workspace_[ng]> g_workspace_[ng+1]);} for (ng=0;ng<n_g;ng++){ flag_g_error += g_workspace_[ng]!=g_val_[g_stride*g_index_[ng]];}
  printf(" %% ndepth %d, flag_g_error %d\n",ndepth,flag_g_error);
  free1(&g_index_);
  free1(&g_workspace_);
  free1(&g_val_);
  GLOBAL_toc(0,1," fquicksort_index_driver_test: ");
}

unsigned int fquicksort_index_index_driver(int n_f,float *f_,int stride,float *f_workspace_,int *index_orig_from_sort_,int *index_sort_from_orig_,int *index_workspace_)
{
  /* finds index listing index_orig_from_sort_ so that f_[stride*index_orig_from_sort_[.]] is sorted. */
  int nf=0;
  unsigned int ndepth0=0,ndepth1=0;
  for (nf=0;nf<n_f;nf++){ f_workspace_[nf] = f_[stride*nf]; if (index_orig_from_sort_!=NULL){ index_orig_from_sort_[nf]=nf;}}
  ndepth0 = fquicksort_index(0,f_workspace_,1,index_orig_from_sort_,0,n_f-1);
  if (index_sort_from_orig_!=NULL){
    for (nf=0;nf<n_f;nf++){ index_workspace_[nf] = index_orig_from_sort_[nf]; index_sort_from_orig_[nf]=nf;}
    ndepth1 = iquicksort_index(0,index_workspace_,1,index_sort_from_orig_,0,n_f-1);
    /* if (index_sort_from_orig_!=NULL){ } */}
  return maximum(ndepth0,ndepth1);
}

void fquicksort_index_index_driver_test()
{
  int n_f = 5;
  int index_orig_from_sort_[n_f];
  int index_sort_from_orig_[n_f];
  float f_[2*n_f];
  float f_workspace_[1*n_f];
  int index_workspace_[1*n_f];
  int index_orig_from_sort_ans_[5] = { 2 , 0 , 3 , 4 , 1 };
  int index_sort_from_orig_ans_[5] = { 1 , 4 , 0 , 2 , 3 };
  memset(f_,0,2*n_f*sizeof(float));
  f_[0] = 15; f_[2] = 25; f_[4] = 13; f_[6] = 18; f_[8] = 20;
  fquicksort_index_index_driver(n_f,f_,2,f_workspace_,index_orig_from_sort_,index_sort_from_orig_,index_workspace_);
  array_printf(f_,"float",1,2*n_f," %% f_: ");
  array_printf(f_workspace_,"float",1,n_f," %% f_workspace_: ");
  array_printf(index_orig_from_sort_,"int",1,n_f," %% index_orig_from_sort_: ");
  array_printf(index_sort_from_orig_,"int",1,n_f," %% index_sort_from_orig_: ");
  printf(" %% index_orig_from_sort_ans_ vs index_orig_from_sort_: relative error %0.16f\n",ifnormn(5,index_orig_from_sort_ans_,index_orig_from_sort_));
  printf(" %% index_sort_from_orig_ans_ vs index_sort_from_orig_: relative error %0.16f\n",ifnormn(5,index_sort_from_orig_ans_,index_sort_from_orig_));
}

/* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */

int dquicksort_partition_index_bkp(double *d_,int stride,int *index_,int l,int r) 
{
  double pivot=0,tmpd=0;
  int i=0,j=0,tmpi=0;
  pivot = d_[stride*l];
  i = l; j = r+1;
  do{
    do{ i++;} while( i<=r && d_[stride*i] <= pivot );
    do{ j--;} while( d_[stride*j] > pivot );
    if( i >= j ) break;
    tmpd = d_[stride*i]; d_[stride*i] = d_[stride*j]; d_[stride*j] = tmpd;
    if (index_!=NULL){ tmpi = index_[i]; index_[i] = index_[j]; index_[j] = tmpi;}
  }while(1);
  tmpd = d_[stride*l]; d_[stride*l] = d_[stride*j]; d_[stride*j] = tmpd;
  if (index_!=NULL){ tmpi = index_[l]; index_[l] = index_[j]; index_[j] = tmpi;}
  return j;
}

int dquicksort_partition_index(double *d_full_,int stride,int *index_full_,int l,int r)
{
  int verbose=0;
  double *d_ = d_full_ + stride*l;
  int *index_ = index_full_ + l;
  double *g_wkspace_=NULL;
  int *index_wkspace_ = NULL;
  int n_d = 1+r-l;
  int m = n_d/2;
  double d_pivot=0;
  int index_pivot=0;
  int g_index_lower=0,g_index_upper=0;
  int nd=0;
  int flag_alt=0;
  double d=0;
  int index=0;
  int m_upd=0;
  if (verbose){ printf(" %% [entering dquicksort_partition_index]\n");}
  if (verbose){ array_printf(d_,"double",1,stride*n_d," %% pre d_: ");}
  if (verbose){ array_printf(index_,"int",1,n_d," %% pre index_: ");}
  g_wkspace_ = (double *) malloc1((unsigned long long int)(2*n_d)*sizeof(double));
  index_wkspace_ = (int *) malloc1((unsigned long long int)(2*n_d)*sizeof(int));
  d_pivot = d_[stride*m];
  index_pivot = index_[m];
  if (verbose){ printf(" %% m %d, d_pivot %f, index_pivot %d\n",m,d_pivot,index_pivot);}
  g_wkspace_[n_d] = d_pivot;
  index_wkspace_[n_d] = index_pivot;
  g_index_lower = n_d-1;
  g_index_upper = n_d+1;
  flag_alt = 0;
  for (nd=0;nd<n_d;nd++){
    if (nd!=m){
      d = d_[stride*nd];
      index = index_[nd];
      if (d< d_pivot){ g_wkspace_[g_index_lower] = d; index_wkspace_[g_index_lower] = index; g_index_lower--;}
      if (d> d_pivot){ g_wkspace_[g_index_upper] = d; index_wkspace_[g_index_upper] = index; g_index_upper++;}
      if (d==d_pivot){
	if (flag_alt==0){ g_wkspace_[g_index_lower] = d; index_wkspace_[g_index_lower] = index; g_index_lower--;}
	if (flag_alt==1){ g_wkspace_[g_index_upper] = d; index_wkspace_[g_index_upper] = index; g_index_upper++;}
	flag_alt = 1-flag_alt;
	/* if (d==d_pivot){ } */}
      /* if (nd!=m){ } */}
    /* for (nd=0;nd<n_d;nd++){ } */}
  for (nd=0;nd<n_d;nd++){ d_[stride*nd] = g_wkspace_[1+g_index_lower+nd]; index_[nd] = index_wkspace_[1+g_index_lower+nd];}
  if (verbose){ array_printf(g_wkspace_,"double",1,2*n_d," %% g_wkspace_: ");}
  if (verbose){ array_printf(index_wkspace_,"int",1,2*n_d," %% index_wkspace_: ");}
  if (verbose){ array_printf(d_,"double",1,stride*n_d," %% pos d_: ");}
  if (verbose){ array_printf(index_,"int",1,n_d," %% pos index_: ");}
  m_upd = n_d - 1 - g_index_lower;
  if (verbose){ printf(" %% local pivot position %d, full pivot position %d\n",m_upd,l+m_upd);}
  free1(&g_wkspace_);
  free1(&index_wkspace_);
  if (verbose){ printf(" %% [finished dquicksort_partition_index]\n");}
  return l+m_upd;
}

unsigned int dquicksort_index(unsigned int recursion_level,double *d_,int stride,int *index_,int l,int r)
{
  /* modified from http://www.comp.dit.ie/rlawlor/Alg_DS/sorting/quicksort.c */
  int j=0; unsigned int n1=recursion_level,n2=recursion_level;
  if( l < r ) { 
    j = dquicksort_partition_index(d_,stride,index_,l,r); 
    if (recursion_level<GLOBAL_quicksort_recursion_limit){ 
      n1 = dquicksort_index(recursion_level+1,d_,stride,index_,l,j-1); 
      n2 = dquicksort_index(recursion_level+1,d_,stride,index_,j+1,r); 
      /* if (recursion_level<GLOBAL_quicksort_recursion_limit){ } */}
    else /* (recursion_level>=GLOBAL_quicksort_recursion_limit) */{ 
      printf(" %% Warning! recursion limit %d reached in dquicksort_index.\n %% This means that many row- and col-scores are identical.\n %% It is likely that your data is mostly zero or constant.\n",recursion_level);
      /* recursion_limit breached */}
    /* if( l < r ) { } */}
  return maximum(n1,n2);
}

unsigned int dquicksort_index_driver(int n_d,double *d_,int stride,double *d_workspace_,int *index_)
{
  /* finds index listing index_ so that d_[stride*index_[.]] is sorted. */
  int nd=0;
  unsigned int ndepth=0;
  for (nd=0;nd<n_d;nd++){ d_workspace_[nd] = d_[stride*nd]; index_[nd]=nd;}
  ndepth = dquicksort_index(0,d_workspace_,1,index_,0,n_d-1);
  return ndepth;
}

void dquicksort_index_driver_test()
{
  int n_f = 5;
  int index_[n_f];
  double f_[2*n_f];
  double f_workspace_[1*n_f];
  int index_ans_[5] = { 2 , 0 , 3 , 4 , 1 };
  int n_g = 1024*8;
  int g_stride = 2;
  int flag_g_error=0;
  int *g_index_=NULL;
  double *g_val_=NULL;
  double *g_workspace_=NULL;
  int ng=0,ndepth=0;
  GLOBAL_tic(0);
  printf(" %% testing length 5\n");
  memset(f_,0,2*n_f*sizeof(double));
  f_[0] = 15; f_[2] = 25; f_[4] = 13; f_[6] = 18; f_[8] = 20;
  array_printf(f_,"double",1,2*n_f," %% f_: ");
  dquicksort_index_driver(n_f,f_,2,f_workspace_,index_);
  array_printf(f_workspace_,"double",1,n_f," %% f_workspace_: ");
  array_printf(index_,"int",1,n_f," %% index_: ");
  printf(" %% index_ans_ vs index_: relative error %0.16f\n",ifnormn(5,index_ans_,index_));
  g_index_ = (int *) malloc1(n_g*sizeof(int));
  g_val_ = (double *) malloc1(g_stride*n_g*sizeof(double));
  g_workspace_ = (double *) malloc1(n_g*sizeof(double));
  printf(" %% testing length %d: constant\n",n_g);
  for (ng=0;ng<n_g;ng++){ g_val_[g_stride*ng] = 0; };
  ndepth = dquicksort_index_driver(n_g,g_val_,g_stride,g_workspace_,g_index_);
  flag_g_error=0; for (ng=0;ng<n_g-1;ng++){ flag_g_error += (g_workspace_[ng]> g_workspace_[ng+1]);} for (ng=0;ng<n_g;ng++){ flag_g_error += g_workspace_[ng]!=g_val_[g_stride*g_index_[ng]];}
  printf(" %% ndepth %d, flag_g_error %d\n",ndepth,flag_g_error);
  printf(" %% testing length %d: mod3\n",n_g);
  for (ng=0;ng<n_g;ng++){ g_val_[g_stride*ng] = ng%3; };
  ndepth = dquicksort_index_driver(n_g,g_val_,g_stride,g_workspace_,g_index_);
  flag_g_error=0; for (ng=0;ng<n_g-1;ng++){ flag_g_error += (g_workspace_[ng]> g_workspace_[ng+1]);} for (ng=0;ng<n_g;ng++){ flag_g_error += g_workspace_[ng]!=g_val_[g_stride*g_index_[ng]];}
  printf(" %% ndepth %d, flag_g_error %d\n",ndepth,flag_g_error);
  printf(" %% testing length %d: sorted\n",n_g);
  for (ng=0;ng<n_g;ng++){ g_val_[g_stride*ng] = ng; };
  ndepth = dquicksort_index_driver(n_g,g_val_,g_stride,g_workspace_,g_index_);
  flag_g_error=0; for (ng=0;ng<n_g-1;ng++){ flag_g_error += (g_workspace_[ng]> g_workspace_[ng+1]);} for (ng=0;ng<n_g;ng++){ flag_g_error += g_workspace_[ng]!=g_val_[g_stride*g_index_[ng]];}
  printf(" %% ndepth %d, flag_g_error %d\n",ndepth,flag_g_error);
  printf(" %% testing length %d: reverse sorted\n",n_g);
  for (ng=0;ng<n_g;ng++){ g_val_[g_stride*ng] = n_g-ng; };
  ndepth = dquicksort_index_driver(n_g,g_val_,g_stride,g_workspace_,g_index_);
  flag_g_error=0; for (ng=0;ng<n_g-1;ng++){ flag_g_error += (g_workspace_[ng]> g_workspace_[ng+1]);} for (ng=0;ng<n_g;ng++){ flag_g_error += g_workspace_[ng]!=g_val_[g_stride*g_index_[ng]];}
  printf(" %% ndepth %d, flag_g_error %d\n",ndepth,flag_g_error);
  free1(&g_index_);
  free1(&g_workspace_);
  free1(&g_val_);
  GLOBAL_toc(0,1," dquicksort_index_driver_test: ");
}

unsigned int dquicksort_index_index_driver(int n_d,double *d_,int stride,double *d_workspace_,int *index_orig_from_sort_,int *index_sort_from_orig_,int *index_workspace_)
{
  /* finds index listing index_orig_from_sort_ so that d_[stride*index_orig_from_sort_[.]] is sorted. */
  int nd=0;
  unsigned int ndepth0=0,ndepth1=0;
  for (nd=0;nd<n_d;nd++){ d_workspace_[nd] = d_[stride*nd]; if (index_orig_from_sort_!=NULL){ index_orig_from_sort_[nd]=nd;}}
  ndepth0 = dquicksort_index(0,d_workspace_,1,index_orig_from_sort_,0,n_d-1);
  if (index_sort_from_orig_!=NULL){
    for (nd=0;nd<n_d;nd++){ index_workspace_[nd] = index_orig_from_sort_[nd]; index_sort_from_orig_[nd]=nd;}
    ndepth1 = iquicksort_index(0,index_workspace_,1,index_sort_from_orig_,0,n_d-1);
    /* if (index_sort_from_orig_!=NULL){ } */}
  return maximum(ndepth0,ndepth1);
}

void dquicksort_index_index_driver_test()
{
  int n_d = 5;
  int index_orig_from_sort_[n_d];
  int index_sort_from_orig_[n_d];
  double d_[2*n_d];
  double d_workspace_[1*n_d];
  int index_workspace_[1*n_d];
  int index_orig_from_sort_ans_[5] = { 2 , 0 , 3 , 4 , 1 };
  int index_sort_from_orig_ans_[5] = { 1 , 4 , 0 , 2 , 3 };
  memset(d_,0,2*n_d*sizeof(double));
  d_[0] = 15; d_[2] = 25; d_[4] = 13; d_[6] = 18; d_[8] = 20;
  dquicksort_index_index_driver(n_d,d_,2,d_workspace_,index_orig_from_sort_,index_sort_from_orig_,index_workspace_);
  array_printf(d_,"double",1,2*n_d," %% d_: ");
  array_printf(d_workspace_,"double",1,n_d," %% d_workspace_: ");
  array_printf(index_orig_from_sort_,"int",1,n_d," %% index_orig_from_sort_: ");
  array_printf(index_sort_from_orig_,"int",1,n_d," %% index_sort_from_orig_: ");
  printf(" %% index_orig_from_sort_ans_ vs index_orig_from_sort_: relative error %0.16f\n",ifnormn(5,index_orig_from_sort_ans_,index_orig_from_sort_));
  printf(" %% index_sort_from_orig_ans_ vs index_sort_from_orig_: relative error %0.16f\n",ifnormn(5,index_sort_from_orig_ans_,index_sort_from_orig_));
}

/* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */

