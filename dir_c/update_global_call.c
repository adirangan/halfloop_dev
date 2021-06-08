#ifndef _MONOLITH
#include "halfloop_header.h"
#endif /* _MONOLITH */

void update_global(char *vname)
{
  /* given a variable name, scan the appropriate format into the appropriate global variable */
  int verbose=GLOBAL_verbose;
  char comma_vs_semicolon[1],tmpchar[FNAMESIZE]; int length=0,nv=0;
  if (0){ /* do nothing */ }
  else if (strcmp(vname,"GLOBAL_verbose")==0){ scanf("%d",&GLOBAL_verbose); if (verbose>0){ printf("%s read to be %d\n",vname,GLOBAL_verbose);}}
  else if (strcmp(vname,"GLOBAL_mode")==0){ scanf("%[^,;]",GLOBAL_mode); if (verbose>0){ printf("%s read to be %s\n",vname,GLOBAL_mode);} /* else if (strcmp(vname,"GLOBAL_mode")==0){ } */}
  else if (strcmp(vname,"GLOBAL_E_base_mda_r4")==0){ scanf("%[^,;]",GLOBAL_E_base_mda_r4); if (verbose>0){ printf("%s read to be %s\n",vname,GLOBAL_E_base_mda_r4);} /* else if (strcmp(vname,"GLOBAL_E_base_mda_r4")==0){ } */}
  else if (strcmp(vname,"GLOBAL_flag_r0drop_vs_rcdrop")==0){ scanf("%d",&GLOBAL_flag_r0drop_vs_rcdrop); if (verbose>0){ printf("%s read to be %d\n",vname,GLOBAL_flag_r0drop_vs_rcdrop);}}
  else if (strcmp(vname,"GLOBAL_gamma")==0){ scanf("%lf",&GLOBAL_gamma); if (verbose>0){ printf("%s read to be %f\n",vname,GLOBAL_gamma);}}
  else if (strcmp(vname,"GLOBAL_n_shuffle")==0){ scanf("%d",&GLOBAL_n_shuffle); if (verbose>0){ printf("%s read to be %d\n",vname,GLOBAL_n_shuffle);}}
  else if (strcmp(vname,"GLOBAL_p_set")==0){ scanf("%lf",&GLOBAL_p_set); if (verbose>0){ printf("%s read to be %f\n",vname,GLOBAL_p_set);}}
  else if (strcmp(vname,"GLOBAL_n_member_lob")==0){ scanf("%d",&GLOBAL_n_member_lob); if (verbose>0){ printf("%s read to be %d\n",vname,GLOBAL_n_member_lob);}}
  else if (strcmp(vname,"GLOBAL_dir_trunk")==0){ scanf("%[^,;]",GLOBAL_dir_trunk); if (verbose>0){ printf("%s read to be %s\n",vname,GLOBAL_dir_trunk);} /* else if (strcmp(vname,"GLOBAL_dir_trunk")==0){ } */}
  else if (strcmp(vname,"GLOBAL_prefix_base")==0){ scanf("%[^,;]",GLOBAL_prefix_base); if (verbose>0){ printf("%s read to be %s\n",vname,GLOBAL_prefix_base);} /* else if (strcmp(vname,"GLOBAL_prefix_base")==0){ } */}
  else if (strcmp(vname,"GLOBAL_flag_force_create")==0){ scanf("%d",&GLOBAL_flag_force_create); if (verbose>0){ printf("%s read to be %d\n",vname,GLOBAL_flag_force_create);}}
  else if (strcmp(vname,"GLOBAL_flag_omp_use")==0){ scanf("%d",&GLOBAL_flag_omp_use); if (verbose>0){ printf("%s read to be %d\n",vname,GLOBAL_flag_omp_use);}}
  else if (strcmp(vname,"END")==0){ /* do nothing */ if (verbose>0){ printf("end of input reached\n");}}
/*   else if (strcmp(vname,"yy")==0){ scanf("%zz",&yy); if (verbose>0){ printf("%s read to be %zz\n",vname,yy);}} */
  else /* if anything else */{ printf(" %% Error! vname %s in update_global\n",vname); exit(EXIT_FAILURE); }
}

void read_input()
{
  /* This reads the piped input file, or standard input.
     the variable names should not be longer than 128 characters */
  int verbose=GLOBAL_verbose;
  char vname[128],equals[128],space[128],semicolon[128];
  do{
    scanf("%[^=]",vname);scanf("%s",equals);scanf("%c",space);update_global(vname);scanf("%c",semicolon);
    if (verbose>1){ printf("At this point variable name is (%s), equals is (%s), semicolon is (%s)\n",vname,equals,semicolon);} 
    scanf("%c",semicolon);}
  while (strcmp(vname,"END")!=0);
}
