#include "client_basic.h"
#include "client_seq2co.h"
#include <stdio.h>
#include <err.h>
#include <errno.h>
#ifdef _OPENMP
  #include <omp.h>
#endif
extern dim_shuffle_t* dim_shuffle;
extern unsigned int hashsize;
extern int component_num;
const char co_dstat[] = "cofiles.stat";
#define LINE_LEN 1024
const char skch_prefix[]="combco";
const char idx_prefix[]="combco.index";

int main (int argc, char**argv){
	printf ("kssd_client version 1.0\nUnit_space_size = %d\n\n",COMPONENT_SZ);
	if (argc < 3){
		printf("USAGE: kssd_client <*.shuf> <outdir> <fastq/fasta file>\n");
		exit(0);
	}
	int has_abund  = 1 ;
	dim_shuffle = read_dim_shuffle_file(argv[1]);
	hashsize = get_hashsz(dim_shuffle);
	seq2co_global_var_initial();
	infile_tab_t *infile_stat = organize_infile_frm_arg( argc-3, argv+3, 1);
	infile_fmt_count_t* qry_fmt_count = infile_fmt_count(infile_stat);
	bool is_valid_fas_fq_in = (infile_stat->infile_num != 0) &&
          (qry_fmt_count->fasta + qry_fmt_count->fastq == infile_stat->infile_num );
	if(!is_valid_fas_fq_in) err(errno,"not valid raw seq format"); 

	llong *CO = (llong *)malloc(hashsize * sizeof(llong) );
	llong all_ctx_ct = 0;
	ctx_obj_ct_t *ctx_ct_list = malloc(sizeof(ctx_obj_ct_t) * infile_stat->infile_num);

	mkdir(argv[2],0700);
	int num_threads = omp_get_num_procs();
	for(int i = 0; i< infile_stat->infile_num; i++){
		char* seqfname = infile_stat->organized_infile_tab[i].fpath;
		char cofname[PATHLEN];
		sprintf(cofname,"%s/%d.co",argv[2],i);
		llong *co;
		
		if(isOK_fmt_infile(seqfname,fastq_fmt,FQ_FMT_SZ)) {
			co = mt_shortreads2koc(seqfname,CO,"", num_threads);
			ctx_ct_list[i] = write_fqkoc2files(cofname,co);
		}
		else{
			  co = fasta2co(seqfname,CO,"");
        ctx_ct_list[i] = wrt_co2cmpn_use_inn_subctx(cofname,co);
				has_abund  = 0 ;
		}

		all_ctx_ct += ctx_ct_list[i] ;
		printf("%d/%d decomposing %s\r",i,infile_stat->infile_num,seqfname) ;
	}
	printf("\n");
	free(CO);
/*****combining all *.co.i into comb.co.i *****/
#pragma omp parallel for schedule(guided)
	for(int c = 0; c < component_num; c++){
 //get cth fsize by cof_index_in_cbdco[j+1] - cof_index_in_cbdco[j]
  	size_t *cof_index_in_cbdco = malloc( (infile_stat->infile_num + 1) * sizeof(size_t) );
    char tmpfname[PATHLEN]; 

    FILE *com_cofp,*indexfp, *com_abund_fp;
    sprintf(tmpfname,"%s/%s.%d",argv[2],skch_prefix,c);
    if( (com_cofp = fopen(tmpfname,"wb")) == NULL) err(errno,"%s",tmpfname);
    sprintf(tmpfname,"%s/%s.%d",argv[2],idx_prefix,c);
    if( (indexfp = fopen(tmpfname,"wb")) == NULL) err(errno,"%s",tmpfname);

    void *tmp_mem = malloc( LD_FCTR * 2 * (1LLU << (4*COMPONENT_SZ - CTX_SPC_USE_L)) *sizeof(unsigned int) );
    struct stat tmpstat;

    FILE *tmpfp;
    cof_index_in_cbdco[0] = 0; //first offset initial to 0

    for(int i = 0; i< infile_stat->infile_num; i++){

      sprintf(tmpfname,"%s/%d.co.%d", argv[2], i, c);
      if( (tmpfp = fopen(tmpfname,"rb")) == NULL) err(errno,"%s",tmpfname);
      stat(tmpfname,&tmpstat);

      int tmpkmerct = tmpstat.st_size/sizeof(unsigned int);
      cof_index_in_cbdco[i+1] = (size_t)cof_index_in_cbdco[i] + tmpkmerct;
      fread(tmp_mem,tmpstat.st_size,1,tmpfp);
      fwrite(tmp_mem,tmpstat.st_size,1,com_cofp);
      fclose(tmpfp);
      remove(tmpfname);
		}
		fclose(com_cofp);
    fwrite(cof_index_in_cbdco,sizeof(size_t), infile_stat->infile_num + 1, indexfp);
		fclose(indexfp);		
		free(cof_index_in_cbdco);

		if(has_abund){
			sprintf(tmpfname,"%s/%s.%d.a",argv[2],skch_prefix,c);
      com_abund_fp = fopen(tmpfname,"wb");
      if( com_abund_fp == NULL) err(errno,"%s",tmpfname);

			for(int i = 0; i< infile_stat->infile_num; i++){
        sprintf(tmpfname,"%s/%d.co.%d.a", argv[2], i, c);
        if( (tmpfp = fopen(tmpfname,"rb")) == NULL) {
					has_abund = 0;
					remove(tmpfname);
					break;
				};
				stat(tmpfname,&tmpstat);

      	fread(tmp_mem, tmpstat.st_size, 1 , tmpfp);
      	fwrite(tmp_mem, tmpstat.st_size, 1 , com_abund_fp );
      	fclose(tmpfp);
      	remove(tmpfname);
      }

      fclose(com_abund_fp) ;

		}
	      free(tmp_mem);	
 	}// component loop

	co_dstat_t co_dstat_wrout;
  co_dstat_wrout.shuf_id = dim_shuffle->dim_shuffle_stat.id ;
  co_dstat_wrout.koc = has_abund? 1: 0 ;
  co_dstat_wrout.kmerlen = dim_shuffle->dim_shuffle_stat.k * 2;
  co_dstat_wrout.dim_rd_len = dim_shuffle->dim_shuffle_stat.drlevel * 2 ;
  co_dstat_wrout.comp_num = component_num ;
  co_dstat_wrout.infile_num = infile_stat->infile_num;
  co_dstat_wrout.all_ctx_ct = all_ctx_ct;

  char *co_dstat_fullname = malloc(PATHLEN*sizeof(char) );
  sprintf(co_dstat_fullname, "%s/%s",argv[2],co_dstat);
  FILE *coutfp;
  if ( ( coutfp = fopen(co_dstat_fullname,"wb")) == NULL ) err(errno,"%s",co_dstat_fullname);
  fwrite(&co_dstat_wrout,sizeof(co_dstat_wrout),1,coutfp);

  //write file ctx_ct and names
  fwrite(ctx_ct_list,sizeof(ctx_obj_ct_t),co_dstat_wrout.infile_num,coutfp);
  free(ctx_ct_list);

  for(int i = 0; i< co_dstat_wrout.infile_num; i++)
    fwrite(infile_stat->organized_infile_tab[i].fpath,PATHLEN,1,coutfp);

  fclose(coutfp);
	
	return 1;
}
