/*used by both global and sumcommand*/ 
#include "client_basic.h"
#include <assert.h>
#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <err.h>
#include <errno.h>
#include <sys/stat.h>
#include <dirent.h>
#include <sys/sysinfo.h>
#include <ctype.h>
#include <math.h>
/*global basemap and mapbase*/

const int Basemap[128] =
{
  [0 ... 127] = DEFAULT,
  ['a'] = 0, ['A'] = 0,
  ['c'] = 1, ['C'] = 1,
  ['g'] = 2, ['G'] = 2,
  ['t'] = 3, ['T'] = 3,
};
const char Mapbase[]={'A','C','G','T'};


//primers for hashtable size
const unsigned int primer[25] =
{ //hashsize nearest larger power of 2:  2^8 ~ 2^32, index 0 ~ 24
  251, 509, 1021, 2039, 4093, 8191, 16381,
  32749, 65521, 131071, 262139, 524287,
  1048573, 2097143, 4194301, 8388593, 16777213,
  33554393, 67108859, 134217689, 268435399,
  536870909, 1073741789, 2147483647, 4294967291
};



/*allowed input file format*/
const char *acpt_infile_fmt[ACPT_FMT_SZ] = {
	"fna",
	"fas",
	"fasta",
	"fq",
	"fastq",
	"fa",
	"co"
};
//stage I fmt
const char *fasta_fmt[FAS_FMT_SZ] = {
	"fasta",
	"fna",
	"fas",
	"fa"	
};
const char *fastq_fmt[FQ_FMT_SZ] = {
	"fq",
	"fastq"
};
// compresssion format
const char *compress_fmt[CMPRESS_FMT_SZ] = {
  ".gz",
  ".bz2"
};

const char *co_fmt[CO_FMT_SZ] = {
  "co"
};



FILE * fpathopen (const char *dpath, const char *fname, const char *mode )
{
  char *fullname = malloc(PATHLEN*sizeof(char));
  sprintf(fullname,"%s/%s",dpath,fname);
	
	struct stat s;
	if(! ( ( stat(dpath, &s) == 0 ) && S_ISDIR(s.st_mode) ) )	
		mkdir(dpath, 0777);

  FILE *fp;
  if( (fp = fopen(fullname, mode) ) == NULL )
    err(errno,"fpathopen()::%s",fullname);

  return fp;
}

/*read file list into organized array********
* detect bad path(too long,empty row, not a file et al.) 
**********************/
//20190910, enable option format check for sra accession format
infile_tab_t * organize_infile_list(char* list_path, int fmt_ck) //20190910, enable option format check
{	
		infile_tab_t *infile_stat = malloc(sizeof(infile_tab_t)); 
		int alloc_usize = 1024;
		infile_stat->organized_infile_tab = malloc( sizeof(infile_entry_t) * alloc_usize );
    struct stat path_stat;
		
    FILE * list;
    list = fopen(list_path,"r");
    if(!list) err(errno,"can't open file %s",list_path);
    char *buf = malloc( LMAX * sizeof(char));
    int file_num = 0;

		if(fmt_ck){ //format check
    	while ( (fgets(buf,LMAX,list))!=NULL){
				while(isspace(*buf)) buf++; //rm spaces before path 
      	buf[strcspn(buf, "\r\n")] = 0; //trim \r\n
      	if( strlen(buf) < 1 )
        	continue;
      	if( strlen(buf) > PATHLEN )
        	err(errno,"the input list: %s\n %dth line:  %s  has %lu characters exceed the maximal allowed length %d",
            list_path,file_num, buf,strlen(buf),PATHLEN);
      //reset path_stat everytime before stat new file
      	memset(&path_stat, 0, sizeof path_stat);
      	stat(buf, &path_stat) ;
      	if(!S_ISREG(path_stat.st_mode))
          err(errno,"%dth line: %s",file_num, buf);
				else if(!isOK_fmt_infile(buf,acpt_infile_fmt,ACPT_FMT_SZ)){				
						printf ("isOK_fmt_infile(): wrong format %dth line: %s\nSupported format are:\n",file_num, buf);
						for(int i=0; acpt_infile_fmt[i]!=NULL;i++)
							printf(".%s ",acpt_infile_fmt[i]);
						printf("\n");	
							err(errno,"program exit");
				}
				else {
					infile_stat->organized_infile_tab[file_num].fsize = path_stat.st_size;
					infile_stat->organized_infile_tab[file_num].fpath = malloc(PATHLEN * sizeof(char) );
					strcpy(infile_stat->organized_infile_tab[file_num].fpath, buf);
					file_num++;
					if( file_num >= alloc_usize){
      			alloc_usize+=alloc_usize;
        		infile_stat->organized_infile_tab
                = realloc(infile_stat->organized_infile_tab, sizeof(infile_entry_t) * alloc_usize );
       		}
				}
    	};// while 
		}//if 
		else{ // no fmt check
			while ( (fgets(buf,LMAX,list))!=NULL){
        while(isspace(*buf)) buf++; //rm spaces before path
        buf[strcspn(buf, "\r\n")] = 0; //trim \r\n
        if( strlen(buf) < 1 )
          continue;
        if( strlen(buf) > PATHLEN )
          err(errno,"the input list: %s\n %dth line:  %s  has %lu characters exceed the maximal allowed length %d",
            list_path,file_num, buf,strlen(buf),PATHLEN);

				infile_stat->organized_infile_tab[file_num].fsize = 0;
				infile_stat->organized_infile_tab[file_num].fpath = malloc(PATHLEN * sizeof(char) );
				strcpy(infile_stat->organized_infile_tab[file_num].fpath, buf);
				file_num++;
				if( file_num >= alloc_usize){
            alloc_usize+=alloc_usize;
            infile_stat->organized_infile_tab
                = realloc(infile_stat->organized_infile_tab, sizeof(infile_entry_t) * alloc_usize );
       }
			
		}//while
	} //else


    infile_stat->infile_num = file_num ;
    fclose(list);
		free(buf);
    return infile_stat;
};
//20190910: enhanced option for format check
infile_tab_t * organize_infile_frm_arg (int num_remaining_args, char ** remaining_args,int fmt_ck)
{
//	printf("organize_infile_frm_arg:%d\t%s\n", num_remaining_args,remaining_args[0]);
	infile_tab_t *infile_stat = malloc(sizeof(infile_tab_t));
	int file_num = 0;	
	struct stat path_stat;
	DIR *dirp;
	struct dirent *dirent;
	char fullpath[PATHLEN];
	int alloc_usize = 1024;
	infile_stat->organized_infile_tab = malloc( sizeof(infile_entry_t) * alloc_usize );

if(fmt_ck) { //normal case: files
	for(int i=0;i<num_remaining_args;i++){
		stat(remaining_args[i],&path_stat);
		if( S_ISDIR(path_stat.st_mode)){
				if (( dirp = opendir(remaining_args[i]) ) == NULL)
					err(errno, "%dth argument: can't open %s",i+1, remaining_args[i] );

				while ((dirent = readdir(dirp)) != NULL){
					if(strlen(remaining_args[i]) + strlen(dirent->d_name) + 1 >PATHLEN) 
						err(errno,"path: %s/%s exceed maximal path lenth %d",remaining_args[i], dirent->d_name,PATHLEN);
					sprintf(fullpath, "%s/%s", remaining_args[i], dirent->d_name);
					stat(fullpath,&path_stat);
			
      		if(isOK_fmt_infile(fullpath,acpt_infile_fmt,ACPT_FMT_SZ)){						
						infile_stat->organized_infile_tab[file_num].fsize = path_stat.st_size;
						infile_stat->organized_infile_tab[file_num].fpath = malloc(PATHLEN * sizeof(char));
						sprintf(infile_stat->organized_infile_tab[file_num].fpath, "%s/%s", remaining_args[i], dirent->d_name);
						file_num++;
						if( file_num >= alloc_usize){
							alloc_usize+=alloc_usize;

							infile_stat->organized_infile_tab 
								= realloc(infile_stat->organized_infile_tab,sizeof(infile_entry_t) * alloc_usize );
						} 						
					}
				}
				closedir(dirp);
		}
		else if(isOK_fmt_infile(remaining_args[i],acpt_infile_fmt,ACPT_FMT_SZ)){
			stat(remaining_args[i],&path_stat);		
			infile_stat->organized_infile_tab[file_num].fsize = path_stat.st_size;
			infile_stat->organized_infile_tab[file_num].fpath = malloc(PATHLEN * sizeof(char));
			strcpy(infile_stat->organized_infile_tab[file_num].fpath, remaining_args[i]);
			file_num++;
			if( file_num >= alloc_usize){
      	alloc_usize+=alloc_usize;

        infile_stat->organized_infile_tab 
					= realloc( infile_stat->organized_infile_tab, sizeof(infile_entry_t) * alloc_usize );
      }
		}
		else{
			printf ("wrong format %dth argument: %s\nSupported format are:\n",i+1,remaining_args[i]);
      for(int i=0; acpt_infile_fmt[i]!=NULL;i++)
      	printf(".%s ",acpt_infile_fmt[i]);
      printf("\n");
      err(errno,"program exit");
		}
	};
}// fmt check mode end
else{ // accession case non-file mode
	for(int i=0;i<num_remaining_args;i++){
		infile_stat->organized_infile_tab[file_num].fsize = 0;
		infile_stat->organized_infile_tab[file_num].fpath = malloc(PATHLEN * sizeof(char));
		strcpy(infile_stat->organized_infile_tab[file_num].fpath, remaining_args[i]);
		file_num++;
		if( file_num >= alloc_usize){
        alloc_usize+=alloc_usize;

        infile_stat->organized_infile_tab
          = realloc( infile_stat->organized_infile_tab, sizeof(infile_entry_t) * alloc_usize );
    }
	}
}

	infile_stat->infile_num = file_num ;
	return infile_stat;
};

/*count file for each fmt type*/
infile_fmt_count_t *infile_fmt_count ( infile_tab_t * infile_tab )
{
	infile_fmt_count_t tmp_fmt_count = {0,0,0,0};
	infile_fmt_count_t* fmt_count = (infile_fmt_count_t*)malloc(sizeof(infile_fmt_count_t));
	*fmt_count = tmp_fmt_count;
	for(int i = 0; i < infile_tab->infile_num; i++ )
	{

		if( isOK_fmt_infile(infile_tab->organized_infile_tab[i].fpath,fasta_fmt,FAS_FMT_SZ) )
			fmt_count->fasta++;
		else if (isOK_fmt_infile (infile_tab->organized_infile_tab[i].fpath, fastq_fmt,FQ_FMT_SZ ) )
			fmt_count->fastq++;
		else if(isOK_fmt_infile (infile_tab->organized_infile_tab[i].fpath, co_fmt,CO_FMT_SZ) )
			fmt_count->co++;
		else if(infile_tab->organized_infile_tab[i].fsize != 0) 
			err(errno,"infile_fmt_count(): %s is not accept format(.fasta,.fastq,.co)",infile_tab->organized_infile_tab[i].fpath);
	}
	return fmt_count;
}


int str_suffix_match(char *str, const char *suf)
{
	int ret = 0;
	if( (strlen(str) > strlen(suf)) && (strcmp( (str + strlen(str) - strlen(suf)),suf) == 0) )
		ret = 1;
	return ret;
};

const char * get_pathname(const char *fullpath, const char *suf)
{
	char *pathcp = malloc( strlen(fullpath) + 1) ;
	strcpy(pathcp,fullpath);
	*(pathcp + strlen(pathcp) - strlen(suf)) = '\0';
	return pathcp; 
}


llong find_lgst_primer_2pow(int w)
{
	if( w < 2 || w > 62 ){
		perror("find_1st_primer_after_2pow: argument should between 8 and 62");
		exit(EXIT_FAILURE);
	}

 	llong n = ( 1llu << w ) ; //memory limit for possible subcontext space
	llong hshsz = (llong) ( (double)n * CTX_SPC_USE_L / LD_FCTR) ; 
	printf("w=%d\tspace_sz=%llu\thashsize=%llu\tkmerlimt=%llu\n",w,n,hshsz,(llong)(hshsz*LD_FCTR) ) ; 
	llong i = 3, c ; llong prime = 0;
	for(i = n - 1 ; i > (n >> 1) ; i--) 
	{
		for ( c = 2  ; c <= (int)pow(i+1,0.5) ; c++ )
		{
			if( i%c == 0 )
				break;
		}
		
		if( c*c > i ){
			prime = i;
			break;
		}
	} 
	printf("nearest prime=%llu\n",prime);
	return prime;
}


int nextPrime(int n){
    int j;
    int tag = 0;

    while (1){
        for(j=2;j<=(int)sqrt(n);j++){
            if(n%j == 0){
                tag = 1;
                break;
            }
        }
        if(tag == 1){
            if(n == 0x7FFFFFFF){
                printf("[ERROR] n exceed 0x7FFFFFFF, Can't find a valid prime\n");
                exit(1);
            }
            n++;
            tag = 0;
        }else{
            return n;
        }
    }
}

const char* test_get_fullpath(const char *parent_path, const char *dstat_f)
{
  struct stat path_stat;
  if( stat(parent_path, &path_stat) < 0 )
    err(errno,"test_get_fullpath()::%s",parent_path);
  if( S_ISDIR(path_stat.st_mode) ){
    char* fullpath = malloc(PATHLEN+1);
    sprintf((char *)fullpath,"%s/%s", parent_path, dstat_f);
    FILE *fp;
    if ( (fp = fopen(fullpath,"rb")) != NULL ){
      fclose(fp);
      return fullpath;
    }
    else{
      free((char*)fullpath);
      return NULL;
    }
  }
  else
    return NULL;
};

//how much subctx should longer than
//drlevel compute according to MIN_SUBCTX_DIM_SMP_SZ
int add_len_drlevel2subk(void)
{
  int min_smp_len = 0, min_subctx_dim_smp_sz;
  min_subctx_dim_smp_sz = MIN_SUBCTX_DIM_SMP_SZ;
  while( min_subctx_dim_smp_sz >>= 1){ min_smp_len++; };
  return ceil((float)min_smp_len/4);
};
dim_shuffle_t* read_dim_shuffle_file(char *dim_shuffle_file)
{
  int basename_len = strlen(dim_shuffle_file) - strlen(".shuf") ;
  if( strcmp( (dim_shuffle_file + basename_len),".shuf") !=0 )
    err(errno,"read_dim_shuffle_file(): input file %s is not .shuf file",dim_shuffle_file);

  FILE *shuf_in;
  if( (shuf_in = fopen(dim_shuffle_file,"rb")) == NULL)
    err(errno,"read_dim_shuffle_file(): open file %s failed",dim_shuffle_file);

  dim_shuffle_t *dim_shuffle = malloc( sizeof(dim_shuffle_t) );
  fread(&(dim_shuffle->dim_shuffle_stat),sizeof(dim_shuffle_stat_t),1,shuf_in);

  int shuf_arr_len = 1<< 4*dim_shuffle->dim_shuffle_stat.subk ;

  dim_shuffle->shuffled_dim = malloc( sizeof(int)*shuf_arr_len );
  fread(dim_shuffle->shuffled_dim,sizeof(int)*shuf_arr_len,1,shuf_in);
  fclose(shuf_in);

  return dim_shuffle;
}

int get_hashsz(dim_shuffle_t *dim_shuffle_in )
{
  //int dim_reduce_rate = 1 << 4*dim_shuffle_in->dim_shuffle_stat.drlevel;
  llong ctx_space_sz = 1LLU << 4*( dim_shuffle_in->dim_shuffle_stat.k - dim_shuffle_in->dim_shuffle_stat.drlevel );
  int primer_ind = 4*( dim_shuffle_in->dim_shuffle_stat.k - dim_shuffle_in->dim_shuffle_stat.drlevel ) - CTX_SPC_USE_L - 7;
  // check primer_ind
  if(primer_ind < 0 || primer_ind > 24 ){
    int k_add = primer_ind < 0 ?  (1 + (0-primer_ind)/4)  :  - (1 + ( primer_ind - 24 )/4)  ;
    err(errno,"get_hashsz(): primer_ind: %d out of range(0 ~ 24), by formula:\n"
    "int primer_ind = 4*(opt_val->k - dim_shuffle->dim_shuffle_stat.drlevel) - CTX_SPC_USE_L - 7\n"
    "this might caused by too small or too large k\n"
    "kmer length = %d\n"
    "dim reduction level = %d\n"
    "ctx_space size = %llu\n"
    "CTX space usage limit = %lf\n\n"
    "try rerun the program with option -k = %d",
    primer_ind, dim_shuffle_in->dim_shuffle_stat.k, dim_shuffle_in->dim_shuffle_stat.drlevel,ctx_space_sz,
       (double)1/(1 << CTX_SPC_USE_L),dim_shuffle_in->dim_shuffle_stat.k + k_add );
  };
  int hashsize_get = primer[primer_ind];
  /*fprintf(logfp,"dimension reduced %d\n"
          "ctx_space size=%llu\n"
          "k=%d\n"
          "drlevel=%d\n"
          "primer_ind=%d\n"
          "hashsize=%u\n",
    dim_reduce_rate,ctx_space_sz,dim_shuffle_in->dim_shuffle_stat.k, dim_shuffle_in->dim_shuffle_stat.drlevel, primer_ind, hashsize_get);
*/
  return hashsize_get ;
}
