/* Concatenated C file */
/* Content of client_basic.h */
#ifndef CLIENT_BASIC
#define CLIENT_BASIC
#include <stdbool.h>
#include <stdio.h>
#include <stdarg.h>
#include <stdlib.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h> // open function
#include <unistd.h> //close function

#define _64MASK 0xffffffffffffffffLLU
#define BIT1MASK 0x0000000000000001LLU
#define DEFAULT (-1)
#define PATHLEN 256 // limit for input file path length
#ifndef COMPONENT_SZ
#define COMPONENT_SZ 7 //6 or 7,default unit component dimension size = 16^(COMPONENT_SZ) or 1<<4*(COMPONENT_SZ)
#endif

#ifndef CTX_SPC_USE_L
#define CTX_SPC_USE_L 8 // 8 for small mem.(< 1g), 4 for larger mem.  //ctx space occupy rate limit = 1/(1<<CTX_SPC_USE_L)  
#endif

#define LD_FCTR 0.6 //hash function load factor

#define LMAX 4096
#define MIN_SUBCTX_DIM_SMP_SZ 4096 //256
typedef unsigned long long int llong;


//*****basic function prototypes ****/

FILE * fpathopen (const char *dpath, const char *fname, const char *mode );
/*completeReverse*/
#define SWAP2  0x3333333333333333ULL
#define SWAP4  0x0F0F0F0F0F0F0F0FULL
#define SWAP8  0x00FF00FF00FF00FFULL
#define SWAP16 0x0000FFFF0000FFFFULL
#define SWAP32 0x00000000FFFFFFFFULL
static inline llong crvs64bits(llong n) {

  n = ((n >> 2 ) & SWAP2 ) | ((n & SWAP2 ) << 2 );
  n = ((n >> 4 ) & SWAP4 ) | ((n & SWAP4 ) << 4 );
  n = ((n >> 8 ) & SWAP8 ) | ((n & SWAP8 ) << 8 );
  n = ((n >> 16) & SWAP16) | ((n & SWAP16) << 16);
  n = ((n >> 32) & SWAP32) | ((n & SWAP32) << 32);
  return ~n;
}

/*basemap for different alphabet */
extern const int Basemap[128];
extern const char Mapbase[];
extern const unsigned int primer[25];
llong find_lgst_primer_2pow(int w);
int nextPrime(int);

/*shuf file type*/
typedef struct dim_shuffle_stat
{
  int id; //random shuffle file id
  int k; //half context len
  int subk; //half suncontex len
  int drlevel; // dimension reduction level
} dim_shuffle_stat_t; //dimension shuffle type

typedef struct dim_shuffle
{
  dim_shuffle_stat_t dim_shuffle_stat;
  int *shuffled_dim;
} dim_shuffle_t;


/* orgnized infile table*/

typedef struct infile_entry { 
	size_t fsize; 
	char* fpath;
} infile_entry_t ;

typedef struct infile_tab {
  int infile_num;
  infile_entry_t* organized_infile_tab; // infile_entry_t arr
} infile_tab_t ;

infile_tab_t * organize_infile_list(char* list_path,int fmt_ck);
infile_tab_t * organize_infile_frm_arg (int num_remaining_args, char ** remaining_args,int fmt_ck);

typedef struct co_dirstat
{
  unsigned int shuf_id;
  bool koc; //kmer occurence or not
  int kmerlen; //2*k
  int dim_rd_len; //2*drlevel
  int comp_num; // components number
  int infile_num; // .co file num, .co file with many components only count as 1
  //int comp_sz[comp_num * infile_num]; // do this if storage all .co per mco bin in one file
  llong all_ctx_ct;
} co_dstat_t;

//********** input file formats test ********************//

/*input file formats array size */
#define ACPT_FMT_SZ 7
#define FAS_FMT_SZ 4
#define FQ_FMT_SZ	2
#define CO_FMT_SZ 1
#define CMPRESS_FMT_SZ 2


/* input file formats */
extern const char
*acpt_infile_fmt[ACPT_FMT_SZ],
*fasta_fmt[FAS_FMT_SZ],
*fastq_fmt[FQ_FMT_SZ],
*co_fmt[CO_FMT_SZ],
*compress_fmt[CMPRESS_FMT_SZ];

/* format test fun.*/
#include <string.h>
static inline int isCompressfile(char *fname)
{
	int ret = 0;
	for(int i=0; i < CMPRESS_FMT_SZ;i++)
	{
		int basename_len = strlen(fname) - strlen(compress_fmt[i]);
		if( strcmp( ( fname + basename_len ),compress_fmt[i]) == 0 )
			return 1;
	}
	return ret;
}

static inline int isOK_fmt_infile (char* fname, const char *test_fmt[], int test_fmt_arr_sz)
{
  int ret = 0 ;
  char suftmp[10];

	for(int i=0; i < CMPRESS_FMT_SZ;i++ ){
  	int basename_len = strlen(fname) - strlen(compress_fmt[i]);
  	// if infile suffix with .gz or other compress format
  	if( strcmp( ( fname + basename_len ),compress_fmt[i]) == 0 ){
    	char cp_fname[PATHLEN];
    	strcpy(cp_fname, fname);
    	*(cp_fname + basename_len) = '\0';
    	fname = cp_fname;
			break;
  	};
	};

  for(int i=0; i< test_fmt_arr_sz; i++){
    sprintf(suftmp,".%s",test_fmt[i]);
    if ( strcmp((char *)(fname+strlen(fname) - strlen(suftmp)), suftmp) == 0  ){
      return 1;
    }
  };
  return ret;
};


static inline void check (int test, const char * message, ...)
{
    if (test) {
        va_list args;
        va_start (args, message);
        vfprintf (stderr, message, args);
        va_end (args);
        fprintf (stderr, "\n");
        exit (EXIT_FAILURE);
    }
}; 
/*from command_shuffle*/
dim_shuffle_t* read_dim_shuffle_file(char *dim_shuffle_file);
int add_len_drlevel2subk(void);
int get_hashsz(dim_shuffle_t *);
int str_suffix_match(char *str, const char *suf); 
const char * get_pathname(const char *fullpath, const char *suf);
const char* test_get_fullpath(const char *parent_path, const char *dstat_f);
// infile fmt count struct
typedef struct
{
  int fasta;
  int fastq;
  int co;
	int mco;
} infile_fmt_count_t ;
// infile fmt count function 
infile_fmt_count_t *infile_fmt_count ( infile_tab_t * infile_tab );

extern const char co_dstat[];
extern const char skch_prefix[];
extern const char idx_prefix[];
typedef unsigned int ctx_obj_ct_t;

#define H1(K,HASH_SZ) ((K)%(HASH_SZ))
#define H2(K,HASH_SZ) ( 1 + (K) % ( (HASH_SZ) - 1 ) )
#define HASH(K,I,HASH_SZ) ( ( H1(K,HASH_SZ) + I * H2(K,HASH_SZ) ) % HASH_SZ )

#endif



/* Content of client_seq2co.h */
#ifndef SEQ2CO_H
#define SEQ2CO_H
#include "client_basic.h"
#define READSEQ_BUFFSZ 65536 
//for koc file //kmer occrence
#define OCCRC_BIT 16 //24 // make sure it is 64 bits machine
#define OCCRC_MAX 0xffffLLU // 0xffffffLLU //must be 1<<OCCRC_BIT - 1

extern void seq2co_global_var_initial(void);

llong * fasta2co(char* seqfname,llong *co,char * pipecmd);
llong * fastq2co(char* seqfname, llong *co, char * pipecmd,int Q, int M );
llong * mt_shortreads2koc (char* seqfname, llong *co, char *pipecmd,int p);

llong wrt_co2cmpn_use_inn_subctx(char* cofilename, llong *co);
unsigned int write_fqkoc2files(char* cofilename, llong *co);


#endif


/* Content of all_in_one.c */


/* Content of client_basic.c */
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


/* Content of client_kssd.c */
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
#ifdef __EMSCRIPTEN__
    int num_threads = 4; // Define a default value for WebAssembly
#else
    int num_threads = omp_get_num_procs();
#endif

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


/* Content of client_seq2co.c */
#include "client_seq2co.h"
#include "client_basic.h"
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <zlib.h>
#include <math.h>
#include <string.h>
#include <err.h>
#include <errno.h>
#include <fcntl.h> // open function
#include <unistd.h> //close function

#define HIBITSET1 0x8000000000000000LLU
#define _64MASK 0xffffffffffffffffLLU

static int rand_id;
static int half_ctx_len;
static int half_subctx_len;
static int half_outctx_len;
static int drlevel;

static int comp_bittl;  //64-BITTL ;
static int crvsaddmove ; //bits the new base should right move before add to crvstuple
static llong tupmask; //TUPMASK;
static int TL ;  //TUPLEN ;
//MASK for Dimension Reduction
static  llong domask ;
static  llong undomask ;
//dim shuffled array arg
dim_shuffle_t* dim_shuffle; //dim stat + dim arr
static int *dim_shuf_arr ;
static int dim_shuf_arr_len;
static int dim_start; //reduction dim start
static int dim_end;  // reduction dim end

unsigned int hashsize;
unsigned int hashlimit; //limit context count in hashtable
int component_num;

void seq2co_global_var_initial(void)
{
	rand_id = dim_shuffle->dim_shuffle_stat.id ;  
	half_ctx_len = dim_shuffle->dim_shuffle_stat.k ;
	half_subctx_len = dim_shuffle->dim_shuffle_stat.subk ;
	half_outctx_len = half_ctx_len - half_subctx_len;
  drlevel = dim_shuffle->dim_shuffle_stat.drlevel;
	hashlimit = hashsize * LD_FCTR ;
	printf("rand_id=%d\thalf_ctx_len=%d\thashsize=%d\thashlimit=%d\n",rand_id,half_ctx_len,hashsize,hashlimit);
//	printf("%d\t%d\t%d\n",COMPONENT_SZ,half_ctx_len,drlevel );
	component_num = half_ctx_len - drlevel > COMPONENT_SZ ? 
								 1LU << 4*(half_ctx_len - drlevel - COMPONENT_SZ ) : 1  ; 
//	printf("compnum=%d\t compsize=%d\n",component_num,1<<4*COMPONENT_SZ) ;	
	comp_bittl = 64-4*half_ctx_len;  //64-BITTL ;
	crvsaddmove = 4*half_ctx_len-2; //4*half_ctx_len, bits the new base should right move before add to crvstuple
	tupmask = _64MASK >> comp_bittl; //TUPMASK;
	TL = 2*half_ctx_len ;  //TUPLEN ;

//MASK for Dimension Reduction
	//( (1LLU << ( half_subctx_len*2 ) ) - 1 ) << (2*half_ctx_len + 2);
	domask = ( (1LLU << ( half_subctx_len*4 ) ) - 1 ) << (2*half_outctx_len); 
	undomask =  ( (1LLU << (half_outctx_len*2)) - 1 ) 
									<< (2*(half_ctx_len + half_subctx_len)); //<< (2*(half_ctx_len + half_subctx_len + 1));	

	dim_shuf_arr = dim_shuffle->shuffled_dim;
	dim_shuf_arr_len = 1LLU << (4*half_subctx_len) ;
	dim_start = 0;
	dim_end = MIN_SUBCTX_DIM_SMP_SZ ;
};

const char gzpipe_cmd[]= "zcat -fc"; //const char gzpipe_cmd[]= "unpigz -fc";

llong * fasta2co(char* seqfname, llong *co, char * pipecmd) //20190910, enhancement: pipecmd 
{
	llong tuple = 0LLU, crvstuple = 0LLU,
    unituple, drtuple, pfilter;

	memset(co,0LLU,hashsize*sizeof(llong));		
	char seqin_buff[ READSEQ_BUFFSZ + 1 ]; // seq infile buffer
	FILE *infp;
	char fas_fname[PATHLEN];
	if(pipecmd[0] != '\0')
		sprintf(fas_fname,"%s %s",pipecmd,seqfname);//other pipecmd except decompress cmd
	else
		sprintf(fas_fname,"%s %s",gzpipe_cmd,seqfname);

	if( (infp=popen(fas_fname,"r")) == NULL ) err(errno,"fasta2co():%s",fas_fname);

	int newLen =  fread(seqin_buff, sizeof(char),READSEQ_BUFFSZ,infp);
	if(! (newLen >0) ) 	err(errno,"fastco():eof or fread error file=%s",seqfname);

	llong base = 1; char ch; int basenum;
	unsigned int keycount = 0;
	// core function begin
	for(int pos = 0; pos <= newLen; pos++) 
	{	
		if(pos == newLen){
				newLen =  fread(seqin_buff, sizeof(char),READSEQ_BUFFSZ,infp);	
			if(newLen > 0)
				pos = 0;
			else break;
		};
		ch = seqin_buff[pos];
		basenum = Basemap[(int)ch];

		if(basenum != DEFAULT) //make sure basenum is not negative
		{
			tuple = ( ( tuple<< 2 ) | (llong)basenum ) & tupmask ;
			crvstuple = ( crvstuple >> 2 ) + (((llong)basenum^3LLU) << crvsaddmove);
			base++;
		}
		else if ( (ch == '\n') || (ch == '\r') ) { continue;}
		else if (isalpha(ch)){ base=1; continue; }
		else if ( ch == '>' )
		{
			while( (pos < newLen ) && ( seqin_buff[pos] != '\n' ) )
			{			
				if (pos < newLen - 1)
					pos++;
				else
				{
					newLen = fread(seqin_buff, sizeof(char),READSEQ_BUFFSZ,infp);
					if(newLen > 0) pos = -1; // pos = 0 is bug ? should be pos = -1 , since continue and pos++
					else err(errno,"fasta2co(): can not find seqences head start from '>' %d",newLen);
				};
			};
      base = 1;
      continue;
		}
		else {
    //  warnx("ignorn illegal charactor%c \n",ch);
      base=1;
      continue;
    };

		if( base > TL ) // if(!(base < TL))
		{
			//make sure unituple == min(tuple,crvstuple);
			unituple = tuple < crvstuple ? tuple:crvstuple;
			//only for 64bits Kmer storage ,
    	//important !!!! make sure is right
			int dim_tup = ((unituple & domask) >> ( (half_outctx_len)*2 ) ) ;
			pfilter = dim_shuf_arr[dim_tup];
			if( ( pfilter >= dim_end) || (pfilter < dim_start ) ) continue;
			pfilter = pfilter - dim_start; //add 190307: prevent pfilter > MIN_SUBCTX_DIM_SMP_SZ when dim_start >0	
			drtuple = ( ( (unituple & undomask) //left half outer ctx
							+ ( ( unituple & ( ( 1LLU<< ( half_outctx_len*2) ) - 1)) << (TL*2 - half_outctx_len*4) ) )										
            	>> ( drlevel*4 ) ) // subctx dim reduced
            	+  pfilter ; //  subctx dim after reduction

			unsigned int i,n ;
    	for(i=0;i<hashsize;i++)
    	{
      	n = HASH(drtuple,i,hashsize);
      	if (co[n] == 0)
      	{
        	co[n] = drtuple;
        	keycount++;
        	if( keycount > hashlimit)
          	err(errno,"the context space is too crowd, try rerun the program using -k%d", half_ctx_len + 1);
        	break;
      	}
      	else if ( co[n] == drtuple )
          break;
			};//end kmer hashing
		};
	};//end file hashing
		pclose(infp);
  return co;
};//end func
// macro for fastq2co
// Q:quality score, M: required least occurence of kmer
// make sure LEN is enough for long read
#define LEN 20000 //4096, 20220720 update for pacbio reads
#define CT_BIT 4 //bits for Kmer count
#define CT_MAX 0xfLLU //make sure smaller than 1LLU<<CT_BTI

llong * fastq2co(char* seqfname, llong *co, char *pipecmd, int Q, int M ) //20190910 enhanced pipecmd
{
	if(M >= CT_MAX) err(errno,"fastq2co(): Occurence num should smaller than %d", (int)CT_MAX);

	llong tuple = 0LLU, crvstuple = 0LLU, unituple, drtuple, pfilter;	
	memset(co,0LLU,hashsize*sizeof(llong));

	FILE *infp;
	char fq_fname[PATHLEN];
	if(pipecmd[0] != '\0')
    sprintf(fq_fname,"%s %s",pipecmd,seqfname);//other pipecmd except decompress cmd
  else
    sprintf(fq_fname,"%s %s",gzpipe_cmd,seqfname);

	if( (infp=popen(fq_fname,"r")) == NULL ) err(errno,"fastq2co():%s",fq_fname);
	//char seq[LEN];
	//char qual[LEN];
	char *seq = malloc(LEN+10);
	char *qual = malloc(LEN+10);
	fgets(seq,LEN,infp); fgets(seq,LEN,infp); 
	fgets(qual,LEN,infp);  fgets(qual,LEN,infp);
//	if(Q != -1) { 		};
	llong base = 1; char ch ; int basenum,line_num = 0 ;
	unsigned int keycount =0 ;
	int sl = strlen(seq); 
	for(int pos = 0; pos < sl; pos++){
		if(seq[pos] == '\n' ){
			fgets(seq,LEN,infp); fgets(seq,LEN,infp);
			fgets(qual,LEN,infp); fgets(qual,LEN,infp);
			sl = strlen(seq);
//			printf("seq=%s\nqual=%s\n==\n",seq,qual); 
			line_num+=4;
			if( !feof(infp) ) {
				base = 1; 
				pos = -1;
				continue ;									
			}
			else break;	
		}		
		else{
			ch = seq[pos];
			basenum = Basemap[(int)ch];
			if( (basenum != DEFAULT) && ( qual[pos] >= Q ) ){
				tuple = ( ( tuple<< 2 ) | (llong)basenum ) & tupmask ;
				crvstuple = ( crvstuple >> 2 ) + (((llong)basenum^3LLU) << crvsaddmove);
				base++;
			} 

			else {
			//	if (qual[pos] < Q)
			//		warnx("in file %s line %d skip low quality charactor '%c[%c]'\n",seqfname, line_num+2 ,ch,qual[pos]);
		 //	else warnx("in file %s line %d ignorn illegal charactor '%c'\n",seqfname,line_num+2,ch);
				base = 1;
				continue;
			};
		};
	
		if( base > TL ){ // serious bug ?!!!:  if( base >= TL ) 

			unituple = tuple < crvstuple ? tuple:crvstuple;
			int dim_tup = ((unituple & domask) >> ( (half_outctx_len)*2 ) ) ;
      pfilter = dim_shuf_arr[dim_tup];
      if( ( pfilter >= dim_end) || (pfilter < dim_start ) ) continue;
			pfilter = pfilter - dim_start;
      drtuple = ( ( (unituple & undomask) //left half outer ctx
              + ( ( unituple & ( ( 1LLU<< ( half_outctx_len*2) ) - 1)) << (TL*2 - half_outctx_len*4) ) )
              >> ( drlevel*4 ) ) // subctx dim reduced
              +  pfilter ;
			
			unsigned int i,n ;
			for(i=0;i<hashsize;i++)	{
				n = HASH(drtuple, i, hashsize);
				if (co[n] == 0LLU){
					//20190910:bug fixed, otherwise occrence==1 k-mer would be skipped
					if( M == 1) co[n] = (drtuple << CT_BIT) | CT_MAX; //--
					else co[n] = (drtuple << CT_BIT) + 1LLU;					
          if( keycount > hashlimit)
            err(errno,"the context space is too crowd, try rerun the program using -k%d", half_ctx_len + 1);
					break;
				} 
				else if ( ( co[n] >> CT_BIT ) == drtuple ) {
					//20190910:bug fixed, test if already to CT_MAX
	          if( (co[n] & CT_MAX) == CT_MAX ) break;
						co[n] += 1LLU;
						if( !((co[n] & CT_MAX) <  M) ) co[n]|= CT_MAX ;
						break ;					
        };
			}; //end kmer hashing
		};
	}// end file hashing 
	printf("%d reads detected\n",line_num);

	pclose(infp);
	free(seq);
  free(qual);
	return co;
}; // end fastq2co()

// fastq2co with occurence
// higher 48 bits for Kmer, lower 16 bits for count
#define OCCRC_BIT 16 // make sure it is 64 bits machine
#define OCCRC_MAX 0xffffLLU //must be 1<<OCCRC_BIT - 1

unsigned int write_fqkoc2files(char* cofilename, llong *co)
{
  int comp_code_bits = half_ctx_len - drlevel > COMPONENT_SZ ? 4*(half_ctx_len - drlevel - COMPONENT_SZ ) : 0  ;
  FILE **outf,**abdf; //sketch and occurence files
  outf = malloc(component_num * sizeof(FILE *));
	abdf = malloc(component_num * sizeof(FILE *));
  char cofilename_with_component[PATHLEN];

  for(int i=0;i<component_num ;i++)
  {
    sprintf(cofilename_with_component,"%s.%d",cofilename,i);
    if ( (outf[i] = fopen(cofilename_with_component,"wb") ) == NULL )
      err(errno,"write_fqkoc2files()") ;
	
		sprintf(cofilename_with_component,"%s.%d.a",cofilename,i);
		if ( (abdf[i] = fopen(cofilename_with_component,"wb") ) == NULL )
      err(errno,"write_fqkoc2files()") ;
  };

  unsigned int count, wr = 0, newid,compi; //compi: component index
	unsigned short abdc; //abundance 
		

  for(count=0;count < hashsize; count++)
  {
    if( co[count] > 0 ) {
			
			compi = (co[count] >> OCCRC_BIT ) % component_num;

      newid = (unsigned int)(co[count] >> (comp_code_bits + OCCRC_BIT));
      fwrite( &newid, sizeof(newid),1,outf[compi] );
			
			abdc = co[count] & OCCRC_MAX;
			fwrite( &abdc, sizeof(abdc),1,abdf[compi] );
		
      wr++;
    }
  }

  for(int i=0;i<component_num ;i++){
    fclose(outf[i]);
		fclose(abdf[i]);
	}
  free(outf);
	free(abdf);
  return wr;
};

//write co to compnent files
//component determeined by modula context,namely use inner subctx
llong wrt_co2cmpn_use_inn_subctx(char* cofilename, llong *co)
{	
	int comp_code_bits = half_ctx_len - drlevel > COMPONENT_SZ ? 4*(half_ctx_len - drlevel - COMPONENT_SZ ) : 0  ;  
	FILE **outf;
	outf = malloc(component_num * sizeof(FILE *));
	char cofilename_with_component[PATHLEN];
  for(int i=0;i<component_num ;i++)
  {
    sprintf(cofilename_with_component,"%s.%d",cofilename,i);
    if ( (outf[i] = fopen(cofilename_with_component,"wb") ) == NULL )
      err(errno,"wrt_co2cmpn_use_inn_subctx()") ;
  };
  unsigned int count, wr = 0, newid;
	for(count=0;count < hashsize; count++)
	{
		if( co[count] != 0 )
		{	
			newid = (unsigned int)( co[count] >> comp_code_bits ) ;
			fwrite( &newid, sizeof(unsigned int),1,outf[(int)( co[count] % component_num )] );
			wr++;
		}
	}
	for(int i=0;i<component_num ;i++)
   fclose(outf[i]);

	free(outf);
  return wr;
};


#define THREAD_MAX 65536
#define FQ_LEN 4096
llong * mt_shortreads2koc (char* seqfname, llong *co, char *pipecmd,int p){
printf("runing mt_shortreads2koc()\n");
char (*fq_buff)[FQ_LEN] = malloc( THREAD_MAX * FQ_LEN );
char tmp[FQ_LEN];		
int l;
unsigned int keycount =0 ;
  memset(co,0LLU,hashsize*sizeof(llong));
	FILE *infp;
  char fq_fname[PATHLEN];
  if(pipecmd[0] != '\0') sprintf(fq_fname,"%s %s",pipecmd,seqfname);
	else sprintf(fq_fname,"%s %s",gzpipe_cmd,seqfname);
//	else sprintf(fq_fname,"%s",seqfname);
  if( (infp=popen(fq_fname,"r")) == NULL ) err(errno,"mtfastq2koc():%s",fq_fname);


	while (!feof(infp)){
		for (l = 0; l < THREAD_MAX && fgets(tmp,FQ_LEN,infp) && fgets(fq_buff[l],FQ_LEN,infp) && fgets(tmp,FQ_LEN,infp) && fgets(tmp,FQ_LEN,infp); l++) ;

#pragma omp parallel for num_threads(p) schedule(guided)
		for (int t = 0 ; t < l ; t++ ){
			int base = 1; char ch;
      llong tuple = 0LLU;
      llong crvstuple = 0LLU;
      llong unituple = 0LLU;
			
			for(int pos = 0; (ch = fq_buff[t][pos]) != '\n'; pos++){
				int basenum = Basemap[(int)ch];
				if(basenum != DEFAULT ){
					tuple = ( ( tuple<< 2 ) | (llong)basenum ) & tupmask ;
					crvstuple = ( crvstuple >> 2 ) + (((llong)basenum^3LLU) << crvsaddmove);
					base++;
				}else{ base = 1;continue;}

				if( base > TL ){
        	unituple = tuple < crvstuple ? tuple:crvstuple;
      		unsigned int dim_tup = ((unituple & domask) >> ( (half_outctx_len)*2 ) ) ;
					llong pfilter = dim_shuf_arr[dim_tup];
					if( ( pfilter >= dim_end) || (pfilter < dim_start ) ) continue;
					pfilter = pfilter - dim_start;
					llong drtuple = ( ( (unituple & undomask)
						+ ( ( unituple & ( ( 1LLU<< ( half_outctx_len*2) ) - 1)) << (TL*2 - half_outctx_len*4) ) )
						>> ( drlevel*4 ) )
						+  pfilter ;

					for(unsigned int i=0;i<hashsize;i++){
						unsigned int n = HASH(drtuple, i, hashsize);
						if (co[n] == 0LLU){
#pragma omp atomic write
						co[n] = (drtuple << OCCRC_BIT) + 1LLU ;
#pragma omp atomic	
							keycount++;
							if( keycount > hashlimit )
								 err(errno,"the context space is too crowd, try rerun the program using -k%d", half_ctx_len + 1);
							break;
						}else if ( ( co[n] >> OCCRC_BIT ) == drtuple ) {
							if( (co[n] & OCCRC_MAX) < OCCRC_MAX )
#pragma omp atomic
								 co[n]+=1LLU;
							 break ;			
						}

					}// end kmer hashing for								
				}
			}// end per read for
		}// end threads block for
	}// end file while

	pclose(infp);
	free(fq_buff);
	return co;
}// end mt_shortreads2koc()




















