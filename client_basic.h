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

