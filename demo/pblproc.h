#ifndef _PBLPROC_
#define _PBLPROC_

#include <ncbi.h>
#include <sequtil.h>
#include <tofasta.h>
#include <netblap2.h>
#include <dust.h>

#define BLAST_DEFAULT   0
#define BLASTN_PROGRAM  1
#define BLASTP_PROGRAM  2
#define BLASTX_PROGRAM  4
#define TBLASTN_PROGRAM 8

#define MIN_BLASTN_LEN 11
#define MIN_BLASTP_LEN 3
#define MAX_BLASTN_LEN 8000       /*if it is greater than MAX, break up */
#define MAX_BLASTX_LEN 3000
#define MAX_TBLASTN_LEN 2000
#define MAX_BLASTP_LEN 4000
#define OVERLAP_SPACE 1000


/*the different possible search databases*/
#define BLAST_NR        ((Uint4)1)
#define BLAST_EST       ((Uint4)2)
#define BLAST_STS       ((Uint4)4)
#define BLAST_MONTH     ((Uint4)8)
#define BLAST_HTGS	((Uint4)16)
#define BLAST_VECTOR    ((Uint4)32)
#define BLAST_MITO      ((Uint4)64)
#define BLAST_KABAT     ((Uint4)128)
#define BLAST_SWISSPROT ((Uint4)256)
#define BLAST_PDB       ((Uint4)512)
#define BLAST_EPD	((Uint4)1024)
#define BLAST_YEAST	((Uint4)2048)
#define BLAST_GSS	((Uint4)4096)
#define BLAST_ALU	((Uint4)8192)
#define BLAST_THC       ((Uint4)16384)

#define DO_BLAST_N 1
#define DO_BLAST_P 2
#define DO_BLAST_X 3
#define DO_T_BLAST_N 4


#define OUTPUT_TEXT     1
#define OUTPUT_SEQALIGN 2
#define OUTPUT_SEQENTRY 4
#define OUTPUT_HTML		8


#define SEQFMT_FASTA	1
#define SEQFMT_ACCESSION	2
#define SEQFMT_GI	3


typedef struct blast_parameter {
	Uint1 blast_program;	/*0 for default */
	CharPtr n_param;	/*parameter for blastn*/
	CharPtr p_param;	/*parameter for blastp */
	CharPtr t_param;	/*parameter for tblastn */
	CharPtr x_param;	/*parameter for blastx */

	Uint4	max_blast_n;
	Uint4	max_blast_p;	/*maximum allowed length for a single search */
	Uint4	max_blast_x;
	Uint4	max_blast_t;

	Uint4   dna_db;		/*the public database */
	CharPtr other_dna;	/*other dna database */
	Uint4	prot_db;	/*the protein database in the public */
	CharPtr other_protein;	/*other protein database */
}BlastParameter, PNTR BlastParamPtr;
void free_blast_param PROTO((BlastParamPtr bpp));

typedef struct pblastoption {
	CharPtr seq_data;
	Uint1 seq_format;
	Char file_name[PATH_MAX];
	BlastParamPtr temp_bpp;
	BlastParamPtr bpp;
	Char repeat_library[PATH_MAX];
	ValNodePtr alu_sep_list;
	ValNodePtr alu_slp_list;
	Uint1 gap_alignment;
	Uint1 filter_org;	/*0=None 1=include 2=exclude */
	CharPtr organism;
	ByteStorePtr aa_gi_bsp;	/*the gis for protein record of the organism */
	ByteStorePtr na_gi_bsp;	/*the gis for nucleotide record of the organism */
	Uint1 output_format;
	Char output_path[PATH_MAX];	/*the path for the output file */
	Boolean dust;		/*dust the sequence? */
	Boolean filter_self;	/* filter the GenBank query itself */
	ValNodePtr output_file;	/*name of the output files */
	Boolean map_align_to_feature;
	Int4 seq_count;
	Boolean monitor;	/*activate the Monitor*/
	Boolean hide_alignment_feature;	/*do NOT show the alignment feature */
	FILE *errfp;
}PBlastOption, PNTR PBlastOptionPtr;

void FreePBlastOption PROTO((PBlastOptionPtr pbop));

/***************************************************************
*
*	some powblast operations, like the repeat filtering 
*	and organism filtering will be operated on all the 
*	sequences. InitPowerBlastOption will load the repeat 
*	library and find if there is any errors in organism 
*	field and decide weather the organism filtering is done 
*	with Eval or the other way 
*
***************************************************************/
Boolean InitPowerBlastOption PROTO((PBlastOptionPtr pbop));


#define POWBLAST_NONE   0       /*reset everything */
#define POWBLAST_HIT    1       /*a successful run */
#define POWBLAST_FATAL  2       /* need to exit the program */
/***********************************************************************
*
*       Run one session of powblast. return POWBLAST_NONE for no hit
*       return POWBLAST_HIT to indicate there is a hit
*       return POWBLAST_FATAL to stop the process and indicate there is
*       fatal error
*
***********************************************************************/
Uint1 RunPowerBlast PROTO((PBlastOptionPtr pbop));


void load_default_blast_param PROTO((PBlastOptionPtr pbop));
/****************************************************************
*
*	write the current parameter to the configure file 
*	.powblastrc
*
*
****************************************************************/
Boolean load_pboption_to_configure PROTO((PBlastOptionPtr pbop));


#endif


