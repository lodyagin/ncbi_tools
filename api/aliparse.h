/*=========================================================================*/
/*                                                                         */
/*  aliparse.h -- Header file for aliparse.c                               */
/*                                                                         */
/*=========================================================================*/

/* Defined constants */

#define ALIGN_FILE_BUFFSIZE   2048

#define SEQUENCE_LINE         1
#define DEFINITION_LINE       2
#define OTHER_LINE            3

/* The alignment file when read in is stored in a   */
/* linked list of ValNodes.  These ValNodes contain */
/* one of the following data types depending on the */
/* choice setting.                                  */
/*                                                  */
/*        choice             type of data.ptrvalue  */
/*        ------             ---------------------  */
/*     0 (ALI_UNKNOWN)           NULL               */
/*     1 (ALI_SEQUENCE)          CharPtr            */
/*     2 (ALI_DEFLINE)           CharPtr            */
/*     3 (ALI_SEQLINE_ID)        CharPtr            */
/*     4 (ALI_DEFLINE_ID)        CharPtr            */
/*     5 (ALI_OWNLINE_ID)        CharPtr            */
/*     6 (ALI_OTHER)             CharPtr            */
/*     7 (ALI_MAYBE_SEQUENCE)    CharPtr            */
/*     8 (ALI_MAYBE_ID)          CharPtr            */

#define ALI_UNKNOWN        0
#define ALI_SEQUENCE       1
#define ALI_DEFLINE        2
#define ALI_SEQLINE_ID     3
#define ALI_DEFLINE_ID     4
#define ALI_OWNLINE_ID     5
#define ALI_OTHER          6
#define ALI_MAYBE_SEQUENCE 7
#define ALI_MAYBE_ID       8

#define ERR_SEQ_BEFORE_ID        1
#define ERR_ID_WITHOUT_SEQ       2
#define ERR_SEQ_WITHOUT_ID       3
#define ERR_DUPLICATE_IDS        4
#define ERR_SEQUENCE_TOO_SHORT   5
#define ERR_SEQUENCE_TOO_LONG    6

#define ERR_OUT_OF_MEMORY        7

/* Data structures */

typedef struct _SeqPart
{
  CharPtr  sequence;
  struct _SeqPart PNTR next;
} SeqPart, PNTR SeqPartPtr;

typedef struct _IdInfo
{
  CharPtr    id;
  Int4       length;
  SeqPartPtr sequence;
  CharPtr    defline;
  struct _IdInfo PNTR next;
} IdInfo, PNTR IdInfoPtr;

typedef struct _ErrInfo
{
  Int2       errNum;
  Int2       level;
  CharPtr    info;
  IdInfoPtr  id;
  SeqPartPtr sequence;
  struct _ErrInfo PNTR next;
} ErrInfo, PNTR ErrInfoPtr;

typedef struct _SequenceList
{
  Int2     dataType;
  CharPtr  data;
  Int4     seqLen;
  CharPtr  seqId;
  struct _SequenceList PNTR next;
} SequenceList, PNTR SequenceListPtr;

typedef struct
{
  IdInfoPtr  sequences;
  IdInfoPtr  maybes;
  ErrInfoPtr errors;
} AlignFileData, PNTR AlignFileDataPtr;


/* Function prototypes */

AlignFileDataPtr Ali_ReadByContent (FILE PNTR alignFilePtr);
void             Ali_FreeFileInfo (AlignFileDataPtr fileInfoPtr);
ErrInfoPtr Ali_AddError (ErrInfoPtr PNTR errorList, Int4 iError);
