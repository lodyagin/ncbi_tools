/* $Id: fdfilter.c,v 6.20 2001/06/22 20:12:52 shavirin Exp $ */
/*****************************************************************************

  
                          PUBLIC DOMAIN NOTICE
              National Center for Biotechnology Information

    This software/database is a "United States Government Work" under the
    terms of the United States Copyright Act.  It was written as part of    
    the author's official duties as a United States Government employee
    and thus cannot be copyrighted.  This software/database is freely
    available to the public for use. The National Library of Medicine and
    the U.S. Government have not placed any restriction on its use or
    reproduction.

    Although all reasonable efforts have been taken to ensure the accuracy
    and reliability of the software and data, the NLM and the U.S.
    Government do not and cannot warrant the performance or results that
    may be obtained by using this software or data. The NLM and the U.S.
    Government disclaim all warranties, express or implied, including
    warranties of performance, merchantability or fitness for any
    particular purpose.

    Please cite the author in any work or product based on this material.

   ***************************************************************************

   File Name:  fdfilter.c

   Author:  Sergei B. Shavirin
   
   Version Creation Date: 05/21/99

   $Revision: 6.20 $

   File Description:  Filter FASTA databases for identical sequences

   $Log: fdfilter.c,v $
   Revision 6.20  2001/06/22 20:12:52  shavirin
   Fixed problem with incorrect sorting of definition lines.

   Revision 6.19  2001/06/15 21:00:32  shavirin
   Added possibility to create sub-database of Taxonomy names for taxonomy
   ids used in Blast database for custom formating of results.

   Revision 6.18  2001/06/13 20:43:29  shavirin
   Fixed function GetPDBRankByGi() for absence of pdbrank file.

   Revision 6.17  2001/05/31 14:42:26  shavirin
   Added possibility to create LocusLink and UniGene affiliations
   for ASN.1 structured definition lines.

   Revision 6.16  2001/05/23 21:16:22  shavirin
   Added population of tax_id and membership fields for sequence to
   database affiliation.

   Revision 6.15  2001/05/22 17:32:12  shavirin
   Fixed few UMR purify errors.

   Revision 6.14  2001/05/14 17:37:28  shavirin
   Added possibility to manipulate with BLAST databases with ASN.1 structured
   deflines.

   Revision 6.13  2001/05/08 22:04:22  shavirin
   Changed usage of function FDBAddSequence() in accordance with it's
   changed definition.

   Revision 6.12  2000/06/07 19:18:42  shavirin
   Added rank sorting of PDB sequences using rank file (optional)

   Revision 6.11  2000/03/14 16:59:32  shavirin
   For equivalent databases added sorting by gi.

   Revision 6.10  2000/02/01 17:19:23  egorov
   Couple functions changed prototypes, so the adjustment made.

   Revision 6.9  1999/09/20 18:35:03  shavirin
   Added dumping of deflines in order defined in tofasta.c file.

   Revision 6.7  1999/09/10 17:46:39  shavirin
   Added parameter for the type of the database and added possibility to
   filter nucleotide databases.

   Revision 6.6  1999/08/26 18:40:58  shavirin
   Added possibility to make subset of database by query string.

   Revision 6.5  1999/08/26 14:11:33  shavirin
   Added parameter -s for creation of sparse indexes in the output.

   Revision 6.4  1999/08/24 14:28:52  shavirin
   Added possibility to create dump info file in filtered database.

   Revision 6.3  1999/05/27 15:53:54  shavirin
   Added protection against non-existent gis.

   Revision 6.2  1999/05/27 14:43:03  shavirin
   Added functionality to filter BLAST database by list of gis.

   Revision 6.1  1999/05/21 21:11:51  shavirin
   Initial revision.

 *
*****************************************************************************/
#include <ncbi.h>
#include <readdb.h>
#include <ncbiwww.h>

#ifdef TAX_CS_LOOKUP
#include <taxblast.h>
#endif

/* We will use regular WWW encoding of request to make specific
   filtering using database information file 

   The most general query string will look like:

   taxid=555,666,777&owner=5,6,7&div=AAA,BBB,CCC

   */

#define TAXID_LABEL  "taxid"
#define OWNER_LABEL  "owner"
#define DIV_LABEL    "div"
#define DBNAME_LABEL "dbname"

#define MAX_SR_ELEMENTS 36

typedef struct _SR_Info
{
    Int4 taxid[MAX_SR_ELEMENTS];
    Int4 owner[MAX_SR_ELEMENTS];
    Char div[MAX_SR_ELEMENTS][4];
    Char dbname[64];
} SR_Info, PNTR SR_InfoPtr;

typedef struct gilist
{
    Int4    count;
    Int4    allocated;
    Int4Ptr seq_num;
} GiList, *GiListPtr;

typedef struct _DefLine 
{
    Int4 type;
    Int4 gi;
    Int4 pdb_rank;
    CharPtr line;
    BlastDefLinePtr bdp;

} DefLine, PNTR DefLinePtr;

typedef struct deflist
{
    Int4       count;
    Int4       allocated;
    DefLinePtr PNTR defs ;
} DefList, *DefListPtr;

#define DEF_ALLOCATE_CHUNK 4

typedef struct hashelm
{
    Int4  hash;
    Int4  seq_num;
    Int4  gi;
    Int4  tax_id;
    Uint1 owner;
    Char  div[3];
    Int4  SequenceLen;
    Int4  date;
} HashElm, *HashElmPtr;

typedef struct HashTable
{
    Int4 count;
    Int4 allocated;
    HashElmPtr hep;
    
} HashTable, *HashTablePtr;

/* Structures used for special PDF sorting of deflines */
typedef struct pdbelm
{
    Int4  gi;
    Int4  group;
    Int4  rank;
} PDBElm, *PDBElmPtr;

typedef struct PDBTable
{
    Int4 count;
    Int4 allocated;
    PDBElmPtr pelms;
} PDBTable, *PDBTablePtr;

/* Structures used for LocuisLink or UniGene identification */
#define IS_UG_ELEMENT 0x1
#define IS_LL_ELEMENT 0x2

typedef struct llugelm
{
    Int4  gi;
    Int4  id;                   /* LL id or UG cluster id */
    Int4  link;                 /* LL or UG ? */
} LLUGElm, *LLUGElmPtr;

typedef struct LLUGTable
{
    Int4 count;
    Int4 allocated;
    LLUGElmPtr pelms;
} LLUGTable, *LLUGTablePtr;

#define NUMARG (sizeof(flt_args)/sizeof(flt_args[0]))

Args flt_args[] = {
    { "Title for output database file",                 /* 0 */
      NULL, NULL, NULL, TRUE, 't', ARG_STRING, 0.0, 0, NULL},
    {"Input file for the filtering (this parameter must be set)", /* 1 */
     NULL, NULL,NULL,FALSE,'i',ARG_FILE_IN, 0.0,0,NULL},
    {"Input file with a list of gis",                   /* 2 */
     NULL, NULL,NULL,TRUE,'g',ARG_FILE_IN, 0.0,0,NULL},
    {"Logfile name:",                                   /* 3 */
     "fdfilter.log", NULL,NULL,TRUE,'l',ARG_FILE_OUT, 0.0,0,NULL},
    {"Create sparse indexes in the filtered database",  /* 4 */
     "F", NULL,NULL,TRUE,'s',ARG_BOOLEAN, 0.0,0,NULL},
    {"Query string for creating database subset",       /* 5 */
     NULL, NULL,NULL,TRUE,'q',ARG_STRING, 0.0,0,NULL},
    {"Input database is proten",                        /* 6 */
     "T", NULL,NULL,TRUE,'p',ARG_BOOLEAN, 0.0,0,NULL},
    {"Reverse query request to negative",               /* 7 */
     "F", NULL,NULL,TRUE,'r',ARG_BOOLEAN, 0.0,0,NULL},
    {"Input file with PDB ranks",                       /* 8 */
     NULL, NULL,NULL,TRUE,'y',ARG_FILE_IN, 0.0,0,NULL},
    {"Input file with UG memberships",                  /* 9 */
     NULL, NULL,NULL,TRUE,'u',ARG_FILE_IN, 0.0,0,NULL},
    {"Input file with LL memberships",                  /* 10 */
     NULL, NULL,NULL,TRUE,'k',ARG_FILE_IN, 0.0,0,NULL},
};

#define FLT_Input    flt_args[1].strvalue
#define IS_Protein   flt_args[6].intvalue
#define ReverseQuery flt_args[7].intvalue

static void FDB_optionsFree(FDB_optionsPtr options)
{
    if(options == NULL)
        return;
    
    MemFree(options->db_title);
    MemFree(options->db_file);
    MemFree(options->LogFileName);
    MemFree(options);
    
    return;
}

static Boolean SRReadCharData(CharPtr buffer, CharPtr PNTR div_in)
{
    CharPtr tmp, ch, ch2;
    Int4 j;
    Char div[MAX_SR_ELEMENTS][4];

    tmp = StringSave(buffer);
    MemSet(div, NULLB, sizeof(div));

    for(ch2 = tmp, j = 0; j < MAX_SR_ELEMENTS; j++) {
        
        if((ch = StringChr(ch2, ',')) == NULL) {
            StringNCpy(div[j], ch2, 3);
            break;
        }
        
        *ch = NULLB;
        ch++;
        StringNCpy(div[j], ch2, 3);
        ch2 = ch;
    }

    MemCpy(div_in, div, sizeof(div));

    MemFree(tmp);
    
    return TRUE;
}

static Boolean SRReadIntData(CharPtr buffer, Int4Ptr id)
{
    CharPtr tmp, ch, ch2;
    Int4 j;
    tmp = StringSave(buffer);
    
    for(ch2 = tmp, j = 0; j < MAX_SR_ELEMENTS; j++) {
        
        if((ch = StringChr(ch2, ',')) == NULL) {
            id[j] = atol(ch2);
            id[j+1] = -1; /* Terminating character */
            break;
        }
        
        *ch = NULLB;
        ch++;
        id[j] = atol(ch2);
        ch2 = ch;
    }

    MemFree(tmp);
    
    return TRUE;
}

static SR_InfoPtr SRReadSRInfo(CharPtr buffer)
{
    WWWInfoPtr info;
    WWWInfoDataPtr info_data;
    SR_InfoPtr srip;
    CharPtr chptr;

    info_data = (WWWInfoDataPtr) MemNew(sizeof(WWWInfoData));
    info_data->query = StringSave(buffer);
    info_data->entries = WWWGetEntries(&info_data->num_entries, 
                                       info_data->query, FALSE);
    info = (VoidPtr) info_data;

    srip = MemNew(sizeof(SR_Info));
    
    if((chptr = WWWGetValueByName(info, TAXID_LABEL)) != NULL)
        SRReadIntData(chptr, srip->taxid);
    if((chptr = WWWGetValueByName(info, OWNER_LABEL)) != NULL)
        SRReadIntData(chptr, srip->owner);
    if((chptr = WWWGetValueByName(info, DIV_LABEL)) != NULL)
        SRReadCharData(chptr, (CharPtr PNTR) srip->div);
    if((chptr = WWWGetValueByName(info, DBNAME_LABEL)) != NULL)
        StringCpy(srip->dbname, chptr);

    /* If no conditions exists - this is an error */
    if(srip->taxid[0] <= 0 && srip->owner[0] <= 0 && **srip->div == NULLB) {
        ErrPostEx(SEV_ERROR, 0, 0, 
                  "No valid conditions exists in query string");
        MemFree(srip);
        srip = NULL;
    }
    WWWInfoFree(info);
    return srip;
}

static FDB_optionsPtr FDB_CreateCLOptions(void)
{
    FDB_optionsPtr options;
    Char buffer[128];

    options = MemNew(sizeof(FDB_options));
    
    if ( !GetArgs ("fdfilter", NUMARG, flt_args) )
        return NULL;
    
    if ( !ErrSetLog (flt_args[3].strvalue) ) { /* Logfile */
        ErrShow();
    } else {
        ErrSetOpts (ERR_CONTINUE, ERR_LOG_ON);
    }
    
    options->db_title = StringSave(flt_args[0].strvalue);

    sprintf(buffer, "%s.flt", FLT_Input);
    options->db_file = StringSave(buffer);

    options->LogFileName = StringSave(flt_args[3].strvalue);
    options->is_protein = IS_Protein;
    options->parse_mode = TRUE;
    options->isASN = FALSE;
    options->asnbin = FALSE;
    options->is_seqentry = FALSE;
    options->base_name = NULL;
    options->dump_info = FALSE;
    options->sparse_idx = flt_args[4].intvalue;

#if 0
    if(flt_args[9].intvalue)
        options->version = FORMATDB_VER;
    else
        options->version = FORMATDB_VER_TEXT;
#endif

    return options;
}
#define HASH_ALLOC_CHUNK 1024
#define GI_ALLOC_CHUNK 1024
#define PDB_ALLOC_CHUNK 1024
#define LLUG_ALLOC_CHUNK 1024

static int GICompare(VoidPtr i, VoidPtr j)
{
    if (*(Int4Ptr)i > *(Int4Ptr)j)
        return (1);
    if (*(Int4Ptr)i < *(Int4Ptr)j)
        return (-1);
    return (0);
}

static void FDBAddNewPDB(PDBTablePtr pdbp, Int4 gi, Int4 group, Int4 rank)
{
    PDBElmPtr pelm;
    
    /* Reallocate if necessary */
    
    if(pdbp->allocated <= pdbp->count) {
        pdbp->allocated += PDB_ALLOC_CHUNK + 1;
        pdbp->pelms = (PDBElmPtr) 
            Realloc (pdbp->pelms, pdbp->allocated * sizeof(PDBElm));
    }
    
    pelm = &pdbp->pelms[pdbp->count];
    
    pelm->gi      = gi;
    pelm->group = group;
    pelm->rank = rank;
    
    pdbp->count++;
    
    return;
}

static PDBTablePtr FDBPDBTableNew(void)
{
    PDBTablePtr pdbp;
    
    pdbp = MemNew(sizeof(PDBTable));
    pdbp->allocated = PDB_ALLOC_CHUNK;
    pdbp->pelms = (PDBElmPtr) MemNew (pdbp->allocated * sizeof(PDBElm));
    
    return pdbp;
}

static void FDBDestroyPDBIndex(PDBTablePtr pdbp)
{
    if(pdbp == NULL)
        return;
    
    MemFree(pdbp->pelms);
    MemFree(pdbp);
    
    return;
}

PDBTablePtr FDBCreatePDBIndex(CharPtr filename)
{
    PDBTablePtr pdbp;
    Int4 length;
    Char buffer[1024];
    FILE *fd;
    Int4 gi, group, rank;
    Char div[32];
    
    pdbp = FDBPDBTableNew();
    length = sizeof(buffer);
    
    if((fd = FileOpen(filename, "r")) == NULL) {
        ErrPostEx(SEV_ERROR, 0,0, "Unable to open input index file");
        return NULL;
    }
    
    while(fgets(buffer, length, fd) != NULL) {
        if(*buffer == '#')
            continue;
        sscanf(buffer+16, "%d %d %d", &gi, &group, &rank);
        FDBAddNewPDB(pdbp, gi, group, rank);
    }
    
    HeapSort(pdbp->pelms, pdbp->count, sizeof(PDBElm), GICompare);
    
    return pdbp;
}

static void FDBAddNewLLUG(LLUGTablePtr llugp, Int4 gi, Int4 id, Int4 link)
{
    LLUGElmPtr pelm;
    
    /* Reallocate if necessary */
    
    if(llugp->allocated <= llugp->count) {
        llugp->allocated += LLUG_ALLOC_CHUNK + 1;
        llugp->pelms = (LLUGElmPtr) 
            Realloc (llugp->pelms, llugp->allocated * sizeof(LLUGElm));
    }
    
    pelm = &llugp->pelms[llugp->count];
    
    pelm->gi   = gi;
    pelm->id   = id;
    pelm->link = link;
    
    llugp->count++;
    
    return;
}
static LLUGTablePtr LLUGTableNew(void)
{
    LLUGTablePtr llugp;
    
    llugp = MemNew(sizeof(LLUGTable));
    llugp->allocated = LLUG_ALLOC_CHUNK;
    llugp->pelms = (LLUGElmPtr) MemNew (llugp->allocated * sizeof(LLUGElm));
    
    return llugp;
}

static void DestroyLLUGIndex(LLUGTablePtr llugp)
{
    if(llugp == NULL)
        return;
    
    MemFree(llugp->pelms);
    MemFree(llugp);
    
    return;
}

LLUGTablePtr CreateLLUGIndex(CharPtr filename)
{
    LLUGTablePtr llugp;
    Int4 length;
    Char buffer[1024];
    FILE *fd;
    Int4 gi, id, link = 0;
    Char div[32];
    
    llugp = LLUGTableNew();
    
    length = sizeof(buffer);
    
    if((fd = FileOpen(filename, "r")) == NULL) {
        ErrPostEx(SEV_ERROR, 0,0, "Unable to open input index file");
        return NULL;
    }
    
    while(fgets(buffer, length, fd) != NULL) {
        if(*buffer == '#')
            continue;
        /* sscanf(buffer, "%d %d %d", &gi, &id, &link); */
        sscanf(buffer, "%d %d", &gi, &id);
        FDBAddNewLLUG(llugp, gi, id, link);
    }
    
    HeapSort(llugp->pelms, llugp->count, sizeof(LLUGElm), GICompare);
    
    return llugp;
}

static HashTablePtr FDBHashTableNew(void)
{
    HashTablePtr htp;

    htp = MemNew(sizeof(HashTable));
    htp->allocated = HASH_ALLOC_CHUNK;
    htp->hep = (HashElmPtr) MemNew (htp->allocated * sizeof(HashElm));

    return htp;
}

static void FDBAddNewHash(HashTablePtr htp, Int4Ptr data, CharPtr div)
{
    HashElmPtr hem;

    /* Reallocate if necessary */
    
    if(htp->allocated <= htp->count) {
        htp->allocated += HASH_ALLOC_CHUNK + 1;
        htp->hep = (HashElmPtr) 
            Realloc (htp->hep, htp->allocated * sizeof(HashElm));
    }
    
    hem = &htp->hep[htp->count];
    
    hem->hash = data[5];
    hem->seq_num = data[0];
    hem->gi      = data[1];
    hem->tax_id  = data[2];
    hem->owner   = data[3];
    StringCpy(hem->div, div);
    hem->SequenceLen = data[4];
    hem->date = data[6];
    
    htp->count++;
    
    return;
}
static int HashCompare(VoidPtr i, VoidPtr j)
{
    if (*(Int4Ptr)i > *(Int4Ptr)j)
        return (1);
    if (*(Int4Ptr)i < *(Int4Ptr)j)
        return (-1);
    return (0);
}

static int DefListCompare(VoidPtr i, VoidPtr j)
{
    DefLinePtr dp, dp1;
    Int4 rank, rank1;

    dp =   *((DefLinePtr *) i);
    dp1 =  *((DefLinePtr *) j);

    if (dp->type > dp1->type)
        return (1);
    if (dp->type < dp1->type)
        return (-1);
    
    if(dp->pdb_rank || dp1->pdb_rank) { /* Both PDB - comaring by rank */
        if (dp->pdb_rank > dp1->pdb_rank)
            return (1);
        if (dp->pdb_rank < dp1->pdb_rank)
            return (-1);
    }
    
    /* If type the same we compare gis */
    
    if (dp->gi > dp1->gi)
        return (1);
    if (dp->gi < dp1->gi)
        return (-1);
    
    return (0);
}

static void FDBDestroyHashIndex(HashTablePtr htp)
{
    if(htp == NULL)
        return;
    
    MemFree(htp->hep);
    MemFree(htp);
    
    return;
}

static HashTablePtr FDBCreateHashIndex(CharPtr filename)
{
    HashTablePtr htp;
    Int4 length;
    Char buffer[1024];
    FILE *fd;
    Int4 data[7];
    Char div[32];
    
    htp = FDBHashTableNew();
    length = sizeof(buffer);

    if((fd = FileOpen(filename, "r")) == NULL) {
        ErrPostEx(SEV_ERROR, 0,0, "Unable to open input index file");
        return NULL;
    }

    while(fgets(buffer, length, fd) != NULL) {
        sscanf(buffer, "%d %d %d %d %s %d %d %d",
               &data[0], &data[1], &data[2], &data[3],
               div, &data[4], &data[5], &data[6]);
        
        FDBAddNewHash(htp, data, div);
        
    }
    HeapSort(htp->hep, htp->count, sizeof(HashElm), HashCompare);
    
    return htp;
}

static CharPtr ConcatDefline(CharPtr dline, CharPtr dline_tmp)
{
    CharPtr buf;
    
    if(dline == NULL) {
        buf = MemNew(StringLen(dline_tmp) + 2);
        sprintf(buf, "%s", dline_tmp);
    } else {
        buf = MemNew(StringLen(dline) + StringLen(dline_tmp) + 2);
        sprintf(buf, "%s%c%s", dline, '\1', dline_tmp);
    }
    
    MemFree(dline); 
    
    return buf;
}
/* --------- Functions shuffleing deflines ----------- */

static DefListPtr DefListNew(void)
{
    DefListPtr dlp;
    Int4 i;

    dlp = MemNew(sizeof(DefList));
    dlp->allocated = DEF_ALLOCATE_CHUNK;
    dlp->defs = MemNew(sizeof(DefLinePtr) * dlp->allocated);

    for(i = 0; i < dlp->allocated; i++) {
        dlp->defs[i] = MemNew(sizeof(DefLine));
    }
    
    return dlp;
}
static Boolean DefListRealloc(DefListPtr dlp)
{
    Int4 i, old_allocated;

    old_allocated = dlp->allocated;
    dlp->allocated += DEF_ALLOCATE_CHUNK;
    if((dlp->defs = Realloc(dlp->defs, 
                            sizeof(DefLinePtr) * dlp->allocated)) == NULL)
        return FALSE;
    
    for(i = old_allocated; i < dlp->allocated; i++) {
        if((dlp->defs[i] = MemNew(sizeof(DefLine))) == NULL)
            return FALSE;
    }
    
    return TRUE;
}
static void DefListFree(DefListPtr dlp)
{
    Int4 i;
    DefLinePtr dp;

    for(i = 0; i < dlp->allocated; i++) {
        dp = dlp->defs[i];
        MemFree(dp->line);
        MemFree(dp);
    }
    MemFree(dlp->defs);
    MemFree(dlp);

    return;
}
/* Note: pointers to Blast defline ASN.1 structures do not belong to
   this def-list - so the should be freed elsewhere */
static Boolean DefListAddLine(DefListPtr dlp, BlastDefLinePtr bdp, 
                              CharPtr line, Int4 type, Int4 gi, Int4 pdb_rank)
{
    DefLinePtr dp;
    
    if(dlp->count >= dlp->allocated) {
        if(!DefListRealloc(dlp))
            return FALSE;
    }
    
    dp = dlp->defs[dlp->count];
    
    dp->line = StringSave(line);
    dp->bdp  = bdp;
    dp->type = type;
    dp->gi = gi;
    dp->pdb_rank = pdb_rank;
    
    dlp->count++;

    return TRUE;
}

static CharPtr FinalDefLineOut(DefListPtr dlp, SeqIdPtr PNTR seqid)
{
    Int4 i, j;
    DefLinePtr dp;
    CharPtr chptr, dline = NULL;
    Char buffer[512];
    
    for(i = 0; i < dlp->count; i++) {
        dp = dlp->defs[i];

        if(i == 0) {
            StringNCpy(buffer, dp->line, sizeof(buffer) -1);
            
            if((chptr = StringChr(buffer, ' ')) != NULL)
                *chptr = NULLB;
            
            *seqid = SeqIdParse(buffer);
            dline = ConcatDefline(NULL, chptr+1);
        } else {
            dline = ConcatDefline(dline, dp->line);
        }
    }
    return dline;
}

static BlastDefLinePtr FinalBdfpOut(DefListPtr dlp)
{
    Int4 i;
    DefLinePtr dp;
    BlastDefLinePtr bdp_head, bdp_tail;
    
    bdp_head = bdp_tail = NULL;
    
    for(i = 0; i < dlp->count; i++) {
        
        dp = dlp->defs[i];
        
        if(i == 0) {
            bdp_head = bdp_tail = dp->bdp;
        } else {
            bdp_tail->next = dp->bdp;
            bdp_tail = dp->bdp;
        }
    }

    return bdp_head;
}

static Int4 GetPDBRankByGi(PDBTablePtr pdbp, Int4 gi)
{

    Int4 m, b, e;
    
    if(pdbp == NULL)
        return 0;
    
    b = 0;
    e = pdbp->count;
    
    while (b < e - 1) {
	m = (b + e) / 2;
	if ((pdbp->pelms[m].gi) > gi)
	    e = m;
	else
	    b = m;
    }
    
    if(pdbp->pelms[b].gi == gi)
        return pdbp->pelms[b].rank;
    else                        /* Gi was not found */
        return 0;
}

static Int4 GetMinimalType(PDBTablePtr pdbp, SeqIdPtr sip, CharPtr defline, 
                           Int4Ptr gip, Int4Ptr pdb_rank)
{
    SeqIdPtr sip_tmp;
    Int4 order, order1 = INT2_MAX;
    Char buffer[512];
    CharPtr chptr;
    Boolean is_pdb = FALSE, sip_allocated = FALSE;

    if(pdbp == NULL || gip == NULL)
        return -1;
    
    *gip = -1;                  /* Default if NOT found */

    if(sip == NULL) {
        StringNCpy(buffer, defline, sizeof(buffer) - 1);
        
        if((chptr = StringChr(buffer, ' ')) != NULL)
            *chptr = NULLB;
        else
            return -1;
        
        sip = SeqIdParse(buffer);
        sip_allocated = TRUE;
    }

    *pdb_rank = 0;
    is_pdb = FALSE;
    *gip = 0;

    for(sip_tmp = sip; sip_tmp != NULL; sip_tmp = sip_tmp->next) {
        if((order = GetOrderBySeqId(sip_tmp->choice, TRUE)) < 0)
            return -1;
        order1 = MIN(order, order1);
        
        /* Extracting gi */
        if(sip_tmp->choice == SEQID_GI)
            *gip = (Int4) sip->data.intvalue;
        
        if(sip_tmp->choice == SEQID_PDB)
            is_pdb = TRUE;

    }
    
    if(is_pdb && *gip != 0)
       *pdb_rank = GetPDBRankByGi(pdbp, *gip);
    
    if(sip_allocated)
        SeqIdSetFree(sip);
    
    return order1;
}

Int4 Check_UGLL_ID(LLUGTablePtr llugp, Int4 gi)
{    
    Int4 m, b, e;
    
    if(llugp == NULL)
        return 0;
    
    b = 0;
    e = llugp->count;
    
    while (b < e - 1) {
	m = (b + e) / 2;
	if ((llugp->pelms[m].gi) > gi)
	    e = m;
	else
	    b = m;
    }
    
    if(llugp->pelms[b].gi == gi)
        return llugp->pelms[b].id;
    else                        /* Gi was not found */
        return 0;
}
    
Boolean FDRAddLinksMembership(BlastDefLinePtr bdp, HashElmPtr hep, 
                              ReadDBFILEPtr rdfp,
                              LLUGTablePtr UGp, LLUGTablePtr LLp)
{
    Int4 mem_int = 0, link_int = 0;
    ValNodePtr vnp, vnp_last;
    Int4 ug_id, ll_id;

    if(bdp == NULL)
        return FALSE;
    
    bdp->taxid = hep->tax_id;

    /* Setting bits for EST_MOUSE and EST_HUMAN databases */
    switch (bdp->taxid) {
    case 9606:
        mem_int += EST_HUMAN_BIT;
        break;
    case 10090:
    case 10091:
    case 10092:
    case 35531:
    case 80274:
    case 57486:
        mem_int += EST_MOUSE_BIT;
        break;
    default:
        break;
    }

    switch(hep->owner) {
    case 6:
        mem_int += SWISSPROT_BIT;
        break;
    case 10:
        mem_int += PDB_BIT;
        break;
    case 20:
        mem_int += REFSEQ_BIT;
        break;
    case 28:
        mem_int += CONTIG_BIT;
        break;
    default:
        break;
    }

    if(mem_int != 0) {
        vnp = ValNodeNew(NULL);
        vnp->data.intvalue = mem_int;
        bdp->memberships = vnp;
    }

    ug_id = Check_UGLL_ID(UGp, hep->gi);
    
    if(ug_id != 0) {
        link_int += IS_UG_ELEMENT;   /* First bit for UNIGENE */
    }
    
    ll_id = Check_UGLL_ID(LLp, hep->gi);
    
    if(ll_id != 0) {
        link_int += IS_LL_ELEMENT;   /* Second bit for Locus link */
    }

    if(link_int > 0) {
        vnp = ValNodeNew(NULL);
        vnp->data.intvalue = link_int;
        bdp->links = vnp;
        vnp_last =  bdp->links;
    }

#ifdef ADD_UGLL_IDS
    if(link_int & IS_UG_ELEMENT) { /* UG ? */
        vnp = ValNodeNew(NULL);
        vnp->data.intvalue = ug_id;
        vnp_last->next = vnp;
        vnp_last = vnp;
    }
    
    if(link_int & IS_LL_ELEMENT) { /* LL ? */
        vnp = ValNodeNew(NULL);
        vnp->data.intvalue = ll_id;
        vnp_last->next = vnp;
        vnp_last = vnp;
    }
#endif

    return TRUE;
}

/* --------------------------------------------------- */

static Int4 NewUniqueFASTA(ReadDBFILEPtr rdfp, HashTablePtr htp, Int4 count,
                           BlastDefLinePtr * bdp_out,
                           ValNodePtr PNTR seqid, CharPtr PNTR defline, 
                           BioseqPtr PNTR bsp, FILE *fd_info, Int4 seq_num,
                           PDBTablePtr pdbp,  
                           LLUGTablePtr UGp, LLUGTablePtr LLp)
{
    Int4 i, hash_val, length, len_seq;
    Int4 first, next_count = 0, type;
    UcharPtr sequence, buffer;
    CharPtr dline = NULL;
    DefListPtr dlp;
    Int4 gi = 0, pdb_rank = 0;
    BlastDefLinePtr bdp;
    
    dlp = DefListNew();
    
    for(i = count, hash_val = htp->hep[count].hash, first = TRUE; 
        htp->hep[i].hash == hash_val; i++) {
        
        if(htp->hep[i].seq_num == -1)
            continue;
        
        length = readdb_get_sequence(rdfp, htp->hep[i].seq_num, &buffer);
        
        if(length <= 0)
            return -1;
        
        bdp = NULL;

        if(first) {
            sequence = buffer;
            len_seq = length;

            *bsp = readdb_get_bioseq(rdfp, htp->hep[i].seq_num);            
            if(*bsp == NULL)
                return -1;

            if(rdfp->formatdb_ver == FORMATDB_VER_TEXT) {            

                readdb_get_defline(rdfp, htp->hep[i].seq_num, &dline);
                
                if(dline == NULL)
                    return -1;

            } else {            
                bdp = FDReadDeflineAsn(rdfp, htp->hep[i].seq_num);
                FDRAddLinksMembership(bdp, &htp->hep[i], rdfp, UGp, LLp);
            }

            
            type = GetMinimalType(pdbp, bdp == NULL ? NULL : bdp->seqid, 
                                  dline, &gi, &pdb_rank);
            DefListAddLine(dlp, bdp, dline, type, gi, pdb_rank);
            
            if(rdfp->formatdb_ver == FORMATDB_VER_TEXT) {            
                MemFree(dline); 
            }
            
            first = FALSE;
        } else {
            /* Comparing sequences - if they are different despite hash: */
            if(length != len_seq || memcmp(buffer, sequence, length)) {
                if(next_count == 0) next_count = i;
                continue;
            }
            
            if(rdfp->formatdb_ver == FORMATDB_VER_TEXT) { 
                readdb_get_defline(rdfp, htp->hep[i].seq_num, &dline);
                if(dline == NULL)
                    return -1; 
            } else {            
                bdp = FDReadDeflineAsn(rdfp, htp->hep[i].seq_num);
                FDRAddLinksMembership(bdp, &htp->hep[i], rdfp, UGp, LLp);
            }
            
            type = GetMinimalType(pdbp, bdp == NULL ? NULL : bdp->seqid, 
                                  dline, &gi, &pdb_rank);
            DefListAddLine(dlp, bdp, dline, type, gi, pdb_rank);

            if(rdfp->formatdb_ver == FORMATDB_VER_TEXT) { 
                MemFree(dline); 
            }

            htp->hep[i].seq_num = -1; /* Label do not pass second time */
        }
        
        fprintf(fd_info, "%d %d %d %d %s %d %d %d\n", 
                seq_num, htp->hep[i].gi, htp->hep[i].tax_id, 
                htp->hep[i].owner, htp->hep[i].div, 
                htp->hep[i].SequenceLen, htp->hep[i].hash, htp->hep[i].date);
        
        fflush(fd_info);
        
        /* DumpInfoFile(&htp->hep[i], fd_info, seq_num); */
    }

    HeapSort(dlp->defs, dlp->count, sizeof(DefLinePtr), DefListCompare);

    if(rdfp->formatdb_ver == FORMATDB_VER_TEXT) {
        *defline = FinalDefLineOut(dlp, seqid);
    } else {
        *bdp_out = FinalBdfpOut(dlp);
    }
    
    DefListFree(dlp);
    
    if(next_count == 0) next_count = i;
    
    return next_count;
}

/* Functions used in filtering by gi number */
static GiListPtr GiListNew(void)
{
    GiListPtr glp;

    glp = MemNew(sizeof(GiList));
    glp->allocated = GI_ALLOC_CHUNK;
    glp->seq_num = MemNew(sizeof(Int4) * glp->allocated);
    glp->count = 0;
    
    return glp;
}
static void GiListFree(GiListPtr glp)
{
    if(glp == NULL)
        return;
    
    MemFree(glp->seq_num);
    MemFree(glp);
    return;
}

static Boolean ReadGiList(ReadDBFILEPtr rdfp, GiListPtr glp, CharPtr filename)
{
    FILE *fd;
    Int4 gi, retvalue, seqnum;

    if((fd = FileOpen(filename, "r")) == NULL)
        return FALSE;

    while((retvalue = fscanf(fd, "%d", &gi)) != EOF) {
        if(retvalue == 0) continue;
        
        if(glp->count >= glp->allocated) {
            glp->allocated += GI_ALLOC_CHUNK;
            glp->seq_num = Realloc(glp->seq_num, 
                                   sizeof(Int4) * glp->allocated);
        }

        seqnum = readdb_gi2seq(rdfp, gi, 0);
            
        if(seqnum < 0) {
            ErrPostEx(SEV_WARNING, 0,0, "Gi %d is not found", gi);
            continue;
        }

        glp->seq_num[glp->count] = seqnum;
        glp->count++;
    }

    FileClose(fd);

    return TRUE;
}
/* Here we will check, that data[2] == tax_id, 
   data[3] == owner, div=div */
static Boolean CheckSRCondition(SR_InfoPtr srip, Int4Ptr data, CharPtr div)
{
    Int4    i;
    CharPtr chptr;
    Boolean  cond_ok = FALSE;

    /* checking tax_id */

    if(srip->taxid[0] > 0) { /* At least one element exists */
        cond_ok = FALSE;
        for(i = 0; srip->taxid[i] > 0 && i < MAX_SR_ELEMENTS; i++) {
            if(data[2] == srip->taxid[i]) {
                cond_ok = TRUE;
                break;
            }
        }
        if(cond_ok == FALSE)
            return FALSE;
    }

    /* checking owner */

    if(srip->owner[0] > 0) { /* At least one element exists */
        cond_ok = FALSE;
        for(i = 0; srip->owner[i] > 0 && i < MAX_SR_ELEMENTS; i++) {
            if(data[3] == srip->owner[i]) {
                cond_ok = TRUE;
                break;
            }
        }
        if(cond_ok == FALSE)
            return FALSE;
    }

    /* checking division */

    if(*srip->div[0] != NULLB) { /* At least one element exists */
        cond_ok = FALSE;
        for(i = 0; *srip->div[i] != NULLB && i < MAX_SR_ELEMENTS; i++) {
            if(!StringCmp(div, srip->div[i])) {
                cond_ok = TRUE;
                break;
            }
        }
        if(cond_ok == FALSE)
            return FALSE;
    }
    
    return TRUE;
}
/* Here we will check, that data[2] == tax_id, 
   data[3] == owner, div=div */
static Boolean CheckSRConditionReverse(SR_InfoPtr srip, 
                                       Int4Ptr data, CharPtr div)
{
    Int4    i;
    CharPtr chptr;
    Boolean  cond_ok = TRUE;
    
    /* checking tax_id */
    
    if(srip->taxid[0] > 0) { /* At least one element exists */
        cond_ok = TRUE;
        for(i = 0; srip->taxid[i] > 0 && i < MAX_SR_ELEMENTS; i++) {
            if(data[2] == srip->taxid[i]) {
                cond_ok = FALSE;
                break;
            }
        }
        if(cond_ok == FALSE)
            return FALSE;
    }

    /* checking owner */

    if(srip->owner[0] > 0) { /* At least one element exists */
        cond_ok = TRUE;
        for(i = 0; srip->owner[i] > 0 && i < MAX_SR_ELEMENTS; i++) {
            if(data[3] == srip->owner[i]) {
                cond_ok = FALSE;
                break;
            }
        }
        if(cond_ok == FALSE)
            return FALSE;
    }

    /* checking division */

    if(*srip->div[0] != NULLB) { /* At least one element exists */
        cond_ok = TRUE;
        for(i = 0; *srip->div[i] != NULLB && i < MAX_SR_ELEMENTS; i++) {
            if(!StringCmp(div, srip->div[i])) {
                cond_ok = FALSE;
                break;
            }
        }
        if(cond_ok == FALSE)
            return FALSE;
    }
    
    return TRUE;
}

static Boolean FDGetGiListByQuery(ReadDBFILEPtr rdfp, SR_InfoPtr srip, 
                                  GiListPtr glp, CharPtr filename, 
                                  Boolean reverse)
{ 
    FILE *fd;
    Int4 length, data[7];
    Char div[32];
    Char buffer[1024];
    Boolean found;

    if((fd = FileOpen(filename, "r")) == NULL) {
        ErrPostEx(SEV_ERROR, 0,0, "Unable to open input info file");
        return FALSE;
    }
    
    length = sizeof(buffer);
    while(fgets(buffer, length, fd) != NULL) {
        sscanf(buffer, "%d %d %d %d %s %d %d %d",
               &data[0], &data[1], &data[2], &data[3],
               div, &data[4], &data[5], &data[6]);
        
        /* If line agree with condition seq_num added to the list */

        if(reverse)
            found = CheckSRConditionReverse(srip, data, div);
        else
            found = CheckSRCondition(srip, data, div);

        if(found) {
            if(glp->count >= glp->allocated) {
                glp->allocated += GI_ALLOC_CHUNK;
                glp->seq_num = Realloc(glp->seq_num, 
                                       sizeof(Int4) * glp->allocated);
            }
            glp->seq_num[glp->count] = data[0];
            glp->count++;
        }
    }
    
    FileClose(fd);
    
    return TRUE;
}

Int2 Main(void)
{
    ReadDBFILEPtr rdfp;
    Int4 i, count, gi, seqnum;
    FDB_optionsPtr options;
    FormatDBPtr	fdbp;
    BioseqPtr bsp;

    HashTablePtr htp = NULL;
    GiListPtr glp = NULL;

    Char buffer[128];
    SeqIdPtr sip = NULL, sip_tmp;
    Int4 next_number = 0;
    Char tmpbuf[128];
    CharPtr defline = NULL;
    FILE *fd_info;
    SR_InfoPtr srip;
    Boolean is_prot;
    PDBTablePtr pdbp;
    BlastDefLinePtr bdp = NULL;
    LLUGTablePtr LLp, UGp;

    /* ---------------------------------------------- */
    /* ----- Initializing formatdb structures ------- */
    /* ---------------------------------------------- */
    
    if((options = FDB_CreateCLOptions()) == NULL)
        return 1;
        
    /* ---------------------------------------------- */
    /* ------ Initializing readdb structures -------- */
    /* ---------------------------------------------- */

    if((rdfp = readdb_new (FLT_Input, IS_Protein)) == NULL) {
        ErrPostEx(SEV_ERROR, 0, 0, 
                  "Failure to intialise database %s", FLT_Input);
        return 1;
    }

    /* New database will have authomatically the same version as
       old one */
    options->version = rdfp->formatdb_ver;
    
    is_prot = (Boolean) (rdfp->parameters & READDB_IS_PROT);
    
    count = readdb_get_num_entries (rdfp);
    
    /* ---------------------------------------------- */
    /* ------ Creating index for hash search -------- */
    /* ---------------------------------------------- */

    if(flt_args[2].strvalue != NULL) {/* This is list of gis */
        glp = GiListNew();
        if(!ReadGiList(rdfp, glp, flt_args[2].strvalue))
            return -1;
    } else if(flt_args[5].strvalue != NULL) {/* This is list of gis */

        if((srip = SRReadSRInfo(flt_args[5].strvalue)) == NULL)
            return -1;

        glp = GiListNew();
        sprintf(buffer, "%s.%cdi", FLT_Input, is_prot? 'p' : 'n');
        if(!FDGetGiListByQuery(rdfp, srip, glp, buffer, ReverseQuery)) 
            return -1;
        
        /* If database name set - replacing default one */
        if(*srip->dbname != NULLB) {
            MemFree(options->db_file);
            options->db_file = StringSave(srip->dbname);
        }
        
        MemFree(srip);

    } else { /* Hash filtering */
        sprintf(buffer, "%s.%cdi", FLT_Input, is_prot? 'p' : 'n');
        
        if((htp = FDBCreateHashIndex(buffer)) == NULL) {
            ErrPostEx(SEV_ERROR, 0,0, "Failure to create hash index");
            return -1;
        }

        if(flt_args[8].strvalue != NULL) 
            pdbp = FDBCreatePDBIndex(flt_args[8].strvalue);
        else
            pdbp = NULL;

        /* Locus Link and UniGene information processing */

        UGp = NULL,  LLp = NULL;
        
        if(rdfp->formatdb_ver > FORMATDB_VER_TEXT) {
            
            if(flt_args[9].strvalue != NULL) {
                UGp = CreateLLUGIndex(flt_args[9].strvalue);
            }
            
            if(flt_args[10].strvalue != NULL) {
                LLp = CreateLLUGIndex(flt_args[10].strvalue);
            }
        }

        sprintf(buffer, "%s.%cdi", options->db_file, is_prot? 'p' : 'n');
        fd_info = FileOpen(buffer, "w");
    }

    /* ---------------------------------------------- */
    /* ----- Initializing formatdb structure  ------- */
    /* ---------------------------------------------- */

    if ((fdbp = FormatDBInit(options)) == NULL)
        return -1;

#ifdef TAX_CS_LOOKUP
    /* These functions will create taxonomy lookup database */
    if(rdfp->formatdb_ver > FORMATDB_VER_TEXT && options->parse_mode) {
        options->tax_lookup = RDTaxLookupInit();
        options->tax_callback = FDBTaxCallback;
    }
#endif    
    
    /* ---------------------------------------------- */
    /* ---------------- Main loop ------------------- */
    /* ---------------------------------------------- */
    
#ifndef TEST_RDB
    if(htp != NULL) { /* filtering by hash value */
        next_number = 0;
        do {
            sip = NULL;
            next_number =  NewUniqueFASTA(rdfp, htp, next_number, &bdp, 
                                          &sip, &defline, &bsp, 
                                          fd_info, fdbp->num_of_seqs, 
                                          pdbp, UGp, LLp);
            
            if(next_number < 0) {
                ErrPostEx(SEV_ERROR, 0, 0, "Failure to get sequence");
                return 1;
            }
            
            SeqIdWrite(sip, tmpbuf, 
                       PRINTID_FASTA_LONG, sizeof(tmpbuf));
            
            FDBAddSequence(fdbp, bdp, bsp->seq_data_type, 
                           &bsp->seq_data, bsp->length, 
                           tmpbuf, defline, 0, 0, 0, 0, 0);

            if(sip)
                sip = SeqIdSetFree(sip);

            bsp = BioseqFree(bsp);        
            
            if(defline)
                defline = MemFree(defline);

            if(bdp)
                BlastDefLineSetFree(bdp);
            
        } while (next_number < count); 
        
        FileClose(fd_info);

        FDBDestroyHashIndex(htp);
        FDBDestroyPDBIndex(pdbp);
        DestroyLLUGIndex(UGp);
        DestroyLLUGIndex(LLp);

    }  else { /* list of gis */
        
        for(i = 0; i < glp->count; i++) {

            if((seqnum = glp->seq_num[i]) == -1)
                continue;
            
            bsp = readdb_get_bioseq(rdfp, seqnum);
            readdb_get_descriptor(rdfp, seqnum, &sip, &defline);
            
            SeqIdWrite(sip, tmpbuf, 
                       PRINTID_FASTA_LONG, sizeof(tmpbuf));
            
            FDBAddSequence(fdbp, NULL, bsp->seq_data_type, 
                           &bsp->seq_data, bsp->length, 
                           tmpbuf, defline, 0, 0, 0, 0, 0);
            
            BioseqFree(bsp);
            MemFree(defline);
            
            do {
                sip_tmp = sip->next;
                SeqIdFree(sip);
                sip = sip_tmp;
            } while(sip != NULL);
        }

        GiListFree(glp);
    }
    
#else
    for(i = 0; i < count; i++) {
        Int4 length;
        Uint1Ptr buffer;
        Char tmpbuf[128];
        CharPtr defline;
        BioseqPtr bsp;

        bsp = readdb_get_bioseq(rdfp, i);
        readdb_get_descriptor(rdfp, i, &sip, &defline);
        
        SeqIdWrite(sip, tmpbuf, PRINTID_FASTA_LONG, sizeof(tmpbuf));
        
        FDBAddSequence(fdbp, NULL, bsp->seq_data_type, 
                       &bsp->seq_data, bsp->length, 
                       tmpbuf, defline, 0, 0, 0, 0, 0);
    
        BioseqFree(bsp);
        MemFree(defline);
        
        do {
            sip_tmp = sip->next;
            SeqIdFree(sip);
            sip = sip_tmp;
        } while(sip != NULL);
    }

#endif

    if(FormatDBClose(fdbp))
        return 3;

#ifdef TAX_CS_LOOKUP
    if(options->tax_lookup != NULL) {
        RDTaxLookupClose(options->tax_lookup);
    }
#endif
    
    readdb_destruct(rdfp);
    FDB_optionsFree(options);
    
    return 0;
}
 
