/* $Id: fdfilter.c,v 6.9 1999/09/20 18:35:03 shavirin Exp $ */
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

   $Revision: 6.9 $

   File Description:  Filter FASTA databases for identical sequences

   $Log: fdfilter.c,v $
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
    CharPtr line;
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

#define NUMARG 8

Args flt_args[NUMARG] = {
    { "Title for output database file", 
      NULL, NULL, NULL, TRUE, 't', ARG_STRING, 0.0, 0, NULL},
    {"Input file for the filtering (this parameter must be set)",
     NULL, NULL,NULL,FALSE,'i',ARG_FILE_IN, 0.0,0,NULL},
    {"Input file with a list of gis",
     NULL, NULL,NULL,TRUE,'g',ARG_FILE_IN, 0.0,0,NULL},
    {"Logfile name:",
     "fdfilter.log", NULL,NULL,TRUE,'l',ARG_FILE_OUT, 0.0,0,NULL},
    {"Create sparse indexes in the filtered database",
     "F", NULL,NULL,TRUE,'s',ARG_BOOLEAN, 0.0,0,NULL},
    {"Query string for creating database subset",
     NULL, NULL,NULL,TRUE,'q',ARG_STRING, 0.0,0,NULL},
    {"Input database is proten",
     "T", NULL,NULL,TRUE,'p',ARG_BOOLEAN, 0.0,0,NULL},
    {"Reverse query request to negative",
     "F", NULL,NULL,TRUE,'r',ARG_BOOLEAN, 0.0,0,NULL},
};

#define FLT_Input    flt_args[1].strvalue
#define IS_Protein   flt_args[6].intvalue
#define ReverseQuery flt_args[7].intvalue

void FDB_optionsFree(FDB_optionsPtr options)
{
    if(options == NULL)
        return;
    
    MemFree(options->db_title);
    MemFree(options->db_file);
    MemFree(options->LogFileName);
    MemFree(options);
    
    return;
}

SeqIdPtr MySeqIdFree(SeqIdPtr sip)
{
    SeqIdPtr sip_tmp;
    do {
        sip_tmp = sip->next;
        SeqIdFree(sip);
        sip = sip_tmp;
    } while(sip != NULL);
    
    return NULL;
}

Boolean SRReadCharData(CharPtr buffer, CharPtr PNTR div_in)
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

Boolean SRReadIntData(CharPtr buffer, Int4Ptr id)
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

SR_InfoPtr SRReadSRInfo(CharPtr buffer)
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

FDB_optionsPtr FDB_CreateCLOptions(void)
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

    return options;
}
#define HASH_ALLOC_CHUNK 1024
#define GI_ALLOC_CHUNK 1024

HashTablePtr FDBHashTableNew(void)
{
    HashTablePtr htp;

    htp = MemNew(sizeof(HashTable));
    htp->allocated = HASH_ALLOC_CHUNK;
    htp->hep = (HashElmPtr) MemNew (htp->allocated * sizeof(HashElm));

    return htp;
}

void FDBAddNewHash(HashTablePtr htp, Int4Ptr data, CharPtr div)
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
    
    dp =   *((DefLinePtr *) i);
    dp1 =  *((DefLinePtr *) j);

    if (dp->type > dp1->type)
        return (1);
    if (dp->type < dp1->type)
        return (-1);
    return (0);
}

void FDBDestroyHashIndex(HashTablePtr htp)
{
    if(htp == NULL)
        return;
    
    MemFree(htp->hep);
    MemFree(htp);
    
    return;
}

HashTablePtr FDBCreateHashIndex(CharPtr filename)
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
/* --------- Functions shuffeling deflines ----------- */

DefListPtr DefListNew(void)
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
Boolean DefListRealloc(DefListPtr dlp)
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
void DefListFree(DefListPtr dlp)
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

Boolean DefListAddLine(DefListPtr dlp, CharPtr line, Int4 type)
{
    DefLinePtr dp;

    if(dlp->count >= dlp->allocated) {
        if(!DefListRealloc(dlp))
            return FALSE;
    }
    
    dp = dlp->defs[dlp->count];
    
    dp->line = StringSave(line);
    dp->type = type;
    
    dlp->count++;
    return TRUE;
}

CharPtr FinalDefLineOut(DefListPtr dlp, SeqIdPtr PNTR seqid)
{
    Int4 i, j;
    DefLinePtr dp;
    CharPtr chptr, dline = NULL;
    Char buffer[512];
    
    HeapSort(dlp->defs, dlp->count, sizeof(DefLinePtr), DefListCompare);

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


Int4 GetMinimalType(CharPtr defline)
{
    SeqIdPtr sip, sip_tmp;
    Int4 order, order1 = INT2_MAX;
    Char buffer[512];
    CharPtr chptr;

    StringNCpy(buffer, defline, sizeof(buffer) - 1);

    if((chptr = StringChr(buffer, ' ')) != NULL)
        *chptr = NULLB;
    else
        return -1;

    sip = SeqIdParse(buffer);

    for(sip_tmp = sip; sip_tmp != NULL; sip_tmp = sip_tmp->next) {
        if((order = GetOrderBySeqId(sip_tmp->choice)) < 0)
            return -1;
        order1 = MIN(order, order1);
    }

    MySeqIdFree(sip);
    return order1;
}

/* --------------------------------------------------- */

Int4 NewUniqueFASTA(ReadDBFILEPtr rdfp, HashTablePtr htp, Int4 count, 
                    ValNodePtr PNTR seqid, CharPtr PNTR defline, 
                    BioseqPtr PNTR bsp, FILE *fd_info, Int4 seq_num)
{
    Int4 i, hash_val, length, len_seq;
    Int4 first, next_count = 0, type;
    UcharPtr sequence, buffer;
    CharPtr dline;
    DefListPtr dlp;

    dlp = DefListNew();
    
    for(i = count, hash_val = htp->hep[count].hash, first = TRUE; 
        htp->hep[i].hash == hash_val; i++) {
        
        if(htp->hep[i].seq_num == -1)
            continue;
        
        length = readdb_get_sequence(rdfp, htp->hep[i].seq_num, &buffer);
        
        if(length <= 0)
            return -1;
        
        if(first) {
            sequence = buffer;
            len_seq = length;
            *bsp = readdb_get_bioseq(rdfp, htp->hep[i].seq_num);
            
            if(*bsp == NULL)
                return -1;
            
            /* readdb_get_descriptor(rdfp, htp->hep[i].seq_num, 
               seqid, &dline); */
            

            readdb_get_defline(rdfp, htp->hep[i].seq_num, &dline);
            type = GetMinimalType(dline);
            
            if(dline == NULL)
                return -1;
            
            DefListAddLine(dlp, dline, type);
            MemFree(dline);          
            
            first = FALSE;
        } else {
            /* Comparing sequences - if they are different despite hash: */
            if(length != len_seq || memcmp(buffer, sequence, length)) {
                if(next_count == 0) next_count = i;
                continue;
            }
            readdb_get_defline(rdfp, htp->hep[i].seq_num, &dline);
            
            if(dline == NULL)
                return -1;

            type = GetMinimalType(dline);
            DefListAddLine(dlp, dline, type);
            MemFree(dline);
                
            htp->hep[i].seq_num = -1; /* Label do not pass second time */
        }

        fprintf(fd_info, "%d %d %d %d %s %d %d %d\n", 
                seq_num, htp->hep[i].gi, htp->hep[i].tax_id, 
                htp->hep[i].owner, htp->hep[i].div, 
                htp->hep[i].SequenceLen, htp->hep[i].hash, htp->hep[i].date);
        
        fflush(fd_info);
        
        /* DumpInfoFile(&htp->hep[i], fd_info, seq_num); */
    }

    *defline = FinalDefLineOut(dlp, seqid);

    DefListFree(dlp);
    
    if(next_count == 0) next_count = i;
    
    return next_count;
}

/* Functions used in filtering by gi number */
GiListPtr GiListNew(void)
{
    GiListPtr glp;

    glp = MemNew(sizeof(GiList));
    glp->allocated = GI_ALLOC_CHUNK;
    glp->seq_num = MemNew(sizeof(Int4) * glp->allocated);
    glp->count = 0;
    
    return glp;
}
void GiListFree(GiListPtr glp)
{
    if(glp == NULL)
        return;
    
    MemFree(glp->seq_num);
    MemFree(glp);
    return;
}

Boolean ReadGiList(ReadDBFILEPtr rdfp, GiListPtr glp, CharPtr filename)
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

        seqnum = readdb_gi2seq(rdfp, gi);
            
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
Boolean CheckSRCondition(SR_InfoPtr srip, Int4Ptr data, CharPtr div)
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
Boolean CheckSRConditionReverse(SR_InfoPtr srip, Int4Ptr data, CharPtr div)
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

Boolean FDGetGiListByQuery(ReadDBFILEPtr rdfp, SR_InfoPtr srip, 
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
    SeqIdPtr sip, sip_tmp;
    Int4 next_number = 0;
    Char tmpbuf[128];
    CharPtr defline;
    FILE *fd_info;
    SR_InfoPtr srip;
   
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
        sprintf(buffer, "%s.%cdi", FLT_Input, rdfp->is_prot? 'p' : 'n');
        if(!FDGetGiListByQuery(rdfp, srip, glp, buffer, ReverseQuery)) 
            return -1;
        
        /* If database name set - replacing default one */
        if(*srip->dbname != NULLB) {
            MemFree(options->db_file);
            options->db_file = StringSave(srip->dbname);
        }
        
        MemFree(srip);

    } else { /* Hash filtering */
        sprintf(buffer, "%s.%cdi", FLT_Input, rdfp->is_prot? 'p' : 'n');
        
        if((htp = FDBCreateHashIndex(buffer)) == NULL) {
            ErrPostEx(SEV_ERROR, 0,0, "Failure to create hash index");
            return -1;
        }
        
        sprintf(buffer, "%s.%cdi", options->db_file, rdfp->is_prot? 'p' : 'n');
        fd_info = FileOpen(buffer, "w");
    }

    /* ---------------------------------------------- */
    /* ----- Initializing formatdb structure  ------- */
    /* ---------------------------------------------- */

    if ((fdbp = FormatDBInit(options)) == NULL)
        return -1;
    
    /* ---------------------------------------------- */
    /* ---------------- Main loop ------------------- */
    /* ---------------------------------------------- */
    
#ifndef TEST_RDB
    if(htp != NULL) { /* filtering by hash value */
        next_number = 0;
        do {
            next_number =  NewUniqueFASTA(rdfp, htp, next_number, &sip, 
                                          &defline, &bsp, fd_info, 
                                          fdbp->num_of_seqs);

            if(next_number < 0) {
                ErrPostEx(SEV_ERROR, 0, 0, "Failure to get sequence");
                return 1;
            }
            
            SeqIdWrite(sip, tmpbuf, 
                       PRINTID_FASTA_LONG, sizeof(tmpbuf));
            
            FDBAddSequence (fdbp, 0, tmpbuf, defline, 0, 0, 0, 
                            bsp->seq_data_type, &bsp->seq_data, 
                            bsp->length, 0);
            sip = MySeqIdFree(sip);
            bsp = BioseqFree(bsp);        
            defline = MemFree(defline);
        } while (next_number < count); 

        FileClose(fd_info);
        FDBDestroyHashIndex(htp);

    }  else { /* list of gis */
        
        for(i = 0; i < glp->count; i++) {

            if((seqnum = glp->seq_num[i]) == -1)
                continue;
            
            bsp = readdb_get_bioseq(rdfp, seqnum);
            readdb_get_descriptor(rdfp, seqnum, &sip, &defline);
            
            SeqIdWrite(sip, tmpbuf, 
                       PRINTID_FASTA_LONG, sizeof(tmpbuf));
            
            FDBAddSequence (fdbp, 0, tmpbuf, defline, 0, 0, 0, 
                            bsp->seq_data_type, &bsp->seq_data, 
                            bsp->length, 0);
            
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
        
        FDBAddSequence (fdbp, 0, tmpbuf, defline, 0, 0, 0, 
                        bsp->seq_data_type,
                        &bsp->seq_data, bsp->length, 0);
    
        /*    SeqIdWrite(bsp->id, tmpbuf, 
              PRINTID_FASTA_LONG, sizeof(tmpbuf));
              FDBAddSequence (fdbp, 0, tmpbuf, BioseqGetTitle(bsp), 0, 0, 0, 
              bsp->seq_data_type,
              &bsp->seq_data, bsp->length, 0); */

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
    
    readdb_destruct(rdfp);
    FDB_optionsFree(options);
    
    return 0;
}
 