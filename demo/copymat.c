/*
* ===========================================================================
*
*                            PUBLIC DOMAIN NOTICE
*               National Center for Biotechnology Information
*
*  This software/database is a "United States Government Work" under the
*  terms of the United States Copyright Act.  It was written as part of
*  the author's official duties as a United States Government employee and
*  thus cannot be copyrighted.  This software/database is freely available
*  to the public for use. The National Library of Medicine and the U.S.
*  Government have not placed any restriction on its use or reproduction.
*
*  Although all reasonable efforts have been taken to ensure the accuracy
*  and reliability of the software and data, the NLM and the U.S.
*  Government do not and cannot warrant the performance or results that
*  may be obtained by using this software or data. The NLM and the U.S.
*  Government disclaim all warranties, express or implied, including
*  warranties of performance, merchantability or fitness for any particular
*  purpose.
*
*  Please cite the author in any work or product based on this material.
*
* ===========================================================================
*/

/*****************************************************************************

File name: copymat.c

Authors: Alejandro Schaffer, Sergei Shavirin

Contents: main routines for copymatrices program to convert
score matrices output by makematrices into a single byte-encoded file.
   
$Log: copymat.c,v $
Revision 6.10  2000/01/13 15:27:10  shavirin
Added concatenation of files into single file (for later formatdb).

Revision 6.9  2000/01/12 14:39:46  shavirin
Added parameter to set cache size in lookup table foe RPS Blast.

Revision 6.8  2000/01/07 22:31:47  shavirin
Lookup table header now has notice, that this is single table.

Revision 6.7  1999/12/30 18:34:20  shavirin
Last row in the matrix for every sequence will be gap-row (-INT2_MAX)

Revision 6.6  1999/12/29 18:49:29  shavirin
Changed a little format of RPS lookup tables file.


*****************************************************************************/


#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <ncbi.h>
#include <objseq.h>
#include <objsset.h>
#include <sequtil.h>
#include <seqport.h>
#include <tofasta.h>
#include <blast.h>
#include <blastpri.h>
#include <txalign.h>
#include <simutil.h>
#include <posit.h>
#include <profiles.h>
#include <gapxdrop.h>

#define NUMARG (sizeof(myargs)/sizeof(myargs[0]))

static Args myargs [] = {
    { "Database for matrix profiles", /* 0 */
      "stdin", NULL, NULL, FALSE, 'P', ARG_FILE_IN, 0.0, 0, NULL},
    { "Print help; overrides all other arguments", /* 1 */
      "F", NULL, NULL, FALSE, 'H', ARG_BOOLEAN, 0.0, 0, NULL},
    { "Create RPS mem map file(s)", /* 2 */
      "F", NULL, NULL, FALSE, 'r', ARG_BOOLEAN, 0.0, 0, NULL},
    { "Cache size for lookup table (RPS Blast)", /* 3 */
      "3", NULL, NULL, FALSE, 's', ARG_INT, 0.0, 0, NULL},
}; 

/*counts the number of items in sequencesFile and matricesFile, assumed to
  be one per line, and checks that the numbers are equal.
  returns the number if equal, 0 if unequal, rewinds the file descriptors
  before returning*/
static Int4 countProfiles(FILE *sequencesFile, FILE *matricesFile)
{
    Int4 sequencesCount = 0; /*count for sequencesFile*/
    Int4 matricesCount = 0; /*count for matricesFile*/
    Char oneFileName[MAXLINELEN]; /*for reading one line per file*/
    
    while (fgets(oneFileName,MAXLINELEN,sequencesFile))
        sequencesCount++;
    while (fgets(oneFileName,MAXLINELEN,matricesFile))
        matricesCount++;
    rewind(matricesFile);
    rewind(sequencesFile);
    if (sequencesCount == matricesCount)
        return(sequencesCount);
    else {
        ErrPostEx(SEV_FATAL, 0, 0, "copymatrices: Sequences file has %d entries; Matrices file has %d entries; these should be equal\n", sequencesCount,matricesCount);
        return(0);
    }
}

/*free the memory associated with the position-specific score matrices*/
static void  freeMatrix(ScoreRow *posMatrix)
{

  MemFree(posMatrix);
}

/*allocate memory for the position-specific score matrices
  enough memory is allocated to hold the largest matrix
  the memory is reused for each different matrix*/
static ScoreRow * allocateMatrix(Int4 maxSequenceLength)
{
  Int4 i; /*row index for matrix*/
  ScoreRow *returnMatrix; /*matrix to return*/

  returnMatrix = (ScoreRow *) MemNew(maxSequenceLength * sizeof(ScoreRow));
  return(returnMatrix);
}

/* read in a position-specific score matrix from thisMatrixFile
   the number of positions is dbSequenceLength
   kbp keeps the Karlin-ALtschul parameters
   returnMatrix is the memory address where the matrix is to be stored*/
static void readNextMatrix(FILE * thisMatrixFile,
              Int4 startPos, Int4 *endPos,
              ScoreRow *bigMatrix)
{
  Int4 i, r; /*row indices for sequence and matrix*/
  Int4 lengthInFile; /*length of query*/
  Char junkChar; /*used to read in useless characters*/
  Nlm_FloatHi junkLambda, junkK, junklogK, junkH; /*used to read in useless
						    Karlin blocks*/
  Char *sequence;  /*sequence to read in*/
  Char rowOfScores[MAXLINELEN]; /*one row of scores to be read in*/

  fscanf(thisMatrixFile, "%d", &lengthInFile);
  sequence = (Char *) MemNew((lengthInFile + 2) * sizeof(Char));
  fscanf(thisMatrixFile,"%s",sequence);
  MemFree(sequence);
  /*read in useless Karlin block*/
  fscanf(thisMatrixFile,"%le", &junkLambda);
  fscanf(thisMatrixFile,"%le", &junkK);
  fscanf(thisMatrixFile,"%le", &junklogK);
  fscanf(thisMatrixFile,"%le", &junkH);
  /*read in useless Karlin block*/
  fscanf(thisMatrixFile,"%le", &junkLambda);
  fscanf(thisMatrixFile,"%le", &junkK);
  fscanf(thisMatrixFile,"%le", &junklogK);
  fscanf(thisMatrixFile,"%le", &junkH);
  /*read in useless Karlin block*/
  fscanf(thisMatrixFile,"%le", &junkLambda);
  fscanf(thisMatrixFile,"%le", &junkK);
  fscanf(thisMatrixFile,"%le", &junklogK);
  fscanf(thisMatrixFile,"%le\n", &junkH);
  for(i = 0, r = startPos; i < lengthInFile; i++, r++) {
    fgets(rowOfScores, MAXLINELEN, thisMatrixFile);
    sscanf(rowOfScores, "%d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d", 
                        &(bigMatrix[r][0]),
                        &(bigMatrix[r][1]),
                        &(bigMatrix[r][2]),
                        &(bigMatrix[r][3]),
                        &(bigMatrix[r][4]),
                        &(bigMatrix[r][5]),
                        &(bigMatrix[r][6]),
                        &(bigMatrix[r][7]),
                        &(bigMatrix[r][8]),
                        &(bigMatrix[r][9]),
                        &(bigMatrix[r][10]),
                        &(bigMatrix[r][11]),
                        &(bigMatrix[r][12]),
                        &(bigMatrix[r][13]),
                        &(bigMatrix[r][14]),
                        &(bigMatrix[r][15]),
                        &(bigMatrix[r][16]),
                        &(bigMatrix[r][17]),
                        &(bigMatrix[r][18]),
                        &(bigMatrix[r][19]),
                        &(bigMatrix[r][20]),
                        &(bigMatrix[r][21]),
                        &(bigMatrix[r][22]),
                        &(bigMatrix[r][23]),
                        &(bigMatrix[r][24]),
                        &(bigMatrix[r][25]));
  }

  if((Boolean) myargs[2].intvalue) {
      /* Last row in the matrix will be gap-row (-INT2_MAX) */
      
      for(i = 0; i < 26; i++) {
          bigMatrix[r][i] = -INT2_MAX;
      }
      r++;
  }

  *endPos = r;
}

/*read each matrix in turn and store its scores in combinedMatrix*/
static void readAllMatrices(FILE *matrixnamefp, ScoreRow *combinedMatrix, 
                            Int4 numProfiles, CharPtr directoryPrefix, 
                            Int4Ptr seqlens)
{
    Int4 i; /*loop index*/
    Char oneMatrixFileName[MAXLINELEN]; /*name of matrix file to read*/
    FILE *thisMatrixFile; /*descriptor for one matrix file*/
    Int4 startPos; /*starting row in big matrix for next small matrix*/
    Int4 endPos; /*ending row + 1 in big matrix for this small matrix*/
    Int4 prefixLength; /*length of directoryPrefix*/
    Int4 c1,c2; /*loop indices over characters*/
    Char relativeMatrixFileName[MAXLINELEN];
    
    startPos = 0;
    endPos = 0;
    if ('\0' != directoryPrefix[0]) {
        strcpy(oneMatrixFileName, directoryPrefix);
        prefixLength = strlen(directoryPrefix);
    }     
    for (i = 0; i < numProfiles; i++) {
        if ('\0' == directoryPrefix[0])
            fscanf(matrixnamefp,"%s", oneMatrixFileName);
        else {
            fscanf(matrixnamefp,"%s", relativeMatrixFileName); 
            for(c1 = prefixLength, c2 = 0; relativeMatrixFileName[c2] != '\0';
                c1++, c2++)
                oneMatrixFileName[c1] = relativeMatrixFileName[c2];
            oneMatrixFileName[c1] = '\0';
        }
        
        if ((thisMatrixFile = FileOpen(oneMatrixFileName, "r")) == NULL)  {
            ErrPostEx(SEV_FATAL, 0, 0, "profiles: Unable to open matrix file %s\n", oneMatrixFileName);
            return;
        }
        readNextMatrix(thisMatrixFile, startPos, &endPos,
                       combinedMatrix);

        if(seqlens != NULL) {
            seqlens[i] = startPos;
        }

        startPos = endPos;
        FileClose(thisMatrixFile);
    }

    if(seqlens != NULL) {   /* Last entry - is the end of last sequence */
        seqlens[i] = startPos;
    }
    
    return;
}

/*findTotalLength scans matrixAuxiliaryFile to find the
  total  number of positions among all the position-specific matrices*/
static Int4 findTotalLength(FILE *matrixAuxiliaryFile, Int4 numProfiles)
{
    Int4 maxLength; /*maximum length of sequence*/
    Int4 thisLength; /*length of next sequence*/
    Int4 totalLength; /*total length to return*/
    Int4 dbLength; /*length of database*/
    Int4 i; /*loop index*/
    Nlm_FloatHi Kungapped, Hungapped; /*two values to read*/
    Char * underlyingMatrixName; /*name of matrix to read*/
    Int4 gap_open, gap_extend; /*gap costs to skip over in reading*/
    Nlm_FloatHi scalingFactor; /*matrix scale to skip over in reading*/
    
    underlyingMatrixName = MemNew(MAXLINELEN * sizeof(Char));
    fscanf(matrixAuxiliaryFile,"%s",underlyingMatrixName);
    fscanf(matrixAuxiliaryFile,"%d\n", &gap_open);
    fscanf(matrixAuxiliaryFile,"%d\n", &gap_extend);
    fscanf(matrixAuxiliaryFile, "%le", &Kungapped);
    fscanf(matrixAuxiliaryFile, "%le", &Hungapped);
    fscanf(matrixAuxiliaryFile, "%d", &maxLength);
    fscanf(matrixAuxiliaryFile, "%d", &dbLength);
    fscanf(matrixAuxiliaryFile, "%lf", &scalingFactor);
    totalLength = 0;
    for (i = 0; i < numProfiles; i++) {
        fscanf(matrixAuxiliaryFile, "%d", &thisLength);
        fscanf(matrixAuxiliaryFile, "%le", &Kungapped);
        totalLength += thisLength;
    }
    rewind(matrixAuxiliaryFile);
    MemFree(underlyingMatrixName);
    return(totalLength);
}
#define RPS_THRESHOLD 11
#define RPS_WORDSIZE  3

Boolean RPSUpdatePointers(LookupTablePtr lookup)
{
    ModLAEntry * mod_lt;
    Uint4 len;
    Int4 index;
    ModLookupPositionPtr start_address;

    len = lookup->array_size;
    mod_lt=lookup->mod_lt;
    start_address = lookup->mod_lookup_table_memory;
    
    /* Walk through table, copying info into mod_lt[] */
    for(index = 0; index < len; index++) {

        if(mod_lt[index].num_used <= 3)
            continue;

        /* Changing value, that lpp now has relative pointer */
        mod_lt[index].entries[1] -= (int) start_address;
        mod_lt[index].entries[2] = 0; /* Not used */
    }
    
    return TRUE;
}

Boolean RPSDumpLookupTable(LookupTablePtr lookup, FILE *fd)
{

    RPSUpdatePointers(lookup);

    FileWrite(lookup->mod_lt, sizeof(ModLAEntry), 1+lookup->array_size, fd);
    if(lookup->mod_lookup_table_size) {
        FileWrite(lookup->mod_lookup_table_memory, 
                  lookup->mod_lookup_table_size, 
                  sizeof(ModLookupPosition), fd);
    }
    
    return TRUE;
}

Boolean RPSCreateLookupFile(ScoreRow *combinedMatrix, Int4 numProfiles,
                            Int4Ptr seqlens, CharPtr filename, 
                            Int4 TheCacheSize)
{
    LookupTablePtr lookup;
    BlastAllWordPtr all_words;
    Int4 start, len, i, header_size, all_length, magicNumber;
    FILE *fd;
    Int4Ptr offsets;
    BLAST_ScorePtr PNTR posMatrix;
    Int4 num_lookups;

    if((fd = FileOpen(filename, "w")) == NULL)
        return FALSE;
    
    num_lookups = 1; /* Single lookup table for all set */
    lookup = lookup_new(PRO_ALPHABET_SIZE, 3, 0, TheCacheSize);

    all_words = BlastPopulateAllWordArrays(3, PRO_ALPHABET_SIZE);

    /* Only Uint4 maximum length for lookup file allowed in current
       implementation */

    header_size = (numProfiles+1)*sizeof(Int4) + 8*sizeof(Int4);
    
    /* Beginning of file will be allocated for lookup offsets */
    fseek(fd, header_size, SEEK_SET);
    
    offsets = MemNew(sizeof(Int4) * (num_lookups + 1));
    
    all_length = seqlens[numProfiles] - seqlens[0];
    
    posMatrix = MemNew((all_length + 1) * sizeof(BLAST_Score));
    for (i = 0; i < all_length; i++) {
        posMatrix[i] = (BLAST_Score *) &(combinedMatrix[i][0]);
    }
    
    /* Last row is necessary */
    posMatrix[all_length] = MemNew(sizeof(BLAST_Score) * PRO_ALPHABET_SIZE);
    for(i = 0; i < PRO_ALPHABET_SIZE; i++) {
        posMatrix[all_length][i] = -INT2_MAX;
    }
#if 0
    for(i = 0; i < numProfiles; i++) {
        
        offsets[i] = ftell(fd);

        start = seqlens[i];
        len = seqlens[i+1] - seqlens[i];

        if(BlastNewFindWordsEx(lookup, &posMatrix[start], 0, len, all_words, 
                               RPS_THRESHOLD, RPS_WORDSIZE, 0) < 0) {
            ErrPostEx(SEV_ERROR, 0,0, "Failure to create llokup table");
            return FALSE;
        }
        
        lookup_position_aux_destruct(lookup);

        RPSDumpLookupTable(lookup, fd);
        
        lookup = lookup_destruct(lookup);
        lookup = lookup_new(PRO_ALPHABET_SIZE, 3, 0);
    }
#else
    offsets[0] = ftell(fd);
    
    start = seqlens[0]; /* 0 */
    
    if(BlastNewFindWordsEx(lookup, &posMatrix[start], 0, all_length, 
                           all_words, RPS_THRESHOLD, RPS_WORDSIZE, 0) < 0) {
        ErrPostEx(SEV_ERROR, 0,0, "Failure to create llokup table");
        return FALSE;
    }
    
    lookup_position_aux_destruct(lookup);
    
    RPSDumpLookupTable(lookup, fd);
    
    /*    lookup = lookup_new(PRO_ALPHABET_SIZE, 3, 0); */
    i = 1;
#endif

    offsets[i] = ftell(fd); /* Last offset also recorded */

    fseek(fd, 0, SEEK_SET);
    magicNumber = RPS_MAGIC_NUMBER;
    FileWrite(&magicNumber, sizeof(Int4), 1, fd); /* header[0] */
    FileWrite(&num_lookups, sizeof(Int4), 1, fd); /* header[1] */
    FileWrite(&lookup->num_pos_added, sizeof(Int4), 1, fd); /* header[2] */
    FileWrite(&lookup->num_unique_pos_added, sizeof(Int4), 1, fd); /* header[3] */
    FileWrite(&lookup->mod_lookup_table_size, sizeof(Int4), 1, fd); /* header[4] */
    
    /* Now writing recorded offsets in the beginning of the file */
    
    fseek(fd, 8*sizeof(Int4), SEEK_SET);
    FileWrite(offsets, sizeof(Int4), num_lookups + 1, fd);
    FileClose(fd);
    
    /* Final memory cleenup */
    
    lookup = lookup_destruct(lookup);

    MemFree(posMatrix[numProfiles]);
    MemFree(posMatrix);

    BlastAllWordDestruct(all_words);
    
    return TRUE;
}

Boolean RPSConcatSequences(FILE *sfp, CharPtr fastaname)
{
    FILE *fasta_fp, *fd;
    Char oneFileName[MAXLINELEN]; /*for reading one line per file*/
    Char buffer[1024];
    Int4 bytes;
    CharPtr chptr;

    if((fasta_fp = FileOpen(fastaname, "w")) == NULL)
        return FALSE;

    rewind(sfp);
    
    while (fgets(oneFileName, MAXLINELEN, sfp)) {

        /* Filtering out '\n' and '\r' characters */
        for(chptr = oneFileName; *chptr != NULLB; chptr++) {
            if(*chptr == '\n' || *chptr == '\r')
                *chptr = NULLB;
        }
        
        if((fd = FileOpen(oneFileName, "r")) == NULL) {
            FileClose(fasta_fp);
            return FALSE;
        }
        
        /* Now concatenating this file into set */
        while((bytes = FileRead(buffer, 1, 1024, fd)) > 0)
            FileWrite(buffer, 1, bytes, fasta_fp);
        FileClose(fd);
    }
    
    FileClose(fasta_fp);
    
    return TRUE;
}

Int2  Main(void)

{
    
    Char *profilesFileName; /*file name for list of profile file names*/
    Char sequencesFileName[MAX_NAME_LENGTH]; /*file anme for list of sequence file names*/
    Char matrixFileName[MAX_NAME_LENGTH]; /*file name for list of matrix file names*/
    Char auxFileName[MAX_NAME_LENGTH]; /*file name for file containing auxiliary information*/
    Char bigFileName[MAX_NAME_LENGTH]; /*file name to store byte-encoded coalesced matrix*/
    Char lookupName[MAX_NAME_LENGTH]; /*file name to store precalculated lookup table */
    FILE *auxiliaryfp; /*file descriptor for matrix auxiliary file*/
    FILE *sequencesfp; /*files descriptor for file containing list of sequences*/
    FILE *matrixnamefp; /*file descriptor for file containing matrix names*/
    FILE *bigmatrixfile; /*file descriptor for file containing single big matrix*/
    Int4  queryLength;  /*length of query sequence*/
    Int4 maxLength; /*maximum length of a sequnce*/
    Int4 numProfiles; /*number of profiles*/
    Int4 totalProfileLength; /*total length of all profiles*/
    ScoreRow *combinedMatrix; /*combined matrix for all profiles*/
    Char *directoryPrefix; /*directory where profile library is kept, used
                             to reach other directories indirectly*/

    Int4Ptr seqlens;

    if (! GetArgs ("copymatrices", NUMARG, myargs)) {
        return (1);
    }
    
    if ((Boolean) myargs[1].intvalue) {
        IMPALAPrintHelp(FALSE, 80, "copymat", stdout);
        return(1);
    }
    profilesFileName = myargs[0].strvalue;
    directoryPrefix = (Char *) MemNew(MAX_NAME_LENGTH *sizeof(char));
    strcpy(directoryPrefix,profilesFileName);
    
    impalaMakeFileNames(profilesFileName, auxFileName, bigFileName,
                        sequencesFileName, matrixFileName, NULL, 
                        directoryPrefix);
    
    if ((matrixnamefp = FileOpen(matrixFileName, "r")) == NULL) {
        ErrPostEx(SEV_FATAL, 0, 0, "copymatrices: Unable to open file with matrix file names %s\n", matrixFileName);
        return (1);
    }
    
    if ((sequencesfp = FileOpen(sequencesFileName, "r")) == NULL) {
        ErrPostEx(SEV_FATAL, 0, 0, "copymatrices: Unable to open file with sequence file names %s\n", sequencesFileName);
        return (1);
    }
    
    if ((auxiliaryfp = FileOpen(auxFileName, "r")) == NULL) {
        ErrPostEx(SEV_FATAL, 0, 0, "profiles: Unable to open auxiliary file %s\n", auxFileName);
        return (1);
    }

    /* Name of matrix file depends on program - RPS or Impala */
    
    if((Boolean) myargs[2].intvalue) {
        sprintf(bigFileName, "%s.rps", profilesFileName);
    }
    
    if ((bigmatrixfile = FileOpen(bigFileName, "w")) == NULL) {
        ErrPostEx(SEV_FATAL, 0, 0, "profiles: Unable to open big matrix file %s\n", bigFileName);
        return (1);
    }
    
    numProfiles =  countProfiles(sequencesfp, matrixnamefp);
    totalProfileLength = findTotalLength(auxiliaryfp, numProfiles);

    /* Additional line in matrix with -INT2_MAX values */
    if((Boolean) myargs[2].intvalue) {
        totalProfileLength += numProfiles;
    }


    combinedMatrix = allocateMatrix(totalProfileLength);
    if (NULL == combinedMatrix) {
        ErrPostEx(SEV_FATAL, 0, 0, "copymatrices: Unable to allocate matrix with%d rows\n", totalProfileLength);
        return (1);
        
    }

    if ((Boolean) myargs[2].intvalue) {
        seqlens = (Int4Ptr) MemNew((numProfiles +1) * sizeof(Int4));
    } else {
        seqlens = NULL;
    }
    
    readAllMatrices(matrixnamefp, combinedMatrix, numProfiles,
                    directoryPrefix, seqlens);
    
    /* For RPS Blast some additional info will be added to the file */
    if ((Boolean) myargs[2].intvalue) {
        Int4 magicNumber = RPS_MAGIC_NUMBER;
        FileWrite(&magicNumber, sizeof(Int4), 1, bigmatrixfile);
        FileWrite(&numProfiles, sizeof(Int4), 1, bigmatrixfile);
        FileWrite(seqlens, sizeof(Int4), numProfiles + 1, bigmatrixfile);
        
        sprintf(lookupName, "%s.loo", profilesFileName);
        RPSCreateLookupFile(combinedMatrix, numProfiles, seqlens, lookupName,
                            myargs[3].intvalue);

        if(!RPSConcatSequences(sequencesfp, profilesFileName)) {
            ErrPostEx(SEV_ERROR, 0,0, "Failure to concatenate sequences");
            return 1;
        }
        
    }
    
    FileWrite((void *) &(combinedMatrix[0]), sizeof(ScoreRow), 
              (size_t) totalProfileLength, bigmatrixfile);
    freeMatrix(combinedMatrix); 
    FileClose(bigmatrixfile);
    FileClose(matrixnamefp);
    FileClose(sequencesfp);
    FileClose(auxiliaryfp);
    return 0;
}
