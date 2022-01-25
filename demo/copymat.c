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

Author: Alejandro Schaffer

Contents: main routines for copymatrices program to convert
score matrices output by makematrices into a single byte-encoded file.
   

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
  *endPos = r;
}

/*read each matrix in turn and store its scores in combinedMatrix*/
static void readAllMatrices(FILE * matrixnamefp, ScoreRow * combinedMatrix, Int4 numProfiles, Char * directoryPrefix)
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
    startPos = endPos;
    FileClose(thisMatrixFile);
  }
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

#define NUMARG 2

static Args myargs [NUMARG] = {
  { "Database for matrix profiles", 
	"stdin", NULL, NULL, FALSE, 'P', ARG_FILE_IN, 0.0, 0, NULL},
  { "Print help; overrides all other arguments",
        "F", NULL, NULL, FALSE, 'H', ARG_BOOLEAN, 0.0, 0, NULL}
}; 

Int2  Main(void)

{

   Char *profilesFileName; /*file name for list of profile file names*/
   Char sequencesFileName[MAX_NAME_LENGTH]; /*file anme for list of sequence file names*/
   Char matrixFileName[MAX_NAME_LENGTH]; /*file name for list of matrix file names*/
   Char auxFileName[MAX_NAME_LENGTH]; /*file name for file containing auxiliary information*/
   Char bigFileName[MAX_NAME_LENGTH]; /*file name to store byte-encoded coalesced matrix*/
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
   
   if (! GetArgs ("copymatrices", NUMARG, myargs))
     {
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

   if ((matrixnamefp = FileOpen(matrixFileName, "r")) == NULL)
     {
       ErrPostEx(SEV_FATAL, 0, 0, "copymatrices: Unable to open file with matrix file names %s\n", matrixFileName);
       return (1);
     }

   if ((sequencesfp = FileOpen(sequencesFileName, "r")) == NULL)
     {
       ErrPostEx(SEV_FATAL, 0, 0, "copymatrices: Unable to open file with sequence file names %s\n", sequencesFileName);
       return (1);
     }

   if ((auxiliaryfp = FileOpen(auxFileName, "r")) == NULL)
     {
       ErrPostEx(SEV_FATAL, 0, 0, "profiles: Unable to open auxiliary file %s\n", auxFileName);
       return (1);
     }

   if ((bigmatrixfile = FileOpen(bigFileName, "w")) == NULL)
     {
       ErrPostEx(SEV_FATAL, 0, 0, "profiles: Unable to open big matrix file %s\n", bigFileName);
       return (1);
     }

   numProfiles =  countProfiles(sequencesfp, matrixnamefp);
   totalProfileLength = findTotalLength(auxiliaryfp, numProfiles);
   combinedMatrix = allocateMatrix(totalProfileLength);
   if (NULL == combinedMatrix) {
       ErrPostEx(SEV_FATAL, 0, 0, "copymatrices: Unable to allocate matrix with%d rows\n", totalProfileLength);
       return (1);

   }
   readAllMatrices(matrixnamefp, combinedMatrix, numProfiles,
         directoryPrefix);
   /*  write(bigmatrixfile, (void *) &(combinedMatrix[0]), sizeof(ScoreRow) * totalProfileLength);*/
   FileWrite((void *) &(combinedMatrix[0]), sizeof(ScoreRow) , (size_t) totalProfileLength, bigmatrixfile);
   freeMatrix(combinedMatrix); 
   FileClose(bigmatrixfile);
   FileClose(matrixnamefp);
   FileClose(sequencesfp);
   FileClose(auxiliaryfp);
   return 0;
}
