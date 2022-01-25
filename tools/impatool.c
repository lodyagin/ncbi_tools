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

File name: impatool.c

Author: Alejandro Schaffer

Contents: utility routines for IMPALA.

*****************************************************************************/



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


/*convert a residue character to its integer representation*/
Char LIBCALL getRes(Char input)
{
    switch(input) 
      {
      case 0: 
	return('-');
      case 1: 
	return('A');
      case 2: 
	return('B');
      case 3: 
	return('C');
      case 4: 
	return('D');
      case 5: 
	return('E');
      case 6: 
	return('F');
      case 7: 
	return('G');
      case 8: 
	return('H');
      case 9: 
	return('I');
      case 10: 
	return('K');
      case 11: 
	return('L');
      case 12: 
	return('M');
      case 13: 
	return('N');
      case 14: 
	return('P');
      case 15: 
	return('Q');
      case 16: 
	return('R');
      case 17: 
	return('S');
      case 18: 
	return('T');
      case 19: 
	return('V');
      case 20: 	
	return('W');
      case 21: 
	return('X');
      case 22: 
	return('Y');
      case 23: 
	return('Z');
      case 24: 
	return('U');
      case 25: 
	return('*');
      default:
        return('?');
    }
} 

/*
	adds the new string to the buffer, separating by a tilde.
	Checks the size of the buffer for FormatBlastParameters and
	allocates longer replacement if needed.
*/

static Boolean 
add_string_to_bufferEx(CharPtr buffer, CharPtr *old, Int2Ptr old_length, Boolean add_tilde)

{
	CharPtr new, ptr;
	Int2 length, new_length;

	length = (StringLen(*old));

	if((StringLen(buffer)+length+3) > *old_length)
	{
		new_length = *old_length + 255;
		new = MemNew(new_length*sizeof(Char));
		if (*old_length > 0 && *old != NULL)
		{
			MemCpy(new, *old, *old_length);
			*old = MemFree(*old);
		}
		*old = new;
		*old_length = new_length;
	}

	ptr = *old;
	ptr += length;
	if (add_tilde)
	{
		*ptr = '~';
		ptr++;

	}

	while (*buffer != NULLB)
	{
		*ptr = *buffer;
		buffer++; ptr++;
	}

	return TRUE;
}


/*get the citation for IMPALA*/

static CharPtr 
IMPALAGetReference(Boolean html)

{
	CharPtr ret_buffer;
	Int2 ret_buffer_length;

	ret_buffer = NULL;
	ret_buffer_length = 0;

	
	if (html) {
          ;
	} else
		add_string_to_bufferEx("Reference: Alejandro A. Schaffer, Yuri I. Wolf, Chris P. Ponting, ", &ret_buffer, &ret_buffer_length, TRUE);
	add_string_to_bufferEx("Eugene V. Koonin, L. Aravind, Stephen F. Altschul (1999), ", &ret_buffer, &ret_buffer_length, TRUE);
	add_string_to_bufferEx("\"IMPALA: Matching a Protein Sequence Against a Collection of ", &ret_buffer, &ret_buffer_length, TRUE);
	add_string_to_bufferEx("\"PSI-BLAST-Constructed Position-Specific Score Matrices\",", &ret_buffer, &ret_buffer_length, TRUE);
	add_string_to_bufferEx("Bioinformatics 15:1000-1011.", &ret_buffer, &ret_buffer_length, TRUE);
	return ret_buffer;
}

/*print the citation for IMPALA*/
Boolean LIBCALL
IMPALAPrintReference(Boolean html, Int4 line_length, FILE *outfp)

{
	CharPtr ret_buffer;
	
	if (outfp == NULL)
		return FALSE;
	
        ret_buffer = IMPALAGetReference(html);
        PrintTildeSepLines(ret_buffer, line_length, outfp);
        ret_buffer = MemFree(ret_buffer);

	return TRUE;
}

static CharPtr 
makematGetHelp(Boolean html)
{

	CharPtr ret_buffer;
	Int2 ret_buffer_length;

	ret_buffer = NULL;
	ret_buffer_length = 0;

	
	if (html) {
          ;
	} else {
	  add_string_to_bufferEx("makemat is the first profile preprocessor in the IMPALA software package.", &ret_buffer, &ret_buffer_length, TRUE);
	  add_string_to_bufferEx("makemat  converts a collection of byte-encoded profiles, created by the -C option of PSI-BLAST, into portable ASCII form.", &ret_buffer, &ret_buffer_length, TRUE);
	  add_string_to_bufferEx("Prepare the following files: ", &ret_buffer, &ret_buffer_length, TRUE);
	  add_string_to_bufferEx("i. a collection of PSI-BLAST-generated profiles with arbitrary   names and suffix .chk  ", &ret_buffer, &ret_buffer_length, TRUE);
	  add_string_to_bufferEx("ii. a collection of master sequences, associated with the profiles, each in a separate file with arbitrary name and a 3 character suffix starting with c", &ret_buffer, &ret_buffer_length, TRUE);
	  add_string_to_bufferEx("iii. a list of profile file names, one per line, named <database_name>.pn", &ret_buffer, &ret_buffer_length, TRUE);
	  add_string_to_bufferEx("iv. a list of master sequence file names, one per line, in the same order as a list of profile names, named <database_name>.sn", &ret_buffer, &ret_buffer_length, TRUE);
	  add_string_to_bufferEx("The following files will be created:", &ret_buffer, &ret_buffer_length, TRUE);
	  add_string_to_bufferEx("i. a collection of ASCII files, corresponding to each of the original profiles, named  <profile_name>.mtx", &ret_buffer, &ret_buffer_length, TRUE);
	  add_string_to_bufferEx("ii. a list of ASCII matrix files, named <database_name>.mn", &ret_buffer, &ret_buffer_length, TRUE);
	  add_string_to_bufferEx("iii. an ASCII file with auxiliary information, named <database_name>.aux", &ret_buffer, &ret_buffer_length, TRUE);
	  add_string_to_bufferEx("For the list of arguments to makemat enter makemat -", &ret_buffer, &ret_buffer_length, TRUE);
        }
	return ret_buffer;
}

static CharPtr 
copymatGetHelp(Boolean html)
{


	CharPtr ret_buffer;
	Int2 ret_buffer_length;

	ret_buffer = NULL;
	ret_buffer_length = 0;

	
	if (html) {
          ;
	} else {
	  add_string_to_bufferEx("copymat is the second profile preprocessor in the IMPALA software package.", &ret_buffer, &ret_buffer_length, TRUE);
	  add_string_to_bufferEx("copymat converts ASCII matrices, produced by the primary preprocessor, makemat into a database that can be read into memory quickly", &ret_buffer, &ret_buffer_length, TRUE);
	  add_string_to_bufferEx("Prepare the following files:", &ret_buffer, &ret_buffer_length, TRUE);
	  add_string_to_bufferEx("i. a collection of ASCII files, correspondingto each of the original profiles, named <profile_name>.mtx (created by makemat)", &ret_buffer, &ret_buffer_length, TRUE);
	  add_string_to_bufferEx("ii. a collection of profile master sequences, associated with the profiles, each in a separate file with arbitrary name and a 3-charecter suffix starting with c", &ret_buffer, &ret_buffer_length, TRUE);
	  add_string_to_bufferEx("iii. a list of ASCII_matrix files, named <database_name>.mn  (created by makemat)", &ret_buffer, &ret_buffer_length, TRUE);
	  add_string_to_bufferEx("iv. a list of master sequence file names, one per  line, in the same order as a list of matrix names, named <database_name>.sn;", &ret_buffer, &ret_buffer_length, TRUE);
	  add_string_to_bufferEx("v. ASCII file with auxiliary information, named  <database_name>.aux (created by makemat)", &ret_buffer, &ret_buffer_length, TRUE);
	  add_string_to_bufferEx("The files input to copymatices are in ASCII format and thus portable between machines with different encodings for machine-readable files", &ret_buffer, &ret_buffer_length, TRUE);
	  add_string_to_bufferEx("The following file will be created:", &ret_buffer, &ret_buffer_length, TRUE);
	  add_string_to_bufferEx("i. a huge binary file, containing all profile matrices, named <database_name>.mat", &ret_buffer, &ret_buffer_length, TRUE);
        }
	return ret_buffer;
}

static CharPtr 
impalaGetHelp(Boolean html)
{
	CharPtr ret_buffer;
	Int2 ret_buffer_length;

	ret_buffer = NULL;
	ret_buffer_length = 0;

	
	if (html) {
          ;
	} else {
	  add_string_to_bufferEx("impala is the main program in the IMPALA software package.", &ret_buffer, &ret_buffer_length, TRUE);
	  add_string_to_bufferEx("impala searches for matches between a query sequence and a library of score matrices, prepared by makemat and copymat.", &ret_buffer, &ret_buffer_length, TRUE);
	  add_string_to_bufferEx("impala produces BLAST-formatted output", &ret_buffer, &ret_buffer_length, TRUE);
	  add_string_to_bufferEx("Before you start searching, check that in the directory with the library you have the following files.", &ret_buffer, &ret_buffer_length, TRUE);
	  add_string_to_bufferEx("If the library has K matrices, you should have:", &ret_buffer, &ret_buffer_length, TRUE);
	  add_string_to_bufferEx("K files with names ending in .mtx", &ret_buffer, &ret_buffer_length, TRUE);
	  add_string_to_bufferEx("K files with names ending in a 3-letter extension starting with c", &ret_buffer, &ret_buffer_length, TRUE);
	  add_string_to_bufferEx("1 file with name ending in .sn", &ret_buffer, &ret_buffer_length, TRUE);
	  add_string_to_bufferEx("1 file with name ending in .aux", &ret_buffer, &ret_buffer_length, TRUE);
	  add_string_to_bufferEx("1 file with name ending in .mn", &ret_buffer, &ret_buffer_length, TRUE);
	  add_string_to_bufferEx("1 file with name ending in .mat", &ret_buffer, &ret_buffer_length, TRUE);
        }
	return ret_buffer;
}

Boolean LIBCALL
IMPALAPrintHelp(Boolean html, Int4 line_length, Char * programName, FILE *outfp)

{
	CharPtr ret_buffer;
	
	if (outfp == NULL)
		return FALSE;
	
        if (0 == strcmp(programName,"makemat"))
	  ret_buffer = makematGetHelp(html);
        if (0 == strcmp(programName,"copymat"))
	  ret_buffer = copymatGetHelp(html);
        if (0 == strcmp(programName,"impala"))
          ret_buffer = impalaGetHelp(html);

        PrintTildeSepLines(ret_buffer, line_length, outfp);
        ret_buffer = MemFree(ret_buffer);

	return TRUE;
}



/*prefix and suffix are strings, returns a string with prefix
  concatenated to suffix. Used to build multiple distinct
  file names with a common prefix*/
Char * LIBCALL addSuffixToName(Char *prefix, Char *suffix)
{
  Char *returnName; /*string to return*/
  Int4 i,j; /*loop indices*/
  Int4 prefixLength, suffixLength, totalLength; /*length of pieces and whole*/


  prefixLength = strlen(prefix);
  suffixLength = strlen(suffix);
  totalLength = prefixLength + suffixLength;
  returnName = (Char *) MemNew((totalLength + 1) * sizeof(Char));
  for(i = 0; i < prefixLength; i++)
    returnName[i] = prefix[i];
  for(j = 0, i = prefixLength; j < suffixLength; i++, j++)
    returnName[i] = suffix[j];
  returnName[totalLength] = '\0';
  return(returnName);
}

/*Use the specified name of the database to make file names
  for the auxiliary information, mmap'd score matrices, list of
  sequence files, and list of matrices, and lists of checkpoints.
  Also extract a directory prefix that is to be prepended
  to all file names for sequences and matrices and checkpoints */
void  LIBCALL impalaMakeFileNames(Char * matrixDbName,
		    Char * auxiliaryFileName, Char * mmapFileName,
		    Char * seqFileName, Char *matrixFileName, 
                    Char * ckptFileName,
		    Char *directoryPrefix)
{
  Int4 commonLength;  /*length of common name prefix*/
  Int4 c; /*loop index*/

  commonLength = strlen(matrixDbName);
  strcpy(directoryPrefix,matrixDbName);
  /*remove file name for profile library to peel back to last '/' */
  for(c = commonLength; c >= 0; c--)
    if ('/' == directoryPrefix[c]) 
      break;
  directoryPrefix[c+1] = '\0';
  if (NULL != auxiliaryFileName)
     strcpy(auxiliaryFileName,matrixDbName);
  if (NULL != mmapFileName)
     strcpy(mmapFileName,matrixDbName);
  if (NULL != seqFileName)
     strcpy(seqFileName,matrixDbName);
  if (NULL != matrixFileName)
     strcpy(matrixFileName,matrixDbName);
  if (NULL != ckptFileName)
     strcpy(ckptFileName,matrixDbName);
  if (NULL != auxiliaryFileName) {
    auxiliaryFileName[commonLength] = '.';
    auxiliaryFileName[commonLength +1] = 'a';
    auxiliaryFileName[commonLength + 2] = 'u';
    auxiliaryFileName[commonLength + 3] = 'x';
    auxiliaryFileName[commonLength + 4] = '\0';
  }
  if (NULL != mmapFileName) {
    mmapFileName[commonLength] = '.';
    mmapFileName[commonLength + 1] = 'm';
    mmapFileName[commonLength + 2] = 'a';
    mmapFileName[commonLength + 3] = 't';
    mmapFileName[commonLength + 4] = '\0';
  }
  if (NULL != seqFileName) {
    seqFileName[commonLength] = '.';
    seqFileName[commonLength + 1] = 's';
    seqFileName[commonLength+2] = 'n';
    seqFileName[commonLength+3] = '\0';
  }
  if (NULL != matrixFileName) {
    matrixFileName[commonLength] = '.';
    matrixFileName[commonLength + 1] = 'm';
    matrixFileName[commonLength+2] = 'n';
    matrixFileName[commonLength+3] = '\0';
  }
  if (NULL != ckptFileName) {
    ckptFileName[commonLength] = '.';
    ckptFileName[commonLength + 1] = 'p';
    ckptFileName[commonLength+2] = 'n';
    ckptFileName[commonLength+3] = '\0';
  }
}

Nlm_FloatHi LIBCALL
IMPALAfindUngappedLambda(Char *matrixName)
{
   if (0 == strcmp(matrixName, "BLOSUM62"))
     return(0.3176);
   if (0 == strcmp(matrixName, "BLOSUM90"))
     return(0.3346);
   if (0 == strcmp(matrixName, "BLOSUM80"))
     return(0.3430);
   if (0 == strcmp(matrixName, "BLOSUM50"))
     return(0.232);
   if (0 == strcmp(matrixName, "BLOSUM45"))
     return(0.2291);
   if (0 == strcmp(matrixName, "PAM30"))
     return(0.340);
   if (0 == strcmp(matrixName, "PAM70"))
     return(0.3345);
   if (0 == strcmp(matrixName, "PAM250"))
     return(0.229);
   return(0);
}

/*Given a sequence of 'length' amino acid residues, compute the
  probability of each residue and put that in the array resProb*/
void LIBCALL
IMPALAfillResidueProbability(Uint1Ptr sequence, Int4 length, Nlm_FloatHi * resProb)
{
  Int4 frequency[PRO_ALPHABET_SIZE]; /*frequency of each letter*/
  Int4 i; /*index*/
  Int4 denominator; /*length not including X's*/

  denominator = length;
  for(i = 0; i < PRO_ALPHABET_SIZE; i++)
    frequency[i] = 0;
  for(i = 0; i < length; i++)
    if (Xchar != sequence[i])
      frequency[sequence[i]]++;
    else
      denominator--;
  for(i = 0; i < PRO_ALPHABET_SIZE; i++) {
    if (frequency[i] == 0)
      resProb[i] = 0.0;
    else
      resProb[i] = ((Nlm_FloatHi) (frequency[i])) /((Nlm_FloatHi) denominator);
  }
}

/*matrix is a position-specific score matrix with matrixLength positions
  queryProbArray is an array containing the probability of occurrence
  of each residue in the query
  scoreArray is an array of probabilities for each score that is
    to be used as a field in return_sfp
  return_sfp is a the structure to be filled in and returned
  range is the size of scoreArray and is an upper bound on the
   difference between maximum score and minimum score in the matrix
  the routine fillSfp computes the probability of each score weighted
   by the probability of each query residue and fills those probabilities
   into scoreArray and puts scoreArray as a field in
   that in the structure that is returned
   for indexing convenience the field storing scoreArray points to the
   entry for score 0, so that referring to the -k index corresponds to
   score -k */
BLAST_ScoreFreqPtr LIBCALL
IMPALAfillSfp(BLAST_Score **matrix, Int4 matrixLength, Nlm_FloatHi *queryProbArray, Nlm_FloatHi *scoreArray,  BLAST_ScoreFreqPtr return_sfp, Int4 range)
{
  Int4 minScore, maxScore; /*observed minimum and maximum scores*/
  Int4 i,j; /* indices */
  Nlm_FloatHi onePosFrac; /*1/matrix length as a double*/

  minScore = maxScore = 0;

  for(i = 0; i < matrixLength; i++) {
    for(j = 0 ; j < PRO_ALPHABET_SIZE; j++) {
      if (Xchar == j)
        continue;
      if ((matrix[i][j] != BLAST_SCORE_MIN) && (matrix[i][j] < minScore))
	minScore = matrix[i][j];
      if (matrix[i][j] > maxScore)
        maxScore = matrix[i][j];
    }
  }
  return_sfp->obs_min = minScore;
  return_sfp->obs_max = maxScore;
  for (i = 0; i < range; i++)
    scoreArray[i] = 0.0;
  return_sfp->sprob = &(scoreArray[-minScore]); /*center around 0*/
  onePosFrac = 1.0/ ((Nlm_FloatHi) matrixLength);
  for(i = 0; i < matrixLength; i++) {
    for (j = 0; j < PRO_ALPHABET_SIZE; j++) {
      if (Xchar == j)
        continue;
      if(matrix[i][j] >= minScore) {
        return_sfp->sprob[matrix[i][j]] += (onePosFrac * queryProbArray[j]);
      }
    }
  }
  return_sfp->score_avg = 0;
  for(i = minScore; i <= maxScore; i++)
    return_sfp->score_avg += i * return_sfp->sprob[i];
  return(return_sfp);
}
