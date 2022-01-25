/* Include files */

#include <ncbi.h>
#include <objalign.h>

#include "aliparse.h"
#include "aliread.h"

/* Defined constants */

#define ALI_CORRUPT_SEQ_THRESHOLD 95

/* Data structures */

typedef struct
{
  CharPtr sequence;
  CharPtr id;
  CharPtr junk;
} SeqLineStruct, PNTR SeqLineStructPtr;

typedef struct
{
  CharPtr definitions;
  CharPtr id;
} DefLineStruct, PNTR DefLineStructPtr;

typedef struct
{
  CharPtr other;
  CharPtr id;
} OtherLine, PNTR OtherLinePtr;

/* Function prototypes */

Boolean     IsNucleotideChar (Char ch);
Boolean     IsProteinChar (Char ch);
Boolean     IsSequenceChar (Char ch, Char gapChar, Char missingChar);
Int2        IsValidIdChar (Char idChar);
Boolean     IsValidId (CharPtr idStr);
CharPtr     ReadAlignFileLine (FILE PNTR alignFilePtr);

static Boolean s_ParseDefLine (CharPtr          lineStr,
			       DefLineStructPtr defLine);
static Boolean s_MightBeCorruptSequence (Int4    seqCharCount,
					 CharPtr seqString,
					 Char    gapChar,
					 Char    missingChar);
static Int2    s_ParseSequenceLine (CharPtr          lineStr,
				    SeqLineStructPtr seqLine);
static Boolean s_ParseOtherLine (CharPtr      lineStr,
				 OtherLinePtr otherLine);

/*=========================================================================*/
/*                                                                         */
/*  IsNucleotideChar ()                                                    */
/*                                                                         */
/*=========================================================================*/

Boolean IsNucleotideChar (Char ch)
{
  if (StringChr("gctaGCTArymkswhbvdnxRYMKSWHBVDNX",ch) != NULL)
    return TRUE;
  else
    return FALSE;
}

/*=========================================================================*/
/*                                                                         */
/*  IsProteinChar ()                                                       */
/*                                                                         */
/*=========================================================================*/

Boolean IsProteinChar (Char ch)
{
  /* !!! TBD -- Need to know allowable protein chars */

  if (StringChr("ABCDEFGHIKLMNPQRSTUVWXYZ*abcdefghiklmnpqrstuvwxyz",ch) != NULL)
    return TRUE;
  else
    return FALSE;
}

/*=========================================================================*/
/*                                                                         */
/*  IsSequenceChar ()                                                      */
/*                                                                         */
/*=========================================================================*/

Boolean IsSequenceChar (Char ch, Char gapChar, Char missingChar)
{
  if (IsNucleotideChar(ch) || 
      IsProteinChar(ch) ||
      (ch == gapChar) ||
      (ch == missingChar))
    return TRUE;
  else
    return FALSE;
}

/*=========================================================================*/
/*                                                                         */
/* IsValidId ()                                                            */
/*                                                                         */
/*=========================================================================*/

#define ID_BAD_CHAR         0
#define ID_GOOD_CHAR_LETTER 1
#define ID_GOOD_CHAR_NUMBER 2
#define ID_GOOD_CHAR_OTHER  3

Int2 IsValidIdChar (Char idChar)
{
  if (StringChr("ABCDEFGHIJKLMNOPQRSTUVWXYZ",idChar) != NULL)
    return ID_GOOD_CHAR_LETTER;

  if (StringChr("abcdefghijklmnopqrstuvwxyz",idChar) != NULL)
    return ID_GOOD_CHAR_LETTER;

  if (StringChr("0123456789",idChar) != NULL)
    return ID_GOOD_CHAR_NUMBER;

  if (StringChr("._-",idChar) != NULL)
    return ID_GOOD_CHAR_OTHER;

  return FALSE;
}

/*=========================================================================*/
/*                                                                         */
/* IsValidId ()                                                            */
/*                                                                         */
/*=========================================================================*/

Boolean IsValidId (CharPtr idStr)
{
  Int4    position;
  Boolean letterFound = FALSE;
  Int2    charType;

  for (position = 0; idStr[position] != '\0'; position++)
    {
      charType = IsValidIdChar(idStr[position]);
      switch (charType)
	{
	case ID_GOOD_CHAR_LETTER :
	  letterFound = TRUE;
	  break;
	case ID_GOOD_CHAR_NUMBER :
	case ID_GOOD_CHAR_OTHER :
	  break;
	default:
	  return FALSE;
	}
    }

  if (!letterFound)
    return FALSE;

  return TRUE;
}

/*=========================================================================*/
/*                                                                         */
/* ReadAlignFileLine() -                                                   */
/*                                                                         */
/*=========================================================================*/

CharPtr ReadAlignFileLine (FILE PNTR alignFilePtr)
{
  CharPtr lineStr = NULL;
  CharPtr tempBuff = NULL;
  Int4    totalLen = 0;
  Int4    segmentLen = 0;
  Int4    segmentCount = 1;
  Boolean done = FALSE;
  Char    ch;

  /* Allocate memory for the line.  More */
  /* can be added later as necessary.    */

  lineStr = (CharPtr) MemNew(sizeof(Char) * ALIGN_FILE_BUFFSIZE);
  if (lineStr == NULL)
    {
      fprintf (stderr, "ERROR: Memory allocation failed.\n");
      return NULL;
    }

  /* If first non-whitespace char is  */
  /* a sequence char then we're good. */

  while (!done && !feof(alignFilePtr))
    {

      /* Process the current character */

      ch = (Char) getc(alignFilePtr);

      if ((ch == '\n') || (ch == '\r'))
	done = TRUE;
      else
	{
	  lineStr[totalLen] = ch;
	  segmentLen++;
	  totalLen++;
	}

      /* Allocate more memory for the */
      /* sequence if needed.          */

      if (segmentLen == ALIGN_FILE_BUFFSIZE)
	{
	  segmentCount++;
	  tempBuff = (CharPtr) MemNew(sizeof(Char) * 
				      segmentCount *
				      ALIGN_FILE_BUFFSIZE);
	  if (tempBuff == NULL)
	    {
	      fprintf (stderr, "ERROR: Memory allocation failed.\n");
	      return FALSE;
	    }
	  MemCpy(tempBuff, lineStr, segmentCount * ALIGN_FILE_BUFFSIZE);
	  MemFree(lineStr);
	  lineStr = tempBuff;
	  segmentLen = 0;
	}

    }

  /* Return successfully */

  lineStr[totalLen] = '\0';
  return lineStr;
}

/*=========================================================================*/
/*                                                                         */
/* s_ParseDefLine () -                                                     */
/*                                                                         */
/*=========================================================================*/

#define DEFLINE_PRE_DATA      0
#define DEFLINE_DEFINITION    1
#define DEFLINE_SEQID         2

static Boolean s_ParseDefLine (CharPtr          lineStr,
			       DefLineStructPtr defLine)
{
  Char    ch;
  CharPtr defStr;
  CharPtr idStr;
  Int4    defPosition;
  Int4    idPosition;
  Int4    position;
  Int2    state;

  defPosition = 0;
  idPosition = 0;

  defStr = (CharPtr) MemNew (StringLen(lineStr));
  idStr = (CharPtr) MemNew (StringLen(lineStr));

  /* Parse the line character by character */

  state       = DEFLINE_PRE_DATA;

  for (position = 0; lineStr[position] != '\0'; position++)
    {
      ch = lineStr[position];

      switch (state)
	{
	case DEFLINE_PRE_DATA :
	  if (IS_WHITESP(ch))
	    continue;
	  else if (ch == '>')
	    state = DEFLINE_SEQID;
	  else
	    return FALSE;  /* Not a defline */
	  break;
	case DEFLINE_SEQID : 
	  if (IS_ALPHANUM(ch))
	    {
	      idStr[idPosition] = ch;
	      idPosition++;
	    }
	  else if (IS_WHITESP(ch))	
	    {
	      if (idPosition > 0)
		{
		  state = DEFLINE_DEFINITION;
		  defStr[defPosition] = ch;
		  defPosition++;
		}
	      else
		continue;
	    }
	  else if (ch == '[')
	    {
	      state = DEFLINE_DEFINITION;
	      defStr[defPosition] = ch;
	      defPosition++;
	    }
	  else
	    {
	      return FALSE;
	    }
	  break;
	case DEFLINE_DEFINITION :
	  defStr[defPosition] = ch;
	  defPosition++;
	  break;
	default:
	  break;
	}
    }

  /* Check for blank line */
  
  if (state == DEFLINE_PRE_DATA)
    {
      MemFree(defStr);
      MemFree(idStr);
      return FALSE;
    }
  
  /* If we made it to here, then */
  /* it's a valid definition line. */

  idStr[idPosition]   = '\0';
  defStr[defPosition] = '\0';

  defLine->definitions = defStr;
  defLine->id          = idStr;

  return TRUE;
}

/*=========================================================================*/
/*                                                                         */
/* s_MightBeCorruptSequence ()                                             */
/*                                                                         */
/*=========================================================================*/

static Boolean s_MightBeCorruptSequence (Int4    seqCharCount,
					 CharPtr seqString,
					 Char    gapChar,
					 Char    missingChar)
{
  Int4    i;
  Int4    badCharCount;
  Int4    seqStrLen;
  FloatLo percentGood;

  seqStrLen = StringLen(seqString);
  badCharCount = 0;

  for (i = 0; i < seqStrLen; i++)
    {
      if (IsSequenceChar(seqString[i], gapChar, missingChar))
	seqCharCount++;
      else
	badCharCount++;
    }

  percentGood = (FloatLo) seqCharCount / ((FloatLo) seqCharCount + 
                                         (FloatLo) badCharCount);

  if ((percentGood * 100) >= ALI_CORRUPT_SEQ_THRESHOLD)
    return TRUE;
  else
    return FALSE;
}

/*=========================================================================*/
/*                                                                         */
/* s_ParseSequenceLine () -                                                */
/*                                                                         */
/*=========================================================================*/

#define PRE_DATA      0
#define FIRST_WORD    1
#define SEQUENCE_DATA 2
#define EOL_JUNK      3
#define POST_JUNK     4

static Int2 s_ParseSequenceLine (CharPtr          lineStr,
				 SeqLineStructPtr seqLine)
{
  CharPtr seqStr;
  Int4    seqPosition = 0;
  CharPtr idStr;
  Int4    idPosition = 0;
  Char    ch;
  Int2    state = PRE_DATA;
  Int4    position;
  Boolean firstWordNotSequence = FALSE;
  Boolean sequenceFound = FALSE;
  Char    gapChar = '-';
  Char    missingChar = '?';
  CharPtr tempStr;
  Boolean corruptSequence = FALSE;

  seqStr = (CharPtr) MemNew (StringLen(lineStr));
  idStr  = (CharPtr) MemNew (StringLen(lineStr));

  for (position = 0; lineStr[position] != '\0'; position++)
    {
      ch = lineStr[position];

      switch (state)
	{
	case PRE_DATA :

	  /* If it's the first non-whitespace char */
	  /* then we've found our first word.      */

	  if (!IS_WHITESP(ch))
	    {
	      state = FIRST_WORD;
	      if (!IsSequenceChar(ch, gapChar, missingChar))
		firstWordNotSequence = TRUE;
	      idStr[idPosition] = ch;
	      idPosition++;
	    }
	  break;
	case FIRST_WORD :
	  if (IS_WHITESP(ch))
	    {
	      state = SEQUENCE_DATA;
	      if ((idPosition > 0)   &&
		  (firstWordNotSequence == FALSE))
		{
		  tempStr = seqStr;
		  seqStr  = idStr;
		  idStr   = tempStr;
		  seqPosition = idPosition;
		  idPosition  = 0;
		  sequenceFound = TRUE;
		}
	    }
	  else
	    {
	      /* If we find a non-sequence char in the */
	      /* first word then it might be an ID,    */
	      /* with the sequence following.          */
	      
	      if (!IsSequenceChar(ch, gapChar, missingChar))
		firstWordNotSequence = TRUE;
	      idStr[idPosition] = ch;
	      idPosition++;
	    }
	  break;
	case SEQUENCE_DATA :
	  if (IS_WHITESP(ch))
	    continue;
	    
	  /* If we're in a sequence, then a non-sequence */
	  /* char invalidates it, although we do allow   */
	  /* 'junk' at the end.                          */
	  
	  if (!IsSequenceChar(ch, gapChar, missingChar))
	    {
	      if ((lineStr[position - 1] == ' ') && sequenceFound)
		state = EOL_JUNK;
	      else if ((corruptSequence == TRUE) ||
		       (s_MightBeCorruptSequence (seqPosition,
						  &(lineStr[position]),
						  gapChar,
						  missingChar)))
		{
		  seqStr[seqPosition] = ch;
		  seqPosition++;
		  sequenceFound = TRUE;
		  corruptSequence = TRUE;
		}
	      else
		{
		  MemFree(seqStr);
		  MemFree(idStr);
		  return ALI_UNKNOWN;
		}
	    }
	  else
	    {
	      seqStr[seqPosition] = ch;
	      seqPosition++;
	      sequenceFound = TRUE;
	    }
	  break;
	case EOL_JUNK :
	  if (IS_WHITESP(ch))
	    state = POST_JUNK;
	  break;
	case POST_JUNK :

	  /* Only one 'word' of junk allowed */

	  if (!IS_WHITESP(ch))
	    {
	      MemFree(seqStr);
	      MemFree(idStr);
	      return ALI_UNKNOWN;
	    }
	  break;
	}
    }

  /* Check for blank line */
  
  if (state == PRE_DATA)
    {
      MemFree(seqStr);
      MemFree(idStr);
      return ALI_UNKNOWN;
    }
  
  if (state == FIRST_WORD)
    {

      /* If there was just one word, and it isn't */
      /* a sequence string, then this isn't a     */
      /* sequence line.                           */
      
      if (firstWordNotSequence == TRUE)
	{
	  MemFree(seqStr);
	  MemFree(idStr);
	  return ALI_UNKNOWN;
	}
      
      /* If there was just one word, and it IS a sequence */
      /* then the idStr is actually the seqStr.           */
      
      else
	{
	  tempStr = seqStr;
	  seqStr  = idStr;
	  idStr   = tempStr;
	  seqPosition = idPosition;
	  idPosition  = 0;
	}
    }

  /* If still no sequence string, */
  /* then not a sequence line.    */

  if (StringLen(seqStr) == 0)
    return ALI_UNKNOWN;

  /* Check to see if the ID is a valid one */

  idStr[idPosition]   = '\0';
  seqStr[seqPosition] = '\0';

  if ((idPosition > 0) && (IsValidId (idStr) == FALSE))
    {
      MemFree(seqStr);
      MemFree(idStr);
      return ALI_UNKNOWN;
    }

  /* If we made it to here, then */
  /* it's a valid sequence line. */

  seqLine->sequence = seqStr;
  seqLine->id       = idStr;

  if (corruptSequence)
    return ALI_MAYBE_SEQUENCE;
  else
    return ALI_SEQUENCE;
}

/*=========================================================================*/
/*                                                                         */
/* s_ParseOtherLine () -                                                   */
/*                                                                         */
/*=========================================================================*/

#define OTHER_PRE_DATA  0
#define OTHER_DATA      1

static Boolean s_ParseOtherLine (CharPtr      lineStr,
				 OtherLinePtr otherLine)
{
  Char    ch;
  CharPtr otherStr;
  Int4    otherPosition;
  Int4    position;
  Int2    state;
  Int4    wordCount;

  /* Parse the line character by character */

  otherStr = (CharPtr) MemNew (StringLen(lineStr));
  otherPosition = 0;
  state     = OTHER_PRE_DATA;
  wordCount = 0;

  for (position = 0; lineStr[position] != '\0'; position++)
    {
      ch = lineStr[position];

      switch (state)
	{
	case OTHER_PRE_DATA :
	  if (IS_WHITESP(ch))
	    continue;
	  else
	    {
	      wordCount = 1;
	      state = OTHER_DATA;
	      otherStr[otherPosition] = ch;
	      otherPosition++;
	    }
	  break;
	case OTHER_DATA : 
	  if (IS_WHITESP(ch))	
	    wordCount++;
	  otherStr[otherPosition] = ch;
	  otherPosition++;
	  break;
	default:
	  break;
	}
    }

  /* Check for blank line */
  
  if (state == OTHER_PRE_DATA)
    {
      MemFree(otherStr);
      return FALSE;
    }
  
  /* If we made it to here, then */
  /* it's a valid definition line. */

  otherStr[otherPosition]   = '\0';

  if ((wordCount == 1) && IsValidId(otherStr))
    {
      otherLine->id    = otherStr;
      otherLine->other = NULL;
    }
  else
    {
      otherLine->id    = NULL;
      otherLine->other = otherStr;
    }

  return TRUE;
}

/*=========================================================================*/
/*                                                                         */
/* Ali_ReadLines ()                                                        */
/*                                                                         */
/*=========================================================================*/

ValNodePtr Ali_ReadLines (FILE PNTR  alignFilePtr,
			  ErrInfoPtr PNTR errorList)
{
  CharPtr        lineStr = NULL;
  ValNodePtr     rowList = NULL;
  ValNodePtr     newRow;
  Boolean        first = TRUE;
  SeqLineStruct  seqLine;
  DefLineStruct  defLine;
  OtherLine      otherLine;
  ValNodePtr     lastId;
  Boolean        nextRowMustBeSeq;
  Boolean        idFound;
  Boolean        lastRowWasOther = FALSE;
  Int2           lineType;

  nextRowMustBeSeq = FALSE;

  while (!feof(alignFilePtr))
    {

      /* Process the line according to its content ... */

      lineStr = ReadAlignFileLine(alignFilePtr);
      if (lineStr == NULL)
	return NULL;

      /* ... DefLine */
      
      if (s_ParseDefLine(lineStr, &defLine))
	{

	  lastRowWasOther = FALSE;

	  if (nextRowMustBeSeq)
	    {
	      lastId->choice = ALI_OTHER;
	      nextRowMustBeSeq = FALSE;
	    }

	  /* Add a record for the ID, if */
	  /* the row contained one.      */

	  if ((defLine.id != NULL) && (StringLen(defLine.id) != 0))
	    {
	      newRow = ValNodeAdd(&rowList);
	      if (NULL == newRow)
		{
		  Ali_AddError (errorList, ERR_OUT_OF_MEMORY);
		  return NULL;
		}

	      nextRowMustBeSeq = TRUE;
	      lastId = newRow;

	      newRow->choice = ALI_DEFLINE_ID;
	      newRow->data.ptrvalue = defLine.id;
	    }
	  else
	    nextRowMustBeSeq = FALSE;


	  /* Add a record for the definitions, */
	  /* if the row contained any.         */

	  if ((defLine.definitions != NULL) &&
	      (StringLen(defLine.definitions) != 0))
	    {
	      newRow = ValNodeAdd(&rowList);
	      if (NULL == newRow)
		{
		  Ali_AddError (errorList, ERR_OUT_OF_MEMORY);
		  return NULL;
		}
	      
	      newRow->choice = ALI_DEFLINE;
	      newRow->data.ptrvalue = defLine.definitions;

	    }
	}
      
      /* ... Sequence Data */

      else if ((lineType = s_ParseSequenceLine(lineStr, &seqLine))
	       != ALI_UNKNOWN)
	{
	  if (lineType == ALI_SEQUENCE)
	    {

	      /* Add a record for the ID, if */
	      /* the row contained one.      */
	      
	      if (StringLen(seqLine.id) != 0)
		{
		  newRow = ValNodeAdd(&rowList);
		  if (NULL == newRow)
		    {
		      Ali_AddError (errorList, ERR_OUT_OF_MEMORY);
		      return NULL;
		    }
	      
		  newRow->choice = ALI_SEQLINE_ID;
		  newRow->data.ptrvalue = seqLine.id;
		  
		  if (nextRowMustBeSeq)
		    lastId->choice = ALI_OTHER;
		  
		  lastId = newRow;
		  
		  lastRowWasOther = FALSE;
		}
	      
	      /* Add a record for the sequence */
	      
	      newRow = ValNodeAdd(&rowList);
	      if (NULL == newRow)
		{
		  Ali_AddError (errorList, ERR_OUT_OF_MEMORY);
		  return NULL;
		}
	      
	      newRow->choice = ALI_SEQUENCE;
	      newRow->data.ptrvalue = seqLine.sequence;
	      
	      if (nextRowMustBeSeq)
		nextRowMustBeSeq = FALSE;
	      
	      /* A sequence must follow either a defline */
	      /* an ID or another sequence.              */
	      
	      if (lastRowWasOther == TRUE)
		newRow->choice = ALI_OTHER;
	      
	      lastRowWasOther = FALSE;
	    }
	  else if (lineType == ALI_MAYBE_SEQUENCE)
	    {
	      /* Add a record for the ID, if */
	      /* the row contained one.      */
	      
	      if (StringLen(seqLine.id) != 0)
		{
		  newRow = ValNodeAdd(&rowList);
		  if (NULL == newRow)
		    {
		      Ali_AddError (errorList, ERR_OUT_OF_MEMORY);
		      return NULL;
		    }
	      
		  newRow->choice = ALI_MAYBE_ID;
		  newRow->data.ptrvalue = seqLine.id;
		}
	      
	      /* Add a record for the sequence */
	      
	      newRow = ValNodeAdd(&rowList);
	      if (NULL == newRow)
		{
		  Ali_AddError (errorList, ERR_OUT_OF_MEMORY);
		  return NULL;
		}
	      
	      newRow->choice = ALI_MAYBE_SEQUENCE;
	      newRow->data.ptrvalue = seqLine.sequence;

	    }
	}      
      /* ... Other */
      
      else
	{
	  if (StringLen(lineStr) > 0)
	    {
	      if (s_ParseOtherLine(lineStr, &otherLine) == TRUE)
		{
		  if (otherLine.id != NULL)
		    {
		      newRow = ValNodeAdd(&rowList);
		      if (NULL == newRow)
			{
			  Ali_AddError (errorList, ERR_OUT_OF_MEMORY);
			  return NULL;
			}

		      lastId = newRow;
		      idFound = TRUE;
		      
		      newRow->choice = ALI_OWNLINE_ID;
		      newRow->data.ptrvalue = otherLine.id;
		      lastRowWasOther = FALSE;
		    }
		  else
		    idFound = FALSE;

		  if (otherLine.other != NULL)
		    {
		      if (otherLine.id == NULL)
			lastRowWasOther = TRUE;
		      newRow = ValNodeAdd(&rowList);
		      if (NULL == newRow)
			{
			  Ali_AddError (errorList, ERR_OUT_OF_MEMORY);
			  return NULL;
			}
		      
		      newRow->choice = ALI_OTHER;
		      newRow->data.ptrvalue = otherLine.other;
		    }

		  /* If the next row needs to be a Sequence, */
		  /* and we're not still on the same row,    */
		  /* then change the previous ID to other.   */
		  
		  if (nextRowMustBeSeq && !idFound)
		    {
		      nextRowMustBeSeq = FALSE;
		      lastId->choice = ALI_OTHER;
		    }
		  
		  if (idFound)
		    nextRowMustBeSeq = TRUE;
		  else
		    nextRowMustBeSeq = FALSE;
		}
	    }
	}
    }

  return rowList;
}
