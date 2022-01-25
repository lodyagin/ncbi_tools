/*=========================================================================*/
/*                                                                         */
/*  aliparse.c                                                             */
/*                                                                         */
/*=========================================================================*/

/* Include files */

#include <ncbi.h>
#include <objalign.h>

#include "aliparse.h"
#include "aliread.h"

/* Defined constants */

#define ALI_USE_MAYBES         FALSE

#define MAX_SEGMENTS           1000
#define NO_SEGMENT             0
#define NUCL_SEGMENT           1
#define GAP_SEGMENT            2

#define ALI_SHOW_SEQUENCES     0x01
#define ALI_SHOW_DEFLINES      0x02
#define ALI_SHOW_IDS           0x04
#define ALI_SHOW_OTHERS        0x08
#define ALI_SHOW_MAYBES        0x10
#define ALI_SHOW_ALL           0xff

/* Function prototypes */

static void             s_FreeErrorList (ErrInfoPtr errorList);
static void             s_FreeSequenceList (SeqPartPtr seqPtr);
static void             s_FreeIdList (IdInfoPtr idList);
static void             s_FreeRowList (ValNodePtr rowList);
static void             s_FreeRowList_Safe (ValNodePtr rowList);
static void             s_DisplayRowList (ValNodePtr rowList,
					  Int2 mask);
static void             s_DisplayAnalyzeError (ValNodePtr errEntry,
					       Int2 error);
static IdInfoPtr        s_ProcessMaybes (ValNodePtr rowList);
static int              s_SegCompare(const void *i,
				     const void *j);
static Boolean          s_IsInterleaved (ValNodePtr rowList,
					 Int2 PNTR idCount);
static AlignFileDataPtr s_AnalyzeInterleaved (ValNodePtr rowList,
					      Int2       idCount);
static AlignFileDataPtr s_AnalyzeContiguous (ValNodePtr rowList);
static AlignFileDataPtr s_AnalyzeContents (ValNodePtr rowList);


Boolean    SegmentMultSequences (SeqAlignPtr sap,
				 SequenceListPtr seqList);


/*=========================================================================*/
/*                                                                         */
/*  s_FreeErrorList () - Free a linked list of error structures and all    */
/*                       the memory that they point to.                    */
/*                                                                         */
/*=========================================================================*/

static void s_FreeErrorList (ErrInfoPtr errorPtr)
{
  ErrInfoPtr currentErr;

  while (errorPtr != NULL)
    {
      MemFree (errorPtr->id);
      currentErr = errorPtr;
      errorPtr = errorPtr->next;
      MemFree (currentErr);
    }
}

/*=========================================================================*/
/*                                                                         */
/*  s_FreeSequenceList () - Free a linked list of SeqPart structures and   */
/*                          all the memory that they point to.             */
/*                                                                         */
/*=========================================================================*/

static void s_FreeSequenceList (SeqPartPtr seqPtr)
{
  SeqPartPtr currentSeq;

  while (seqPtr != NULL)
    {
      MemFree (seqPtr->sequence);
      currentSeq = seqPtr;
      seqPtr = seqPtr->next;
      MemFree (currentSeq);
    }
}

/*=========================================================================*/
/*                                                                         */
/*  s_FreeIdList () - Free a linked list of ID structures and all the      */
/*                    memory that they point to.                           */
/*                                                                         */
/*=========================================================================*/

static void s_FreeIdList (IdInfoPtr idPtr)
{
  IdInfoPtr currentId;

  while (idPtr != NULL)
    {
      MemFree (idPtr->id);
      s_FreeSequenceList (idPtr->sequence);
      MemFree (idPtr->defline);
      currentId = idPtr;
      idPtr = idPtr->next;
      MemFree (currentId);
    }
}

/*=========================================================================*/
/*                                                                         */
/*  Ali_FreeFileInfo () - Free a AlignFileData structure and all the       */
/*                        memory that it points to.                        */
/*                                                                         */
/*=========================================================================*/

void Ali_FreeFileInfo (AlignFileDataPtr fileInfoPtr)
{

  s_FreeIdList (fileInfoPtr->sequences);
  s_FreeIdList (fileInfoPtr->maybes);
  s_FreeErrorList (fileInfoPtr->errors);

  MemFree (fileInfoPtr);

  return;
}

/*=========================================================================*/
/*                                                                         */
/*  s_FreeRowList () - Free all row data structures and the strings that   */
/*                     they point to.                                      */
/*                                                                         */
/*         NOTE: The actual data strings in the row list may be pointed    */
/*               to by other structures, in which case                     */
/*               s_FreeRowList_Safe () should be used instead.             */
/*                                                                         */
/*                                                                         */
/*=========================================================================*/

static void s_FreeRowList (ValNodePtr rowList)
{
  CharPtr    deleteStr;
  ValNodePtr currentRow;

  while (rowList != NULL)
    {
      deleteStr = (CharPtr) rowList->data.ptrvalue;
      if (deleteStr != NULL)
	MemFree (deleteStr);
      currentRow = rowList;
      rowList = rowList->next;
      MemFree (currentRow);
    }
}

/*=========================================================================*/
/*                                                                         */
/*  s_FreeRowList_Safe () - Free all row data structures, but don't free   */
/*                          the strings that they point to unless we're    */
/*                          sure that they aren't pointed to elsewhere.    */
/*                          i.e. -- Only free OTHER lines and UNKNOWN      */
/*                          lines.                                         */
/*                                                                         */
/*=========================================================================*/

static void s_FreeRowList_Safe (ValNodePtr rowList)
{
  CharPtr    deleteStr;
  ValNodePtr currentRow;

  while (rowList != NULL)
    {
      if ((rowList->choice == ALI_OTHER) ||
	  (rowList->choice == ALI_UNKNOWN))
	{
	  deleteStr = (CharPtr) rowList->data.ptrvalue;
	  if (deleteStr != NULL)
	    MemFree (deleteStr);
	}
      currentRow = rowList;
      rowList = rowList->next;
      MemFree (currentRow);
    }
}

/*=========================================================================*/
/*                                                                         */
/*  s_ProcessMaybes ()                                                     */
/*                                                                         */
/*=========================================================================*/

static IdInfoPtr s_ProcessMaybes (ValNodePtr rowList)
{
  ValNodePtr    currentEntry;
  IdInfoPtr     badIdList = NULL;
  IdInfoPtr     existingId = NULL;
  IdInfoPtr     currentId = NULL;
  IdInfoPtr     lastId = NULL;
  CharPtr       currentIdStr;
  SeqPartPtr    newSeqPart;
  SeqPartPtr    lastSeqPart;

  currentEntry = rowList;

  while (currentEntry != NULL)
    {
      if (((currentEntry->choice == ALI_SEQLINE_ID) ||
	   (currentEntry->choice == ALI_DEFLINE_ID) ||
	   (currentEntry->choice == ALI_OWNLINE_ID)) ||
	   (currentEntry->choice == ALI_MAYBE_ID))
	{
	  currentIdStr = (CharPtr) currentEntry->data.ptrvalue;
	}
      else if (currentEntry->choice == ALI_MAYBE_SEQUENCE)
	{
	  existingId = badIdList;
	  while (existingId != NULL)
	    {
	      if (StringCmp(existingId->id,currentIdStr) == 0)
		break;
	      existingId = existingId->next;
	    }

	  if (existingId != NULL)
	    currentId = existingId;
	  else
	    {
	      currentId = (IdInfoPtr) MemNew (sizeof(IdInfo));
	      if (currentId == NULL)
		return NULL;
	      
	      currentId->sequence = NULL;
	      currentId->id       = currentIdStr;
	      currentId->length   = 0;
	      currentId->next     = NULL;
	      
	      if (badIdList == NULL)
		badIdList = currentId;
	      else
		{
		  lastId = badIdList;
		  while (lastId->next != NULL)
		    lastId = lastId->next;
		  lastId->next = currentId;
		}
	    }

	  /* Add the sequence to the current ID */

	  newSeqPart = (SeqPartPtr) MemNew(sizeof(SeqPartPtr));
	  if (newSeqPart == NULL)
	    return NULL;

	  newSeqPart->sequence = (CharPtr) currentEntry->data.ptrvalue;
	  newSeqPart->next     = NULL;

	  if (currentId->sequence == NULL)
	    currentId->sequence = newSeqPart;
	  else
	    lastSeqPart->next = newSeqPart;

	  currentId->length += StringLen (newSeqPart->sequence);
	  lastSeqPart = newSeqPart;

	}
      currentEntry = currentEntry->next;
    }

  return badIdList;
}


/*=========================================================================*/
/*                                                                         */
/*  s_DisplayAnalyzeError ()                                               */
/*                                                                         */
/*=========================================================================*/

static void s_DisplayAnalyzeError (ValNodePtr errEntry, Int2 error)
{
  if (errEntry != NULL)
    {
      if ((errEntry->choice == ALI_SEQUENCE) ||
	  (errEntry->choice == ALI_MAYBE_SEQUENCE))
	fprintf(stderr, "\nInvalid sequence : %s\n",
		(CharPtr) errEntry->data.ptrvalue);
      else
	fprintf(stderr, "\nInvalid ID : %s\n",
		(CharPtr) errEntry->data.ptrvalue);
    }

  switch (error)
    {
    case ERR_SEQ_BEFORE_ID :
      fprintf (stderr, "ERROR: The first sequence appears before the first ID.\n");
      break;
    case ERR_ID_WITHOUT_SEQ :
      fprintf (stderr, "ERROR: An ID occurs that cannot be matched to a sequence.\n");
      break;
    case ERR_SEQ_WITHOUT_ID :
      fprintf (stderr, "ERROR: A sequence occurs that cannot be matched to an ID.\n");
      break;
    case ERR_DUPLICATE_IDS :
      fprintf (stderr, "ERROR: Duplicate IDs.\n");
      break;
    case ERR_SEQUENCE_TOO_SHORT :
      fprintf (stderr, "ERROR: The sequence is too short.\n");
      break;
    case ERR_SEQUENCE_TOO_LONG :
      fprintf (stderr, "ERROR: The sequence is too long.\n");
      break;
    default:
      fprintf (stderr, "ERROR: Unknown error.\n");
      break;
    }
}

/*=========================================================================*/
/*                                                                         */
/*  DisplayRowList() - Prints to stderr the linked list of ValNodes that   */
/*                     contain the data read in from the alignment file.   */
/*                                                                         */
/*     NOTE: This function is for debugging purposes only and should       */
/*           not be called in production code.                             */
/*                                                                         */
/*=========================================================================*/

static void s_DisplayRowList (ValNodePtr rowList,
			      Int2       mask)
{
  ValNodePtr currRow;

  currRow = rowList;
  while (currRow != NULL)
    {
      if ((currRow->choice == ALI_SEQUENCE) &&
	  ((mask & ALI_SHOW_SEQUENCES) ||
	   (mask == ALI_SHOW_ALL)))
	fprintf(stderr,"SEQUENCE : %s\n", currRow->data.ptrvalue);
      else if ((currRow->choice == ALI_DEFLINE) &&
	       ((mask & ALI_SHOW_DEFLINES) ||
		(mask == ALI_SHOW_ALL)))
	fprintf(stderr,"DEFLINE  : %s\n", currRow->data.ptrvalue);
      else if ((currRow->choice == ALI_SEQLINE_ID) &&
	       ((mask & ALI_SHOW_IDS) ||
		(mask == ALI_SHOW_ALL)))
	fprintf(stderr,"INLINE ID: %s\n", currRow->data.ptrvalue);
      else if ((currRow->choice == ALI_DEFLINE_ID) &&
	       ((mask & ALI_SHOW_IDS) ||
		(mask == ALI_SHOW_ALL)))
	fprintf(stderr,"DEFLINE ID: %s\n", currRow->data.ptrvalue);
      else if ((currRow->choice == ALI_OWNLINE_ID) &&
	       ((mask & ALI_SHOW_IDS) ||
		(mask == ALI_SHOW_ALL)))
	fprintf(stderr,"OWNLINE ID: %s\n", currRow->data.ptrvalue);
      else if ((currRow->choice == ALI_MAYBE_ID) &&
	       ((mask & ALI_SHOW_MAYBES) ||
		(mask == ALI_SHOW_ALL)))
	fprintf(stderr,"MAYBE ID: %s\n", currRow->data.ptrvalue);
      else if ((currRow->choice == ALI_MAYBE_SEQUENCE) &&
	       ((mask & ALI_SHOW_MAYBES) ||
		(mask == ALI_SHOW_ALL)))
	fprintf(stderr,"MAYBE SEQUENCE: %s\n", currRow->data.ptrvalue);
      else if ((currRow->choice == ALI_OTHER) &&
	       ((mask & ALI_SHOW_OTHERS) ||
		(mask == ALI_SHOW_ALL)))
	fprintf(stderr,"OTHER    : %s\n", currRow->data.ptrvalue);

      currRow = currRow->next;
    }

  return;
}


/*=========================================================================*/
/*                                                                         */
/* s_isInterleaved ()                                                      */
/*                                                                         */
/*=========================================================================*/

static Boolean s_IsInterleaved (ValNodePtr rowList,
				Int2 PNTR idCount)
{
  ValNodePtr    currentEntry;
  CharPtr       newIdStr;
  IdInfoPtr     idList = NULL;
  IdInfoPtr     lastId = NULL;
  IdInfoPtr     currentId = NULL;
  IdInfoPtr     existingId = NULL;
  Boolean       isInterleaved;
  Int4          patternRowCount;
  Int4          patternCharCount;
  Int4          currentRowCount;
  Int4          currentCharCount;
  Boolean       isFirstId;
  CharPtr       seqStr;

  isInterleaved = FALSE;
  currentEntry  = rowList;

  patternRowCount  = 0;
  patternCharCount = 0;
  currentRowCount  = 0;
  currentCharCount = 0;
  isFirstId        = TRUE;
  *idCount         = 0;

  /* Search the row list for IDs */

  while (currentEntry != NULL)
    {
      if (((currentEntry->choice == ALI_SEQLINE_ID) ||
	   (currentEntry->choice == ALI_DEFLINE_ID) ||
	   (currentEntry->choice == ALI_OWNLINE_ID)) ||
	  ((ALI_USE_MAYBES == TRUE) &&
	   (currentEntry->choice == ALI_MAYBE_ID)))
	{

	  /* Check to see if we've already found */
	  /* this particular ID.                 */

	  newIdStr = (CharPtr) currentEntry->data.ptrvalue;
	  existingId = idList;
	  while (existingId != NULL)
	    {
	      if (StringCmp(existingId->id,newIdStr) == 0)
		break;
	      existingId = existingId->next;
	    }

	  /* If we have, break and return TRUE */

	  if (existingId != NULL)
	    {
	      isInterleaved = TRUE;
	      break;
	    }

	  /* Otherwise, add the ID to the list */

	  if (idList != NULL)
	    isFirstId = FALSE;

	  (*idCount)++;
	  
	  currentId = (IdInfoPtr) MemNew (sizeof(IdInfo));
	  if (currentId == NULL)
	    return FALSE;
	  
	  currentId->sequence = NULL;
	  currentId->id       = newIdStr;
	  currentId->length   = 0;
	  currentId->next     = NULL;
	  
	  if (idList == NULL)
	    idList = currentId;
	  else
	    {
	      lastId = idList;
	      while (lastId->next != NULL)
		lastId = lastId->next;
	      lastId->next = currentId;
	    }
	}

      /* Process the row as a sequence */

      else if ((currentEntry->choice == ALI_SEQUENCE) ||
	       ((ALI_USE_MAYBES == TRUE) &&
		(currentEntry->choice == ALI_MAYBE_SEQUENCE)))
	{

	  /* There must be an ID before the first sequence */

	  if (currentId == NULL)
	    {
	      isInterleaved = FALSE;
	      break;
	    }

	  /* Look for sequences that probably */
	  /* have no ID assigned to them.     */

	  seqStr = (CharPtr) currentEntry->data.ptrvalue;
	  if (isFirstId)
	    {
	      patternRowCount++;
	      patternCharCount += StringLen (seqStr);
	    }
	  else
	    {
	      currentRowCount++;
	      currentCharCount++;
	      if ((currentRowCount > patternRowCount) &&
		  (currentCharCount > patternCharCount))
		{
		  isInterleaved = TRUE;
		  break;
		}
	    }

	}

      /* Go to next row */

      currentEntry = currentEntry->next;
    }

  /* Delete the ID records that we created */
  /*  NOTE -- The ID strings themselves    */
  /*          are stored elsewhere and     */
  /*          only pointed to here, so     */
  /*          DON"T delete them.           */

  while (idList != NULL)
    {
      lastId = idList;
      idList = idList->next;
      MemFree(lastId);
    }

  /* Return result of search */

  return isInterleaved;
}

/*=========================================================================*/
/*                                                                         */
/* s_AnalyzeInterleaved ()                                                 */
/*                                                                         */
/*=========================================================================*/

static AlignFileDataPtr s_AnalyzeInterleaved (ValNodePtr rowList,
					      Int2       idCount)
{
  ValNodePtr       currentEntry;
  ValNodePtr       lastEntry = NULL;
  Boolean          isFirstGroup;
  Boolean          isValidPattern;
  IdInfoPtr        currentId = NULL;
  IdInfoPtr        lastId = NULL;
  IdInfoPtr        existingId = NULL;
  IdInfoPtr        idList = NULL;
  IdInfoPtr        maybeIdList = NULL;
  SeqPartPtr       newSeqPart = NULL;
  SeqPartPtr       lastSeqPart = NULL;
  Boolean          isFirstId;
  Int4             previousLength;
  CharPtr          newIdStr;
  Int2             error;
  Boolean          maybesFound;
  Boolean          gotAllIds;
  Int2             currentIdCount;
  AlignFileDataPtr fileInfoPtr;

  /* Match the sequences up with the IDs */

  currentEntry    = rowList;
  isFirstId       = TRUE;
  isValidPattern  = TRUE;
  isFirstGroup    = TRUE;
  maybesFound     = FALSE;
  gotAllIds       = FALSE;
  currentIdCount  = 0;

  while (currentEntry != NULL)
    {

      /* Process the row as an ID */

      if (((currentEntry->choice == ALI_SEQLINE_ID) ||
	   (currentEntry->choice == ALI_DEFLINE_ID) ||
	   (currentEntry->choice == ALI_OWNLINE_ID)) ||
	  ((ALI_USE_MAYBES == TRUE) &&
	   (currentEntry->choice == ALI_MAYBE_ID)))
	{

	  /* If we've already got all our IDs */
	  /* then ignore any more.            */

	  if (gotAllIds == TRUE)
	    {
	      currentEntry = currentEntry->next;
	      continue;
	    }

	  /* All ID's, except for the first one, should */
	  /* immediately follow a sequence line.        */

	  if (isFirstId)
	    {
	      isFirstId = FALSE;
	      if (lastEntry != NULL)
		{
		  isValidPattern = FALSE;
		  error = ERR_SEQ_BEFORE_ID;
		  break;
		}
	    }
	  else
	    {
	      if (ALI_USE_MAYBES == FALSE)
		{
		  if (lastEntry->choice != ALI_SEQUENCE)
		    {
		      isValidPattern = FALSE;
		      error = ERR_ID_WITHOUT_SEQ;
		      break;
		    }
		  else
		    isFirstGroup = FALSE;
		}
	      else /* ALI_USE_MAYBES == TRUE */
		{
		  if ((lastEntry->choice != ALI_SEQUENCE) &&
		      (lastEntry->choice != ALI_MAYBE_SEQUENCE))
		    {
		      isValidPattern = FALSE;
		      error = ERR_ID_WITHOUT_SEQ;
		      break;
		    }
		  else
		    isFirstGroup = FALSE;
		}
	    }

	  /* If this id already exists, */
	  /* make it the current ID.    */

	  newIdStr = (CharPtr) currentEntry->data.ptrvalue;
	  existingId = idList;
	  while (existingId != NULL)
	    {
	      if (StringCmp(existingId->id,newIdStr) == 0)
		break;
	      existingId = existingId->next;
	    }

	  if (existingId != NULL)
	    currentId = existingId;

	  /* Otherwise create a new Id record */
	  /* and add it to the end of list.   */
	  
	  else
	    {
	      currentId = (IdInfoPtr) MemNew (sizeof(IdInfo));
	      if (currentId == NULL)
		return NULL;
	      
	      currentId->sequence = NULL;
	      currentId->id       = newIdStr;
	      currentId->length   = 0;
	      currentId->next     = NULL;
	      
	      if (idList == NULL)
		idList = currentId;
	      else
		{
		  lastId = idList;
		  while (lastId->next != NULL)
		    lastId = lastId->next;
		  lastId->next = currentId;
		}

	      currentIdCount++;
	      if (currentIdCount == idCount)
		gotAllIds = TRUE;
	    }

	  /* Mark this as the last row processed */

	  lastEntry = currentEntry;
	}

      /* Process the row as a sequence */

      else if ((currentEntry->choice == ALI_SEQUENCE) ||
	       ((ALI_USE_MAYBES == TRUE) &&
		(currentEntry->choice == ALI_MAYBE_SEQUENCE)))
	{

	  /* There must be an ID before the first sequence */

	  if (currentId == NULL)
	    {
	      isValidPattern = FALSE;
	      error = ERR_SEQ_WITHOUT_ID;
	      break;
	    }

	  /* Add the sequence to the current ID */

	  newSeqPart = (SeqPartPtr) MemNew(sizeof(SeqPartPtr));
	  if (newSeqPart == NULL)
	    return NULL;

	  newSeqPart->sequence = (CharPtr) currentEntry->data.ptrvalue;
	  newSeqPart->next     = NULL;

	  if (currentId->sequence == NULL)
	    currentId->sequence = newSeqPart;
	  else
	    lastSeqPart->next = newSeqPart;

	  currentId->length += StringLen (newSeqPart->sequence);
	  lastSeqPart = newSeqPart;

	  /* If we've started repeating IDs then */
	  /* rotate through the id list.         */

	  if (gotAllIds == TRUE)
	    {
	      if (currentId->next == NULL)
		currentId = idList;
	      else
		currentId = currentId->next;
	    }

	  /* Mark this as the last row processed */

	  lastEntry = currentEntry;
	}
      else if ((ALI_USE_MAYBES == FALSE) &&
	    (currentEntry->choice == ALI_MAYBE_SEQUENCE))
	maybesFound = TRUE;

      currentEntry = currentEntry->next;
    }

  /* Make sure that we don't have a dangling ID */

  if ((lastEntry == NULL) || (lastEntry->choice != ALI_SEQUENCE))
    isValidPattern = FALSE;

  /* If pattern still not found, return failure */

  if (!isValidPattern)
    {
      s_DisplayAnalyzeError(currentEntry, error);
      return NULL;
    }

  /* Sequences should all be the same length. */

  currentId = idList;
  isFirstId = TRUE;

  while (currentId != NULL)
    {
      if (isFirstId)
	isFirstId = FALSE;
      else
	if (previousLength != currentId->length)
	  break;
      previousLength = currentId->length;
      currentId = currentId->next;
    }

  if (currentId != NULL)
    fprintf (stderr, "\nERROR: Sequences are not of matched length.\n");

  /* Process the maybes if they weren't used already */

  if ((ALI_USE_MAYBES == FALSE) && (maybesFound == TRUE))
    maybeIdList = s_ProcessMaybes (rowList);

  /* Return successfully */

  if (currentId == NULL)
    {
      fileInfoPtr = (AlignFileDataPtr) MemNew (sizeof(AlignFileDataPtr));
      fileInfoPtr->sequences = idList;
      fileInfoPtr->maybes    = maybeIdList;
      fileInfoPtr->errors    = NULL;
      return fileInfoPtr;
    }
  else
    return NULL;
} 

/*=========================================================================*/
/*                                                                         */
/* s_AnalyzeContiguous ()                                                  */
/*                                                                         */
/*=========================================================================*/

static AlignFileDataPtr s_AnalyzeContiguous (ValNodePtr rowList)
{
  ValNodePtr       currentEntry;
  ValNodePtr       lastEntry = NULL;
  Boolean          isFirstGroup;
  Int4             curPatternCount;
  Int4             totPatternCount;
  Boolean          isValidPattern;
  IdInfoPtr        currentId = NULL;
  IdInfoPtr        lastId = NULL;
  IdInfoPtr        nextToLastId = NULL;
  IdInfoPtr        existingId = NULL;
  IdInfoPtr        idList = NULL;
  IdInfoPtr        maybeIdList = NULL;
  SeqPartPtr       newSeqPart = NULL;
  SeqPartPtr       lastSeqPart = NULL;
  Boolean          isFirstId;
  CharPtr          newIdStr;
  Int2             error;
  Boolean          maybesFound;
  AlignFileDataPtr fileInfoPtr;

  /* Match the sequences up with the IDS */

  currentEntry    = rowList;
  isFirstId       = TRUE;
  isValidPattern  = TRUE;
  isFirstGroup    = TRUE;
  maybesFound     = FALSE;
  curPatternCount = 0;
  totPatternCount = 0;

  while (currentEntry != NULL)
    {

      /* Process the row as an ID */

      if (((currentEntry->choice == ALI_SEQLINE_ID) ||
	   (currentEntry->choice == ALI_DEFLINE_ID) ||
	   (currentEntry->choice == ALI_OWNLINE_ID)) ||
	  ((ALI_USE_MAYBES == TRUE) &&
	   (currentEntry->choice == ALI_MAYBE_ID)))
	{

	  /* All ID's, except for the first one, should */
	  /* immediately follow a sequence line.        */

	  if (isFirstId)
	    {
	      isFirstId = FALSE;
	      if (lastEntry != NULL)
		{
		  isValidPattern = FALSE;
		  error = ERR_SEQ_BEFORE_ID;
		  break;
		}
	    }
	  else
	    {
	      if (ALI_USE_MAYBES == FALSE)
		{
		  if (lastEntry->choice != ALI_SEQUENCE)
		    {
		      isValidPattern = FALSE;
		      error = ERR_ID_WITHOUT_SEQ;
		      break;
		    }
		  else
		    isFirstGroup = FALSE;
		}
	      else /* ALI_USE_MAYBES == TRUE */
		{
		  if ((lastEntry->choice != ALI_SEQUENCE) &&
		      (lastEntry->choice != ALI_MAYBE_SEQUENCE))
		    {
		      isValidPattern = FALSE;
		      error = ERR_ID_WITHOUT_SEQ;
		      break;
		    }
		  else
		    isFirstGroup = FALSE;
		}
	    }

	  /* The length of the last pattern must match */
	  /* the length of previous ones.              */

	  if (curPatternCount < totPatternCount)
	    {
	      isValidPattern = FALSE;
	      error = ERR_SEQUENCE_TOO_SHORT;
	      break;
	    }

	  curPatternCount = 0;

	  /* If this id already exists, */
	  /* make it the current ID.    */

	  newIdStr = (CharPtr) currentEntry->data.ptrvalue;
	  existingId = idList;
	  while (existingId != NULL)
	    {
	      if (StringCmp(existingId->id,newIdStr) == 0)
		{
		  error = ERR_DUPLICATE_IDS;
		  break;
		}
	      existingId = existingId->next;
	    }

	  if (existingId != NULL)
	    currentId = existingId;

	  /* Otherwise create a new Id record */
	  /* and add it to the end of list.   */
	  
	  else
	    {
	      currentId = (IdInfoPtr) MemNew (sizeof(IdInfo));
	      if (currentId == NULL)
		return NULL;
	      
	      currentId->sequence = NULL;
	      currentId->id       = newIdStr;
	      currentId->length   = 0;
	      currentId->next     = NULL;
	      
	      if (idList == NULL)
		idList = currentId;
	      else
		{
		  lastId = idList;
		  while (lastId->next != NULL)
		    lastId = lastId->next;
		  lastId->next = currentId;
		}
	    }

	  /* Mark this as the last row processed */

	  lastEntry = currentEntry;
	}

      /* Process the row as a sequence */

      else if ((currentEntry->choice == ALI_SEQUENCE) ||
	       ((ALI_USE_MAYBES == TRUE) &&
		(currentEntry->choice == ALI_MAYBE_SEQUENCE)))
	{

	  /* There must be an ID before we get a sequence */

	  if (currentId == NULL)
	    {
	      isValidPattern = FALSE;
	      error = ERR_SEQ_WITHOUT_ID;
	      break;
	    }

	  /* Add the sequence to the current ID */

	  newSeqPart = (SeqPartPtr) MemNew(sizeof(SeqPartPtr));
	  if (newSeqPart == NULL)
	    return NULL;

	  newSeqPart->sequence = (CharPtr) currentEntry->data.ptrvalue;
	  newSeqPart->next     = NULL;

	  if (currentId->sequence == NULL)
	    currentId->sequence = newSeqPart;
	  else
	    lastSeqPart->next = newSeqPart;

	  currentId->length += StringLen (newSeqPart->sequence);
	  lastSeqPart = newSeqPart;

	  if (isFirstGroup)
	    {
	      totPatternCount += StringLen (newSeqPart->sequence);
	      curPatternCount += StringLen (newSeqPart->sequence);
	    }
	  else
	    {
	      curPatternCount += StringLen (newSeqPart->sequence);

	      if (curPatternCount > totPatternCount)
		{
		  isValidPattern = FALSE;
		  error = ERR_SEQUENCE_TOO_LONG;
		  break;
		}
	    }

	  /* Mark this as the last row processed */

	  lastEntry = currentEntry;
	}
      else if ((ALI_USE_MAYBES == FALSE) &&
	    (currentEntry->choice == ALI_MAYBE_SEQUENCE))
	maybesFound = TRUE;

      currentEntry = currentEntry->next;
    }

  /* Make sure that we don't have a dangling ID */

  if ((lastEntry == NULL) || (lastEntry->choice != ALI_SEQUENCE))
    {
      error = ERR_ID_WITHOUT_SEQ;
      isValidPattern = FALSE;
    }

  /* If the last sequence is too short, mark */
  /* it as a maybe.                          */
  
  else if (lastEntry->choice == ALI_SEQUENCE)
    {
      if (ALI_USE_MAYBES == FALSE)
	{
	  maybesFound = TRUE;
	  if (curPatternCount < totPatternCount)
	    {
	      lastEntry->choice = ALI_MAYBE_SEQUENCE;
	      nextToLastId = NULL;
	      lastId = idList;
	      while (lastId->next != NULL)
		{
		  nextToLastId = lastId;
		  lastId = lastId->next;
		}
	      MemFree(lastId);
	      if (nextToLastId == NULL)
		idList = NULL;
	      else
		nextToLastId->next = NULL;
	    }
	}
      else
	{
	  error = ERR_SEQUENCE_TOO_SHORT;
	  isValidPattern = FALSE;
	}
    }

  /* If pattern not found, return failure */

  if (!isValidPattern)
    {
      s_DisplayAnalyzeError(currentEntry, error);
      return NULL;
    }

  /* If we have some possibly bad sequences that */
  /* weren't used, process them seperately.      */

  if ((ALI_USE_MAYBES == FALSE) && (maybesFound == TRUE))
    maybeIdList = s_ProcessMaybes (rowList);

  /* Return successfully */

  if (currentId != NULL)
    return NULL;
  else
    {
      fileInfoPtr = (AlignFileDataPtr) MemNew (sizeof(AlignFileDataPtr));
      fileInfoPtr->sequences = idList;
      fileInfoPtr->maybes    = maybeIdList;
      fileInfoPtr->errors    = NULL;
      return fileInfoPtr;
    }
}

/*=========================================================================*/
/*                                                                         */
/* Ali_AddError ()                                                         */
/*                                                                         */
/*=========================================================================*/

ErrInfoPtr Ali_AddError (ErrInfoPtr PNTR errorList,
			 Int4            iError)
{
  return NULL;
}

/*=========================================================================*/
/*                                                                         */
/* s_AnalyzeContents () -                                                  */
/*                                                                         */
/*=========================================================================*/

static AlignFileDataPtr s_AnalyzeContents (ValNodePtr rowList)
{
  AlignFileDataPtr  fileInfoPtr;
  Int2              idCount;

  if (s_IsInterleaved (rowList, &idCount))
    fileInfoPtr = s_AnalyzeInterleaved (rowList, idCount);
  else 
    fileInfoPtr = s_AnalyzeContiguous (rowList);

  return fileInfoPtr;
}

/*=========================================================================*/
/*                                                                         */
/* Ali_ReadByContent ()                                                    */
/*                                                                         */
/*=========================================================================*/

AlignFileDataPtr Ali_ReadByContent (FILE PNTR alignFilePtr)
{
  ValNodePtr        rowList = NULL;
  AlignFileDataPtr  fileInfoPtr;
  IdInfoPtr         idList = NULL;
  IdInfoPtr         maybeIdList = NULL;
  ErrInfoPtr        errorList = NULL;

  /* Check parameters */

  if (alignFilePtr == NULL)
    return FALSE;

  /* Read in and parse each row until */
  /* we reach the end of file.        */

  rowList = Ali_ReadLines (alignFilePtr, &errorList);
  if (rowList == NULL)
    return FALSE;

  /*
  s_DisplayRowList (rowList, ALI_SHOW_ALL);
  */

  /* Analyze the IDs and sequences for consistancy */

  fileInfoPtr = s_AnalyzeContents (rowList);

  /* Clean up and return successfully */

  if (fileInfoPtr == NULL)
    s_FreeRowList (rowList);
  else
    s_FreeRowList_Safe(rowList);

  return fileInfoPtr;
}

/*=========================================================================*/
/*                                                                         */
/*  s_SegCompare () -- Comparison function for qsort in the function       */
/*                     SegmentMultSequences().                             */
/*                                                                         */
/*=========================================================================*/

static int s_SegCompare(const void *i, const void *j)
{
  Int4 PNTR first;
  Int4 PNTR second;

  first  = (Int4 *) i;
  second = (Int4 *) j;

  if (*first > *second)
    return (1);
  if (*first < *second)
    return (-1);

  return (0);
}

/*=========================================================================*/
/*                                                                         */
/*  SegmentMultSequences ()                                                */
/*                                                                         */
/*=========================================================================*/

Boolean SegmentMultSequences (SeqAlignPtr sap,
			      SequenceListPtr seqList)
{
  Int2            sequenceCount = 0;
  Int4            segmentPoints[MAX_SEGMENTS];
  Int4            i;
  Int4            sequenceLen;
  SequenceListPtr currentSeq = NULL;
  Int4            segmentCount = 0;
  Int2            currSegType;
  Int2            currCharType;

  /* Check parameters */

  if (sap == NULL)
    return FALSE;

  if (seqList == NULL)
    return FALSE;

  /* Determine where each sequence would */
  /* individually segment.               */

  currentSeq = seqList;
  while (currentSeq != NULL)
    {
      if (currentSeq->data == NULL)
	return FALSE;
      sequenceLen = StringLen(currentSeq->data);
      currSegType = NO_SEGMENT;
      for (i = 0; i < sequenceLen; i++)
	{
	  if ((currentSeq->data[i] == 'A') ||
	      (currentSeq->data[i] == 'T') ||
	      (currentSeq->data[i] == 'C') ||
	      (currentSeq->data[i] == 'G') ||
	      (currentSeq->data[i] == 'a') ||
	      (currentSeq->data[i] == 't') ||
	      (currentSeq->data[i] == 'c') ||
	      (currentSeq->data[i] == 'g'))
	    currCharType = NUCL_SEGMENT;
	  else
	    currCharType = GAP_SEGMENT;

	  if (currCharType != currSegType)
	    {
	      segmentPoints[segmentCount] = i;
	      currSegType = currCharType;
	      segmentCount++;
	    }
	}

      currentSeq = currentSeq->next;
      sequenceCount++;
    }

  /* Sort the list of segment breaks to */
  /* come up with an overall list.      */

  qsort(segmentPoints,  segmentCount,  sizeof(Int4), s_SegCompare);

  /* Break each sequence up into segments and */
  /* store those segments in a SeqAlignPtr.   */

  currentSeq = seqList;
  while (currentSeq != NULL)
    {
      currentSeq = currentSeq->next;
    }

  /* Return successfully */

  return TRUE;
}
