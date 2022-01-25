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

File name: impalatool.c

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
		add_string_to_bufferEx("Reference: Alejandro A. Schaffer, Yuri I. Wolf, Eugene V. Koonin, ", &ret_buffer, &ret_buffer_length, TRUE);
	add_string_to_bufferEx("L. Aravind, Stephen F. Altschul (1999), ", &ret_buffer, &ret_buffer_length, TRUE);
	add_string_to_bufferEx("\"IMPALA: Integrating Matrix Profiles and Local Alignments\",", &ret_buffer, &ret_buffer_length, TRUE);
	add_string_to_bufferEx("submitted.", &ret_buffer, &ret_buffer_length, TRUE);
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
