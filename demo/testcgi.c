/*
* ===========================================================================
*
*                            PUBLIC DOMAIN NOTICE
*            National Center for Biotechnology Information (NCBI)
*
*  This software/database is a "United States Government Work" under the
*  terms of the United States Copyright Act.  It was written as part of
*  the author's official duties as a United States Government employee and
*  thus cannot be copyrighted.  This software/database is freely available
*  to the public for use. The National Library of Medicine and the U.S.
*  Government do not place any restriction on its use or reproduction.
*  We would, however, appreciate having the NCBI and the author cited in
*  any work or product based on this material
*
*  Although all reasonable efforts have been taken to ensure the accuracy
*  and reliability of the software and data, the NLM and the U.S.
*  Government do not and cannot warrant the performance or results that
*  may be obtained by using this software or data. The NLM and the U.S.
*  Government disclaim all warranties, express or implied, including
*  warranties of performance, merchantability or fitness for any particular
*  purpose.
*
* ===========================================================================
*
* File Name:  testcgi.c
*
* Author:  Jonathan Kans
*
* Version Creation Date:   4/24/98
*
* $Revision: 6.17 $
*
* File Description: 
*
* Modifications:  
* --------------------------------------------------------------------------
* Date     Name        Description of modification
* -------  ----------  -----------------------------------------------------
*
*
* ==========================================================================
*/

/*
*  To compile on Solaris:
*
*   cc -xildoff -o testcgi.cgi testcgi.c -lgen -lm
*
*  To compile on SGI:
*
*   cc -mips1 -o testcgi.cgi testcgi.c -lm -lPW -lsun
*/

#include <stddef.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>

/* convenient defines and typedefs from NCBI toolkit */

#ifndef NULL
#define NULL ((void *)0)
#endif

#ifndef Pointer
typedef void * Pointer;
#endif

#ifndef Char
typedef char Char, * CharPtr;
#endif

#ifndef Bool
typedef unsigned char Bool, * BoolPtr;
#endif

#ifndef TRUE
#define TRUE ((Bool)1)
#endif

#ifndef FALSE
#define FALSE ((Bool)0)
#endif

#ifndef Int2
typedef short Int2, * Int2Ptr;
#endif

#ifndef Int4
typedef long Int4, * Int4Ptr;
#endif

/* useful portable character macros from NCBI toolkit (assumes ASCII) */

#define IS_DIGIT(c)	('0'<=(c) && (c)<='9')
#define IS_UPPER(c)	('A'<=(c) && (c)<='Z')
#define IS_LOWER(c)	('a'<=(c) && (c)<='z')
#define IS_ALPHA(c)	(IS_UPPER(c) || IS_LOWER(c))
#define TO_LOWER(c)	((Char)(IS_UPPER(c) ? (c)+' ' : (c)))
#define TO_UPPER(c)	((Char)(IS_LOWER(c) ? (c)-' ' : (c)))
#define IS_WHITESP(c) (((c) == ' ') || ((c) == '\n') || ((c) == '\r') || ((c) == '\t'))
#define IS_ALPHANUM(c) (IS_ALPHA(c) || IS_DIGIT(c))
#define IS_PRINT(c)	(' '<=(c) && (c)<='~')

/* copy of QUERY_STRING is parsed with strtok into tag and val paired arrays */

#define MAX_ENTRIES  32

static CharPtr  query = NULL;
static Int2     num_tags = 0;
static CharPtr  tag [MAX_ENTRIES];
static CharPtr  val [MAX_ENTRIES];

static Bool ParseQuery (Bool queryRequired)

{
  Char     ch;
  size_t   len;
  CharPtr  ptr;
  CharPtr  to;

/*
*  given a sample URL:
*
*   http://myserver:80/cgi-bin/mydigest.cgi?enzyme=EcoRI&pattern=GAATTC
*
*  the QUERY_STRING variable contains enzyme=EcoRI&pattern=GAATTC
*/

/* The >Message prefix causes Sequin to display the message to the user */

  ptr = getenv ("QUERY_STRING");
  if (ptr == NULL) {
    printf ("Content-type: text/html\n\n");
    printf (">Message\nFAILURE - Query not detected\n");
    fflush (stdout);
    return FALSE;
  }

/* allocates a copy of query string that can be modified during parsing */

  len = strlen (ptr);
  query = malloc (len + 3);
  memset (query, 0, len + 2);
  /* strcpy (query, ptr); */

/* may need to convert %20 encoding in URL to a space */

  to = query;
  ch = *ptr;
  while (ch != '\0') {
    if (ch == '%' && *(ptr + 1) == '2' && *(ptr + 2) == '0') {
      *to = ' ';
      to++;
      ptr += 3;
    } else {
      *to = ch;
      to++;
      ptr++;
    }
    ch = *ptr;
  }
  *to = '\0';

/* parse tag=value&tag=value query into arrays for easier interpretation */

  memset (tag, 0, sizeof (tag));
  memset (val, 0, sizeof (val));

  ptr = strtok (query, "=");
  for (num_tags = 0; num_tags < MAX_ENTRIES && ptr != NULL; num_tags++) {
    tag [num_tags] = ptr;
    ptr = strtok (NULL, "&");
    if (ptr != NULL) {
      val [num_tags] = ptr;
      ptr = strtok (NULL, "=");
    }
  }

/*
*  given the above query example, the tag and val arrays are now:
*
*   tag [0] = "enzyme"    val [0] = "EcoRI"
*   tag [1] = "pattern"   val [1] = "GAATTC"
*   tag [2] = NULL        val [2] = NULL
*
*  and num_tags is 2
*/

/* verify that any required query string is present in the URL */

  if (queryRequired && num_tags == 0) {
    printf ("Content-type: text/html\n\n");
    printf (">Message\nFAILURE - Unable to parse tokens from query:\n'%s'\n", query);
    fflush (stdout);
    return FALSE;
  }

  return TRUE;
}

static CharPtr FindByName (CharPtr find)

{
  Int2  i;

/* search the tag array for the desired name, returning the associated value */

  if (find == NULL) return NULL;
  for (i = 0; i < num_tags; i++) {
    if (tag [i] != NULL && strcmp (tag [i], find) == 0) {
      return val [i];
    }
  }
  return NULL;
}

static Int2 ListHasString (CharPtr list [], CharPtr str)

{
  Int2  i;

/* search a null-terminated array of strings, returning an integer index */

  if (str == NULL || list == NULL) return -1;
  for (i = 0; list [i] != NULL; i++) {
    if (strcmp (str, list [i]) == 0) return i;
  }
  return -1;
}

#define SEND_FILENAME_ARGS_BEFORE  1
#define SEND_FILENAME_ARGS_AFTER   2
#define SEND_STDIN                 3

static void RunCustom (CharPtr tempfile)

{
  CharPtr   arguments;
  Char      buf [256];
  Char      cmmd [256];
  size_t    ct;
  FILE*     fp;
  Int2      meth = SEND_STDIN;
  CharPtr   method;
  CharPtr   port;
  CharPtr   program;
  long int  val;

/* protect against custom requests on public server */

  port = getenv ("SERVER_PORT");
  if (port != NULL) {
    if (sscanf (port, "%ld", &val) == 1 && val == 80) {
      return;
    }
  }

/* program argument would be actual path on cgi server machine */

  program = FindByName ("program");
  if (program == NULL) {
    printf (">Message\nFAILURE - no program path specified\n");
    fflush (stdout);
    return;
  }

/* method defaults to sending data to program via stdin */

  method = FindByName ("method");
  if (method == NULL || strstr (method, "stdin") != NULL) {
    meth = SEND_STDIN;
  } else if (strstr (method, "filename") != NULL) {
    meth = SEND_FILENAME_ARGS_AFTER;
  } else if (strstr (method, "arguments") != NULL) {
    meth = SEND_FILENAME_ARGS_BEFORE;
  } else {
    meth = SEND_STDIN;
  }

/* note that for arguments, %20 in a query string is converted to a space */

  arguments = FindByName ("arguments");
  if (arguments == NULL) {
    arguments = "";
  }

/* launch program sending arguments and filename or data in appropriate order */

  switch (meth) {
    case SEND_FILENAME_ARGS_BEFORE :
      sprintf (cmmd, "%s %s %s", program, arguments, tempfile);
      break;
    case SEND_FILENAME_ARGS_AFTER :
      sprintf (cmmd, "%s %s %s", program, tempfile, arguments);
      break;
    case SEND_STDIN :
    default :
      sprintf (cmmd, "%s %s < %s", program, arguments, tempfile);
      break;
  }
  fp = popen (cmmd, "r");
  if (fp == NULL) return;

/* assumes program sends output to stdout */

  while ((ct = fread (buf, 1, sizeof (buf), fp)) > 0) {
    fwrite (buf, 1, ct, stdout);
    fflush (stdout);
  }
  pclose (fp);
}

static void RunEcho (CharPtr tempfile)

{
  CharPtr  ptr;

/* simply echo URL query string */

  ptr = getenv ("QUERY_STRING");
  if (ptr != NULL) {
    printf ("%s\n", ptr);
    fflush (stdout);
  }
}

static void RunSeg (CharPtr tempfile)

{
  Char     buf [256];
  Char     cmmd [256];
  size_t   ct;
  FILE*    fp;
  float    hi;
  CharPtr  hicut;
  CharPtr  lowcut;
  Char     tmp [16];
  CharPtr  window;

/* launch seg with -x parameter, arguments, and name of data file */

  window = FindByName ("window");
  if (window == NULL) {
    window = "12";
  }
  lowcut = FindByName ("lowcut");
  if (lowcut == NULL) {
    lowcut = "2.2";
  }
  hicut = FindByName ("hicut");
  if (hicut == NULL || hicut [0] == '\0') {
    if (sscanf (lowcut, "%f", &hi) == 1) {
      hi += 0.3;
      sprintf (tmp, "%4.2f", hi);
      hicut = tmp;
    } else {
      hicut = "2.5";
    }
  }

  sprintf (cmmd, "/usr/ncbi/bin/seg %s %s %s %s -x", tempfile, window, lowcut, hicut);
  fp = popen (cmmd, "r");
  if (fp == NULL) return;

/* send processed FASTA data from seg directly to stdout and calling program */

  while ((ct = fread (buf, 1, sizeof (buf), fp)) > 0) {
    fwrite (buf, 1, ct, stdout);
    fflush (stdout);
  }
  pclose (fp);
}

#define MAX_FIELDS  9

static void RunTrnaScan (CharPtr tempfile)

{
  CharPtr   aa;
  CharPtr   beg;
  Char      buf [256];
  Char      cmmd [256];
  CharPtr   end;
  CharPtr   field [MAX_FIELDS];
  FILE*     fp;
  CharPtr   id;
  Int2      idNotSent = TRUE;
  Int2      inBody = FALSE;
  CharPtr   intronBeg;
  CharPtr   intronEnd;
  long int  intronStart;
  long int  intronStop;
  Int2      numFields = 0;
  CharPtr   ptr;
  long int  start;
  long int  stop;

/* launch tRNAscan-SE with -q parameter and name of data file */

  sprintf (cmmd, "/am/MolBio/trnascan-SE/bin/tRNAscan-SE -q %s", tempfile);
  fp = popen (cmmd, "r");
  if (fp == NULL) return;

/* line by line processing of tRNAscan-SE output table */

  while (fgets (buf, sizeof (buf), fp) != NULL) {

    if (inBody) {
      memset (field, 0, sizeof (field));

/*
*  parse tab-delimited output line into array of fields, avoiding use of
*  strtok so that empty columns (adjacent tabs) are properly assigned to
*  field array
*/

      ptr = buf;
      for (numFields = 0; numFields < MAX_FIELDS && ptr != NULL; numFields++) {
        field [numFields] = ptr;
        ptr = strchr (ptr, '\t');
        if (ptr != NULL) {
          *ptr = '\0';
          ptr++;
        }
      }

/* interested in ID, start, stop, amino acid, and intron start and stop */

      id = field [0];
      beg = field [2];
      end = field [3];
      aa = field [4];
      intronBeg = field [6];
      intronEnd = field [7];

      if (numFields > 7 &&
          sscanf (beg, "%ld", &start) == 1 &&
          sscanf (end, "%ld", &stop) == 1 &&
          sscanf (intronBeg, "%ld", &intronStart) == 1 &&
          sscanf (intronEnd, "%ld", &intronStop) == 1) {

/* first line of output gives SeqId from FASTA definition line */

        if (idNotSent) {
          printf (">Features %s tRNAscan-SE\n", id);
          fflush (stdout);
          idNotSent = FALSE;
        }

/* first line of feature has start (tab) stop (tab) feature key */
/* multiple intervals would have lines of start (tab) stop */

        if (intronStart == 0 && intronStop == 0) {
          printf ("%ld\t%ld\ttRNA\n", (long) start, (long) stop);
        } else {
          printf ("%ld\t%ld\ttRNA\n", (long) start, (long) (intronStart - 1));
          printf ("%ld\t%ld\t\n", (long) (intronStop + 1), (long) stop);
        }

/* qualifier lines are (tab) (tab) (tab) qualifier key (tab) value */

        if (strstr (aa, "Pseudo") != NULL) {
          printf ("\t\t\tnote\ttRNA-Pseudo\n");
        } else {
          printf ("\t\t\tproduct\t%s\n", aa);
        }
        fflush (stdout);
      }
    }

/* detect last line of table header, ignoring everything before data section */

    if (strstr (buf, "-----") != NULL) {
      inBody = TRUE;
    }
  }
  pclose (fp);
}

/* this program can provide several services, specified in URL query string */

static CharPtr services [] = {
  "", "custom", "echo", "seg", "trnascan", NULL
};

/* main function of cgi program, called by HTTPD server */

main (int argc, char *argv[])

{
  Char     buf [256];
  size_t   ct;
  FILE*    fp;
  CharPtr  method;
  CharPtr  request;
  Int2     service;
  Char     tempfile [L_tmpnam];

/* at startup, first verify environment */

  method = getenv ("REQUEST_METHOD");
  if (method == NULL) {
    printf ("Content-type: text/html\n\n");
    printf (">Message\nFAILURE - Program launched from command line\n");
    fflush (stdout);
    return 1;
  }

/* ensure that the POST method is being sent from the HTTPD server */

  if (strcmp (method, "POST") != 0) {
    printf ("Content-type: text/html\n\n");
    printf (">Message\nFAILURE - Method (%s) must be POST\n", method);
    fflush (stdout);
    return 1;
  }

/* parse POST query into tag and val arrays */

  if (! ParseQuery (TRUE)) {
    fflush (stdout);
    return 1;
  }

/* expect request=custom, request=echo, request=seg, or request=trnascan */

  request = FindByName ("request");
  if (request == NULL) {
    printf ("Content-type: text/html\n\n");
    printf (">Message\nFAILURE - No service request\n");
    fflush (stdout);
    return 1;
  }

/* compare request value against list of available services */

  service = ListHasString (services, request);
  if (service < 1) {
    printf ("Content-type: text/html\n\n");
    printf (">Message\nFAILURE - Unable to match request '%s'\n", request);
    fflush (stdout);
    return 1;
  }

/*
*  copy all of stdin to temporary file before sending anything to stdout,
*  to get around an apparent limitation in the HTTPD implementation that
*  can cause socket deadlock
*/

  tmpnam (tempfile);
  fp = fopen (tempfile, "w");
  if (fp == NULL) {
    printf ("Content-type: text/html\n\n");
    printf (">Message\nFAILURE - file open failed\n");
    fflush (stdout);
    return 1;
  }

  while ((ct = fread (buf, 1, sizeof (buf), stdin)) > 0) {
    fwrite (buf, 1, ct, fp);
  }
  fflush (fp);
  fclose (fp);

/* now send required first header information to stdout */

  printf ("Content-type: text/html\n\n");
  fflush (stdout);

/* call appropriate external analysis program */

  switch (service) {

#ifdef ALLOW_CUSTOM_PROGRAM

/* for security, custom programs not allowed without symbol at compile time */

    case 1 :
      RunCustom (tempfile);
      break;

#endif /* ALLOW_CUSTOM_PROGRAM */

    case 2 :
      RunEcho (tempfile);
      break;
    case 3 :
      RunSeg (tempfile);
      break;
    case 4 :
      RunTrnaScan (tempfile);
      break;
    default :
      break;
  }

/* flush buffer, cleanup temporary files and allocated memory, and exit */

  fflush (stdout);
  remove (tempfile);
  free (query);
  return 0;
}
















