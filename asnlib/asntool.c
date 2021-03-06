/*  asntool.c
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
*
* File Name: asntool.c
*
* Author:  James Ostell
*
* Version Creation Date: 1/1/91
*
* $Revision: 6.21 $
*
* File Description:
*   Main routine for asntool.  Uses the ASN.1 library routines to perform
*   analysis and display of ASN.1 module specifications and associated
*   value files.  Can also produce the ASN.1 parse tree as a C header file.
*
* Modifications:  
* --------------------------------------------------------------------------
* Date     Name        Description of modification
* -------  ----------  -----------------------------------------------------
*
* ==========================================================================
*/
/*****************************************************************************
*
*   asntool.c
*   	-m ModuleFile  	(an asn.1 specification)
*   	-v ValueFile 	(print values for input)
*   	-o OutputFile	(for module structs header)
*   	-e EncodeFile	(for output of BER encode value)
*   	-d DecodeFile	(for input of BER encoded value)
*   	-p PrintFile    (print value to file)
*       -f ModoutFile   (print modules to file)
*       -l Loadable File + header file
*       -b Size         (set AsnIo buffer size)
*       -w Length       (Word length for #defines in output header file)
*   
*   	read, validate, and display an asn1 module specification
*
*****************************************************************************/

#include "asnbuild.h"
#include "asntool.h"

extern void AsnTxtReadValFile PROTO((AsnModulePtr amp, AsnIoPtr aip, AsnIoPtr aipout,
				     AsnIoPtr encode, AsnIoPtr xaipout));
extern void AsnBinReadValFile PROTO((AsnTypePtr atp, AsnIoPtr aip, AsnIoPtr aipout,
				     AsnIoPtr encode, AsnIoPtr xaipout));

#define NUMARGS 30

Args asnargs[NUMARGS] = {
	{"ASN.1 Module File",NULL,NULL,NULL,FALSE,'m',ARG_FILE_IN,0.0,0,NULL},
	{"ASN.1 Module File", NULL,NULL,NULL,TRUE,'f',ARG_FILE_OUT,0.0,0,NULL},
	{"XML DTD File\n\t(\"m\" to print each module to a separate file)", NULL,NULL,NULL,TRUE,'X',ARG_FILE_OUT,0.0,0,NULL},
	{"ASN.1 Tree Dump File", NULL,NULL,NULL,TRUE,'T',ARG_FILE_OUT,0.0,0,NULL},
	{"Print Value File",NULL,NULL,NULL,TRUE,'v',ARG_FILE_IN,0.0,0,NULL},
	{"Print Value File",NULL,NULL,NULL,TRUE,'p',ARG_FILE_OUT,0.0,0,NULL},
	{"XML Data File", NULL,NULL,NULL,TRUE,'x',ARG_FILE_OUT,0.0,0,NULL},
	{"Binary Value File (type required)",NULL,NULL,NULL,TRUE,'d',ARG_FILE_IN,0.0,0,NULL},
	{"Binary Value Type",NULL,NULL,NULL,TRUE,'t',ARG_STRING,0.0,0,NULL},
	{"Binary Value File",NULL,NULL,NULL,TRUE,'e',ARG_FILE_OUT,0.0,0,NULL},
	{"Header File",NULL,NULL,NULL,TRUE,'o',ARG_FILE_OUT,0.0,0,NULL},
	{"Loader File",NULL,NULL,NULL,TRUE,'l',ARG_FILE_OUT,0.0,0,NULL},
    {"Buffer Size","1024","512","10000",TRUE,'b',ARG_INT,0.0,0,NULL},
    {"Word length maximum for #defines","31","31","128",TRUE,'w',ARG_INT,0.0,0,NULL},

    {"Fix Non-Printing Characters\n"
     "     0 - Replace with #, post ERROR\n"
     "     1 - Replace with # silently\n"
     "     2 - Pass through silently\n"
     "     3 - Replace with #, post FATAL\n","0","0","3",TRUE,'F',ARG_INT,0.0,0,NULL},
    {"UTF8 Input Conversion\n"
     "     0 - Convert silently\n"
     "     1 - Convert, post WARNING first time\n"
     "     2 - Convert, post WARNING each time\n"
     "     3 - Do not convert\n","0","0","3",TRUE,'N',ARG_INT,0.0,0,NULL},
    {"UTF8 Output Conversion\n"
     "     0 - Convert silently\n"
     "     1 - Convert, post WARNING first time\n"
     "     2 - Convert, post WARNING each time\n"
     "     3 - Do not convert\n","0","0","3",TRUE,'U',ARG_INT,0.0,0,NULL},

   {"Generate object loader .c and .h files, \n\tif used, see below parameters:","F",NULL,NULL,TRUE,'G',ARG_BOOLEAN,0.0,0,NULL},


   {"ASN.1 module filenames, comma separated used for external refs from the \'m\',\n     but no other action taken", NULL, NULL, NULL, TRUE, 'M', ARG_FILE_IN, 0.0, 0, NULL},
   {"Base for filename, without extensions, for generated objects and code", NULL, NULL, NULL, TRUE, 'B', ARG_FILE_OUT, 0.0, 0, NULL},
   {"During code generation, debugging level\n"
    "     0 - No debugging\n"
    "     1 - Shallow debugging\n"
    "     2 - Deep\n", "0", "0", "9", TRUE, 'D', ARG_INT, 0.0, 0, NULL},
   {"Debugging filename", "stderr", NULL, NULL, TRUE, 'S', ARG_FILE_OUT, 0.0, 0, NULL},
   {"In generated .c, add #include to this filename", NULL, NULL, NULL, TRUE, 'I', ARG_STRING, 0.0, 0, NULL},
   {"Bit twiddle for optional zero value base slots","F",NULL,NULL,TRUE,'Z',ARG_BOOLEAN,0.0,0,NULL},
   {"In generated .c, forces name of #included asn header", NULL, NULL, NULL, TRUE, 'K', ARG_STRING, 0.0, 0, NULL},
   {"Register type with object manager", NULL, NULL, NULL, TRUE, 'J', ARG_STRING, 0.0, 0, NULL},
   {"Label for registered type", NULL, NULL, NULL, TRUE, 'L', ARG_STRING, 0.0, 0, NULL},
   {"XML module prefix for DOCTYPE", NULL, NULL, NULL, TRUE, 'P', ARG_STRING, 0.0, 0, NULL},
   {"Force choice to use structure instead of ValNodePtr","F",NULL,NULL,TRUE,'V',ARG_BOOLEAN,0.0,0,NULL},
   /*-- not used now {"Asnload directory [for parsetrees]", NULL, NULL, NULL, TRUE, 'd', ARG_FILE_IN, 0.0, 0, NULL}, -----------*/

   {"Use quoted syntax form for generated include files","F",NULL,NULL,TRUE,'Q',ARG_BOOLEAN,0.0,0,NULL},

};

typedef enum {
	m_argModuleIn = 0,
	f_arg_moduleOut,
	X_argDTDModuleOut,
	T_argTreeDumpOut,
	v_argPrintFileIn,
	p_argPrintFileOut,
	x_argXMLDataOut,
	d_argBinaryFileIn,
	t_argAsnTypeName,
	e_argBinaryFileOut,
	o_argHeadFile,
	l_argLoadFile,
	b_argBufferSize,
	w_argTokenMax,
	F_argFixPrint,
	N_argUtfInput,
	U_argUtfOutput,
	/*--- args below here are capitilized and for code generation, only--*/
	/*--  Except for the 'M' arg, which will also affect normal use ---*/
	G_argGenerateCode,
	M_argMoreModuleFiles,
	B_argCodeFileName,
	D_argCodeGenDebugLevel,
	S_argDebugFileName,
	I_argExtraIncludeName,
	Z_argBitTwiddle,
	K_argLoadName,
	J_objMgrEntry,
	L_objMgrLabel,
	P_argXMLmodulePrefix,
	V_argChoiceStruct,
	Q_argQuoted,
} Arguments;

Int2 Main (void)
{
	AsnIoPtr aip = NULL,
		aipout = NULL,
		xaipout = NULL,
		aipencode = NULL;
	FILE * fp;
	AsnModulePtr amp = NULL,
		currentmod = NULL,
		nextmod, thisamp;
	AsnTypePtr atp;
	Boolean print_each_module = FALSE;
	AsnCodeInfoPtr acip = (AsnCodeInfoPtr)MemNew(sizeof(AsnCodeInfo));
	CharPtr filename = NULL, p, last_comma;
	int len;


    /* never abort on error, but show it */
    ErrSetFatalLevel(SEV_MAX);
    ErrSetMessageLevel(SEV_MIN);
    asnargs[P_argXMLmodulePrefix].defaultvalue = (const char *)AsnGetXMLmodulePrefix();

    if (! GetArgs("AsnTool 7", NUMARGS, asnargs)) {
	  return 1;
    }

    ErrClear();

    AsnSetXMLmodulePrefix((CharPtr)(asnargs[P_argXMLmodulePrefix].strvalue));

    if (! AsnIoSetBufsize(NULL, (Int2)asnargs[b_argBufferSize].intvalue))
        return 1;

    if ((aip = AsnIoOpen(asnargs[m_argModuleIn].strvalue, "r")) == NULL)
	{
	    ErrShow();
	    return 1;
	}

    acip -> loadname = asnargs[m_argModuleIn].strvalue;

    if (asnargs[K_argLoadName].strvalue != NULL) {
	acip -> loadname = asnargs[K_argLoadName].strvalue; /* overrides m_argModuleIn, if set */
    }

    if (asnargs[e_argBinaryFileOut].strvalue != NULL) {   /* output a binary value file */
	    if ((aipencode = AsnIoOpen(asnargs[e_argBinaryFileOut].strvalue, "wb")) == NULL)
	    {
		    ErrShow();
		    return 1;
	    }
        aipencode->fix_non_print = asnargs[F_argFixPrint].intvalue;
        aipencode->fix_utf8_out = asnargs[U_argUtfOutput].intvalue;
    }

				/**  parse the module(s)  ***/
	
    if (asnargs[f_arg_moduleOut].strvalue != NULL)
	{
	    if ((aipout = AsnIoOpen(asnargs[f_arg_moduleOut].strvalue, "w")) == NULL)
		{
		    ErrShow();
		    return 1;
		}
	}

    thisamp = NULL;
    while ((nextmod = AsnLexTReadModule(aip)) != NULL )
	{
	    if (thisamp == NULL)
		thisamp = nextmod;

	    if (amp == NULL)
		amp = nextmod;
	    else
		currentmod->next = nextmod;
	    currentmod = nextmod;
	}
    acip ->  last_amp = currentmod; /* last module of main file */
    AsnStoreTree(acip->loadname, thisamp); /* store and link tree */

    /*--- read additional module files that will be used for everything
      but code generation.
      ---*/

    if (asnargs[M_argMoreModuleFiles].strvalue != NULL)
	for (p = asnargs[M_argMoreModuleFiles].strvalue; *p; p = last_comma + 1) {
	    /*--- extract next filename for reading ASN.1 definitions ---*/
	    for (last_comma = p; *last_comma; last_comma++) {
	        if (*last_comma == ',')
		    break;
	    }
	    len = last_comma - p;
	    filename = (char *)MemFree (filename);
	    filename = (char *)MemNew (len + 1);
	    StringNCpy (filename, p, len);
	    filename[len] = '\0';
#ifdef WIN_DUMB
	    printf ("Loading %s \n", filename);
#endif

	    if ((aip = AsnIoOpen(filename, "r")) == NULL)
		{
		    ErrShow();
		    return 1;
		}
	    /*--- read the modules in this current file ---*/
	    thisamp = NULL;
	    while ((nextmod = AsnLexTReadModule(aip)) != NULL )
		{
		    if (thisamp == NULL)
			thisamp = nextmod;
		    if (amp == NULL)
			amp = nextmod;
		    else
			currentmod->next = nextmod;
		    currentmod = nextmod;
		}
	    AsnStoreTree(filename, thisamp); /* store and link tree */
	    if (!*last_comma)
	        break;
	    aip = AsnIoClose(aip);
	}

    aip = AsnIoClose(aip);

    if (amp == NULL)
	{
	    ErrPostEx(SEV_FATAL,0,0, "Unable to continue due to bad ASN.1 module");
	    ErrShow();
	    return 1;
	}


    if (asnargs[f_arg_moduleOut].strvalue != NULL)
	{
	    if ((aipout = AsnIoOpen(asnargs[f_arg_moduleOut].strvalue, "w")) == NULL)
		{
		    ErrShow();
		    return 1;
		}

	    currentmod = amp;
	    do
		{
		    AsnPrintModule(currentmod, aipout);
		    if (currentmod == acip->last_amp)  /* last main module */
			currentmod = NULL;
		    else
			currentmod = currentmod->next;
		} while (currentmod != NULL);

	    aipout = AsnIoClose(aipout);
	}
	
    if (asnargs[X_argDTDModuleOut].strvalue != NULL)
	{
	    Char tbuf[250];
	    CharPtr ptr;

	    if (! StringCmp(asnargs[X_argDTDModuleOut].strvalue, "m"))
		{
		    print_each_module = TRUE;
		}
	    else
		{
		    if ((aipout = AsnIoOpen(asnargs[X_argDTDModuleOut].strvalue, "wx")) == NULL)
			{
			    ErrShow();
			    return 1;
			}
		}

	    currentmod = amp;
	    do
		{
		    if (print_each_module)
			{
			    StringMove(tbuf, currentmod->modulename);
			    for (ptr = tbuf; *ptr != '\0'; ptr++)
				{
				    if (*ptr == '-')
					*ptr = '_';
				}
			    StringMove(ptr, ".dtd");

			    AsnPrintModuleXMLInc(currentmod, tbuf);

			    StringMove(ptr, ".mod");
				
			    aipout = AsnIoOpen(tbuf, "wx");
			}
			

		    AsnPrintModuleXML(currentmod, aipout);

		    if (print_each_module)
			aipout = AsnIoClose(aipout);

		    if (currentmod == acip->last_amp)  /* last main module */
			currentmod = NULL;
		    else
			currentmod = currentmod->next;
		} while (currentmod != NULL);

	    if (! print_each_module)
		aipout = AsnIoClose(aipout);
	}
	
    if (asnargs[T_argTreeDumpOut].strvalue != NULL)
	{
	    if ((fp = FileOpen(asnargs[T_argTreeDumpOut].strvalue, "w")) == NULL)
		{
		    ErrShow();
		    return 1;
		}

	    currentmod = amp;
	    do
		{
		    AsnPrintTreeModule(currentmod, fp);
		    if (currentmod == acip->last_amp)  /* last main module */
			currentmod = NULL;
		    else
			currentmod = currentmod->next;
		} while (currentmod != NULL);

	    FileClose(fp);
	}


    acip -> amp = amp;
	
    /* print a value file */

    if (asnargs[p_argPrintFileOut].strvalue != NULL)
	{
	    if ((aipout = AsnIoOpen(asnargs[p_argPrintFileOut].strvalue, "w")) == NULL)
		{
		    ErrShow();
		    return 1;
		}
		aipout->fix_non_print = asnargs[F_argFixPrint].intvalue;
		aipout->fix_utf8_out = asnargs[U_argUtfOutput].intvalue;
	}
    /* print an XML file */

    if (asnargs[x_argXMLDataOut].strvalue != NULL)
	{
	    if ((xaipout = AsnIoOpen(asnargs[x_argXMLDataOut].strvalue, "wx")) == NULL)
		{
		    ErrShow();
		    return 1;
		}
		xaipout->fix_non_print = asnargs[F_argFixPrint].intvalue;
		xaipout->fix_utf8_out = asnargs[U_argUtfOutput].intvalue;
	}

    if (asnargs[v_argPrintFileIn].strvalue != NULL)        /* read a printvalue file */
	{
	    if ((aip = AsnIoOpen(asnargs[v_argPrintFileIn].strvalue, "r")) == NULL)
		{
		    ErrShow();
		    return 1;
		}
		aip->fix_non_print = asnargs[F_argFixPrint].intvalue;
		aip->fix_utf8_in = asnargs[N_argUtfInput].intvalue;

	    AsnTxtReadValFile(amp, aip, aipout, aipencode, xaipout);
	    ErrShow();
	}

    aip = AsnIoClose(aip);

    if (asnargs[d_argBinaryFileIn].strvalue != NULL)        /* read a ber file */
	{
	    if ((asnargs[t_argAsnTypeName].strvalue == NULL) || (! TO_UPPER(*asnargs[t_argAsnTypeName].strvalue)))
                {
		    ErrPostEx(SEV_FATAL,0,0, "Must use -t Type to define contents of decode file");
		    ErrShow();
		    return 1;
		}

	    atp = AsnTypeFind(amp, asnargs[t_argAsnTypeName].strvalue);
	    if (atp == NULL)
                {
#ifdef WIN_MSWIN
		    ErrPostEx(SEV_FATAL,0,0, "Couldn't find Type %Fs", asnargs[t_argAsnTypeName].strvalue);
#else
		    ErrPostEx(SEV_FATAL,0,0, "Couldn't find Type %s", asnargs[t_argAsnTypeName].strvalue);
#endif
		    ErrShow();
		    return 1;
		}

	    if ((aip = AsnIoOpen(asnargs[d_argBinaryFileIn].strvalue, "rb")) == NULL)
		{
		    ErrShow();
		    return 1;
		}
        aip->fix_non_print = asnargs[F_argFixPrint].intvalue;
        aip->fix_utf8_in = asnargs[N_argUtfInput].intvalue;


	    AsnBinReadValFile(atp, aip, aipout, aipencode, xaipout);
	    ErrShow();
	}

    AsnIoClose(xaipout);
    AsnIoClose(aipout);
    AsnIoClose(aip);
    AsnIoClose(aipencode);

    if (asnargs[o_argHeadFile].strvalue != NULL)         /* produce header file */
	AsnOutput(asnargs[o_argHeadFile].strvalue, amp, FALSE, (Int2)asnargs[w_argTokenMax].intvalue);

    if (asnargs[l_argLoadFile].strvalue != NULL)        /* produce loader file */
        AsnOutput(asnargs[l_argLoadFile].strvalue, amp, TRUE, (Int2)asnargs[w_argTokenMax].intvalue);

    if (asnargs[G_argGenerateCode ].intvalue != 0)
	{
	    acip ->  filename = asnargs[B_argCodeFileName].strvalue;
	    acip ->  do_bit_twiddle = asnargs[Z_argBitTwiddle].intvalue;
	    acip ->  force_choice_struct = asnargs[V_argChoiceStruct].intvalue;
	    acip ->  include_filename = asnargs[I_argExtraIncludeName].strvalue;
	    acip ->  maxDefineLength = asnargs[w_argTokenMax].intvalue;
	    acip -> debug_level = asnargs[D_argCodeGenDebugLevel ].intvalue;
	    acip -> object_manager_entry = asnargs[J_objMgrEntry].strvalue;
	    acip -> object_label = asnargs[L_objMgrLabel].strvalue;
	    if (asnargs[S_argDebugFileName].strvalue != NULL) {
		(acip -> bug_fp) = FileOpen (asnargs[S_argDebugFileName].strvalue, "w");
	    } else {
		(acip -> bug_fp) = FileOpen ("stderr", "w");
	    }
	    acip -> use_quoted_include = asnargs[Q_argQuoted ].intvalue != 0;
	    AsnCode(acip);
	}
	
	MemFree(acip);
	MemFree(filename);

    return 0;
}
/*****************************************************************************
*
*   void AsnTxtReadValFile(amp, aip, aipout, aipencode)
*   	reads a file of values
*   	prints to aipout if aipout != NULL
*
*****************************************************************************/
void AsnTxtReadValFile (AsnModulePtr amp, AsnIoPtr aip, AsnIoPtr aipout,
			AsnIoPtr aipencode, AsnIoPtr xaipout)

{
	AsnTypePtr atp;
	DataVal value;
	Boolean read_value, print_value, encode_value, restart, xmlvalue;

	if (aipout != NULL)
		print_value = TRUE;
	else
		print_value = FALSE;

	if (xaipout != NULL)
		xmlvalue = TRUE;
	else
		xmlvalue = FALSE;

	if (aipencode != NULL)
		encode_value = TRUE;
	else
		encode_value = FALSE;

	if (print_value || encode_value || xmlvalue)
		read_value = TRUE;
	else
		read_value = FALSE;


	atp = NULL;
	restart = FALSE;

	while ((atp = AsnTxtReadId(aip, amp, atp)) != NULL)
	{
		if (restart == TRUE)
		{
			if (encode_value)
			{
				ErrPostEx(SEV_FATAL,0,0, "Cannot binary encode multiple values in one message");
				return;
			}
			if (print_value)         /* new line */
			{
				AsnPrintNewLine(aipout);
				AsnPrintNewLine(aipout);
			}
			if (xmlvalue)
			{
				AsnPrintNewLine(xaipout);
				AsnPrintNewLine(xaipout);
			}
			restart = FALSE;
		}

		if (read_value)
		{
			if (! AsnTxtReadVal(aip, atp, &value))
                        {
				return;
			}
			if (print_value)
			{
				if (! AsnTxtWrite(aipout, atp, &value))
					return;
			}
			if (xmlvalue)
			{
				if (! AsnTxtWrite(xaipout, atp, &value))
					return;
			}
			if (encode_value)
			{
				if (! AsnBinWrite(aipencode, atp, &value))
					return;
			}
			AsnKillValue(atp, &value);
		}
		else
		{
			if (! AsnTxtReadVal(aip, atp, NULL))
				return;
		}

		if (! aip->type_indent)      /* finished reading an object */
		{
			atp = NULL;              /* restart */
			restart = TRUE;
		}
	}
	return;
}

/*****************************************************************************
*
*   void AsnBinReadValFile(aip, aipout, aipencode)
*   	reads a file of binary values
*       MUST start knowing the Type of the value
*   	prints to aipout if aipout != NULL
*
*****************************************************************************/
void AsnBinReadValFile (AsnTypePtr atp, AsnIoPtr aip, AsnIoPtr aipout,
			AsnIoPtr aipencode, AsnIoPtr xaipout)        /* type of message */

{
	DataVal value;
	Boolean read_value, print_value, encode_value, xmlvalue;

	if (aipout != NULL)
		print_value = TRUE;
	else
		print_value = FALSE;

	if (xaipout != NULL)
		xmlvalue = TRUE;
	else
		xmlvalue = FALSE;

	if (aipencode != NULL)
		encode_value = TRUE;
	else
		encode_value = FALSE;

	if (print_value || encode_value || xmlvalue)
		read_value = TRUE;
	else
		read_value = FALSE;


	while ((atp = AsnBinReadId(aip, atp)) != NULL)
	{
		if (read_value)
		{
			if (! AsnBinReadVal(aip, atp, &value))
				return;
			if (print_value)
			{
				if (! AsnTxtWrite(aipout, atp, &value))
					return;
			}
			if (xmlvalue)
			{
				if (! AsnTxtWrite(xaipout, atp, &value))
					return;
			}
			if (encode_value)
			{
				if (! AsnBinWrite(aipencode, atp, &value))
					return;
			}
			AsnKillValue(atp, &value);
		}
		else
			AsnBinReadVal(aip, atp, NULL);
	}
	return;
}



