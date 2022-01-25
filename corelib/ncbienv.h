/*   ncbimain.h
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
* File Name:  ncbimain.h
*
* Author:  Gish, Kans, Ostell, Schuler
*
* Version Creation Date:   7/7/91
*
* $Revision: 6.0 $
*
* File Description:
*   	protokeys for portable string routines
*
* Modifications:
* --------------------------------------------------------------------------
* Date      Name        Description of modification
* --------  ----------  -----------------------------------------------------
* 06-14-94  Schuler     Created this file.  These definitions previously
*                       resided in ncbimain.h
* 06-14-94  Schuler     Added some new functions:
*                       GetAppParamBoolean() , SetAppParamBoolean(),
*                       GetAppParamInt2(), SetAppParamInt2(),
*                       GetAppParamInt4(), SetAppParamInt4()
* 06-14-94  Schuler     Added LIBCALL to ProgramPath
*
* $Log: ncbienv.h,v $
* Revision 6.0  1997/08/25 18:15:21  madden
* Revision changed to 6.0
*
* Revision 5.1  1996/12/03 21:48:33  vakatov
* Adopted for 32-bit MS-Windows DLLs
*
 * Revision 5.0  1996/05/28  13:18:57  ostell
 * Set to revision 5.0
 *
 * Revision 4.1  1995/10/06  15:54:00  epstein
 * add CacheAppParam() and FlushAppParam()
 *
 * Revision 4.0  1995/07/26  13:46:50  ostell
 * force revision to 4.0
 *
 * Revision 1.2  1995/05/15  18:45:58  ostell
 * added Log line
 *
*
* ==========================================================================
*/

#ifndef __NCBIenv_h__
#define __NCBIenv_h__

#undef NLM_EXTERN
#ifdef NLM_IMPORT
#define NLM_EXTERN NLM_IMPORT
#else
#define NLM_EXTERN extern
#endif


#ifdef __cplusplus
extern "C" {
#endif

NLM_EXTERN Nlm_Boolean LIBCALL Nlm_FindPath PROTO((Nlm_CharPtr file, Nlm_CharPtr section, Nlm_CharPtr type, Nlm_CharPtr buf, Nlm_Int2 buflen));
NLM_EXTERN Nlm_Boolean LIBCALL Nlm_TransientSetAppParam PROTO((Nlm_CharPtr file, Nlm_CharPtr section, Nlm_CharPtr type, Nlm_CharPtr value));
NLM_EXTERN void LIBCALL Nlm_FreeConfigStruct PROTO((void));
NLM_EXTERN void LIBCALL Nlm_ProgramPath PROTO((Nlm_CharPtr buf, size_t maxsize));

NLM_EXTERN void LIBCALL Nlm_FlushAppParam PROTO((void));
NLM_EXTERN Nlm_Boolean LIBCALL Nlm_CacheAppParam PROTO((Nlm_Boolean value));

NLM_EXTERN Nlm_Int2 LIBCALL Nlm_GetAppParam PROTO((Nlm_CharPtr file, Nlm_CharPtr section, Nlm_CharPtr type, Nlm_CharPtr dflt, Nlm_CharPtr buf, Nlm_Int2 buflen));
NLM_EXTERN Nlm_Boolean LIBCALL Nlm_SetAppParam PROTO((Nlm_CharPtr file, Nlm_CharPtr section, Nlm_CharPtr type, Nlm_CharPtr value));
NLM_EXTERN Nlm_Boolean LIBCALL GetAppParamBoolean PROTO((const char *filebase, const char *sect, const char *key, Nlm_Boolean dflt));
NLM_EXTERN Nlm_Boolean LIBCALL SetAppParamBoolean PROTO((const char *filebase, const char *sect, const char *key, Nlm_Boolean value));
NLM_EXTERN Nlm_Boolean LIBCALL SetAppParamLong PROTO((const char *filebase, const char *sect, const char *key, long value));
NLM_EXTERN long LIBCALL GetAppParamLong PROTO((const char *filebase, const char *sect, const char *key, long dflt));
#define GetAppParamInt(a,b,c,d)   (int)GetAppParamLong(a,b,c,(long)(d))
#define SetAppParamInt(a,b,c,d)   SetAppParamLong(a,b,c,(long)(d))
#define GetAppParamShort(a,b,c,d) (short)GetAppParamLong(a,b,c,(long)(d))
#define SetAppParamShort(a,b,c,d) SetAppParamLong(a,b,c,(long)(d))

#define FindPath Nlm_FindPath
#define FlushAppParam Nlm_FlushAppParam
#define CacheAppParam Nlm_CacheAppParam
#define GetAppParam Nlm_GetAppParam
#define SetAppParam Nlm_SetAppParam
#define TransientSetAppParam Nlm_TransientSetAppParam
#define FreeConfigStruct Nlm_FreeConfigStruct
#define ProgramPath Nlm_ProgramPath
#define GetAppParamInt2 GetAppParamShort
#define GetAppParamInt4 GetAppParamLong

#ifdef __cplusplus
}
#endif


#undef NLM_EXTERN
#ifdef NLM_EXPORT
#define NLM_EXTERN NLM_EXPORT
#else
#define NLM_EXTERN
#endif

#endif

