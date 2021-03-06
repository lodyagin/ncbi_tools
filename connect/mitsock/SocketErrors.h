/* $Copyright:
 *
 * Copyright ? 1998-1999 by the Massachusetts Institute of Technology.
 * 
 * All rights reserved.
 * 
 * Permission to use, copy, modify, and distribute this software and its
 * documentation for any purpose and without fee is hereby granted,
 * provided that the above copyright notice appear in all copies and that
 * both that copyright notice and this permission notice appear in
 * supporting documentation, and that the name of M.I.T. not be used in
 * advertising or publicity pertaining to distribution of the software
 * without specific, written prior permission.  Furthermore if you modify
 * this software you must label your software as modified software and not
 * distribute it in such a fashion that it might be confused with the
 * original MIT software. M.I.T. makes no representations about the
 * suitability of this software for any purpose.  It is provided "as is"
 * without express or implied warranty.
 * 
 * THIS SOFTWARE IS PROVIDED ``AS IS'' AND WITHOUT ANY EXPRESS OR IMPLIED
 * WARRANTIES, INCLUDING, WITHOUT LIMITATION, THE IMPLIED WARRANTIES OF
 * MERCHANTIBILITY AND FITNESS FOR A PARTICULAR PURPOSE.
 * 
 * Individual source code files are copyright MIT, Cygnus Support,
 * OpenVision, Oracle, Sun Soft, FundsXpress, and others.
 * 
 * Project Athena, Athena, Athena MUSE, Discuss, Hesiod, Kerberos, Moira,
 * and Zephyr are trademarks of the Massachusetts Institute of Technology
 * (MIT).  No commercial use of these trademarks may be made without prior
 * written permission of MIT.
 * 
 * "Commercial use" means use of a name in a product or other for-profit
 * manner.  It does NOT prevent a commercial firm from referring to the MIT
 * trademarks in order to convey information (although in doing so,
 * recognition of their trademark status should be given).
 * $
 */

/* $Header: /src/NCBI/vault.ncbi/distrib/connect/mitsock/SocketErrors.h,v 1.1 2001/04/03 20:35:40 juran Exp $ */

/* 
 *
 * SocketErrors.h -- Error codes for socket errors.
 *
 */

/* NOTE: Before you add a new error code, check the other libraries to make sure that you
   have not taken an error code range designated for another library. */


#ifndef _SOCKET_ERRORS_
#define _SOCKET_ERRORS_


/* New error definitions */
#define kSocketsFirstErr	200						/*  The beginning of the sockets errors */
#define kENOTDIRErr			201						/*	Not a directory						*/
#define kEISDIRErr			202						/*	Is a directory						*/
#define kEPFNOSUPPORTErr	203						/*	Protocol family not supported		*/
#define kEAFNOSUPPORTErr	204						/*	Address family not supported		*/
#define kSockNotInitErr		205						/*	Sockets lib is not initialized		*/
#define kSockAlreadyInitErr	206						/*	Sockets lib is already initialized	*/
#define kNoOTErr			207						/*	Open Transport is unavailable		*/
#define kSocketIsIdleErr	208						/*  No operation in progress on socket	*/
#define kEMFILEErr			209						/*  Too many sockets open				*/
#define kENOSPCErr			210						/*  Not enough space to write output	*/
#define kEBADDOMAINNAMEErr	211						/*  Bad domain name in TCP/IP settings	*/
#define kEBADNAMESERVERSErr	212						/*  Bad name servers in TCP/IP settings	*/
#define kENONETWORKErr		213						/*  No network: check TCP/IP settings	*/
#define kSocketsLastErr		299						/*  The last sockets error	*/

#ifndef rez  /* This part cannot be included by a rez file */

#include <OpenTptInternet.h>
//#include <errno.h>

#undef ERANGE  /* defined by errno.h -- but we want the OT definition */

/* Mappings from errno.h errors to OT LibC errors */
enum OTCompatErrors
{
  EPERM					= kEPERMErr,				/*  Permission denied					*/
  ENOENT				= kENOENTErr,				/*  No such file or directory			*/
  ESRCH					= kESRCHErr,				/* 	No such process						*/
  ENORSRC				= kENORSRCErr,				/*  No such resource					*/
  EINTR					= kEINTRErr,				/*  Interrupted system service			*/
  EIO					= kEIOErr,					/*  I/O error							*/
  ENXIO					= kENXIOErr,				/*  No such device or address			*/
  EBADF					= kEBADFErr,				/*  Bad file number						*/
  EAGAIN				= kEAGAINErr,				/*  Try operation again later			*/
  ENOMEM				= kENOMEMErr,				/*  Not enough space					*/
  EACCES				= kEACCESErr,				/*  Permission denied					*/
  EFAULT				= kEFAULTErr,				/*  Bad address							*/
  EBUSY					= kEBUSYErr,				/*  Device or resource busy				*/
  EEXIST				= kEEXISTErr,				/*  File exists							*/
  ENODEV				= kENODEVErr,				/*  No such device						*/
  ENOTDIR				= kENOTDIRErr,				/*	Not a directory						*/
  EISDIR				= kEISDIRErr,				/*	Is a directory						*/
  EMFILE				= kEMFILEErr,				/*  Too many sockets open				*/
  EINVAL				= kEINVALErr,				/*  Invalid argument					*/
  ENOTTY				= kENOTTYErr,				/*  Not a character device				*/
  EPIPE					= kEPIPEErr,				/*  Broken pipe							*/
  ENOSPC				= kENOSPCErr,				/*  Not enough space to write output	*/

  ERANGE				= kERANGEErr,				/*  Message size too large for STREAM	*/
  EDEADLK				= kEDEADLKErr,				/*  or a deadlock would occur			*/

  EPROTO				= kEPROTOErr,				/* 	Protocol error						*/
  EBADMSG				= kEBADMSGErr,				/* 	Trying to read unreadable message	*/
  ECANCEL				= kECANCELErr,				/* 	Operation cancelled					*/
  ENOMSG				= kENOMSGErr,				/* 	No message of desired type			*/

  ENOSTR				= kENOSTRErr,				/* 	Device not a stream					*/
  ENODATA				= kENODATAErr,				/* 	No data (for no delay I/O)			*/
  ETIME					= kETIMEErr,				/* 	Timer expired						*/
  ENOSR					= kENOSRErr,				/* 	Out of streams resources			*/

  ENOTSOCK				= kENOTSOCKErr,				/*  Socket operation on non-socket		*/
  EDESTADDRREQ			= kEDESTADDRREQErr,			/*  Destination address required		*/
  EMSGSIZE				= kEMSGSIZEErr,				/*  Message too long					*/
  EPROTOTYPE			= kEPROTOTYPEErr,			/*  Protocol wrong type for socket		*/
  ENOPROTOOPT			= kENOPROTOOPTErr,			/*  Protocol not available				*/
  EPROTONOSUPPORT		= kEPROTONOSUPPORTErr, 		/*  Protocol not supported				*/
  ESOCKTNOSUPPORT		= kESOCKTNOSUPPORTErr, 		/*  Socket type not supported			*/
  EOPNOTSUPP			= kEOPNOTSUPPErr,			/*  Operation not supported on socket	*/
  EPFNOSUPPORT			= kEPFNOSUPPORTErr,			/*	Protocol family not supported		*/
  EAFNOSUPPORT			= kEAFNOSUPPORTErr,			/*	Address family not supported		*/

  EADDRINUSE			= kEADDRINUSEErr,			/*  Address already in use				*/
  EADDRNOTAVAIL			= kEADDRNOTAVAILErr,		/*  Can't assign requested address		*/
  ENETDOWN				= kENETDOWNErr,				/*  No network, check TCP/IP settings	*/
  ENETUNREACH			= kENETUNREACHErr,			/*  Network is unreachable				*/
  ENETRESET				= kENETRESETErr,			/*  Network dropped connection on reset	*/

  ECONNABORTED			= kECONNABORTEDErr,			/*  Software caused connection abort	*/
  ECONNRESET			= kECONNRESETErr,			/*  Connection reset by peer			*/
  ENOBUFS				= kENOBUFSErr,				/*  No buffer space available			*/
  EISCONN				= kEISCONNErr,				/*  Socket is already connected			*/
  ENOTCONN				= kENOTCONNErr,				/*  Socket is not connected				*/
  ESHUTDOWN				= kESHUTDOWNErr,			/*  Can't send after socket shutdown	*/
  ETOOMANYREFS			= kETOOMANYREFSErr,			/*  Too many references: can't splice	*/
  ETIMEDOUT				= kETIMEDOUTErr,			/*  Connection timed out				*/
  ECONNREFUSED			= kECONNREFUSEDErr,			/*  Connection refused					*/
  EHOSTDOWN				= kEHOSTDOWNErr,			/*  Host is down						*/
  EHOSTUNREACH			= kEHOSTUNREACHErr,			/*  No route to host					*/
  EWOULDBLOCK			= kEWOULDBLOCKErr,			/*  Call would block, so was aborted	*/
  EALREADY				= kEALREADYErr,				/* 	Operation already in progress		*/
  EINPROGRESS			= kEINPROGRESSErr,			/* 	Operation now in progress			*/

  EBADDOMAINNAME		= kEBADDOMAINNAMEErr,		/*  Bad domain name in TCP/IP settings	*/
  EBADNAMESERVERS		= kEBADNAMESERVERSErr,		/*  Bad name servers in TCP/IP settings	*/
  ENONETWORK			= kENONETWORKErr			/*  No network: check TCP/IP settings	*/
};

#endif /* !rez */

#endif /* _SOCKET_ERRORS_ */
