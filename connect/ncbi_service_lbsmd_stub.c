/*  $Id: ncbi_service_lbsmd_stub.c,v 6.3 2002/04/13 06:40:44 lavr Exp $
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
 * Author:  Anton Lavrentiev
 *
 * File Description:
 *   Dummy LBSMD mapper for non-UNIX platforms.
 *
 * --------------------------------------------------------------------------
 * $Log: ncbi_service_lbsmd_stub.c,v $
 * Revision 6.3  2002/04/13 06:40:44  lavr
 * Few tweaks to reduce the number of syscalls made
 *
 * Revision 6.2  2001/09/10 21:25:35  lavr
 * Unimportant code style compliance change
 *
 * Revision 6.1  2000/10/06 18:06:03  lavr
 * Initial revision
 *
 * ==========================================================================
 */

#include "ncbi_servicep_lbsmd.h"


const SSERV_VTable* SERV_LBSMD_Open(SERV_ITER iter,
                                    SSERV_Info** info, char** env)
{
    return 0;
}
