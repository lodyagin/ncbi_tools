/*   $Id: PubStructAsn.c,v 6.14 1998/09/03 21:30:11 kimelman Exp $
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
 * Author:  Michael Kimelman
 *
 * File Description: PubStruct DB Asn (down)loader.
 *                   
 * Modifications:  
 * --------------------------------------------------------------------------
 * $Log: PubStructAsn.c,v $
 * Revision 6.14  1998/09/03 21:30:11  kimelman
 * added number of retries to connect to DB server
 * added softer processing to server connection failure
 * added transmission timeouts.
 *
 * Revision 6.13  1998/08/08 04:52:26  kimelman
 * bugfixes:
 * 1. enforced loading mode.
 * 2. memory  leaks
 * 3. 'pos' value in "check on load" processing
 *
 * Revision 6.12  1998/08/05 21:12:33  kimelman
 * skip/load bugfix
 *
 * Revision 6.11  1998/07/16 20:02:51  kimelman
 * enforce option added
 *
 * Revision 6.10  1998/07/14 20:24:45  kimelman
 * FT schema & smart load
 *
 * Revision 6.9  1998/06/25 19:05:49  kimelman
 * move context statics from ctlibutils to PubStructAsn
 *
 * Revision 6.8  1998/06/17 15:33:52  kimelman
 * added commit transaction at the end of removeasn.
 *
 * Revision 6.7  1998/05/28 17:33:39  kimelman
 * throw away obsolete code
 *
 * Revision 6.6  1998/05/27 18:09:16  kimelman
 * put compression stuff into production
 *
 * Revision 6.5  1998/05/15 20:20:01  kimelman
 * compr -> nlmzip
 *
 * Revision 6.4  1998/05/14 16:11:12  kimelman
 * Compression stuff added in debug mode.
 * few bugs fixed, related to OpenServer/SQL Server switching
 *
 * Revision 6.3  1998/05/08 03:03:41  kimelman
 * Open Server fix
 *
 * Revision 6.2  1998/04/15 14:53:54  kimelman
 * 1. Make retrieval unifirm from open server and sql server.
 * 2. mmdbsrv retrival performance tuning:
 * 	- size of read-in buffers
 * 	- direct pdb2mmdb DB lookup instead of full loading pdb to mmdb translataion
 * 	  table
 * 3. cleaning mmdblocl.* and make.mmdbsrv in order to remove obsolete code
 *
 * Revision 6.1  1998/04/03 20:25:15  kimelman
 * PubStruct access code added to mmdbsrv
 *
 *
 * ==========================================================================
 */


#include <PubStructAsn.h>
#include <ctlibutils.h>
#include <mmdbapi.h>
#include <nlmzip.h>

#include <assert.h>

#define DEF_SRV "BACH10:PubStruct=anyone,allowed"


/*-------------------------------------------
 * SQL server return code processing macros *
 -------------------------------------------*/

#ifdef DEBUG_MODE
#define RRC(sb_cmd)  { ErrPostEx(SEV_INFO, ERR_SYBASE,0,"%s running", #sb_cmd);   \
                       retcode = sb_cmd; }
#else
#define RRC(sb_cmd)  { retcode = sb_cmd; }
#endif

#define RC1(sb_cmd,level) {                                                             \
  RRC(sb_cmd);                                                                          \
  if(typecode(retcode,0)<level) {                                                       \
    ErrPostEx(SEV_ERROR, ERR_SYBASE,0,"\n%s:%d:'%s' failed",__FILE__,__LINE__,#sb_cmd); \
    typecode(retcode,1);                                                                \
    goto errexit;                                                                       \
} }

#define RC(sb_cmd) RC1(sb_cmd,4)

typedef enum {
  PS_NEW,
  PS_READ,
  PS_UPDATE,
  PS_STORE,
  PS_DELETE
} ps_action_t;

struct ps_chunk {
  struct ps_chunk  *next;
  int                 len;
  int                 size;
  int                 start;
  char               *data;
};

typedef struct {
  CTLibUtils       clu;
  CS_IODESC        iodesc;
  struct ps_chunk *top;
  struct ps_chunk *bottom;
  ps_action_t      action;
  char            *srv;
  int              acc;
  int              eos;
  int              cache_size;
  int              open_server;
} DB_stream_t ;


/*****************************************************************************
*****************************************************************************/

static int
typecode(CS_RETCODE retcode,int enforce)
{
    switch(retcode)
      {
#ifdef DEBUG_MODE
#define QQ(code,rc) case code: \
        if (enforce>rc) { \
          ErrPostEx(SEV_ERROR, ERR_SYBASE,0,"retcode = %d(%s)\n",retcode,#code); return (rc); \
        } else { \
          ErrPostEx(SEV_INFO , ERR_SYBASE,0,"retcode = %d(%s)\n",retcode,#code); return (rc); \
        }
#else
#define QQ(code,rc) case code: \
        if (enforce>rc) ErrPostEx(SEV_ERROR, ERR_SYBASE,0,"retcode = %d(%s)\n",retcode,#code); return (rc);
#endif
      QQ(CS_SUCCEED,4);
      QQ(CS_END_ITEM,3);
      QQ(CS_END_DATA,2);
      QQ(CS_END_RESULTS,1);
      QQ(CS_FAIL,0);    
      QQ(CS_BUSY,0);    
      QQ(CS_PENDING,0); 
      QQ(CS_CANCELED,0);
      default:
        ErrPostEx(SEV_INFO, ERR_SYBASE,0,"retcode = %d(default)\n",retcode);
      }
#undef QQ
    return 0;
}

#ifdef DEBUG_MODE
static int
typeres(CS_INT restype)
{
    switch(restype)
      {
#define QQ(code) case code: ErrPostEx(SEV_INFO, ERR_SYBASE,0,"res_type = %d(%s)\n",restype,#code); break;
        QQ(CS_CMD_DONE);
        QQ(CS_CMD_FAIL);
        QQ(CS_CMD_SUCCEED);
        QQ(CS_COMPUTE_RESULT);
        QQ(CS_CURSOR_RESULT);
        QQ(CS_PARAM_RESULT);
        QQ(CS_ROW_RESULT);
        QQ(CS_STATUS_RESULT);
        QQ(CS_COMPUTEFMT_RESULT);
        QQ(CS_ROWFMT_RESULT);
        QQ(CS_MSG_RESULT);
        QQ(CS_DESCRIBE_RESULT);
      default:
        ErrPostEx(SEV_INFO, ERR_SYBASE,0,"retcode = %d(default)\n",restype);
      }
#undef QQ
    return 0;
}
#else
#define typeres(rc)
#endif


static struct ps_chunk *
new_piece(int min_size, int max_size)
{
  struct ps_chunk *piece;
  piece = MemNew(sizeof(*piece));
  assert(piece);
  piece->start = 0;
  piece->len = 0;
  piece->next = NULL;
  piece->size = max_size;
  while ((piece->data = MemNew(piece->size)) == NULL)
    {
      if (piece->size == min_size)
        {
          MemFree(piece);
          ErrPostEx(SEV_ERROR, ERR_SYBASE,0,"\n%s:%d: memory exhausted '%d' ",
                    __FILE__,__LINE__,min_size);
          return NULL;
        }
      else
        {
          piece->size /= 2;
          if (piece->size < min_size)
            piece->size = min_size;
        }
    } 
  return piece;
}

/*****************************************************************************
*   db handler interface  
*****************************************************************************/
static void
pubstruct_db_close(DB_stream_t *db)
{
  CS_INT restype;
  if (!db)
    return; /*???? kind of assert */
  if (db->clu.ctcmd && db->action == PS_READ)
    ct_cancel(db->clu.connection,NULL, CS_CANCEL_ALL);
  db->clu.context = NULL; /* avoid cleaning context */
  CTLibDrop(&db->clu);
  while (db->top)
    {
      struct ps_chunk *piece = db->top ;
      db->top = db->top->next;
      MemFree(piece->data);
      MemFree(piece);
    }
  if (db->srv)
    MemFree(db->srv);
  MemFree(db);
  ErrShow();
}

static DB_stream_t *
pubstruct_db_open(char *server,ps_action_t action)
{
  static int               done    = 0   ;
  static int               retries = 0   ;
  static CS_CONTEXT PNTR   context = NULL;
         DB_stream_t      *db ;
  
  if (!done)
    {
      if (getenv("DB_DEBUG"))
        {
          ErrSetLogLevel(0);
          ErrSetOptFlags(EO_LOGTO_STDERR);
        }
      done = 1;
    }
  
  assert ( action == PS_READ || action == PS_NEW || action == PS_UPDATE);
  
  db = MemNew(sizeof(*db));
  db->top = db->bottom = NULL;
  if (server == NULL)
    server = DEF_SRV;
  db->srv = MemNew(strlen(server)+1);
  strcpy(db->srv,server);
  {
    char *os = strstr(db->srv,"_OS");
    if (os)
      {
        db->open_server = 1;
        db->clu.used4os = TRUE;
      }
  }
  db->clu.context = context;
  if(!CTLibInit(&db->clu,db->srv,NULL,NULL,NULL,0,NULL))
    goto errexit;
  
  if(context==NULL)
    {
      CS_INT to = 60;
      context=db->clu.context;
      ct_config(context,CS_SET,CS_TIMEOUT,(CS_VOID*)&to,CS_UNUSED,NULL);
    }
  
  if(action != PS_READ)
    if(!CTLibSimpleSQL_Ex(db->clu.ctcmd,"begin transaction"))
      goto errexit;
  db->action = action;
  return db;
errexit:
  pubstruct_db_close(db);
  db = NULL;
  ErrPostEx(SEV_ERROR,  ERR_SYBASE, 0,"Connection to %s - failed",server);
  ErrShow();
  retries++;
  if (retries < 20)
    {
      sleep(60);
      db = pubstruct_db_open(server,action);
    }
  retries--;
  return db;
}

/*
 *  DB Asn callbacks
 */

static Int4 LIBCALLBACK
dbio_read(Pointer ptr, CharPtr obuf, Int4 count)
{
  DB_stream_t     *db = (DB_stream_t*)ptr;
  CS_COMMAND PNTR  cmd = db->clu.ctcmd;
  CS_INT           bytes = 0;
  CS_RETCODE       retcode;

  assert(db->action == PS_READ);

  if (db->eos)
    return 0;
  RRC(ct_get_data(cmd, 1,(CS_TEXT*)obuf,(CS_INT)count, &bytes));
  typecode(retcode,2); /* print result code in debug mode */
  if (bytes < count )
    db->eos = 1 ;
  return bytes;
}

static Int4 LIBCALLBACK
dbio_write(Pointer ptr, CharPtr buf, Int4 count)
{
  DB_stream_t     *db    = (DB_stream_t*)ptr;
  CS_COMMAND PNTR  cmd   = db->clu.ctcmd;
  CS_RETCODE       retcode;
  struct ps_chunk *piece;
  int              desired_size = 100 * 1024;
  
  assert(db->action != PS_READ);
  assert(count >= 0);
  if(count==0)
    return 0;
  /*
   * because sybase required to say it the full length of data before
   * providing it - we can only collect data in memory until we got a
   * 'close' command. At that moment we will be able to calculate and store a whole
   * bunch of data -  a bit stupid activity - but...
   */
  if (db->bottom)
    piece = db->bottom;
  else
    db->bottom = db->top = piece = new_piece(count, desired_size);
  if (piece->size - piece->len < count)
    db->bottom = db->bottom->next = piece = new_piece(count, desired_size);
  assert (piece);
  assert (piece->size - piece->len >= count );
  memcpy (piece->data + piece->len,buf,count);
  piece->len += count;
  db->iodesc.total_txtlen += count;
  
  return count;
}

#define dbio_close pubstruct_closeasn

#define dbio_open(db) fci_open(db,(db->action==PS_READ?dbio_read:dbio_write),dbio_close)

/*****************************************************************************
*
*   AsnIoPtr AsnIo_PubStruct_Close (file_name, mode)
*
*****************************************************************************/

static int
pubstruct_parseasn(char *postupdate_cmd, Int4 acc, char* buf, int buflen)
{
  static AsnTypePtr atp, biostruc,mmdb_tp,bs;
  static AsnIoPtr aip;
  static AsnIoMemPtr aimp ;
  static AsnModulePtr amp = NULL;
  Int4  mmdb_id = 0;
  BiostrucSourcePtr bsp = NULL;
  
  if (amp == NULL)
    {
      /* initialization */
      
      if (! (objmmdb1AsnLoad() &&
             objmmdb2AsnLoad() &&
             objmmdb3AsnLoad()))   /* load Biostruc defintions */
        {
          ErrShow();
          return 0;
        }
      
      amp = AsnAllModPtr();
      biostruc = AsnFind("Biostruc");
      bs       = AsnFind("Biostruc-history.data-source");
      mmdb_tp  = AsnFind("Biostruc.id.E.mmdb-id");
      assert(bs);
    }
  
  aimp = AsnIoMemOpen("rb", (UcharPtr)buf, buflen);
  aip = aimp->aip;
  atp = biostruc;

  while (1)
    {
      atp = AsnReadId(aip, amp, atp);
      if (atp == bs)
        {
          assert (bsp == NULL);
          assert (mmdb_id != 0 );
          bsp = BiostrucSourceAsnRead(aip, atp);
          break;
        }
      else if (atp == mmdb_tp)
        {
          DataVal dv;
          AsnReadVal(aip, atp, &dv);
          if (mmdb_id == 0)
            mmdb_id = dv.intvalue;
          else
            if ( mmdb_id != dv.intvalue )
              {
                /* seems to be in replacement branch */
                ErrPostEx(SEV_INFO, 0, 0," %d replaced %d",
                          mmdb_id,dv.intvalue);
              }
            else
              ErrPostEx(SEV_ERROR, 0, 0,"PubStruct : second occurance of mmdb_id ");
          if(postupdate_cmd==NULL)
            goto exit;
        }
      else if (atp == NULL)
        goto err;
      else
        AsnReadVal(aip, atp, NULL);  /* skip it */
    }
  assert (bsp != NULL);
  assert (mmdb_id != 0 );

  sprintf(postupdate_cmd,"exec new_struct1 %d,%d,",acc,mmdb_id);

  {  /*date*/
    DatePtr    d  = bsp->database_entry_date;
    ValNodePtr vn = bsp->VersionOfDatabase_version_of_database;

    if (vn->choice == VersionOfDatabase_version_of_database_release_date)
        d = (DatePtr)(vn->data.ptrvalue);
    
    if (d->data[0]==0)
      sprintf(postupdate_cmd+strlen(postupdate_cmd),"\'%s\',",d->str);
    else
      sprintf(postupdate_cmd+strlen(postupdate_cmd),"\'%d/%d/%d\',",
              (int)d->data[2],(int)d->data[3],1900+(int)d->data[1]);
  }
  {  /*pdb*/
    ValNodePtr vn =  bsp->database_entry_id;
    DbtagPtr   dbtag = (DbtagPtr)(vn->data.ptrvalue);
    assert(vn->choice == BiostrucId_other_database);
    assert( strcmp(dbtag->db,"PDB")==0);
    sprintf(postupdate_cmd+strlen(postupdate_cmd),"\'%s\'",dbtag->tag->str);
  }
#ifdef DEBUG_MODE
  ErrPostEx(SEV_INFO, 0, 0,"PubStruct : cmd generated \"%s\"",
            postupdate_cmd);
#endif
exit:
  BiostrucSourceFree(bsp);
  AsnIoMemClose(aimp);
  return mmdb_id;
err:
  ErrPostEx(SEV_FATAL, CTX_NCBIASN1, 81,
            "PubStruct : asn data error");
  return 0;
}

static Int4 LIBCALLBACK
pubstruct_closeasn(Pointer ptr,int commit)
{
  DB_stream_t *db = (DB_stream_t *)ptr;
  CS_COMMAND PNTR cmd = db->clu.ctcmd;

  switch(db->action)
    { 
    case PS_READ:
      break;
    case PS_STORE:
      return commit;
    case PS_NEW:
    case PS_UPDATE:
      {
        CS_INT          count,restype;
        Int4            status;
        CS_RETCODE      retcode;
        struct ps_chunk *piece;
        Int4            expectation;
        
        char            postupdate_cmd[1024];
        CS_INT          mmdb;
        CS_CHAR         pdb[4];
        
        if(!commit)
          goto errexit;

        postupdate_cmd[0]=0;
        if (db->action == PS_NEW)
          if(pubstruct_parseasn(postupdate_cmd,db->acc,db->top->data,db->top->len)==0)
            goto errexit;

        db->action = PS_STORE;

        assert(strlen(postupdate_cmd)< sizeof(postupdate_cmd));
        {
          struct ps_chunk *top; 
          fci_t            compr;
          int              cache_size       = 1024;
          
          top = db->top;
          db->top = db->bottom = NULL;
          
          expectation=db->iodesc.total_txtlen;
          db->iodesc.total_txtlen = 0;
          
          compr = compressor_open(dbio_open(db),30*1024,0);
          
          cache_size = 1024;
          while (top)
            {
              Int4   len;
              Int4   len1;
              piece = top ;
              expectation -= piece->len;
              if ( expectation<0 )
                goto errexit;
              
              for (piece->start = 0; piece->start < piece->len; piece->start += len1)
                {
                  len1 = piece->len - piece->start;
                  if (len1>cache_size)
                    len1 = cache_size;
                  len = compr->proc_buf(compr->data, (CharPtr)piece->data+piece->start,len1);
                  if (len!= len1)
                    {
                      compr->close(compr->data,0);
		      MemFree(compr);
                      goto errexit;
                    }
                  cache_size *=2;
                  if (cache_size > piece->len)
                    cache_size = piece->len;
                }
              top = top->next;
              MemFree(piece->data);
              MemFree(piece);
            }
          compr->close(compr->data,1);
	  MemFree(compr);
        }
        /* flush data to DB */
        expectation=db->iodesc.total_txtlen;
        RC ( ct_data_info(cmd,CS_SET,CS_UNUSED,&db->iodesc) );
        
        while (db->top)
          {
            piece = db->top ;
            expectation -= piece->len;
            if ( expectation<0 )
              goto errexit;
            db->top = db->top->next;
            RC(ct_send_data(cmd, (CS_VOID*)piece->data,(CS_INT)piece->len));
            MemFree(piece->data);
            MemFree(piece);
          }
        RC ( ct_send(cmd) );
        RC ( ct_results(cmd,&restype));
        if (restype == CS_CMD_FAIL)
          goto errexit;
        assert(restype == CS_PARAM_RESULT);
        RC (ct_cancel(NULL,cmd, CS_CANCEL_ALL));
        if (strlen(postupdate_cmd)>0)
          if(!CTLibSimpleSQL_Ex(cmd,postupdate_cmd))
            goto errexit;
        if(!CTLibSimpleSQL_Ex(cmd,"commit transaction"))
          goto errexit;
        break;
      default:
        ErrPostEx(SEV_FATAL, 0, 0,"Internal error at %s:%d: action = %d",
                  __FILE__,__LINE__,db->action);
      }
    }
  pubstruct_db_close(db);
  return 1;
  
errexit:
  ErrPostEx(SEV_FATAL, 0, 0,"PubStruct update unsuccessfull");
  CTLibSimpleSQL_Ex(cmd,"roolback transaction");
  pubstruct_db_close(db);
  return 0;
}

/*****************************************************************************
*****************************************************************************/
static AsnIoPtr LIBCALL
pubstruct_openasn (DB_stream_t *db, CS_INT *accp,int *expectation,int mmdb_id)
{
  AsnIoPtr        aip     = NULL;
  
  CS_COMMAND PNTR cmd;
  CS_INT          count,restype;
  Int4            status;
  CS_RETCODE      retcode;

  char            buf[1024];

  if (accp)
    {
      assert(*accp> 0);
      assert (!db->open_server);
      db->acc = *accp;
      if (db->action == PS_READ)
        sprintf(buf,"exec id_get_asn 0,%d,10,0,0 ",(int)db->acc);
      else
        sprintf(buf,"select blob from Struct where acc = %d",(int)db->acc);
    }
  else
    {
      assert (mmdb_id > 0);
      assert (db->action == PS_READ);
      db->acc = 0;
      sprintf(buf,"exec id_get_asn %d,0,10,0,0 ",mmdb_id);
    }
  cmd = db->clu.ctcmd;
#ifdef DEBUG_MODE
  ErrPostEx(SEV_INFO,  ERR_SYBASE, 0,"execute(%s)",buf);
#endif
  RC ( ct_command(cmd,CS_LANG_CMD,(Pointer)buf,CS_NULLTERM,CS_UNUSED) );
  RC ( ct_send(cmd) );
  RC ( ct_results(cmd,&restype));
  typeres(restype);

  /* skip 'exec' status line */
  if (db->action == PS_READ)
    {
      count = 0;
      while(1)
        {
          if (restype == CS_ROW_RESULT || restype == CS_PARAM_RESULT)
            {
              count ++;
              switch (count)
                {
                case 1: /* get_asnprop */ break ;
                case 2: /* asn */
                  goto loopexit;
                }
            }
          if (restype == CS_STATUS_RESULT)
            goto errexit;
          ct_cancel(NULL,db->clu.ctcmd, CS_CANCEL_CURRENT);
          RC ( ct_results(cmd,&restype));
          typeres(restype);
        }
    }
loopexit:

  /* we are ready to read the blob */
  assert(restype == CS_ROW_RESULT || restype == CS_PARAM_RESULT);
  RC ( ct_fetch(cmd,CS_UNUSED,CS_UNUSED,CS_UNUSED,&count)) ;

  retcode = ct_get_data(cmd,1,buf,(CS_INT)0,NULL);
  if(!typecode(retcode,0))
    goto errexit;

  RC ( ct_data_info(cmd,CS_GET,1,&db->iodesc) );

  if ( db->action == PS_READ )
    {
      if (expectation)
        *expectation = db->iodesc.total_txtlen;
      aip = asnio2fci_open(1,compressor_open(cacher_open(dbio_open(db),100*1024,1
                                                         ),100*1024,1
                                             )
                           );
    }
  else /* write case */
    {
      do  { ct_cancel(NULL,db->clu.ctcmd, CS_CANCEL_CURRENT); }
      while ( ct_results(db->clu.ctcmd,&restype) == CS_SUCCEED);
      RC ( ct_command(cmd,CS_SEND_DATA_CMD,NULL,CS_UNUSED,CS_COLUMN_DATA) );
      if (expectation)
        {
          assert(*expectation > 0 );
          db->iodesc.total_txtlen = *expectation;
          RC ( ct_data_info(cmd,CS_SET,CS_UNUSED,&db->iodesc) );
        }
      else
        db->iodesc.total_txtlen = 0;
      aip = asnio2fci_open(0,dbio_open(db));
    }
  return aip;
errexit:
  pubstruct_db_close(db);
  return NULL;
}

/**
 * PubStruct_closeasn closes AsnIO stream, which was open by some of functions
 * above and does some "termination procedure" which determined by thw function,
 * which opens connection.
 */

int      LIBCALL
PubStruct_closeasn(AsnIoPtr aip,int commit)
{
      return asnio2fci_close(aip,commit);
}

/**
 * PubStruct_newasn opens AsnIo stream to created new entry in the database. When
 * this stream is closed, database table's fields will be populated by data,
 * extracted from written asn. state is the only data - which is absent in asn.
 * (beside DB accession)
 * accp argument is optional - side effect return of new accession number
 */
static AsnIoPtr
pubstruct_newasn (DB_stream_t *db, int state, Int4 *accp)
{
  CS_COMMAND PNTR cmd = db->clu.ctcmd;
  char buffer[1024];
  CS_INT acc;

  if (!accp)
    accp = &acc;
  sprintf(buffer,"exec new_struct %d",state);
  
  if (CTlibSingleValueSelect(cmd,buffer,accp,sizeof(Int4)))
    return pubstruct_openasn (db, accp,NULL,0) ;  /* SUCCESSFULL exit */

  /* FAILURE exit */
  ErrPostEx(SEV_FATAL, 0, 0,"PubStruct insert unsuccessfull");
  CTLibSimpleSQL_Ex(cmd,"roolback transaction");
  pubstruct_db_close(db);
  ErrShow();
  return NULL;
}


AsnIoPtr LIBCALL
PubStruct_newasn (char *server,int state, Int4 *accp)
{
  DB_stream_t     *db = pubstruct_db_open(server,PS_NEW);
  return (db?pubstruct_newasn (db,state,accp):NULL);
}

/**
 * PubStruct_readasn opens AsnIo stream for reading asn found by accession number
 */

AsnIoPtr LIBCALL
PubStruct_readasn    (char *server,Int4 acc)
{
  DB_stream_t     *db = pubstruct_db_open(server,PS_READ);
  if (db)
    return pubstruct_openasn (db, &acc,NULL,0);
  return NULL;
}

/**
 * PubStruct_viewasn opens AsnIo stream for reading indexed asn found by mmdb_id
 */

AsnIoPtr LIBCALL
PubStruct_viewasn    (char *server,Int4 mmdbid)
{
  DB_stream_t     *db = pubstruct_db_open(server,PS_READ);
  if (db)
    return pubstruct_openasn (db, NULL,NULL,mmdbid);
  return NULL;
}

/**
 * PubStruct_updateasn opens AsnIo stream for updating existing asn. asn identified
 * by accession number. after updating the 'state' of the data become 'newstate'. 
 */
AsnIoPtr LIBCALL
PubStruct_updateasn  (char *server,Int4 acc, int newstate)
{
  DB_stream_t     *db; 
  CS_COMMAND PNTR cmd; 
  char buf[1024];

  db  = pubstruct_db_open(server,PS_UPDATE);
  if(!db)
    return NULL;
  
  cmd = db->clu.ctcmd;
  sprintf(buf,"exec push_struct_by_ticket %d,%d",acc,newstate);
  if(CTLibSimpleSQL_Ex(cmd,buf))
    return pubstruct_openasn (db, &acc,NULL,0) ;

  /* FAIL way */
  pubstruct_db_close(db);
  ErrShow();
  return NULL;
}

/**
 * PubStruct_removeasn suppress given asn.
 */
int      LIBCALL
PubStruct_removeasn  (char *server,Int4 acc)
{
  DB_stream_t     *db;
  CS_COMMAND PNTR cmd;
  char buf[1024];
  int rc = 0;
  
  db = pubstruct_db_open(server,PS_UPDATE);
  if(!db)
    return 0;
  cmd= db->clu.ctcmd;
  sprintf(buf,"exec rm_struct %d",acc);
  if (CTLibSimpleSQL_Ex(cmd,buf))
    if(CTLibSimpleSQL_Ex(cmd,"commit transaction"))
      rc = 1;
  pubstruct_db_close(db);
  return rc;
}

/**
 * File reading wrappers. ( a bit optimized )
 *
 * PubStruct_load read given asn from file stream and put it in database
 * (using ..._newasn) returns accession number if everything ok. or 0 in
 * case of fail
 */

int LIBCALL
PubStruct_load(FILE *infile, int state_out, char *server)
{
  AsnIoPtr aip;
  Int4 acc;
  DB_stream_t *db;
  int    old_asn = 1 ;

  if(!infile)
    return 0;
  db = pubstruct_db_open(server,PS_NEW);
  if(!db)
    return 0;
  assert(db->bottom==NULL);

  while (!feof(infile)) /* read data */
    {
      struct ps_chunk *piece;
      int len;

      if (db->bottom)
        piece = db->bottom;
      else
        db->bottom = db->top = piece = new_piece(0x4000, 1024L*1024L);
      if (piece->size == piece->len)
        db->bottom =  db->bottom->next = piece = new_piece(0x4000, 1024L*1024L);
      
      assert(piece);
      
      len = fread(piece->data + piece->len,1,piece->size - piece->len,infile);
      piece->len += len;
      db->iodesc.total_txtlen += len;
    }
  if ( state_out >=0 )
    { /* check on load */
      int  mmdb_id;
      Int4 acc1;
      
      mmdb_id = pubstruct_parseasn(NULL,0,db->top->data,db->top->len);
      if ( mmdb_id <= 0 )
        goto err;
      acc = PubStruct_lookup (server,mmdb_id,-state_out-1);
      if ( acc <= 0 )
        old_asn = 0;
      else
        {
          struct ps_chunk *piece = db->top;
          int    pos = 0;
          Int2   len;
          Int2   len1;
          
          aip = PubStruct_readasn ( server, acc );
          
          if(!aip)
            goto err;
          while(piece)
            {
              char buf[0x4000];
              
              assert(pos>=0 && pos <= piece->len);
              len1 = sizeof(buf);
              if ( piece->len - pos < len1)
                len1 = piece->len - pos;
              len = aip->readfunc(aip->iostruct, buf, len1);
              if (len <=0 )
                break;
              /* compare data */
              assert(pos>=0 && pos <= piece->len);
              assert(len <= piece->len - pos);
              if(memcmp(buf,piece->data+pos, len)!=0)
                {
                  old_asn = 0 ;
                  break;
                }
              pos+= len;
              if (piece->len == pos)
                {
                  piece = piece -> next ;
                  pos = 0;
                }
              if (len < len1)
                {
                  old_asn = 0;
                  break;
                }
            }
          PubStruct_closeasn (aip,0);
          aip = 0;
        }
    }
  if(old_asn && state_out >=0 )
    {
      acc = -acc ;
      printf("skipped...      ");
      pubstruct_db_close(db);
    }
  else
    {
      int txt = db->iodesc.total_txtlen;
      aip = pubstruct_newasn (db,(state_out>=0?state_out:-state_out-1),&acc);
      db->iodesc.total_txtlen = txt;
      if (!PubStruct_closeasn (aip,1))
        return 0;
    }
  return acc;
err:
  if(aip)
    PubStruct_closeasn (aip,0);
  pubstruct_db_close(db);
  return 0;
}

/**
 * PubStruct_download read given asn from DB and dump it to file stream
 */
int LIBCALL
PubStruct_download(char *server, Int4 acc, Int4 mmdb, FILE *outfile)
{
  char buf[0x4000];
  AsnIoPtr aip;

#ifdef PURIFY
  purify_new_inuse();
  purify_new_leaks();
#endif
  if(!outfile)
    return 0;
  if (acc>0)
    aip = PubStruct_readasn (server, acc);
  else
    aip = PubStruct_viewasn (server, mmdb);
  
  if(!aip)
    return 0;
  while(1)
    {
      int len;
      int len1;
      len = aip->readfunc(aip->iostruct, buf, sizeof(buf));
      if (len <=0 )
        break;
      len1 = fwrite(buf,1,len,outfile);
      assert(len == len1);
      if (len < sizeof(buf))
        break;
    }
  return PubStruct_closeasn (aip,1);
}

/**
 * PubStruct_lookup transforms mmdb and state into accession number (return value)
 * the meaning of state is as follows.
 * state = 0 : Production data
 * state > 0 : intermediate stages of asn assembling. "up to user"
 * state < 0 : request for Struct.state <= abs('state')
 */
NLM_EXTERN Int4     LIBCALL
PubStruct_lookup(char *server,Int4 mmdb,int state)
{
  DB_stream_t     *db;
  CS_COMMAND PNTR cmd;
  CS_INT          count,restype;
  CS_RETCODE      retcode;
  Int4            status;
  char buf[1024];
  CS_INT acc;
  
  db = pubstruct_db_open(server,PS_READ);
  if(!db)
    return 0;
  cmd = db->clu.ctcmd;
  sprintf(buf,"exec id_find_gi %d,%d",(int)mmdb,(int)state);
  
  ErrPostEx(SEV_INFO, 0, 0,buf);
#ifdef DEBUG_MODE
  ErrPostEx(SEV_INFO,  ERR_SYBASE, 0,"execute(%s)",buf);
#endif
  RC ( ct_command(cmd,CS_LANG_CMD,(Pointer)buf,CS_NULLTERM,CS_UNUSED) );
  RC ( ct_send(cmd) );
  RC ( ct_results(cmd,&restype));

  /* skip 'exec' status line */
  typeres(restype);
  if (restype == CS_STATUS_RESULT)
    {
      ct_cancel(NULL,db->clu.ctcmd, CS_CANCEL_CURRENT);
      RC ( ct_results(cmd,&restype));
      typeres(restype);
    }

  if (!(restype == CS_ROW_RESULT || restype == CS_PARAM_RESULT))
    goto errexit;
  count = 0;
  RC1 ( ct_fetch(cmd,CS_UNUSED,CS_UNUSED,CS_UNUSED,&count),2) ;
  acc = 0;
  if (count==1)
    {    
      retcode = ct_get_data(cmd,1,buf,(CS_INT)sizeof(buf),NULL);
      retcode = ct_get_data(cmd,2,&acc,(CS_INT)sizeof(acc),NULL);
      if(!typecode(retcode,0))
        goto errexit;
    }
  pubstruct_db_close(db);
  return acc;  /* SUCCESSFULL exit */

errexit:
  /* FAILURE exit */
  ErrPostEx(SEV_FATAL, 0, 0,"PubStruct lookup unsuccessfull");
  CTLibSimpleSQL_Ex(cmd,"roolback transaction");
  pubstruct_db_close(db);
  ErrShow();
  return NULL;
}

/**
 * PubStruct_pdb2mmdb make a fast lookup for given pdb code
 */

Int4     LIBCALL
PubStruct_pdb2mmdb(char *server,CharPtr pdb)
{
  char buf[1024];
  DB_stream_t     *db;
  CS_INT          mmdb_id = 0  ;
  
  db = pubstruct_db_open(server,PS_READ);
  if(!db)
    return 0;
  sprintf(buf,"exec pdb2mmdb '%s'",pdb);
  if (!CTlibSingleValueSelect(db->clu.ctcmd,buf,&mmdb_id,sizeof(mmdb_id)))
    mmdb_id = 0;
  
  pubstruct_db_close(db);
  return mmdb_id;
}

/**********************************************************************/
