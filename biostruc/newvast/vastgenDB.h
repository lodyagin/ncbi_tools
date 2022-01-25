#ifndef _VASTGEN_
#define _VASTGEN_

#include "/netopt/structure/include/dart.h"
#include "/netopt/structure/include/dartutil.h"
#include "/home/he/proj/VastSrv/PubVastApi.h" 
#include "/home/he/proj/MmdbSrv/PubStructApi.h"
#include "../../../bin/utilisCJ.h"

#define MaxSeqImgSize   620
#define GraphWidth      770

#define MAX_TBUFF	8192
#define LOG10_500       2.69897    /* -log10(500); database size correction */
#define LOG_10          2.302585        /* log(10.0) */
#define ASP_SCALE_FACTOR        10000
#define ParURL 		"\"%s%s?sdid=%d&sort=%d&"
#define PageSubsetURL	"doclistpage=%d&subset=%d&presubset=%d"
#define	VSParURL	"\"%s%s?chaindom=%d&sort=%d&"

#endif
