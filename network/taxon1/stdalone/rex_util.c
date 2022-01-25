/*******************************************
 * File: rex_util.c
 * Description: regular expression searcher
 */

#include "rex_util.h"
#include <string.h>
#include <stdlib.h>

rex_handler rex_setExpr(char* re_in, int mode)
{
  struct re_pattern_buffer *buf;
  int i;
  char* trn;
  char* re;


  /* replace "*" with ".*" , "+" with ".+" and "?" with ".?" */
  if((re= malloc(strlen(re_in)*2 + 1)) == NULL) return NULL;
  
  for(i= 0; *re_in != '\0'; i++, re_in++) {
    if((*re_in == '*') || (*re_in == '+') || (*re_in == '?')) {
      if((i == 0) || (re[i-1] != '\\')) {
	re[i++]= '.';
      }
    }
    re[i]= *re_in;
  }
  re[i]= '\0';

  if((buf= malloc(sizeof(struct re_pattern_buffer))) == NULL) return NULL;

  memset(buf, 0, sizeof(struct re_pattern_buffer)); 

  /* prepare re_syntax_options */
#if 1
  re_syntax_options= RE_BACKSLASH_ESCAPE_IN_LISTS | 
                     RE_CONTEXT_INDEP_ANCHORS |
                     RE_CONTEXT_INDEP_OPS |
                     RE_INTERVALS |
                     RE_NO_BK_BRACES |
                     RE_NO_BK_VBAR;
#else
  re_syntax_options= RE_SYNTAX_AWK;
#endif

  /* initialize some fields in pattern buffer */
  buf->buffer= NULL;
  buf->allocated= 0;
  buf->fastmap= NULL;
  if(mode == 0) {
    buf->translate= NULL;
  }
  else {
    if((trn= malloc(256)) != NULL) {
      for(i= 0; i != 256; i++) trn[i]= i;
      if((mode & REX_NO_CASE) != 0) {
	for(i= 1; i != 256; i++) trn[i]= toupper(i);
      }
      if((mode & REX_NO_SPACE) != 0) {
	for(i= 1; i != 256; i++) {
	  if(isspace(i)) trn[i]= ' ';
	}
      }
    }
    buf->translate= trn;
  }

  /* compile the regular expression */
  if((trn= (char *) re_compile_pattern(re, strlen(re), buf)) != NULL) {
    fprintf(stderr,"Error in regular expression: %s\n", trn);
    free(buf);
    free(re);
    return NULL;
  }
  
  free(re);
  return buf;
}
  
int rex_cmp(rex_handler re, char* str)
{
  int res;

  res= re_match(re, str, strlen(str), 0, NULL);
  return (res < 0)? 0 : 1;
}
  

