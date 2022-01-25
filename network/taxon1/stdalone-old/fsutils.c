/****************************
 * File: fsutils.c
 * Description: utility functions
 */

#include "fsutils.h"

CharPtr fs_getString(FILE* f, CharPtr marker)
{
  CharPtr buff;
  int i, buffSize= 16;
  int l;

  l= strlen(marker);

  if((buff= MemNew(buffSize + 2)) == NULL) {
    return NULL;
  }

  for(buff[i= 0]= fgetc(f); feof(f) == 0; buff[++i]= fgetc(f)) {
    if(i == buffSize) {
      buffSize+= buffSize/2;
      if((buff= MemMore(buff, buffSize + 2)) == NULL) return NULL;
    }
    
    /* check for marker*/
    if(i >= (l-1)) {
      switch (l) {
      case 0: 
	continue;
      case 1:
	if(buff[i] == *marker) {
	  buff[i]= '\0';
	  break;
	}
	continue;
      case 2:
	if((buff[i-1] == marker[0]) && (buff[i] == marker[1])) {
	  buff[i-1]= '\0';
	  break;
	}
	continue;
      case 3:
	if((buff[i-2] == marker[0]) && (buff[i-1] == marker[1]) && (buff[i] == marker[2])) {
	  buff[i-2]= '\0';
	  break;
	}
	continue;
      default:
	if((buff[i-l+1] == *marker) && (strncmp(&buff[i-l+1], marker, l) == 0)) {
	  buff[i-l+1]= '\0';
	  break;
	}
	continue;
      }
      /* match marker */
      buff= MemMore(buff, strlen(buff)+1);
      return buff;
    }
  }
  /* eof matched */
  buff[i]= '\0';
  buff= MemMore(buff, i+1);
  return buff;
}

Int2 fs_getInt(FILE* f, CharPtr marker)
{
  CharPtr buff;
  int res;

  if((buff= fs_getString(f, marker)) == NULL) return 0;

  if(buff[0] == '\0') {
    res= 0;
  }
  else {
    res= atoi(buff);
  }

  MemFree(buff);
  return (Int2)res;
}

Int4 fs_getLong(FILE* f, CharPtr marker)
{
  CharPtr buff;
  long res;

  if((buff= fs_getString(f, marker)) == NULL) return 0;

  if(buff[0] == '\0') {
    res= 0;
  }
  else {
    res= atol(buff);
  }

  MemFree(buff);
  return (Int4) res;
}

Int4 fs_putString(FILE* f, CharPtr str, CharPtr marker)
{
  Int4 l1= 0, l2= 0;

  if((str != NULL) && ((l1= strlen(str)) > 0)) {
    if(FilePuts(str, f) == EOF) return -1;
  }
  if((marker != NULL) && ((l2= strlen(marker)) > 0)) {
    if(FilePuts(marker, f)) return -1;
  }
  return l1+l2;
}

Int4 fs_putInt2(FILE* f, Int2 val, CharPtr marker)
{
  char buff[32];

  sprintf(buff, "%d", val);
  return fs_putString(f, buff, marker);
}

Int4 fs_putInt4(FILE* f, Int4 val, CharPtr marker)
{
  char buff[32];

  sprintf(buff, "%d", (int)val);
  return fs_putString(f, buff, marker);
}

Int2 fsb_getInt2(FILE* f)
{
  Int2 res;

  if(fread(&res, 2, 1, f) > 0) return res;
  else return 0;
}

Int4 fsb_getInt4(FILE* f)
{
  Int4 res;

  if(fread(&res, 4, 1, f) > 0) return res;
  else return 0;
}

CharPtr fsb_getString(FILE* f)
{
  Int2 len;
  CharPtr buff;

  if((len= fsb_getInt2(f)) > 0) {
    if((buff= MemNew(len)) != NULL) {
      fread(buff, 1, len, f);
      return buff;
    }
  }
  
  return NULL;
}

Int4 fsb_putInt2(FILE* f, Int2 v)
{
  return fwrite(&v, 2, 1, f)*2;
}

Int4 fsb_putInt4(FILE* f, Int4 v)
{
  return fwrite(&v, 4, 1, f)*4;
}

Int4 fsb_putString(FILE* f, CharPtr str)
{
  Int4 len, res;
  Int2 l;

  len= strlen(str)+1;
  if(len > (32*1024 - 1)) l= 32*1024 - 1;
  else l= (Int2) len;

  if(l == 1) l= 0;

  len= l;

  res= fwrite(&l, 2, 1, f)*2;
  
  return  res + ((len > 0)? (fwrite(str, len, 1, f)*len) : 0);
}
  
static char nameBuff[1024];

CharPtr fs_fileName(CharPtr path, CharPtr name)
{
  int n;
#define PATH_SEPARATOR '/'

  if((path == NULL) || (path[0] == '\0')) return name;

  n= StringLen(path);

  StringCpy(nameBuff, path);
  if(path[n-1] != PATH_SEPARATOR) {
    nameBuff[n]= PATH_SEPARATOR;
    nameBuff[n+1]= '\0';
  }

  StringCat(nameBuff, name);

  return nameBuff;
}
