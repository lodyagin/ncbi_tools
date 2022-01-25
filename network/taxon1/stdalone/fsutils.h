/***************************
 * File: fsutils.h
 * Description: header for fsutils.c
 */

#ifndef FSUTILS_H_DONE
#define FSUTILS_H_DONE

#include <stdlib.h>
#include <stdio.h>
#include <ncbi.h>

CharPtr fs_getString(FILE* f, CharPtr marker);
Int2 fs_getInt2(FILE* f, CharPtr marker);
Int4 fs_getInt4(FILE* f, CharPtr marker);
Int4 fs_putString(FILE* f, CharPtr str, CharPtr marker);
Int4 fs_putInt2(FILE* f, Int2 val, CharPtr marker);
Int4 fs_putInt4(FILE* f, Int4 val, CharPtr marker);

#define fs_getInt(f, m) fs_getInt2(f, m)
#define fs_getLong(f, m) fs_getInt4(f, m)
#define fs_putInt(f, v, m) fs_putInt2(f, v, m)
#define fs_putLong(f, v, m) fs_putInt4(f, v, m)

#define fsb_getBytes(f, n, b) fread(b, n, 1, f)
#define fsb_getInt1(f) getc(f)

Int2 fsb_getInt2(FILE* f);
Int4 fsb_getInt4(FILE* f);
CharPtr fsb_getString(FILE* f);

#define fsb_putBytes(f, n, b) (fwrite(b, n, 1, f)*(n))
#define fsb_putInt1(f, v) ((putc(v, f) != EOF) ? 1 : 0)

Int4 fsb_putInt2(FILE* f, Int2 v);
Int4 fsb_putInt4(FILE* f, Int4 v);
Int4 fsb_putString(FILE* f, CharPtr str);

#define fsb_getInt(f) fsb_getInt2(f)
#define fsb_getLong(f) fsb_getInt4(f);

CharPtr fs_fileName(CharPtr path, CharPtr name);

#endif
