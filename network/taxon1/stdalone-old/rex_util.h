/*******************************************
 * File: rex_util.h
 * Description: regular expression searcher
 */

#ifndef REX_UTIL_H_DONE
#define REX_UTIL_H_DONE

#include <string.h>
#include <stdio.h>
#include <regex.h>

typedef struct re_pattern_buffer * rex_handler;

#define REX_NO_CASE 0x1
#define REX_NO_SPACE 0x2

rex_handler rex_setExpr(char* expr, int mode);
int rex_cmp(rex_handler re, char* str);

#endif
