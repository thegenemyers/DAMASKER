/*******************************************************************************************
 *
 *  Filter interface for the datandem.
 *
 *  Author:  Gene Myers
 *  Date  :  March 2016
 *
 ********************************************************************************************/

#ifndef _TANDEM

#define _TANDEM

#include "DB.h"
#include "align.h"

extern int    VERBOSE;
extern int    MINOVER;
extern char  *SORT_PATH;

int Set_Filter_Params(int kmer, int binshift, int hitmin, int nthreads); 

void Match_Self(char *aname, HITS_DB *ablock, Align_Spec *settings);

#endif
