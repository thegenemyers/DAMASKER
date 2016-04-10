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

#define NTHREADS  4    //  Must be a power of 2
#define NSHIFT    2    //  log_2 NTHREADS

int Set_Filter_Params(int kmer, int binshift, int hitmin); 

void Match_Self(char *aname, HITS_DB *ablock, Align_Spec *settings);

#endif
