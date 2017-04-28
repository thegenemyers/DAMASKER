/*******************************************************************************************
 *
 *  MASKrep takes a .las file and builds a .rep mask of every interval of a read that
 *    is covered by -c or more LA's.
 *
 *  Author:  Gene Myers
 *  Date  :  March 27 2016
 *
 *******************************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>

#include "DB.h"
#include "align.h"

#ifdef HIDE_FILES
#define PATHSEP "/."
#else
#define PATHSEP "/"
#endif

static char *Usage = "[-v] [-m<track(rep)] -c<int> <source:db> <overlaps:las> ...";

#undef  DEBUG_BLOCKS
#undef  DEBUG_GAP_MERGE


//  Partition Constants

#define TINY_HILL   3    //  Minimum peak to trough coverage between two repeat intervals
#define MIN_HILL   10    //  Minimum peak to trough coverage between two overlapping rep.int.s 
#define MAX_OVL    20    //  Maximum overlap between two repeat intervals
#define PEEL_BACK 300    //  First and last PEEL_BACK bases don't count for coverage


//  Global Data Structures

static int VERBOSE;
static int MIN_COVER;

static HITS_DB _DB, *DB = &_DB;   //  Data base

static int DB_PART;               //  Input is an Overlap block
static int DB_FIRST;              //     for reads DB_FIRST to DB_LAST-1
static int DB_LAST;

static int TRACE_SPACING;         //  Trace spacing (from .las file)
static int TBYTES;                //  Bytes per trace segment (from .las file)

static FILE   *MSK_AFILE;         //  .rep.anno
static FILE   *MSK_DFILE;         //  .rep.data
static int64   MSK_INDEX;         //  Current index into .rep.data file as it is being written

//  Statistics

static int64 nreads, totlen;
static int64 nmasks, masked;


/*******************************************************************************************
 *
 *  HIGH_COVERAGE INTERVAL PARTITIONER:
 *    Find intervals of high coverage separated by low or no coverage interavls.
 *    Alignments are peeled back by PEEL_BACK base pairs so that pile ups whose boundaries
 *    overlap by this much are still found.  Small dips below the coverage threshold are
 *    merged.
 *
 *******************************************************************************************/

typedef struct
  { int add;
    int pos;
  } Event;

static int EVENT_SORT(const void *l, const void *r)
{ Event *x = ((Event *) l);
  Event *y = ((Event *) r);

  if (x->pos != y->pos)
    return (x->pos - y->pos);
  return (x->add - y->add);
}

static int *blocks(Overlap *ovls, int novl, int *ptrim)
{ static int    nmax = 0;
  static Event *ev = NULL;
  static int   *trim = NULL;
  static int   *flim;

  int ecnt, ntrim;

  if (novl > nmax)
    { nmax = novl*1.2 + 1000; 
      ev    = (Event *) Realloc(ev,sizeof(Event)*2*nmax,"Reallocating event vector");
      trim  = (int *) Realloc(trim,sizeof(int)*4*nmax,"Reallocating trim vector");
      flim  = trim + 2*nmax;
    }
 
  //  Set up and sort event queue

  { int i, ab, ae;

    ecnt = 0;
    for (i = 0; i < novl; i++)
      { ab = ovls[i].path.abpos+PEEL_BACK;
        ae = ovls[i].path.aepos-PEEL_BACK;

        if (ae < ab) 
          ab = ae = (ovls[i].path.abpos + ovls[i].path.aepos)/2;

        ev[ecnt].add = 1;
        ev[ecnt].pos = ab;
        ecnt += 1;
        ev[ecnt].add = 0;
        ev[ecnt].pos = ae;
        ecnt += 1;
      }

    qsort(ev,ecnt,sizeof(Event),EVENT_SORT);
  }

  //  Compute in - out intervals (over the contracted alignment intervals) with respect to
  //    coverage depth threshold MIN_COVER and also the max (in interavls) or min (out intervals)
  //    coverage.

  { int i;
    int cov, min, max;

    ntrim = 0;
    cov   = 0;
    min   = 0;
    max   = 0;
    for (i = 0; i < ecnt; i++)
      if (ev[i].add)
        { cov += 1;
          if (cov > max)
            max += 1;
          if (cov == MIN_COVER) 
            { trim[ntrim]   = ev[i].pos-PEEL_BACK;
              flim[ntrim++] = min;
              max = MIN_COVER;
#ifdef DEBUG_BLOCKS
              printf("    In %4d\n",ev[i].pos-PEEL_BACK);
#endif
            }
#ifdef DEBUG_BLOCKS
          printf("  Add %4d (%3d)\n",ev[i].pos,cov);
#endif
        }
      else
        { if (cov == MIN_COVER)
            { trim[ntrim]   = ev[i].pos+PEEL_BACK;
              flim[ntrim++] = max;
              min = MIN_COVER;
#ifdef DEBUG_BLOCKS
              printf("    Out %4d\n",ev[i].pos+PEEL_BACK);
#endif
            }
          cov -= 1;
          if (cov < min)
            min = cov;
#ifdef DEBUG_BLOCKS
          printf("  Del %4d (%3d)\n",ev[i].pos,cov);
#endif
        }
  }

  //  Merge intervals that overlap by more than 20bp and have a small "hill" to either
  //    the left or right.

  { int i, j;

    j = 2;
    for (i = 2; i < ntrim; i += 2)
      { int deep, over;

        deep = flim[i-1] - flim[i];
        if (deep > flim[i+1] - flim[i])
          deep = flim[i+1] - flim[i];
        over = trim[i-1] - trim[i];
#ifdef DEBUG_BLOCKS
        printf("  Trim [%5d,%5d] [%3d,%3d]",trim[i],trim[i+1],flim[i],flim[i+1]);
        printf("  DEEP %4d OVER %4d",deep,over);
#endif
        if (deep <= TINY_HILL || (deep < MIN_HILL && over > MAX_OVL))
          { trim[j-1] = trim[i+1];
            if (flim[i+1] > flim[j-1])
              flim[j-1] = flim[i+1];
#ifdef DEBUG_BLOCKS
            printf(" FOLD");
#endif
          }
        else
          { trim[j]   = trim[i];
            trim[j+1] = trim[i+1];
            flim[j]   = flim[i];
            flim[j+1] = flim[i+1];
            j += 2;
          }
#ifdef DEBUG_BLOCKS
        printf("\n");
#endif
      }
    if (ntrim > 0)
      ntrim = j;
  }

  *ptrim = ntrim;
  return (trim);
}


/*******************************************************************************************
 *
 *  FORMULATE POSSIBLE REPEAT INTERVALS
 *
 *******************************************************************************************/

static void PARTITION(int aread, Overlap *ovls, int novl)
{ int   ntrim, *trim;

#if defined(DEBUG_BLOCKS) || defined(DEBUG_GAP_MERGE)
  printf("\nAREAD %d (%d)\n",aread,DB->reads[aread].rlen);
#endif

  if (novl <= 0)
    { fwrite(&MSK_INDEX,sizeof(int64),1,MSK_AFILE);
      return;
    }

  //  Merge overlap pairs that appear to have a low-quality induced gap between them

  { int   i, j, k;
    int   bread, agap, bgap;
    Path *ipath, *jpath;

#ifdef DEBUG_GAP_MERGE
    printf("\nGAPS\n");
    if (novl > 0)
      { ipath = &(ovls[0].path);
        bread = ovls[0].bread;
        printf("    %5d %5d  [%5d,%5d] [%5d,%5d]\n",
               aread,bread,ipath->abpos,ipath->aepos,ipath->bbpos,ipath->bepos);
      }
#endif

    k = 0;
    for (i = 1; i < novl; i++)
      { ipath = &(ovls[i].path);
        bread = ovls[i].bread;
#ifdef DEBUG_GAP_MERGE
        printf("    %5d %5d  [%5d,%5d] [%5d,%5d]\n",
               aread,bread,ipath->abpos,ipath->aepos,ipath->bbpos,ipath->bepos);
#endif
        for (j = k; j >= 0; j--)
          if (bread == ovls[j].bread && COMP(ovls[i].flags) == COMP(ovls[j].flags))
            { jpath = &(ovls[j].path);
              if (jpath->aepos < ipath->abpos+PEEL_BACK)
                { agap = ipath->abpos - jpath->aepos;
                  bgap = ipath->bbpos - jpath->bepos;
                  if (abs(agap-bgap) < .2*(agap+bgap) + 200)
                    { if (ipath->aepos > jpath->aepos)
                        jpath->aepos = ipath->aepos;
                      if (ipath->bepos > jpath->bepos)
                        jpath->bepos = ipath->bepos;
                      jpath->tlen += agap;
#ifdef DEBUG_GAP_MERGE
                      printf("  Span off = %d ->%4d %4d",k-j,agap,bgap);
                      printf("    %5d %5d  [%5d,%5d] [%5d,%5d]\n",
                             aread,bread,jpath->abpos,jpath->aepos,jpath->bbpos,jpath->bepos);
#endif
                      break;
                    }
                }
            }
          else
            { j = -1;
              break;
            }
        if (j < 0)
          { k += 1;
            ovls[k] = ovls[i];
          }
      }
    novl = k+1;
  }

  //  Find the high-coverage intervals over the pair-merged alignment intervals

  trim = blocks(ovls,novl,&ntrim);

  if (VERBOSE)
    { int i;

      for (i = 0; i < ntrim; i += 2)
        masked += trim[i+1]-trim[i];
      nmasks += ntrim/2;
      nreads += 1;
      totlen += DB->reads[aread].rlen;
    }

  //  Write out the trim intervals for the read

  { int i;

    for (i = 0; i < ntrim; i += 2)
      { fwrite(trim+i,sizeof(int),1,MSK_DFILE);
        fwrite(trim+(i+1),sizeof(int),1,MSK_DFILE);
      }
    MSK_INDEX += ntrim*sizeof(int);
    fwrite(&MSK_INDEX,sizeof(int64),1,MSK_AFILE);
  }
}


  //  Read in each successive pile and call ACTION on it.  Read in the traces only if
  //   "trace" is nonzero

static int make_a_pass(FILE *input, void (*ACTION)(int, Overlap *, int), int trace)
{ static Overlap *ovls = NULL;
  static int      omax = 500;
  static uint16  *paths = NULL;
  static int      pmax = 100000;

  int64 i, j, novl;
  int   n, a;
  int   pcur;
  int   max;

  if (ovls == NULL)
    { ovls = (Overlap *) Malloc(sizeof(Overlap)*omax,"Allocating overlap buffer");
      if (ovls == NULL)
        exit (1);
    }
  if (trace && paths == NULL)
    { paths = (uint16 *) Malloc(sizeof(uint16)*pmax,"Allocating path buffer");
      if (paths == NULL)
        exit (1);
    }

  rewind(input);
  fread(&novl,sizeof(int64),1,input);
  fread(&TRACE_SPACING,sizeof(int),1,input);
  if (TRACE_SPACING <= TRACE_XOVR)
    TBYTES = sizeof(uint8);
  else
    TBYTES = sizeof(uint16);

  Read_Overlap(input,ovls);
  if (trace)
    { if (ovls[0].path.tlen > pmax)
        { pmax  = 1.2*(ovls[0].path.tlen)+10000;
          paths = (uint16 *) Realloc(paths,sizeof(uint16)*pmax,"Expanding path buffer");
          if (paths == NULL) exit (1);
        }
      fread(paths,TBYTES,ovls[0].path.tlen,input);
      if (TBYTES == 1)
        { ovls[0].path.trace = paths;
          Decompress_TraceTo16(ovls);
        }
    }
  else
    fseek(input,TBYTES*ovls[0].path.tlen,SEEK_CUR);

  if (ovls[0].aread < DB_FIRST)
    { fprintf(stderr,"%s: .las file overlaps don't correspond to reads in block %d of DB\n",
                     Prog_Name,DB_PART);
      exit (1);
    }

  pcur = 0;
  n = max = 0;
  for (j = DB_FIRST; j < DB_LAST; j++)
    { ovls[0] = ovls[n];
      a = ovls[0].aread;
      if (a != j)
        n = 0;
      else
        { if (trace)
            memmove(paths,paths+pcur,sizeof(uint16)*ovls[0].path.tlen);
          n = 1;
          pcur = ovls[0].path.tlen;
          while (1)
            { if (Read_Overlap(input,ovls+n) != 0)
                { ovls[n].aread = INT32_MAX;
                  break;
                }
              if (trace)
                { if (pcur + ovls[n].path.tlen > pmax)
                    { pmax = 1.2*(pcur+ovls[n].path.tlen)+10000;
                      paths = (uint16 *) Realloc(paths,sizeof(uint16)*pmax,"Expanding path buffer");
                      if (paths == NULL) exit (1);
                    }
                  fread(paths+pcur,TBYTES,ovls[n].path.tlen,input);
                  if (TBYTES == 1)
                    { ovls[n].path.trace = paths+pcur;
                      Decompress_TraceTo16(ovls+n);
                    }
                }
              else
                fseek(input,TBYTES*ovls[n].path.tlen,SEEK_CUR);
              if (ovls[n].aread != a)
                break;
              pcur += ovls[n].path.tlen;
              n    += 1;
              if (n >= omax)
                { omax = 1.2*n + 100;
                  ovls = (Overlap *) Realloc(ovls,sizeof(Overlap)*omax,"Expanding overlap buffer");
                  if (ovls == NULL) exit (1);
                }
            }

          if (n >= max)
            max = n;
          pcur = 0;
          for (i = 0; i < n; i++)
            { ovls[i].path.trace = paths+pcur;
              pcur += ovls[i].path.tlen;
            }
        }
      ACTION(j,ovls,n);
    }

  return (max);
}


int main(int argc, char *argv[])
{ FILE  *input;
  char  *root, *dpwd;
  char  *las, *lpwd;
  int64  novl;
  int    c;
  char  *MASK_NAME;

  //  Process arguments

  { int  i, j, k;
    int  flags[128];
    char *eptr;

    ARG_INIT("REPmask")

    MIN_COVER = -1;
    MASK_NAME = "rep";

    j = 1;
    for (i = 1; i < argc; i++)
      if (argv[i][0] == '-')
        switch (argv[i][1])
        { default:
            ARG_FLAGS("v")
            break;
          case 'c':
            ARG_POSITIVE(MIN_COVER,"Repeat coverage threshold")
            break;
          case 'm':
            MASK_NAME = argv[i]+2;
            break;
        }
      else
        argv[j++] = argv[i];
    argc = j;

    VERBOSE = flags['v'];

    if (argc < 3)
      { fprintf(stderr,"Usage: %s %s\n",Prog_Name,Usage);
        exit (1);
      }
    if (MIN_COVER <= 0)
      { fprintf(stderr,"%s: Must supply -c parameter for repeat threshold\n",Prog_Name);
        exit (1);
      }
  }

  //  Open trimmed DB

  { int status;

    status = Open_DB(argv[1],DB);
    if (status < 0)
      exit (1);
    if (status == 1)
      { fprintf(stderr,"%s: Cannot be called on a .dam index: %s\n",Prog_Name,argv[1]);
        exit (1);
      }
    if (DB->part)
      { fprintf(stderr,"%s: Cannot be called on a block: %s\n",Prog_Name,argv[1]);
        exit (1);
      }
    Trim_DB(DB);
  }

  //  Initialize statistics gathering

  if (VERBOSE)
    { int i;

      nreads = 0;
      totlen = 0;
      masked = 0;
      nmasks = 0;

      printf("\nREPmask -c%d -m%s %s",MIN_COVER,MASK_NAME,argv[1]);
      for (i = 2; i < argc; i++)
        printf(" %s",argv[i]);
      printf("\n");
    }

  //  Determine if overlap block is being processed and if so get first and last read
  //    from .db file

  dpwd = PathTo(argv[1]);
  root = Root(argv[1],".db");

  for (c = 2; c < argc; c++)
    { las  = Root(argv[c],".las");

      DB_PART  = 0;
      DB_FIRST = 0;
      DB_LAST  = DB->nreads;

      { FILE *dbfile;
        char  buffer[2*MAX_NAME+100];
        char *p, *eptr;
        int   i, part, nfiles, nblocks, cutoff, all, oindx;
        int64 size;

        p = rindex(las,'.');
        if (p != NULL)
          { part = strtol(p+1,&eptr,10);
            if (*eptr == '\0' && eptr != p+1)
              { dbfile = Fopen(Catenate(dpwd,"/",root,".db"),"r");
                if (dbfile == NULL)
			exit (1);
                if (fscanf(dbfile,DB_NFILE,&nfiles) != 1)
                  SYSTEM_ERROR
                for (i = 0; i < nfiles; i++)
                  if (fgets(buffer,2*MAX_NAME+100,dbfile) == NULL)
                    SYSTEM_ERROR
                if (fscanf(dbfile,DB_NBLOCK,&nblocks) != 1)
			SYSTEM_ERROR
				if (fscanf(dbfile,DB_PARAMS,&size,&cutoff,&all) != 3)
                  SYSTEM_ERROR
                for (i = 1; i <= part; i++)
                  if (fscanf(dbfile,DB_BDATA,&oindx,&DB_FIRST) != 2)
                    SYSTEM_ERROR
                if (fscanf(dbfile,DB_BDATA,&oindx,&DB_LAST) != 2)
                  SYSTEM_ERROR
                fclose(dbfile);
                DB_PART = part;
                *p = '\0';
              }
          }
      }

      //   Set up preliminary trimming track

      { int   len, size;
        char  ans[strlen(MASK_NAME)+7];
        char  dts[strlen(MASK_NAME)+7];

        strcpy(ans,Catenate(".",MASK_NAME,".","anno"));
        strcpy(dts,Catenate(".",MASK_NAME,".","data"));
        if (DB_PART > 0)
          { MSK_AFILE = Fopen(Catenate(dpwd,PATHSEP,root,
                                       Numbered_Suffix(".",DB_PART,ans)),"w");
            MSK_DFILE = Fopen(Catenate(dpwd,PATHSEP,root,
                                       Numbered_Suffix(".",DB_PART,dts)),"w");
          }
        else
          { MSK_AFILE = Fopen(Catenate(dpwd,PATHSEP,root,ans),"w");
            MSK_DFILE = Fopen(Catenate(dpwd,PATHSEP,root,dts),"w");
          }
        if (MSK_AFILE == NULL || MSK_DFILE == NULL)
          exit (1);

        len  = DB_LAST - DB_FIRST;
        size = 0;
        fwrite(&len,sizeof(int),1,MSK_AFILE);
        fwrite(&size,sizeof(int),1,MSK_AFILE);
        MSK_INDEX = 0;
        fwrite(&MSK_INDEX,sizeof(int64),1,MSK_AFILE);
      }

      //  Open overlap file

      lpwd = PathTo(argv[c]);
      if (DB_PART > 0)
        input = Fopen(Catenate(lpwd,"/",las,Numbered_Suffix(".",DB_PART,".las")),"r");
      else
        input = Fopen(Catenate(lpwd,"/",las,".las"),"r");
      if (input == NULL)
        exit (1);

      free(lpwd);
      free(las);

      //  Get trace point spacing information

      fread(&novl,sizeof(int64),1,input);
      fread(&TRACE_SPACING,sizeof(int),1,input);

      //  Process each read pile

      make_a_pass(input,PARTITION,1);

      fclose(MSK_AFILE);
      fclose(MSK_DFILE);
    }

  if (VERBOSE)
    { printf("\nInput:    ");
      Print_Number((int64) nreads,7,stdout);
      printf(" (100.0%%) reads     ");
      Print_Number(totlen,12,stdout);
      printf(" (100.0%%) bases\n");

      printf("Masks:    ");
      Print_Number(nmasks,7,stdout);
      printf(" (%5.1f%%) masks     ",(100.*nmasks)/nreads);
      Print_Number(masked,12,stdout);
      printf(" (%5.1f%%) bases\n",(100.*masked)/totlen);
    }

  free(dpwd);
  free(root);

  Close_DB(DB);
  free(Prog_Name);

  exit (0);
}
