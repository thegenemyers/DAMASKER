/*******************************************************************************************
 *
 *  MASKtan takes as input a .las file of self alignments (produced by datandem) and
 *    builds a .tan mask track that encodes the union of all self-overlapping LA's (signature)
 *    of a tandem repeat) of length greater than MIN_LEN (set by -l parameter).
 *
 *  Author:  Gene Myers
 *  Date  :  March 27 2016
 *
 *******************************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <math.h>

#include "DB.h"
#include "align.h"

#ifdef HIDE_FILES
#define PATHSEP "/."
#else
#define PATHSEP "/"
#endif

static char *Usage = "[-v] [-m<track(tan)>] [-l<int(500)>] <source:db> <overlaps:las> ...";

//  Partition Constants

#define SEP_FUZZ 20


//  Global Data Structures

static int VERBOSE;
static int MIN_LEN;

static HITS_DB _DB, *DB = &_DB;   //  Data base

static int DB_PART;               //  Input is an Overlap block
static int DB_FIRST;              //     for reads DB_FIRST to DB_LAST-1
static int DB_LAST;

static int TRACE_SPACING;         //  Trace spacing (from .las file)
static int TBYTES;                //  Bytes per trace segment (from .las file)

static FILE   *TN_AFILE;          //  .tan.anno
static FILE   *TN_DFILE;          //  .tan.data
static int64   TN_INDEX;          //  Current index into .tan.data file as it is being written

//  Statistics

static int64 nreads, totlen;
static int64 nmasks, masked;


static int ISORT(const void *l, const void *r)
{ int x = *((int *) l);
  int y = *((int *) r);
  return (x-y);
}

static void TANDEM(int aread, Overlap *ovls, int novl)
{ static int  nmax = -1;
  static int *add = NULL;
  static int *del;

  int evnum;

  (void) aread;

  if (VERBOSE)
    { nreads += 1;
      totlen += DB->reads[aread].rlen;
    }

  if (novl == 0)
    { fwrite(&TN_INDEX,sizeof(int64),1,TN_AFILE);
      return;
    }

#ifdef DEBUG
  printf("\nAREAD %d:\n",aread);
#endif

  if (novl > nmax)
    { nmax = 1.2*novl + 500;
      add  = (int *) Realloc(add,4*sizeof(int)*nmax,"Allocating sort vector");
      if (add == NULL)
        exit (1);
      del = add + nmax;
    }

  //  For each overlapping LA record mask interval as an add and del event
  //    that are then sorted

  { int   i;
    Path *ipath;

    evnum = 0;
    for (i = 0; i < novl; i++)
      { ipath = &(ovls[i].path);
        if (ipath->abpos - ipath->bepos <= SEP_FUZZ)
          { if (ipath->aepos - ipath->bbpos > MIN_LEN)
              { add[evnum] = ipath->bbpos;
                del[evnum] = ipath->aepos;
                evnum += 1;
              }
          }
      }
    qsort(add,evnum,sizeof(int),ISORT);
    qsort(del,evnum,sizeof(int),ISORT);
  }

  //  Output the union of the mask intervals to the .tan track

  { int i, j, x, a;

    x = a = 0;
    i = j = 0;
    while (j < evnum)
      if (i < evnum && add[i] <= del[j])
        { if (x == 0)
            { fwrite(add+i,sizeof(int),1,TN_DFILE);
              TN_INDEX += sizeof(int);
#ifdef DEBUG
              printf("  + %5d: %3d\n",add[i],x);
#endif
              a = add[i];
            }
          x += 1;
          i += 1;
        }
      else
        { x -= 1;
          if (x == 0)
            { fwrite(del+j,sizeof(int),1,TN_DFILE);
              TN_INDEX += sizeof(int);
#ifdef DEBUG
              printf("  - %5d: %3d\n",del[j],x);
#endif
              if (VERBOSE)
                { masked += del[j]-a;
                  nmasks += 1;
                }
            }
          j += 1;
        }
#ifdef DEBUG
    printf("   %lld\n",TN_INDEX);
#endif
    fwrite(&TN_INDEX,sizeof(int64),1,TN_AFILE);
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

  int rv = Read_Overlap(input,ovls);
  if (rv)
    {
      return 0;
    }

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
            memcpy(paths,paths+pcur,sizeof(uint16)*ovls[0].path.tlen);
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
{ FILE       *input;
  char       *root, *dpwd;
  char       *las, *lpwd;
  int         c;
  char       *MASK_NAME;

  //  Process arguments

  { int  i, j, k;
    int  flags[128];
    char *eptr;

    ARG_INIT("TANmask")

    MIN_LEN   = 500;
    MASK_NAME = "tan";

    j = 1;
    for (i = 1; i < argc; i++)
      if (argv[i][0] == '-')
        switch (argv[i][1])
        { default:
            ARG_FLAGS("v")
            break;
          case 'm':
            MASK_NAME = argv[i]+2;
            break;
          case 'l':
            ARG_POSITIVE(MIN_LEN,"Minimum retained segment length")
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
  }

  //  Open trimmed DB

  { int status;

    status = Open_DB(argv[1],DB);
    if (status < 0)
      exit (1);
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

      printf("\nTANmask -l%d -m%s %s",MIN_LEN,MASK_NAME,argv[1]);
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

      //  Set up mask track

      { int len, size;
        char  ans[strlen(MASK_NAME)+7];
        char  dts[strlen(MASK_NAME)+7];

        strcpy(ans,Catenate(".",MASK_NAME,".","anno"));
        strcpy(dts,Catenate(".",MASK_NAME,".","data"));
        if (DB_PART > 0)
          { TN_AFILE = Fopen(Catenate(dpwd,PATHSEP,root,
                                      Numbered_Suffix(".",DB_PART,ans)),"w");
            TN_DFILE = Fopen(Catenate(dpwd,PATHSEP,root,
                                      Numbered_Suffix(".",DB_PART,dts)),"w");
          }
        else
          { TN_AFILE = Fopen(Catenate(dpwd,PATHSEP,root,ans),"w");
            TN_DFILE = Fopen(Catenate(dpwd,PATHSEP,root,dts),"w");
          }
        if (TN_AFILE == NULL || TN_DFILE == NULL)
          exit (1);

        len  = DB_LAST - DB_FIRST;
        size = 0;
        fwrite(&len,sizeof(int),1,TN_AFILE);
        fwrite(&size,sizeof(int),1,TN_AFILE);
        TN_INDEX = 0;
        fwrite(&TN_INDEX,sizeof(int64),1,TN_AFILE);
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

      //  Process each read pile

      make_a_pass(input,TANDEM,0);

      fclose(TN_AFILE);
      fclose(TN_DFILE);
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
