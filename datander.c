/*********************************************************************************************\
 *
 *  Find all local self-alignment between long, noisy DNA reads:
 *    Compare sequences in each supplied blocks against themselves search for local alignments
 *    of MIN_OVERLAP or more above the diagonal (A start coord > B start coord).  An output
 *    stream of 'Overlap' records (see align.h) is written in binary to a set of files, one
 *    per thread, each encoding a given found local alignment.  The -v option turns on a verbose
 *    reporting mode that reports progress and gives statistics on each major stage.
 *
 *    The filter operates by looking for a pair of diagonal bands of width 2^'s' that contain
 *    a collection of exact matching 'k'-mers between positions of a sequence, such that the total
 *    number of bases covered by 'k'-mer hits is 'h'.  k cannot be larger than 15 in the
 *    current implementation.
 *
 *    For each subject, say XXX, the program outputs a file containing LAs of the form
 *    XXX.XXX.T#.las where # is the thread that detected and wrote out the collection of LAs.
 *    For example, if NTHREAD in the program is 4, then 4 files are output for each subject block.
 *
 *  Author:  Gene Myers
 *  Date  :  March 27, 2016
 *
 *********************************************************************************************/
 
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <ctype.h>
#include <unistd.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <dirent.h>

#include <sys/param.h>
#if defined(BSD)
#include <sys/sysctl.h>
#endif

#include "DB.h"
#include "tandem.h"

static char *Usage[] =
  { "[-v] [-k<int(12)>] [-w<int(4)>] [-h<int(35)>] [-T<int(4)>] [-P<dir(/tmp)>]",
    "     [-e<double(.70)] [-l<int(500)>] [-s<int(100)>] <subject:db|dam> ...",
  };

int     VERBOSE;   //   Globally visible to tandem.c
char   *SORT_PATH;
int     MINOVER;

static int read_DB(HITS_DB *block, char *name, int kmer)
{ int i, isdam;

  isdam = Open_DB(name,block);
  if (isdam < 0)
    exit (1);

  Trim_DB(block);

  if (block->cutoff < kmer)
    { for (i = 0; i < block->nreads; i++)
        if (block->reads[i].rlen < kmer)
          { fprintf(stderr,"%s: Block %s contains reads < %dbp long !  Run DBsplit.\n",
                           Prog_Name,name,kmer);
            exit (1);
          }
    }

  Read_All_Sequences(block,0);

  return (isdam);
}

static char *CommandBuffer(char *bname)
{ static char *cat = NULL;
  static int   max = -1;
  int len;

  len = 2*strlen(bname) + 200;
  if (len > max)
    { max = ((int) (1.2*len)) + 100;
      if ((cat = (char *) realloc(cat,max+1)) == NULL)
        { fprintf(stderr,"%s: Out of memory (Making path name)\n",Prog_Name);
          exit (1);
        }
    }
  return (cat);
}

int main(int argc, char *argv[])
{ HITS_DB    _bblock;
  HITS_DB    *bblock = &_bblock;
  char       *bfile;
  char       *broot;
  Align_Spec *settings;
  int         isdam;
  char       *command;

  int    KMER_LEN;
  int    BIN_SHIFT;
  int    HIT_MIN;
  double AVE_ERROR;
  int    SPACING;
  int    NTHREADS;

  { int    i, j, k;
    int    flags[128];
    char  *eptr;
    DIR   *dirp;

    ARG_INIT("datander")

    KMER_LEN  = 12;
    HIT_MIN   = 35;
    BIN_SHIFT = 4;
    AVE_ERROR = .70;
    SPACING   = 100;
    MINOVER   = 500;    //   Globally visible to filter.c
    NTHREADS  = 4;
    SORT_PATH = "/tmp";

    j    = 1;
    for (i = 1; i < argc; i++)
      if (argv[i][0] == '-')
        switch (argv[i][1])
        { default:
            ARG_FLAGS("v")
            break;
          case 'k':
            ARG_POSITIVE(KMER_LEN,"K-mer length")
            if (KMER_LEN > 32)
              { fprintf(stderr,"%s: K-mer length must be 32 or less\n",Prog_Name);
                exit (1);
              }
            break;
          case 'w':
            ARG_POSITIVE(BIN_SHIFT,"Log of bin width")
            break;
          case 'h':
            ARG_POSITIVE(HIT_MIN,"Hit threshold (in bp.s)")
            break;
          case 'e':
            ARG_REAL(AVE_ERROR)
            if (AVE_ERROR < .6 || AVE_ERROR >= 1.)
              { fprintf(stderr,"%s: Average correlation must be in [.6,1.) (%g)\n",
                               Prog_Name,AVE_ERROR);
                exit (1);
              }
            break;
          case 'l':
            ARG_POSITIVE(MINOVER,"Minimum alignment length")
            break;
          case 's':
            ARG_POSITIVE(SPACING,"Trace spacing")
            break;
          case 'P':
            SORT_PATH = argv[i]+2;
            if ((dirp = opendir(SORT_PATH)) == NULL)
              { fprintf(stderr,"%s: -P option: cannot open directory %s\n",Prog_Name,SORT_PATH);
                exit (1);
              }
            closedir(dirp);
            break;
          case 'T':
            ARG_POSITIVE(NTHREADS,"Number of threads")
            break;
        }
      else
        argv[j++] = argv[i];
    argc = j;

    VERBOSE = flags['v'];   //  Globally declared in filter.h

    if (argc <= 1)
      { fprintf(stderr,"Usage: %s %s\n",Prog_Name,Usage[0]);
        fprintf(stderr,"       %*s %s\n",(int) strlen(Prog_Name),"",Usage[1]);
        exit (1);
      }
  }

  MINOVER *= 2;
  if (Set_Filter_Params(KMER_LEN,BIN_SHIFT,HIT_MIN,NTHREADS))
    { fprintf(stderr,"Illegal combination of filter parameters\n");
      exit (1);
    }

  /* Compare each block against itself */

  { int i;

    broot = NULL;
    for (i = 1; i < argc; i++)
      { bfile = argv[i];
        isdam = read_DB(bblock,bfile,KMER_LEN);
        if (isdam)
          broot = Root(bfile,".dam");
        else
          broot = Root(bfile,".db");

        settings = New_Align_Spec( AVE_ERROR, SPACING, bblock->freq);

        Match_Self(broot,bblock,settings);

        Free_Align_Spec(settings);
        Close_DB(bblock);

        command = CommandBuffer(broot);

        sprintf(command,"LAsort %s/%s.T*.las",SORT_PATH,broot);
        if (VERBOSE)
          printf("\n%s\n",command);
        system(command);
        sprintf(command,"LAmerge TAN.%s.las %s/%s.T*.S.las",broot,SORT_PATH,broot);
        if (VERBOSE)
          printf("%s\n",command);
        system(command);
        sprintf(command,"rm %s/%s.T*.las",SORT_PATH,broot);
        if (VERBOSE)
          printf("%s\n",command);
        system(command);
      }
  }

  exit (0);
}
