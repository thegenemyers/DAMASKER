/*********************************************************************************************\
 *
 *  Produce a script to compute tandem repeat masks for all blocks or range of blocks of a DB
 *
 *  Author:  Gene Myers
 *  Date  :  March 24, 2016
 *
 *********************************************************************************************/
 
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <ctype.h>
#include <unistd.h>
#include <sys/types.h>
#include <sys/stat.h>

#include "DB.h"
#include "tandem.h"

#undef  LSF  //  define if want a directly executable LSF script

static char *Usage[] =
  { "[-v] [-k<int(12)>] [-w<int(4)>] [-h<int(35)>]",
    "     [-e<double(.70)] [-l<int(500)>] [-s<int(100)]",
    "     <reads:db|dam> [<first:int>[-<last:int>]"
  };

#define LSF_TAND "bsub -q medium -n 4 -o TANDEM.out -e TANDEM.err -R span[hosts=1] -J tandem#%d"
#define LSF_SORT "bsub -q short -n 12 -o SORT.TAN.out -e SORT.TAN.err -R span[hosts=1] -J sort#%d"
#define LSF_MERGE \
    "bsub -q short -n 12 -o MERGE%d.REP%d.out -e MERGE%d.REP%d.err -R span[hosts=1] -J merge#%d"
#define LSF_MASK \
          "bsub -q short -n 12 -o MASK.TAN.out -e MASK.TAN.err -R span[hosts=1] -J masktan#%d"

int main(int argc, char *argv[])
{ int   nblocks;
  int   useblock;
  int   fblock, lblock;
#ifdef LSF
  int   jobid;
#endif

  char *pwd, *root;

#define BUNIT  4

  int    VON;
  int    WINT, HINT, KINT, SINT, LINT;
  double EREL;

  { int    i, j, k;         //  Process options
    int    flags[128];
    char  *eptr;

    ARG_INIT("HPC.TANmask")

    KINT  = 12;
    WINT  = 4;
    HINT  = 35;
    EREL  = 0.;
    LINT  = 500;
    SINT  = 100;

    j = 1;
    for (i = 1; i < argc; i++)
      if (argv[i][0] == '-')
        switch (argv[i][1])
        { default:
            ARG_FLAGS("v");
            break;
          case 'e':
            ARG_REAL(EREL)
            if (EREL < .7 || EREL >= 1.)
              { fprintf(stderr,"%s: Average correlation must be in [.7,1.) (%g)\n",Prog_Name,EREL);
                exit (1);
              }
            break;
          case 'h':
            ARG_POSITIVE(HINT,"Hit threshold (in bp.s)")
            break;
          case 'k':
            ARG_POSITIVE(KINT,"K-mer length")
            break;
          case 'l':
            ARG_POSITIVE(LINT,"Minimum ovlerap length")
            break;
          case 's':
            ARG_POSITIVE(SINT,"Trace spacing")
            break;
          case 'w':
            ARG_POSITIVE(WINT,"Log of bin width")
            break;
        }
      else
        argv[j++] = argv[i];
    argc = j;

    VON = flags['v'];

    if (argc < 2 || argc > 3)
      { fprintf(stderr,"Usage: %s %s\n",Prog_Name,Usage[0]);
        fprintf(stderr,"       %*s %s\n",(int) strlen(Prog_Name),"",Usage[1]);
        fprintf(stderr,"       %*s %s\n",(int) strlen(Prog_Name),"",Usage[2]);
        exit (1);
      }
  }

  //  Make sure DB exists and is partitioned, get number of blocks in partition

  pwd = PathTo(argv[1]);
  if (strcmp(argv[1]+(strlen(argv[1])-4),".dam") == 0)
    root = Root(argv[1],".dam");
  else
    root = Root(argv[1],".db");

  { int i, nfiles;
    FILE *dbvis;

    dbvis = fopen(Catenate(pwd,"/",root,".dam"),"r");
    if (dbvis == NULL)
      { dbvis = Fopen(Catenate(pwd,"/",root,".db"),"r");
        if (dbvis == NULL)
          exit (1);
      }

    if (fscanf(dbvis,"files = %d\n",&nfiles) != 1)
      SYSTEM_ERROR
    for (i = 0; i < nfiles; i++)
      { char buffer[30001];

        if (fgets(buffer,30000,dbvis) == NULL)
          SYSTEM_ERROR
      }

    useblock = 1;
    if (fscanf(dbvis,"blocks = %d\n",&nblocks) != 1 || nblocks == 1)
      { useblock = 0;
        nblocks  = 1;
      }
  }

  //  Set range fblock-lblock checking that DB.<fblock-1>.las exists & DB.<fblock>.las does not

  { char *eptr, *fptr;
    FILE *file;

    if (argc == 3)
      { fblock = strtol(argv[2],&eptr,10);
        if ((*eptr != '\0' && *eptr != '-') || eptr <= argv[2])
          { fprintf(stderr,"%s: final argument '%s' does not start with an integer\n",
                           Prog_Name,argv[2]);
            exit (1);
          }
        useblock = 1;
        if (*eptr == '-')
          { lblock = strtol(eptr+1,&fptr,10);
            if (*fptr != '\0' || fptr <= eptr+1)
              { fprintf(stderr,"%s: second part of range '%s' is not an integer\n",
                               Prog_Name,eptr+1);
                exit (1);
              }
          }
        else
          lblock = fblock;
        if (fblock < 1 || lblock > nblocks || fblock > lblock)
          { fprintf(stderr,"%s: range %d-%d is empty or out of bounds\n",Prog_Name,fblock,lblock);
            exit (1);
          }
      }
    else
      { fblock = 1;
        lblock = nblocks;
      }

    if (fblock > 1)
      { file = fopen(Catenate(pwd,"/.",root,Numbered_Suffix(".",fblock-1,".tan.anno")),"r");
        if (file == NULL)
          { fprintf(stderr,"%s: File %s/.%s.%d.tan.anno should already be present!\n",
                           Prog_Name,pwd,root,fblock-1);
            exit (1);
          }
        else
          fclose(file);
      }
    file = fopen(Catenate(pwd,"/.",root,Numbered_Suffix(".",fblock,".tan.anno")),"r");
    if (file != NULL)
      { fprintf(stderr,"%s: File %s/%s.%d.tan.anno should not yet exist!\n",
                       Prog_Name,pwd,root,fblock);
        exit (1);
      }
  }

  { int njobs;
    int i, j, k;
    int usepath;

    //  Produce all necessary datandem jobs ...

    usepath = (strcmp(pwd,".") != 0);

    njobs = (lblock - fblock)/BUNIT + 1;

    printf("# Datander jobs (%d)\n",njobs);

#ifdef LSF
    jobid = 1;
#endif
    for (i = fblock; i <= lblock; i += BUNIT)
      {
#ifdef LSF
        printf(LSF_TAND,jobid++);
        printf(" \"");
#endif
        printf("datander");
        if (VON)
          printf(" -v");
        if (KINT != 12)
          printf(" -k%d",KINT);
        if (WINT != 4)
          printf(" -w%d",WINT);
        if (HINT != 35)
          printf(" -h%d",HINT);
        if (EREL > .1)
          printf(" -e%g",EREL);
        if (LINT != 500)
          printf(" -l%d",LINT);
        if (SINT != 100)
          printf(" -s%d",SINT);
        j = i+BUNIT;
        if (j > lblock+1)
          j = lblock+1;
        for (k = i; k < j; k++)
          if (useblock)
            if (usepath)
              printf(" %s/%s.%d",pwd,root,k);
            else
              printf(" %s.%d",root,k);
          else
            if (usepath)
              printf(" %s/%s",pwd,root);
            else
              printf(" %s",root);
#ifdef LSF
        printf("\"");
#endif
        printf("\n");
      }

    //  ... and then all the sort & merge jobs for each block

    printf("# Sort & merge jobs (%d)\n", lblock - fblock + 1);

#ifdef LSF
    jobid = 1;
#endif
    for (i = fblock; i <= lblock; i++)
      {
#ifdef LSF
        printf(LSF_SORT,jobid++);
        printf(" \"");
#endif
        printf("LAsort");
        if (VON)
          printf(" -v");
        for (k = 0; k < NTHREADS; k++)
          if (useblock)
            printf(" %s.%d.T%d",root,i,k);
          else
            printf(" %s.T%d",root,k);
        printf(" && LAmerge");
        if (VON)
          printf(" -v");
        if (useblock)
          printf(" %s.T.%d",root,i);
        else
          printf(" %s.T",root);
        for (k = 0; k < NTHREADS; k++)
          if (useblock)
            printf(" %s.%d.T%d.S",root,i,k);
          else
            printf(" %s.T%d.S",root,k);
#ifdef LSF
          printf("\"");
#endif
          printf("\n");
        }

    //  Check .las (option)

    printf("# Check all .las files (optional but recommended)\n");

    for (i = fblock; i <= lblock; i++)
      { printf("LAcheck -vS");
        if (usepath)
          printf(" %s/%s",pwd,root);
        else
          printf(" %s",root);
        if (useblock)
          printf(" %s.T.%d\n",root,i);
        else
          printf(" %s.T\n",root);
        }

    //  Clean up (optional)

    printf("# Cleanup all intermediate .las files (optional)\n");

    for (i = fblock; i <= lblock; i++)
      { printf("rm");
        for (k = 0; k < NTHREADS; k++)
          if (useblock)
            printf(" %s.%d.T%d.S.las %s.%d.T%d.las",root,i,k,root,i,k);
          else
            printf(" %s.T%d.S.las %s.T%d.las",root,k,root,k);
          printf("\n");
        }

    //  Finish with MASKtan

    printf("# TANmask jobs (%d)\n",njobs);

#ifdef LSF
    jobid = 1;
#endif
    for (i = fblock; i <= lblock; i += BUNIT)
      {
#ifdef LSF
        printf(LSF_MASK,jobid++);
        printf(" \"");
#endif
        printf("TANmask");
        if (VON)
          printf(" -v");
        if (LINT != 500)
          printf(" -l%d",LINT);
        if (usepath)
          printf(" %s/%s",pwd,root);
        else
          printf(" %s",root);
        j = i+BUNIT;
        if (j > lblock+1)
          j = lblock+1;
        for (k = i; k < j; k++)
          if (useblock)
            printf(" %s.T.%d",root,k);
          else
            printf(" %s.T",root);
#ifdef LSF
        printf("\"");
#endif
        printf("\n");
      }

    //  Clean up (optional)

    printf("# Cleanup all T.las files (optional)\n");

    for (i = fblock; i <= lblock; i += BUNIT)
      { printf("rm");
        j = i+BUNIT;
        if (j > lblock+1)
          j = lblock+1;
        for (k = i; k < j; k++)
          if (useblock)
            printf(" %s.T.%d.las",root,k);
          else
            printf(" %s.T.las",root);
        printf("\n");
      }
  }

  printf("# Once all the .tan masks have been computed for every block\n");
  printf("#   you should call 'Catrack' to merge them, and then you should\n");
  printf("#   remove the individual block tracks\n");

  free(root);
  free(pwd);

  exit (0);
}
