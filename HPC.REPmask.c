/*********************************************************************************************\
 *
 *  Produce a script to compute overlaps for all block pairs of a DB, and then sort and merge
 *    them into as many .las files as their are blocks.
 *
 *  Author:  Gene Myers
 *  Date  :  June 1, 2014
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

#define NTHREADS 4

#undef  LSF  //  define if want a directly executable LSF script

static char *Usage[] =
  { "[-vb] [-t<int>] [-w<int(6)>] [-l<int(1000)>] [-s<int(100)>]",
    "      [-M<int>] [-B<int(4)>] [-D<int( 250)>] [-m<track>]+",
    "      [-k<int(14)>] [-h<int(35)>] [-e<double(.70)>]",
    "      -g<int> -c<int> <reads:db|dam> [<block:int>[-<range:int>]"
  };

static int power(int base, int exp)
{ int i, pow;

  pow = 1;
  for (i = 0; i < exp; i++)
    pow *= base;
  return (pow);
}

#define LSF_ALIGN \
            "bsub -q medium -n 4 -o DAL.REP%d.out -e DAL.REP%d.err -R span[hosts=1] -J align#%d"
#define LSF_SORT  \
            "bsub -q short -n 12 -o SORT.REP%d.out -e SORT.REP%d.err -R span[hosts=1] -J sort#%d"
#define LSF_MERGE \
    "bsub -q short -n 12 -o MERGE%d.REP%d.out -e MERGE%d.REP%d.err -R span[hosts=1] -J merge#%d"
#define LSF_MASK \
    "bsub -q short -n 12 -o MASK.REP%d.out -e MASK.REP%d.err -R span[hosts=1] -J maskrep#%d"

int main(int argc, char *argv[])
{ int   nblocks;
  int   usepath;
  int   useblock;
  int   fblock, lblock;
#ifdef LSF
  int   jobid;
#endif

  char *pwd, *root;

  int    CINT, SPAN;
  int    DUNIT, BUNIT;
  int    VON, BON, AON, ION;
  int    WINT, TINT, HINT, KINT, SINT, LINT, MINT;
  double EREL;
  int    MMAX, MTOP;
  char **MASK;

  { int    i, j, k;         //  Process options
    int    flags[128];
    char  *eptr;

    ARG_INIT("HPC.REPmask")

    BUNIT = 4;
    DUNIT = 250;
    KINT  = 14;
    WINT  = 6;
    HINT  = 35;
    TINT  = 0;
    EREL  = 0.;
    LINT  = 1000;
    SINT  = 100;
    MINT  = -1;

    MTOP = 0;
    MMAX = 10;
    MASK = (char **) Malloc(MMAX*sizeof(char *),"Allocating mask track array");
    if (MASK == NULL)
      exit (1);

    SPAN = -1;
    CINT = -1;

    j = 1;
    for (i = 1; i < argc; i++)
      if (argv[i][0] == '-')
        switch (argv[i][1])
        { default:
            ARG_FLAGS("vbAI");
            break;
          case 'c':
            ARG_POSITIVE(CINT,"Repeat coverage threshold")
            break;
          case 'e':
            ARG_REAL(EREL)
            if (EREL < .7 || EREL >= 1.)
              { fprintf(stderr,"%s: Average correlation must be in [.7,1.) (%g)\n",Prog_Name,EREL);
                exit (1);
              }
            break;
          case 'g':
            ARG_POSITIVE(SPAN,"Block span")
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
          case 't':
            ARG_POSITIVE(TINT,"Tuple suppression frequency")
            break;
          case 's':
            ARG_POSITIVE(SINT,"Trace spacing")
            break;
          case 'w':
            ARG_POSITIVE(WINT,"Log of bin width")
            break;
          case 'm':
            if (MTOP >= MMAX)
              { MMAX = 1.2*MTOP + 10;
                MASK = (char **) Realloc(MASK,MMAX*sizeof(char *),"Reallocating mask track array");
                if (MASK == NULL)
                  exit (1);
              }
            MASK[MTOP++] = argv[i]+2;
            break;
          case 'B':
            ARG_NON_NEGATIVE(BUNIT,"Blocks per command")
            break;
          case 'D':
            ARG_NON_NEGATIVE(DUNIT,"Files per merge")
            if (DUNIT < 3)
              { fprintf(stderr,"%s: Files per merge must be at least 3 (%d)\n",
                               Prog_Name,DUNIT);
                exit (1);
              }
            break;
          case 'M':
            ARG_NON_NEGATIVE(MINT,"Memory allocation (in Gb)")
            break;
        }
      else
        argv[j++] = argv[i];
    argc = j;

    VON = flags['v'];
    BON = flags['b'];
    AON = flags['A'];
    ION = flags['I'];

    if (argc < 2 || argc > 3)
      { fprintf(stderr,"Usage: %s %s\n",Prog_Name,Usage[0]);
        fprintf(stderr,"       %*s %s\n",(int) strlen(Prog_Name),"",Usage[1]);
        fprintf(stderr,"       %*s %s\n",(int) strlen(Prog_Name),"",Usage[2]);
        fprintf(stderr,"       %*s %s\n",(int) strlen(Prog_Name),"",Usage[3]);
        exit (1);
      }
    if (SPAN < 0)
      { fprintf(stderr,"%s: Must supply -g parameter\n",Prog_Name);
        exit (1);
      }
    if (CINT < 0)
      { fprintf(stderr,"%s: Must supply -c parameter\n",Prog_Name);
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

    usepath = (strcmp(pwd,".") != 0);
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
      { file = fopen(Catenate(pwd,"/",root,Numbered_Suffix(".",fblock-1,".las")),"r");
        if (file == NULL)
          { if (usepath)
              fprintf(stderr,"%s: File %s/%s.%d.las should already be present!\n",
                             Prog_Name,pwd,root,fblock-1);
            else
              fprintf(stderr,"%s: File %s.%d.las should already be present!\n",
                             Prog_Name,root,fblock-1);
            exit (1);
          }
        else
          fclose(file);
      }
    if (useblock)
      file = fopen(Catenate(pwd,"/",root,Numbered_Suffix(".",fblock,".las")),"r");
    else
      file = fopen(Catenate(pwd,"/",root,".las"),"r");
    if (file != NULL)
      { if (usepath)
          if (useblock)
            fprintf(stderr,"%s: File %s/%s.%d.las should not yet exist!\n",
                           Prog_Name,pwd,root,fblock);
          else
            fprintf(stderr,"%s: File %s/%s.las should not yet exist!\n",Prog_Name,pwd,root);
        else
          if (useblock)
            fprintf(stderr,"%s: File %s.%d.las should not yet exist!\n",Prog_Name,root,fblock);
          else
            fprintf(stderr,"%s: File %s.las should not yet exist!\n",Prog_Name,root);
        exit (1);
      }
  }

  { int level, njobs;
    int i, j, k;
    int beg, end;

    //  Produce all necessary daligner jobs ...

    njobs = 0;
    for (i = fblock; i <= lblock; i++)
      { beg = ((i-1)/SPAN)*SPAN + 1;
        njobs += (i-beg)/BUNIT+1;
      }

    printf("# Daligner jobs (%d)\n",njobs);

#ifdef LSF
    jobid = 1;
#endif
    for (i = fblock; i <= lblock; i++)
      { int bits, base;
        int low, hgh;

        base = ((i-1)/SPAN)*SPAN + 1;
        bits = (i-base)/BUNIT+1;
        low  = base;
        for (j = 1; j <= bits; j++)
          {
#ifdef LSF
            printf(LSF_ALIGN,SPAN,SPAN,jobid++);
            printf(" \"");
#endif
            printf("daligner");
            if (VON)
              printf(" -v");
            if (BON)
              printf(" -b");
            if (AON)
              printf(" -A");
            if (ION)
              printf(" -I");
            if (KINT != 14)
              printf(" -k%d",KINT);
            if (WINT != 6)
              printf(" -w%d",WINT);
            if (HINT != 35)
              printf(" -h%d",HINT);
            if (TINT > 0)
              printf(" -t%d",TINT);
            if (EREL > .1)
              printf(" -e%g",EREL);
            if (LINT != 1000)
              printf(" -l%d",LINT);
            if (SINT != 100)
              printf(" -s%d",SINT);
            if (MINT >= 0)
              printf(" -M%d",MINT);
            for (k = 0; k < MTOP; k++)
              printf(" -m%s",MASK[k]);
            if (useblock)
              if (usepath)
                printf(" %s/%s.%d",pwd,root,i);
              else
                printf(" %s.%d",root,i);
            else
              if (usepath)
                printf(" %s/%s",pwd,root);
              else
                printf(" %s",root);
            hgh = ((i-base+1)*j)/bits + base;
            for (k = low; k < hgh; k++)
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
            low = hgh;
          }
      }

    //  ... and then all the initial sort & merge jobs for each block pair

    beg = ((fblock-1)/SPAN)*SPAN+1;
    end = ((lblock)/SPAN)*SPAN;

    njobs = ((end-beg)+1)*SPAN - (fblock-beg)*(fblock-beg) + (lblock-end)*(lblock-end);
    printf("# Initial sort jobs (%d)\n",njobs);

#ifdef LSF
    jobid = 1;
#endif
    for (i = beg; i <= lblock; i++)
      { int low, hgh;

        low = ((i-1)/SPAN)*SPAN + 1;
        hgh = low + (SPAN-1);
        if (hgh > lblock)
          hgh = lblock;
        if (i < fblock)
          low = fblock;
        for (j = low; j <= hgh; j++)
          {
#ifdef LSF
            printf(LSF_SORT,SPAN,SPAN,jobid++);
            printf(" \"");
#endif
            printf("LAsort");
            if (VON)
              printf(" -v");
            for (k = 0; k < NTHREADS; k++)
              if (useblock)
                { printf(" %s.%d.%s.%d.C%d",root,i,root,j,k);
                  printf(" %s.%d.%s.%d.N%d",root,i,root,j,k);
                }
              else
                { printf(" %s.%s.C%d",root,root,k);
                  printf(" %s.%s.N%d",root,root,k);
                }
            printf(" && LAmerge");
            if (VON)
              printf(" -v");
            if (useblock)
              if (hgh-low <= 1)
                printf(" %s.R%d.%d",root,SPAN,j);
              else
                printf(" L1.%d.%d",i,(j-low)+1);
            else
              printf(" %s.R%d",root,SPAN);
            for (k = 0; k < NTHREADS; k++)
              if (useblock)
                { printf(" %s.%d.%s.%d.C%d.S",root,i,root,j,k);
                  printf(" %s.%d.%s.%d.N%d.S",root,i,root,j,k);
                }
              else
                { printf(" %s.%s.C%d.S",root,root,k);
                  printf(" %s.%s.N%d.S",root,root,k);
                }
#ifdef LSF
            printf("\"");
#endif
            printf("\n");
          }
      }

    //  Check .las files (optional)

    printf("# Check all level 1 .las files (optional but recommended)\n");

    for (i = beg; i <= lblock; i++)
      { int low, hgh;

        low = ((i-1)/SPAN)*SPAN + 1;
        hgh = low + (SPAN-1);
        if (hgh > lblock)
          hgh = lblock;
        if (i < fblock)
          low = fblock;
        for (j = low; j <= hgh; j++)
          { printf("LAcheck -vS");
            if (usepath)
              printf(" %s/%s",pwd,root);
            else
              printf(" %s",root);
            if (useblock)
              if (hgh-low <= 1)
                printf(" %s.R%d.%d",root,SPAN,j);
              else
                printf(" L1.%d.%d",i,(j-low)+1);
            else
              printf(" %s.R%d",root,SPAN);
            printf("\n");
          }
      }

    //  Clean up (optional)

    printf("# Remove initial .las files (optional)\n");

    for (i = beg; i <= lblock; i++)
      { int low, hgh;

        low = ((i-1)/SPAN)*SPAN + 1;
        hgh = low + (SPAN-1);
        if (hgh > lblock)
          hgh = lblock;
        if (i < fblock)
          low = fblock;
        for (j = low; j <= hgh; j++)
          { printf("rm");
            for (k = 0; k < NTHREADS; k++)
              if (useblock)
                { printf(" %s.%d.%s.%d.C%d.las",root,i,root,j,k);
                  printf(" %s.%d.%s.%d.N%d.las",root,i,root,j,k);
                }
              else
                { printf(" %s.%s.C%d.las",root,root,k);
                  printf(" %s.%s.N%d.las",root,root,k);
                }
            printf("\n");
            printf("rm");
            for (k = 0; k < NTHREADS; k++)
              if (useblock)
                { printf(" %s.%d.%s.%d.C%d.S.las",root,i,root,j,k);
                  printf(" %s.%d.%s.%d.N%d.S.las",root,i,root,j,k);
                }
              else
                { printf(" %s.%s.C%d.S.las",root,root,k);
                  printf(" %s.%s.N%d.S.las",root,root,k);
                }
            printf("\n");
          }
      }

    //  Higher level merges (if lblock > 1)

    if (lblock > 1)
      { int pow, mway;

        //  Determine most balanced mway for merging in ceil(log_mrg lblock) levels

        pow = 1;
        for (level = 0; pow < SPAN; level++)
          pow *= DUNIT;

        for (mway = DUNIT; mway >= 3; mway--)
          if (power(mway,level) < SPAN)
            break;
        mway += 1;

        //  Issue the commands for each merge level

        { int  p, cnt, dnt, ent;

          cnt = SPAN;
          dnt = SPAN - (fblock-beg);
          ent = lblock - end;
          for (i = 1; i <= level; i++)
            { int cits, dits, eits;
              int low, hgh;

              cits = (cnt-1)/mway+1;
              dits = (dnt-1)/mway+1;
              eits = (ent-1)/mway+1;

              //  Incremental update merges

              njobs = ((end-fblock)+1)*cits;
              if (dnt > 1)
                njobs += (fblock-beg)*dits;
              if (ent > 1)
                njobs += (lblock-end)*eits;
              printf("# Level %d jobs (%d)\n",i,njobs);

#ifdef LSF
              jobid = 1;
#endif
              if (dnt > 1)
                for (j = beg; j < fblock; j++)
                  {
#ifdef LSF
                    printf(LSF_MERGE,i,SPAN,i,SPAN,jobid++);
                    printf(" \"");
#endif
                    if (dits == 1)
                      printf("mv %s.R%d.%d.las L%d.%d.0.las && ",root,SPAN,j,i,j);
                    low = 1;
                    for (p = 1; p <= dits; p++)
                      { hgh = (dnt*p)/dits;
#ifdef LSF
                        if (p > 1)
                          { printf(LSF_MERGE,i,SPAN,i,SPAN,jobid++);
                            printf(" \"");
                          }
#endif
                        printf("LAmerge");
                        if (VON)
                          printf(" -v");
                        if (dits == 1)
                          printf(" %s.R%d.%d L%d.%d.0",root,SPAN,j,i,j);
                        else
                          printf(" L%d.%d.%d",i+1,j,p);
                        for (k = low; k <= hgh; k++)
                          printf(" L%d.%d.%d",i,j,k);
#ifdef LSF
                        printf("\"");
#endif
                        printf("\n");
                        low = hgh+1;
                      }
                  }

              //  New block merges

              for (j = fblock; j <= end; j++) 
                { low = 1;
                  for (p = 1; p <= cits; p++)
                    { hgh = (cnt*p)/cits;
#ifdef LSF
                      printf(LSF_MERGE,i,SPAN,i,SPAN,jobid++);
                      printf(" \"");
#endif
                      printf("LAmerge");
                      if (VON)
                        printf(" -v");
                      if (cits == 1)
                        printf(" %s.R%d.%d",root,SPAN,j);
                      else
                        printf(" L%d.%d.%d",i+1,j,p);
                      for (k = low; k <= hgh; k++)
                        printf(" L%d.%d.%d",i,j,k);
#ifdef LSF
                      printf("\"");
#endif
                      printf("\n");
                      low = hgh+1;
                    }
                }

              if (ent > 1)
                for (j = end+1; j <= lblock; j++)
                  { low = 1;
                    for (p = 1; p <= eits; p++)
                      { hgh = (ent*p)/eits;
#ifdef LSF
                        printf(LSF_MERGE,i,SPAN,i,SPAN,jobid++);
                        printf(" \"");
#endif
                        printf("LAmerge");
                        if (VON)
                          printf(" -v");
                        if (eits == 1)
                          printf(" %s.R%d.%d",root,SPAN,j);
                        else
                          printf(" L%d.%d.%d",i+1,j,p);
                        for (k = low; k <= hgh; k++)
                          printf(" L%d.%d.%d",i,j,k);
#ifdef LSF
                        printf("\"");
#endif
                        printf("\n");
                        low = hgh+1;
                      }
                  }

              //  Check new .las (optional)

              printf("# Check all level %d .las files (optional but recommended)\n",i+1);

              if (dnt > 1)
                for (j = beg; j < fblock; j++)
                  { low = 1;
                    for (p = 1; p <= dits; p++)
                      { hgh = (dnt*p)/dits;
                        printf("LAcheck -vS");
                        if (usepath)
                          printf(" %s/%s",pwd,root);
                        else
                          printf(" %s",root);
                        if (dits == 1)
                          printf(" %s.R%d.%d",root,SPAN,j);
                        else
                          printf(" L%d.%d.%d",i+1,j,p);
                        printf("\n");
                        low = hgh+1;
                      }
                  }

              for (j = fblock; j <= end; j++) 
                { low = 1;
                  for (p = 1; p <= cits; p++)
                    { hgh = (cnt*p)/cits;
                      printf("LAcheck -vS");
                      if (usepath)
                        printf(" %s/%s",pwd,root);
                      else
                        printf(" %s",root);
                      if (cits == 1)
                        printf(" %s.R%d.%d",root,SPAN,j);
                      else
                        printf(" L%d.%d.%d",i+1,j,p);
                      printf("\n");
                      low = hgh+1;
                    }
                }

              if (ent > 1)
                for (j = end+1; j <= lblock; j++)
                  { low = 1;
                    for (p = 1; p <= eits; p++)
                      { hgh = (ent*p)/eits;
                        printf("LAcheck -vS");
                        if (usepath)
                          printf(" %s/%s",pwd,root);
                        else
                          printf(" %s",root);
                        if (eits == 1)
                          printf(" %s.R%d.%d",root,SPAN,j);
                        else
                          printf(" L%d.%d.%d",i+1,j,p);
                        printf("\n");
                        low = hgh+1;
                      }
                  }

              //  Cleanup (optional)

              printf("# Remove level %d .las files (optional)\n",i);

              if (dnt > 1)
                for (j = beg; j < fblock; j++)
                  { low = 1;
                    for (p = 1; p <= dits; p++)
                      { hgh = (dnt*p)/dits;
                        printf("rm");
                        if (dits == 1)
                          printf(" L%d.%d.0.las",i,j);
                        for (k = low; k <= hgh; k++)
                          printf(" L%d.%d.%d.las",i,j,k);
                        printf("\n");
                        low = hgh+1;
                      }
                  }

              for (j = fblock; j <= end; j++) 
                { low = 1;
                  for (p = 1; p <= cits; p++)
                    { hgh = (cnt*p)/cits;
                      printf("rm");
                      for (k = low; k <= hgh; k++)
                        printf(" L%d.%d.%d.las",i,j,k);
                      printf("\n");
                      low = hgh+1;
                    }
                }

              if (ent > 1)
                for (j = end+1; j <= lblock; j++)
                  { low = 1;
                    for (p = 1; p <= eits; p++)
                      { hgh = (ent*p)/eits;
                        printf("rm");
                        for (k = low; k <= hgh; k++)
                          printf(" L%d.%d.%d.las",i,j,k);
                        printf("\n");
                        low = hgh+1;
                      }
                  }

              dnt = dits;
              cnt = cits;
              ent = eits;
            }
        }
      }

    //  Finish with MASKrep

    printf("# REPmask jobs (%d)\n",(end-beg)/BUNIT+1);

#ifdef LSF
    jobid = 1;
#endif
    for (i = beg; i <= end; i += BUNIT)
      {
#ifdef LSF
        printf(LSF_MASK,SPAN,SPAN,jobid++);
        printf(" \"");
#endif
        printf("REPmask");
        if (VON)
          printf(" -v");
        printf(" -c%d -mrep%d",CINT,SPAN);
        if (usepath)
          printf(" %s/%s",pwd,root);
        else
          printf(" %s",root);
        j = i+BUNIT;
        if (j > end)
          j = end+1;
        for (k = i; k < j; k++)
          if (useblock)
            printf(" %s.R%d.%d",root,SPAN,k);
          else
            printf(" %s.R%d",root,SPAN);
#ifdef LSF
        printf("\"");
#endif
        printf("\n");
      }

    //  Cleanup (optinoal)

    printf("# Remove final R%d.las files (optional)\n",SPAN);

    for (i = beg; i <= end; i += BUNIT)
      { printf("rm");
        j = i+BUNIT;
        if (j > end)
          j = end+1;
        for (k = i; k < j; k++)
          if (useblock)
            printf(" %s.R%d.%d.las",root,SPAN,k);
          else
            printf(" %s.R%d.las",root,SPAN);
        printf("\n");
      }
  }

  printf("# Once all the .rep masks have been computed for every block\n");
  printf("#   you should call 'Catrack' to merge them, and then you should\n");
  printf("#   remove the block tracks\n");

  free(root);
  free(pwd);

  exit (0);
}
