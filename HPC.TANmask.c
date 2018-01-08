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

#undef  LSF  //  define if want a directly executable LSF script

static char *Usage[] =
  { "[-v] [-k<int(12)>] [-w<int(4)>] [-h<int(35)>] [-T<int(4)>] [-P<dir(/tmp)>]",
    "     [-e<double(.70)] [-l<int(500)>] [-s<int(100)] [-f<name>]",
    "     <reads:db|dam> [<first:int>[-<last:int>]"
  };

#define LSF_TAND "bsub -q medium -n 4 -o TANDEM.out -e TANDEM.err -R span[hosts=1] -J tandem#%d"
#define LSF_CHECK \
          "bsub -q short -n 12 -o CHECK%d.DAL.out -e CHECK%d.DAL.err -R span[hosts=1] -J check#%d"
#define LSF_MASK \
          "bsub -q short -n 12 -o MASK.TAN.out -e MASK.TAN.err -R span[hosts=1] -J masktan#%d"

int main(int argc, char *argv[])
{ int   nblocks;
  int   usepath;
  int   useblock;
  int   fblock, lblock;
#ifdef LSF
  int   jobid;
#endif

  FILE *out;
  char  name[100];
  char *pwd, *root;

#define BUNIT  4

  int    VON;
  int    WINT, HINT, KINT, SINT, LINT;
  int    NTHREADS;
  double EREL;
  char  *ONAME;
  char  *PDIR;

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
    ONAME = NULL;
    PDIR  = NULL;
    out   = stdout;

    NTHREADS = 4;

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
          case 'f':
            ONAME = argv[i]+2;
            break;
          case 'h':
            ARG_POSITIVE(HINT,"Hit threshold (in bp.s)")
            break;
          case 'k':
            ARG_POSITIVE(KINT,"K-mer length")
            if (KINT > 32)
              { fprintf(stderr,"%s: K-mer length must be 32 or less\n",Prog_Name);
                exit (1);
              }
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
          case 'P':
            PDIR = argv[i]+2;
            break;
          case 'T':
            ARG_POSITIVE(NTHREADS,"Number of threads")
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

    for (j = 1; 2*j <= NTHREADS; j *= 2)
      ;
    NTHREADS = j;
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
      SYSTEM_READ_ERROR
    for (i = 0; i < nfiles; i++)
      { char buffer[30001];

        if (fgets(buffer,30000,dbvis) == NULL)
          SYSTEM_READ_ERROR
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
      { file = fopen(Catenate(pwd,"/.",root,Numbered_Suffix(".",fblock-1,".tan.anno")),"r");
        if (file == NULL)
          { if (usepath)
              fprintf(stderr,"%s: File %s/.%s.%d.tan.anno should already be present!\n",
                             Prog_Name,pwd,root,fblock-1);
            else
              fprintf(stderr,"%s: File .%s.%d.tan.anno should already be present!\n",
                             Prog_Name,root,fblock-1);
            exit (1);
          }
        else
          fclose(file);
      }
    if (useblock)
      file = fopen(Catenate(pwd,"/.",root,Numbered_Suffix(".",fblock,".tan.anno")),"r");
    else
      file = fopen(Catenate(pwd,"/.",root,".tan.anno"),"r");
    if (file != NULL)
      { if (useblock)
          { if (usepath)
              fprintf(stderr,"%s: File %s/.%s.%d.tan.anno should not yet exist!\n",
                             Prog_Name,pwd,root,fblock);
            else
              fprintf(stderr,"%s: File .%s.%d.tan.anno should not yet exist!\n",
                             Prog_Name,root,fblock);
            exit (1);
          }
        else
          { if (usepath)
              fprintf(stderr,"%s: File %s/.%s.tan.anno should not yet exist!\n",Prog_Name,pwd,root);
            else
              fprintf(stderr,"%s: File .%s.tan.anno should not yet exist!\n",Prog_Name,root);
            exit (1);
          }
      }
  }

  { int njobs;
    int i, j, k;

    //  Produce all necessary datandem jobs ...

    if (ONAME != NULL)
      { sprintf(name,"%s.01.OVL",ONAME);
        out = fopen(name,"w");
      }

    njobs = (lblock - fblock)/BUNIT + 1;

    fprintf(out,"# Datander jobs (%d)\n",njobs);

#ifdef LSF
    jobid = 1;
#endif
    for (i = fblock; i <= lblock; i += BUNIT)
      {
#ifdef LSF
        fprintf(out,LSF_TAND,jobid++);
        fprintf(out," \"");
#endif
        fprintf(out,"datander");
        if (VON)
          fprintf(out," -v");
        if (KINT != 12)
          fprintf(out," -k%d",KINT);
        if (WINT != 4)
          fprintf(out," -w%d",WINT);
        if (HINT != 35)
          fprintf(out," -h%d",HINT);
        if (EREL > .1)
          fprintf(out," -e%g",EREL);
        if (LINT != 500)
          fprintf(out," -l%d",LINT);
        if (SINT != 100)
          fprintf(out," -s%d",SINT);
        if (PDIR != NULL)
          fprintf(out," -P%s",PDIR);
        if (NTHREADS != 4)
          fprintf(out," -T%d",NTHREADS);
        j = i+BUNIT;
        if (j > lblock+1)
          j = lblock+1;
        for (k = i; k < j; k++)
          if (useblock)
            if (usepath)
              fprintf(out," %s/%s.%d",pwd,root,k);
            else
              fprintf(out," %s.%d",root,k);
          else
            if (usepath)
              fprintf(out," %s/%s",pwd,root);
            else
              fprintf(out," %s",root);

#ifdef LSF
        fprintf(out,"\"");
#endif
        fprintf(out,"\n");
      }

    //  Check .las (option)

    if (ONAME != NULL)
      { fclose(out);
        sprintf(name,"%s.02.CHECK.OPT",ONAME);
        out = fopen(name,"w");
      }

    fprintf(out,"# Check all .las files jobs (%d) (optional but recommended)\n",
                (lblock-fblock)/(BUNIT+1)+1);

#ifdef LSF
    jobid = 1;
#endif
    for (i = fblock; i <= lblock; )
      { k = i+BUNIT;
        if (k > lblock)
          k = lblock;
#ifdef LSF
        fprintf(out,LSF_CHECK,0,0,jobid++);
        fprintf(out," \"");
#endif
        fprintf(out,"LAcheck -vS");
        if (usepath)
          fprintf(out," %s/%s",pwd,root);
        else
          fprintf(out," %s",root);
        while (i <= k)
          { if (useblock)
              fprintf(out," TAN.%s.%d",root,i);
            else
              fprintf(out," TAN.%s",root);
            i += 1;
          }
#ifdef LSF
        fprintf(out,"\"");
#endif
        fprintf(out,"\n");
      }

    //  Finish with MASKtan

    if (ONAME != NULL)
      { fclose(out);
        sprintf(name,"%s.03.MASK",ONAME);
        out = fopen(name,"w");
      }

    fprintf(out,"# TANmask jobs (%d)\n",njobs);

#ifdef LSF
    jobid = 1;
#endif
    for (i = fblock; i <= lblock; i += BUNIT)
      {
#ifdef LSF
        fprintf(out,LSF_MASK,jobid++);
        fprintf(out," \"");
#endif
        fprintf(out,"TANmask");
        if (VON)
          fprintf(out," -v");
        if (LINT != 500)
          fprintf(out," -l%d",LINT);
        if (usepath)
          fprintf(out," %s/%s",pwd,root);
        else
          fprintf(out," %s",root);
        j = i+BUNIT;
        if (j > lblock+1)
          j = lblock+1;
        for (k = i; k < j; k++)
          if (useblock)
            fprintf(out," TAN.%s.%d",root,k);
          else
            fprintf(out," TAN.%s",root);
#ifdef LSF
        fprintf(out,"\"");
#endif
        fprintf(out,"\n");
      }

    //  Clean up

    if (ONAME != NULL)
      { fclose(out);
        sprintf(name,"%s.04.RM",ONAME);
        out = fopen(name,"w");
      }

    fprintf(out,"# Cleanup all T.las files\n");

    for (i = fblock; i <= lblock; i += BUNIT)
      { fprintf(out,"rm");
        j = i+BUNIT;
        if (j > lblock+1)
          j = lblock+1;
        for (k = i; k < j; k++)
          if (useblock)
            fprintf(out," TAN.%s.%d.las",root,k);
          else
            fprintf(out," TAN.%s.las",root);
        fprintf(out,"\n");
      }

    if (ONAME != NULL)
      fclose(out);
  }

  printf("# Once all the .tan masks have been computed for every block\n");
  printf("#   you should call 'Catrack' to merge them, and then you should\n");
  printf("#   remove the individual block tracks, e.g.:\n");
  if (usepath)
    { printf("#      Catrack -v %s/%s tan\n",pwd,root);
      printf("#      rm %s/.%s.*.tan.*\n",pwd,root);
    }
  else
    { printf("#      Catrack -v %s tan\n",root);
      printf("#      rm .%s.*.tan.*\n",root);
    }

  free(root);
  free(pwd);

  exit (0);
}
