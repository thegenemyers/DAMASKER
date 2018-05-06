/*********************************************************************************************\
 *
 *  Produce a script to compute overlaps for all block pairs of a DB, and then sort and merge
 *    them into as many .las files as their are blocks.
 *
 *  Author:  Gene Myers
 *  Date  :  June 1, 2014
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

#undef  LSF    //  define if want a directly executable LSF script
#undef  SLURM  //  define if want a directly executable SLURM script

static char *Usage[] =
  { "[-vbd] [-t<int>] [-w<int(6)>] [-l<int(1000)>] [-s<int(100)>] [-M<int>]",
    "       [-n<name(rep-g)>] [-P<dir(/tmp)>] [-B<int(4)>] [T<int(4)>] [-f<name>]",
    "       [-k<int(14)>] [-h<int(35)>] [-e<double(.70)>] [-m<track>]+",
    "       -g<int> -c<int> <reads:db|dam> [<block:int>[-<range:int>]"
  };

#ifdef LSF

#define HPC

#define HPC_ALIGN \
    "bsub -q medium -n %d -o DAL.REP%d.out -e DAL.REP%d.err -R span[hosts=1] -J ralign#%d"
#define HPC_MERGE \
    "bsub -q short -n 12 -o MERGE.REP%d.out -e MERGE.REP%d.err -R span[hosts=1] -J rmerge#%d"
#define HPC_CHECK \
    "bsub -q short -n 12 -o CHECK.REP%d.out -e CHECK.REP%d.err -R span[hosts=1] -J rcheck#%d"
#define HPC_MASK \
    "bsub -q short -n 12 -o MASK.REP%d.out -e MASK.REP%d.err -R span[hosts=1] -J maskrep#%d"

#endif

#ifdef SLURM

#define HPC

#define HPC_ALIGN \
    "srun -p batch -n 1 -c %d --mem_per_cpu=%d -o DAL.REP%d.out -e DAL.REP%d.err -J ralign#%d"
#define HPC_MERGE \
    "srun -p batch -n 1 -c 12 -t 00:05:00 -o MERGE.REP%d.DAL.out -e MERGE.REP%d.err -J rmerge#%d"
#define HPC_CHECK \
    "srun -p batch -n 1 -c 12 -t 00:05:00 -o CHECK.REP%d.out -e CHECK.REP%d.err -J rcheck#%d"
#define HPC_MASK \
    "srun -p batch -n 1 -c 12 -t 00:05:00 -o MASK.REP%d.out -e MASK.REP%d.err -J maskrep#%d"

#endif

int main(int argc, char *argv[])
{ int   nblocks;
  int   usepath;
  int   useblock;
  int   fblock, lblock;
#ifdef HPC
  int   jobid;
#endif

  FILE *out;
  char  name[100];
  char *pwd, *root;

  int    CINT, SPAN;
  int    BUNIT;
  int    VON, BON, DON;
  int    WINT, TINT, HINT, KINT, SINT, LINT, MINT;
  int    NTHREADS;
  char  *MASK_NAME, defname[25];
  double EREL;
  int    MMAX, MTOP;
  char **MASK;
  char  *ONAME;
  char  *PDIR;

  { int    i, j, k;         //  Process options
    int    flags[128];
    char  *eptr;

    ARG_INIT("HPC.REPmask")

    BUNIT = 4;
    KINT  = 14;
    WINT  = 6;
    HINT  = 35;
    TINT  = 0;
    EREL  = 0.;
    LINT  = 1000;
    SINT  = 100;
    MINT  = -1;
    PDIR  = NULL;

    MTOP = 0;
    MMAX = 10;
    MASK = (char **) Malloc(MMAX*sizeof(char *),"Allocating mask track array");
    if (MASK == NULL)
      exit (1);
    ONAME = NULL;
    out   = stdout;

    NTHREADS  = 4;

    MASK_NAME = NULL;
    SPAN = -1;
    CINT = -1;

    j = 1;
    for (i = 1; i < argc; i++)
      if (argv[i][0] == '-')
        switch (argv[i][1])
        { default:
            ARG_FLAGS("vbd");
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
          case 'f':
            ONAME = argv[i]+2;
            break;
          case 'g':
            ARG_POSITIVE(SPAN,"Block span")
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
          case 'n':
            MASK_NAME = argv[i]+2;
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
          case 'M':
            ARG_NON_NEGATIVE(MINT,"Memory allocation (in Gb)")
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
    BON = flags['b'];
    DON = flags['d'];

    if (argc < 2 || argc > 3)
      { fprintf(stderr,"Usage: %s %s\n",Prog_Name,Usage[0]);
        fprintf(stderr,"       %*s %s\n",(int) strlen(Prog_Name),"",Usage[1]);
        fprintf(stderr,"       %*s %s\n",(int) strlen(Prog_Name),"",Usage[2]);
        fprintf(stderr,"       %*s %s\n",(int) strlen(Prog_Name),"",Usage[3]);
        fprintf(stderr,"\n");
        fprintf(stderr,"     Passed through to daligner.\n");
        fprintf(stderr,"      -k: k-mer size (must be <= 32).\n");
        fprintf(stderr,"      -w: Look for k-mers in averlapping bands of size 2^-w.\n");
        fprintf(stderr,"      -h: A seed hit if the k-mers in band cover >= -h bps in the");
        fprintf(stderr," targest read.\n");
        fprintf(stderr,"      -t: Ignore k-mers that occur >= -t times in a block.\n");
        fprintf(stderr,"      -M: Use only -M GB of memory by ignoring most frequent k-mers.\n");
        fprintf(stderr,"\n");
        fprintf(stderr,"      -e: Look for alignments with -e percent similarity.\n");
        fprintf(stderr,"      -l: Look for alignments of length >= -l.\n");
        fprintf(stderr,"      -s: Use -s as the trace point spacing for encoding alignments.\n");
        fprintf(stderr,"\n");
        fprintf(stderr,"      -T: Use -T threads.\n");
        fprintf(stderr,"      -P: Do first level sort and merge in directory -P.\n");
        fprintf(stderr,"      -m: Soft mask the blocks with the specified mask.\n");
        fprintf(stderr,"      -b: For AT/GC biased data, compensate k-mer counts (deprecated).\n");
        fprintf(stderr,"\n");
        fprintf(stderr,"     Passed through to REPmask.\n");
        fprintf(stderr,"      -c: coverage threshold for repeat intervals.\n");
        fprintf(stderr,"      -n: use this name for the repeat mask track.\n");
        fprintf(stderr,"\n");
        fprintf(stderr,"     Script control.\n");
        fprintf(stderr,"      -v: Run all commands in script in verbose mode.\n");
        fprintf(stderr,"      -g: # of blocks per comparison group.\n");
        fprintf(stderr,"      -d: Put .las files for each target block in a sub-directory\n");
        fprintf(stderr,"      -B: # of block compares per daligner job\n");
        fprintf(stderr,"      -f: Place script bundles in separate files with prefix <name>\n");
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

    if (MASK_NAME == NULL)
      { sprintf(defname,"rep%d",SPAN);
        MASK_NAME = defname;
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

    if (nblocks < SPAN)
      { fprintf(stderr,"%s: There a fewer than -g = %d blocks in the DB!\n",Prog_Name,SPAN);
        exit (1);
      }

    usepath = (strcmp(pwd,".") != 0);

    fclose(dbvis);
  }

  //  Set range fblock-lblock checking that DB.<fblock-1>.las exists & DB.<fblock>.las does not

  { char *eptr, *fptr, *sfx;
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
            if (lblock < SPAN)
              { fprintf(stderr,"%s: End of range %d must be >= group span -g = %d\n",
                               Prog_Name,lblock,SPAN);
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

    sfx = Strdup(Numbered_Suffix(".rep",SPAN,".anno"),"Allocating small string!");
    if (sfx == NULL)
      exit (1);

    if (fblock > 1)
      { file = fopen(Catenate(pwd,"/.",root,Numbered_Suffix(".",fblock-1,sfx)),"r");
        if (file == NULL)
          { if (usepath)
              fprintf(stderr,"%s: File %s/%s.%d%s should already be present!\n",
                             Prog_Name,pwd,root,fblock-1,sfx);
            else
              fprintf(stderr,"%s: File %s.%d%s should already be present!\n",
                             Prog_Name,root,fblock-1,sfx);
            exit (1);
          }
        else
          fclose(file);
      }
    if (useblock)
      file = fopen(Catenate(pwd,"/.",root,Numbered_Suffix(".",fblock,sfx)),"r");
    else
      file = fopen(Catenate(pwd,"/.",root,sfx),"r");
    if (file != NULL)
      { if (usepath)
          if (useblock)
            fprintf(stderr,"%s: File %s/.%s.%d%s should not yet exist!\n",
                           Prog_Name,pwd,root,fblock,sfx);
          else
            fprintf(stderr,"%s: File %s/.%s%s as should not yet exist!\n",Prog_Name,pwd,root,sfx);
        else
          if (useblock)
            fprintf(stderr,"%s: File .%s.%d%s should not yet exist!\n",Prog_Name,root,fblock,sfx);
          else
            fprintf(stderr,"%s: File .%s%s should not yet exist!\n",Prog_Name,root,sfx);
        exit (1);
      }

    DON = (DON && (SPAN > 1));
  }

  { int njobs;
    int i, j, k;

    if (DON)
      { if (ONAME != NULL)
          { sprintf(name,"%s.00.MKDIR",ONAME);
            out = fopen(name,"w");
          }

        fprintf(out,"# Create work subdirectories\n");
        for (i = fblock; i <= lblock; i++)
          fprintf(out,"mkdir temp%d\n",i);

        if (ONAME != NULL)
          fclose(out);
      }

    //  Produce all necessary daligner jobs ...

    if (ONAME != NULL)
      { sprintf(name,"%s.01.OVL",ONAME);
        out = fopen(name,"w");
      }

    njobs = 0;
    for (i = fblock; i <= lblock; i++)
      { int base;

        base = fblock + ((i-fblock)/SPAN)*SPAN;
        if (base + SPAN > lblock+1)
          base = (lblock+1) - SPAN;
        njobs += (i-base)/BUNIT+1;
      }

    fprintf(out,"# Daligner jobs (%d)\n",njobs);

#ifdef HPC
    jobid = 1;
#endif
    for (i = fblock; i <= lblock; i++)
      { int bits, base;
        int low, hgh;

        base = fblock + ((i-fblock)/SPAN)*SPAN;
        if (base + SPAN > lblock+1)
          base = (lblock+1) - SPAN;
        bits = (i-base)/BUNIT+1;
        low  = base;
        for (j = 1; j <= bits; j++)
          {
#ifdef LSF
            fprintf(out,HPC_ALIGN,NTHREADS,SPAN,SPAN,jobid++);
            fprintf(out," \"");
#endif
#ifdef SLURM
            if (MINT >= 0)
              fprintf(out,HPC_ALIGN,NTHREADS,(MINT*1024)/NTHREADS,SPAN,SPAN,jobid++);
            else
              fprintf(out,HPC_ALIGN,NTHREADS,(16*1024)/NTHREADS,SPAN,SPAN,jobid++);
            fprintf(out," \"");
#endif
            fprintf(out,"daligner");
            if (VON)
              fprintf(out," -v");
            if (BON)
              fprintf(out," -b");
            if (KINT != 14)
              fprintf(out," -k%d",KINT);
            if (WINT != 6)
              fprintf(out," -w%d",WINT);
            if (HINT != 35)
              fprintf(out," -h%d",HINT);
            if (TINT > 0)
              fprintf(out," -t%d",TINT);
            if (EREL > .1)
              fprintf(out," -e%g",EREL);
            if (LINT != 1000)
              fprintf(out," -l%d",LINT);
            if (SINT != 100)
              fprintf(out," -s%d",SINT);
            if (MINT >= 0)
              fprintf(out," -M%d",MINT);
            if (PDIR != NULL)
              fprintf(out," -P%s",PDIR);
            if (NTHREADS != 4)
              fprintf(out," -T%d",NTHREADS);
            for (k = 0; k < MTOP; k++)
              fprintf(out," -m%s",MASK[k]);
            if (useblock)
              if (usepath)
                fprintf(out," %s/%s.%d",pwd,root,i);
              else
                fprintf(out," %s.%d",root,i);
            else
              if (usepath)
                fprintf(out," %s/%s",pwd,root);
              else
                fprintf(out," %s",root);
            hgh = ((i-base+1)*j)/bits + base;
            for (k = low; k < hgh; k++)
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

            if (SPAN == 1)   // ==> [low,hgh) = [i,i+1)
              if (useblock)
                fprintf(out," && mv %s.%d.%s.%d.las %s.R1.%d.las",root,i,root,i,root,i);
              else
                fprintf(out," && mv %s.%s.las %s.R1.%d.las",root,root,root,i);
            else if (DON)
              { fprintf(out," && mv");
                for (k = low; k < hgh; k++)
                  fprintf(out," %s.%d.%s.%d.las",root,i,root,k);
                fprintf(out," temp%d",i);
                for (k = low; k < hgh; k++)
                  if (k != i)
                    fprintf(out," && mv %s.%d.%s.%d.las temp%d",root,k,root,i,k);
              }

#ifdef HPC
            fprintf(out,"\"");
#endif
            fprintf(out,"\n");
            low = hgh;
          }
      }

    //  Check .las files (optional)

    if (ONAME != NULL)
      { fclose(out);
        sprintf(name,"%s.02.CHECK.OPT",ONAME);
        out = fopen(name,"w");
      }

    njobs = (lblock-fblock)+1;

    fprintf(out,"# Check inital .las jobs (%d) (optional but recommended)\n",njobs);

#ifdef HPC
    jobid = 1;
#endif
    for (i = fblock; i <= lblock; i++)
      { int base;

        base = fblock + ((i-fblock)/SPAN)*SPAN;
        if (base + SPAN > lblock+1)
          base = (lblock+1)-SPAN;
#ifdef HPC
        fprintf(out,HPC_CHECK,SPAN,SPAN,jobid++);
        fprintf(out," \"");
#endif
        fprintf(out,"LAcheck -vS");
        if (usepath)
          fprintf(out," %s/%s",pwd,root);
        else
          fprintf(out," %s",root);
        if (SPAN == 1)
          fprintf(out," %s.R1.%d",root,i);
        else if (DON)
          fprintf(out," temp%d/%s.%d.%s.%c%d-%d",i,root,i,root,BLOCK_SYMBOL,base,base+(SPAN-1));
        else
          fprintf(out," %s.%d.%s.%c%d-%d",root,i,root,BLOCK_SYMBOL,base,base+(SPAN-1));
#ifdef HPC
        fprintf(out,"\"");
#endif
        fprintf(out,"\n");
      }

    if (ONAME != NULL)
      fclose(out);

    //  Merges required if SPAN > 1

    if (SPAN > 1)
      { if (ONAME != NULL)
          { sprintf(name,"%s.03.MERGE",ONAME);
            out = fopen(name,"w");
          }

        fprintf(out,"# Merge jobs (%d)\n",njobs);

#ifdef HPC
        jobid = 1;
#endif
        for (i = fblock; i <= lblock; i++)
          { int base;

            base = fblock + ((i-fblock)/SPAN)*SPAN;
            if (base + SPAN > lblock+1)
              base = (lblock+1)-SPAN;
#ifdef HPC
            fprintf(out,HPC_MERGE,SPAN,SPAN,jobid++);
            fprintf(out," \"");
#endif
            fprintf(out,"LAmerge ");
            if (VON)
              fprintf(out," -v");
            if (DON)
              fprintf(out," temp%d/%s.R%d.%d",i,root,SPAN,i);
            else
              fprintf(out," %s.R%d.%d",root,SPAN,i);
            if (usepath)
              fprintf(out," %s/%s",pwd,root);
            else
              fprintf(out," %s",root);
            if (DON)
              fprintf(out," temp%d/%s.%d.%s.%c%d-%d",i,root,i,root,BLOCK_SYMBOL,base,base+(SPAN-1));
            else
              fprintf(out," %s.%d.%s.%c%d-%d",root,i,root,BLOCK_SYMBOL,base,base+(SPAN-1));

            fprintf(out," && LAcheck -vS");
            if (DON)
              fprintf(out," temp%d/%s.R%d.%d",i,root,SPAN,i);
            else
              fprintf(out," %s.R%d.%d",root,SPAN,j);
#ifdef HPC
            fprintf(out,"\"");
#endif
            fprintf(out,"\n");
          }

        if (ONAME != NULL)
          { fclose(out);
            sprintf(name,"%s.04.RM",ONAME);
            out = fopen(name,"w");
          }

        fprintf(out,"# Remove block-pair .las files\n");

        for (j = fblock; j <= lblock; j++) 
          if (DON)
            fprintf(out,"rm -r temp%d\n",j);
          else
            fprintf(out,"rm %s.%d.%s.*.las\n",root,j,root);

        if (ONAME != NULL)
          fclose(out);
      }

    //  Finish with MASKrep

    if (ONAME != NULL)
      { sprintf(name,"%s.05.MASK",ONAME);
        out = fopen(name,"w");
      }
 
#ifdef HPC
    jobid = 1;
#endif
    njobs  = (lblock-fblock) / BUNIT + 1;

    fprintf(out,"# REPmask jobs (%d)\n",njobs);

    { int low, hgh;

      low = fblock;
      for (j = 1; j <= njobs; j++)
        { hgh = (fblock-1) + (((lblock-fblock)+1)*j)/njobs;
          
#ifdef HPC
          fprintf(out,HPC_MASK,SPAN,SPAN,jobid++);
          fprintf(out," \"");
#endif
          fprintf(out,"REPmask");
          if (VON)
            fprintf(out," -v");
          fprintf(out," -c%d -n%s",CINT,MASK_NAME);
          if (usepath)
            fprintf(out," %s/%s",pwd,root);
          else
            fprintf(out," %s",root);
          if (DON)
            fprintf(out," temp%d/%s.R%d.%c%d-%d",k,root,SPAN,BLOCK_SYMBOL,low,hgh);
          else
            fprintf(out," %s.R%d.%c%d-%d",root,SPAN,BLOCK_SYMBOL,low,hgh);
          low = hgh+1;
#ifdef HPC
          fprintf(out,"\"");
#endif
          fprintf(out,"\n");
        }
    }

    if (ONAME != NULL)
      { sprintf(name,"%s.06.RM",ONAME);
        out = fopen(name,"w");
      }

    if (DON)
      fprintf(out,"# Cleanup all temporary directories\n");
    else
      fprintf(out,"# Cleanup all R%d.las files\n",SPAN);

    if (DON)
      fprintf(out,"rm -r temp*\n");
    else
      fprintf(out,"rm %s.R%d.*.las\n",root,SPAN);
  
    if (ONAME != NULL)
      fclose(out);
  }

  printf("# Once all the .rep masks have been computed for every block\n");
  printf("#   you should call 'Catrack' to merge them, and then you should\n");
  printf("#   remove the individual block tracks, e.g.:\n");
  if (usepath)
    { printf("#      Catrack -v %s/%s %s\n",pwd,root,MASK_NAME);
      printf("#      rm %s/.%s.*.%s.*\n",pwd,root,MASK_NAME);
    }
  else
    { printf("#      Catrack -v %s %s\n",root,MASK_NAME);
      printf("#      rm .%s.*.%s.*\n",root,MASK_NAME);
    }


  free(root);
  free(pwd);

  exit (0);
}
