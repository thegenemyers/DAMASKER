/*********************************************************************************************\
 *
 *  HPC script to update coverage statistics upon an update to the database
 *
 *  Produce a script to
 *  1. create low complexity mask for the given DB (DBdust)
 *  2. create tandem repeat mask for the given DB (datander + TANmask)
 *  3. create overlap based repeat mask for the given DB,
 *       based on diagonal blocks N vs N (daligner + REPmask)
 *  4. compute overlaps for block 1 of pairs of a DB, and then sort and merge them
 *      into as many .las files as their are blocks. Run DAScover on the resulting .las files
 *
 *
 *  Author:  Martin Pippel, Gene Myers
 *  Date  :  March 19, 2018
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

static char *Usage[] =
  { "[-vlF] [-U(w<int(64)> |t<double(2.)> |m<int(10)> |b) ]",
    "   [-S(k<int(12)> |w<int(4)> |h<int(35)> |e<double(.7)> |l<int(500)>) ]",
    "   [-L(k<int(14)> |w<int(6)> |h<int(35)> |e<double(.7)> |l<int(1000)> |t<int>) ]",
    "   [-c<int(10)>] [-s<int(100)] [-M<int>] [-P<dir(/tmp)>] [-T<int(4)>]",
    "   [-B<int(4)>] [-f<name>] <reads:db|dam> [<target:int(1)>]"
  };

int main(int argc, char *argv[])
{ int   nblocks;
  int   usepath;
  int   fblock, lblock, tblock;
  char *pwd, *root;

  int   VERBOSE;     //  Option parameters
  int   ALLBLKS;
  int   FORCE;

  int   BLOCKID;
  char *PROJECT;
  int   BUNIT;

  int   COVERAGE;
  int   TRACE;
  char *TMPDIR;
  int   NTHREADS;
  int   MEMORY;
 
  int   D_BIAS;
  int   D_WINDOW;
  float D_THRESH;
  int   D_MINBP;

  int   S_KMER;
  int   S_BAND;
  int   S_HITS;
  float S_ERATE;
  int   S_OLEN;
    
  int   D_KMER;
  int   D_BAND;
  int   D_HITS;
  float D_ERATE;
  int   D_OLEN;
  int   D_FREQ;

  { int i, j, k;         //  Process options
    int flags[128];
    char *eptr;

    ARG_INIT("HPC.DAScover")

    BLOCKID = 1;
    PROJECT = NULL;
    BUNIT   = 4;

    COVERAGE = 10;
    TRACE    = 100;
    TMPDIR   = "/tmp";
    NTHREADS = 4;
    MEMORY   = 16;

    D_BIAS   = 0;
    D_WINDOW = 64;
    D_THRESH = -1.;
    D_MINBP  = 10;

    S_KMER  = 12;
    S_BAND  = 4;
    S_HITS  = 35;
    S_ERATE = 0.;
    S_OLEN  = 500;

    D_FREQ  = -1;
    D_KMER  = 14;
    D_BAND  = 6;
    D_HITS  = 35;
    D_ERATE = 0.;
    D_OLEN  = 1000;

    j = 1;
    for (i = 1; i < argc; i++)
      if (argv[i][0] == '-')
        switch (argv[i][1])
        { default:
            ARG_FLAGS("vlF")
            break;
          case 'c':
            ARG_POSITIVE(COVERAGE, "Repeat coverage threshold")
            break;
          case 'f':
            PROJECT = argv[i]+2;
            break;
          case 's':
            ARG_POSITIVE(TRACE, "Trace spacing")
            break;
          case 'B':
            ARG_POSITIVE(BUNIT, "Blocks per command")
            break;
          case 'M':
            ARG_POSITIVE(MEMORY, "Daligner memory usage")
            break;
          case 'P':
            TMPDIR = argv[i]+2;
            break;
          case 'T':
            ARG_POSITIVE(NTHREADS, "Number of threads (daligner+datander)")
            break;
          case 'U':
            argv[i] += 1;
            switch (argv[i][1])
            { case 'b':
                D_BIAS = 1;
                break;
              case 'm':
                ARG_POSITIVE(D_MINBP, "DBdust minimum bases")
                break;
              case 'w':
                ARG_POSITIVE(D_WINDOW, "DBdust window size")
                break;
              case 't':
                ARG_REAL(D_THRESH)
                break;
              default:
                fprintf(stderr, "%s: option '-%s' not supported for DBdust\n",
                               Prog_Name, argv[i]);
                exit(1);
            }
            break;
          case 'S':
            argv[i] += 1;
            switch (argv[i][1])
            { case 'e':
                ARG_REAL(S_ERATE)
                if (S_ERATE < .7 || S_ERATE >= 1.)
                  { fprintf(stderr, "%s: Average correlation must be in [.7,1.) (%g)\n",
                                   Prog_Name, S_ERATE);
                    exit(1);
                  }
                break;
              case 'h':
                ARG_POSITIVE(S_HITS, "Datander hit threshold")
                break;
              case 'k':
                ARG_POSITIVE(S_KMER, "Datander kmer")
                break;
              case 'l':
                ARG_POSITIVE(S_OLEN, "Datander minimum alignment length")
                break;
              case 'w':
                ARG_POSITIVE(S_BAND, "Datander band width")
                break;
              default:
                fprintf(stderr, "%s: option '-%s' not supported for datander\n",
                               Prog_Name, argv[i]);
                exit(1);
            }
            break;
          case 'L':
            argv[i] += 1;
            switch (argv[i][1])
            { case 'e':
                ARG_REAL(D_ERATE)
                if (D_ERATE < .7 || D_ERATE >= 1.)
                  { fprintf(stderr, "%s: Average correlation must be in [.7,1.) (%g)\n",
                                   Prog_Name, D_ERATE);
                    exit(1);
                  }
                break;
              case 'h':
                ARG_POSITIVE(D_HITS, "Daligner hit threshold")
                break;
              case 'k':
                ARG_POSITIVE(D_KMER, "Daligner kmer")
                break;
              case 'l':
                ARG_POSITIVE(D_OLEN, "Daligner minimum alignment length")
                break;
              case 't':
                ARG_POSITIVE(D_FREQ, "Daligner tuple suppresion frequency")
                break;
              case 'w':
                ARG_POSITIVE(D_BAND, "Daligner band width")
                break;
              default:
                fprintf(stderr, "%s: option '-%s' not supported for daligner\n",
                               Prog_Name, argv[i]);
                exit(1);
            }
            break;
        }
      else
        argv[j++] = argv[i];
    argc = j;

    VERBOSE = flags['v'];
    ALLBLKS = flags['l'];
    FORCE   = flags['F'];

    if (argc < 2 || argc > 3)
      { fprintf(stderr,"Usage: %s %s\n",Prog_Name,Usage[0]);
        fprintf(stderr,"       %*s %s\n",(int) strlen(Prog_Name),"",Usage[1]);
        fprintf(stderr,"       %*s %s\n",(int) strlen(Prog_Name),"",Usage[2]);
        fprintf(stderr,"       %*s %s\n",(int) strlen(Prog_Name),"",Usage[3]);
        fprintf(stderr,"       %*s %s\n",(int) strlen(Prog_Name),"",Usage[4]);
        fprintf(stderr,"\n");
        fprintf(stderr, "  Create an HPC workflow to estimate the coverage of a given DB:\n");
        fprintf(stderr, "  1. Create low complexity mask for the given DB (DBdust)\n");
        fprintf(stderr, "  2. Create tandem masks for each DB block (datander + TANmask)\n");
        fprintf(stderr, "  3. Create repeat mask for each DB block by self-comparison");
        fprintf(stderr,     " (daligner + REPmask)\n");
        fprintf(stderr, "  4. Compute all soft-masked overlaps for <target> block (default 1) vs");
        fprintf(stderr,     " all other DB blocks.\n");
        fprintf(stderr, "  5. Run DAScover on the target block .las file\n");
        fprintf(stderr,"\n");
        fprintf(stderr,"     DBdust parameters.\n");
        fprintf(stderr,"      -Uw: DUST algorithm window size.\n");
        fprintf(stderr,"      -Ut: DUST algorithm threshold.\n");
        fprintf(stderr,"      -Um: Record only low-complexity intervals >= this size.\n");
        fprintf(stderr,"      -Ub: Take into account base composition bias.\n");
        fprintf(stderr,"\n");
        fprintf(stderr,"     Datander parameters.\n");
        fprintf(stderr,"      -Sk: k-mer size (must be <= 32).\n");
        fprintf(stderr,"      -Sw: Look for k-mers in averlapping bands of size 2^-w.\n");
        fprintf(stderr,"      -Sh: A seed hit if the k-mers in band cover >= -h bps in the");
        fprintf(stderr,          " targest read.\n");
        fprintf(stderr,"      -Se: Look for alignments with -e percent similarity.\n");
        fprintf(stderr,"      -Sl: Look for alignments of length >= -l.\n");
        fprintf(stderr,"\n");
        fprintf(stderr,"     Daligner parameters.\n");
        fprintf(stderr,"      -Lk: k-mer size (must be <= 32).\n");
        fprintf(stderr,"      -Lw: Look for k-mers in averlapping bands of size 2^-w.\n");
        fprintf(stderr,"      -Lh: A seed hit if the k-mers in band cover >= -h bps in the");
        fprintf(stderr,          " targest read.\n");
        fprintf(stderr,"      -Le: Look for alignments with -e percent similarity.\n");
        fprintf(stderr,"      -Ll: Look for alignments of length >= -l.\n");
        fprintf(stderr,"      -Lt: Ignore k-mers that occur >= -t times in a block.\n");
        fprintf(stderr,"\n");
        fprintf(stderr,"     Common parameters.\n");
        fprintf(stderr,"      -s: Use -s as the trace point spacing for encoding alignments");
        fprintf(stderr,         " (daligner, datander).\n");
        fprintf(stderr,"      -M: Use only -M GB of memory by ignoring most frequent k-mers");
        fprintf(stderr,         " (daligner.\n");
        fprintf(stderr,"      -P: Use this directory for all scratch files");
        fprintf(stderr,         " (daligner, datander, LAsort, LAmerge).\n");
        fprintf(stderr,"      -T: Use -T threads (daligner, datander).\n");
        fprintf(stderr,"\n");
        fprintf(stderr,"     Script control.\n");
        fprintf(stderr,"      -v: Run all commands in script in verbose mode.\n");
        fprintf(stderr,"      -f: Place script bundles in separate files with prefix <name>\n");
        fprintf(stderr,"      -l: Include the last block in the calculation\n");
        fprintf(stderr,"      -B: # of block compares per command job\n");
        fprintf(stderr,"      -F: Restart computation from the beginnning\n");
        exit (1);
      }
  }

  //  Make sure DB exists and is partitioned, get number of blocks in partition

  pwd = PathTo(argv[1]);
  if (strcmp(argv[1] + (strlen(argv[1]) - 4), ".dam") == 0)
    { fprintf(stderr, "\n%s: A mapper database (.dam) is not supported!\n\n", Prog_Name);
      exit(1);
    }
  root = Root(argv[1], ".db");

  { int i, nfiles;
    FILE *dbvis;

    dbvis = Fopen(Catenate(pwd, "/", root, ".db"), "r");
    if (dbvis == NULL)
      exit(1);

    if (fscanf(dbvis, "files = %d\n", &nfiles) != 1)
      SYSTEM_READ_ERROR
    for (i = 0; i < nfiles; i++)
      { char buffer[30001];

        if (fgets(buffer, 30000, dbvis) == NULL)
          SYSTEM_READ_ERROR
      }

    if (fscanf(dbvis, "blocks = %d\n", &nblocks) != 1)
      { fprintf(stderr,"%s: Database must be split!\n",Prog_Name);
        exit (1);
      }

    usepath = (strcmp(pwd, ".") != 0);

    fclose(dbvis);
  }

  if (argc == 3)
    { char *eptr;

      BLOCKID = strtol(argv[2],&eptr,10);
      if (*eptr != '\0')
        { fprintf(stderr,"%s: target block '%s' should be an integer\n",Prog_Name,argv[2]);
          exit (1);
        }
    }

  if (BLOCKID > (nblocks-1) + ALLBLKS)
    { fprintf(stderr,"%s: Target block %d is greater than # of%sblocks in the db (%d)!\n",
                     Prog_Name, BLOCKID, ALLBLKS ?" ":" complete ", (nblocks-1)+ALLBLKS);
      exit(1);
    }

  // determine which blocks are new

  { char *tail;
    int   i;

    if (ALLBLKS)
      lblock = nblocks;
    else
      lblock = nblocks-1;

    if (lblock <= 0)
      { fprintf(stderr,"%s: DB has 1 block, if its complete use -l\n",Prog_Name);
        exit (1);
      }

    tail = Strdup(Numbered_Suffix(".",BLOCKID,".las"),"Creating .las string");
    for (i = nblocks; i >= 1; i--)
      if (access(Numbered_Suffix(Catenate(pwd,"/",root,"."),i,tail),F_OK) == 0)
        break;
    tblock = i;
    free(tail);

    if ( ! FORCE && tblock >= lblock)
      { fprintf(stderr,"%s: Nothing to do!\n",Prog_Name);
        exit (1);
      }

    fblock = 1;
    if ( ! FORCE)
      { for (i = 1; i <= nblocks; i++)
          { if (access(Numbered_Suffix(Catenate(pwd, "/.", root, "."), i, ".dust.anno"),F_OK) < 0)
              break;
            if (access(Numbered_Suffix(Catenate(pwd, "/.", root, "."), i, ".tan.anno"),F_OK) < 0)
              break;
            if (access(Numbered_Suffix(Catenate(pwd, "/.", root, "."), i, ".rep.anno"),F_OK) < 0)
              break;
          }
        fblock = i;
      }

    if (VERBOSE)
      { fprintf(stderr,"# Masking blocks %d-%d, add blocks from %d to target .las\n",
                       fblock,lblock,tblock+1);
      }
  }

  { FILE *out;
    char  name[1000];
    int   njobs;
    int   i, j, k;

    //  CREATE BLOCK MASKS

    if (fblock <= lblock)
      {
        // 01 create low complexity (dust) track

        if (PROJECT != NULL)
          { sprintf(name, "%s.01.DUST", PROJECT);
            out = fopen(name, "w");
          }
        else
          out = stdout;
      
        njobs = (lblock - fblock) + 1;
        if (VERBOSE)
          fprintf(stderr, "# 01 DBdust - jobs (%d)\n", njobs);
      
        for (i = fblock; i <= lblock; i++)
          { fprintf(out, "DBdust");
            if (D_BIAS)
              fprintf(out, " -b");
            if (D_WINDOW != 64)
              fprintf(out, " -w%d", D_WINDOW);
            if (D_THRESH >= 0.)
              fprintf(out, " -t%.2f", D_THRESH);
            if (D_MINBP != 10)
              fprintf(out, " -m%d", D_MINBP);
      
            if (usepath)
              fprintf(out, " %s/%s.%d", pwd, root, i);
            else
              fprintf(out, " %s.%d", root, i);
            fprintf(out, "\n");
          }
      
        // 02 create tandem repeat (tan) track - step 1 datander
      
        if (PROJECT != NULL)
          { fclose(out);
            sprintf(name, "%s.02.TANDEM", PROJECT);
            out = fopen(name, "w");
          }
      
        njobs = (lblock - fblock) / BUNIT + 1;
        if (VERBOSE)
          fprintf(stderr, "# 02 Datander - jobs (%d)\n", njobs);
      
        for (i = fblock; i <= lblock; i = j)
          { j = i + BUNIT;
            if (j > lblock + 1)
              j = lblock + 1;
      
            fprintf(out, "datander");
            if (VERBOSE)
              fprintf(out, " -v");
            if (S_KMER != 12)
              fprintf(out, " -k%d", S_KMER);
            if (S_BAND != 4)
              fprintf(out, " -w%d", S_BAND);
            if (S_HITS != 35)
              fprintf(out, " -h%d", S_HITS);
            if (S_ERATE > 0.)
              fprintf(out, " -e%g", S_ERATE);
            if (S_OLEN != 500)
              fprintf(out, " -l%d", S_OLEN);
            if (TRACE != 100)
              fprintf(out, " -s%d", TRACE);
            if (strcmp(TMPDIR, "/tmp"))
              fprintf(out, " -P%s", TMPDIR);
            if (NTHREADS != 4)
              fprintf(out, " -T%d", NTHREADS);
      
            for (k = i; k < j; k++)
              if (usepath)
                fprintf(out, " %s/%s.%d", pwd, root, k);
              else
                fprintf(out, " %s.%d", root, k);
            fprintf(out, "\n");
          }
      
        // 03 create tandem repeat (tan) track - step 2 check datander
      
        if (VERBOSE)
          fprintf(stderr, "# 03 Check all TAN.*.las - jobs (%d)\n", njobs);
      
        if (PROJECT != NULL)
          { fclose(out);
            sprintf(name, "%s.03.CHECK.TAN", PROJECT);
            out = fopen(name, "w");
          }
      
        for (i = fblock; i <= lblock; i = j)
          { j = i + BUNIT;
            if (j > lblock + 1)
              j = lblock + 1;
      
            fprintf(out, "LAcheck -S");
            if (VERBOSE)
              fprintf(out, " -v");
      
            if (usepath)
              fprintf(out, " %s/%s", pwd, root);
            else
              fprintf(out, " %s", root);
      
            fprintf(out, " TAN.%s.@%d-%d", root, i, j - 1);
            fprintf(out, "\n");
          }
      
        // 04 create tandem repeat (tan) track - step 3 TANmask
      
        if (VERBOSE)
          fprintf(stderr, "# 04 TANmask - jobs (%d)\n", njobs);
      
        if (PROJECT != NULL)
          { fclose(out);
            sprintf(name, "%s.04.TANMASK", PROJECT);
            out = fopen(name, "w");
          }
      
        for (i = fblock; i <= lblock; i = j)
          { j = i + BUNIT;
            if (j > lblock + 1)
              j = lblock + 1;
      
            fprintf(out, "TANmask");
            if (VERBOSE)
              fprintf(out, " -v");
            if (S_OLEN != 500)
              fprintf(out, " -l%d", S_OLEN);
            if (usepath)
              fprintf(out, " %s/%s", pwd, root);
            else
              fprintf(out, " %s", root);
      
            fprintf(out, " TAN.%s.@%d-%d", root, i, j - 1);
            fprintf(out, "\n");
          }
      
        // 05 create tandem repeat (tan) track - step 4 remove .las files
      
        if (VERBOSE)
          fprintf(stderr, "# 05 Remove all TAN.*.las jobs (%d)\n", njobs);
      
        if (PROJECT != NULL)
          { fclose(out);
            sprintf(name, "%s.05.RM.TAN", PROJECT);
            out = fopen(name, "w");
          }
      
        for (i = fblock; i <= lblock; i = j)
          { j = i + BUNIT;
            if (j > lblock + 1)
              j = lblock + 1;
      
            fprintf(out, "rm");
            for (k = i; k < j; k++)
              fprintf(out, " TAN.%s.%d.las", root, k);
            fprintf(out, "\n");
          }
      
        // REPEAT MASKING
        // 06 create overlap based repeat (rep) track - step 1 daligner
      
        njobs = (lblock - fblock) + 1;
        if (VERBOSE)
          fprintf(stderr, "# 06 Repeat Masking: Daligner - jobs (%d)\n", njobs);
      
        if (PROJECT != NULL)
          { fclose(out);
            sprintf(name, "%s.06.REPEAT", PROJECT);
            out = fopen(name, "w");
          }
      
        for (i = fblock; i <= lblock; i++)
          { fprintf(out, "daligner");
            if (VERBOSE)
              fprintf(out, " -v");
            if (D_KMER != 14)
              fprintf(out, " -k%d", D_KMER);
            if (D_BAND != 6)
              fprintf(out, " -w%d", D_BAND);
            if (D_HITS != 35)
              fprintf(out, " -h%d", D_HITS);
            if (D_ERATE > 0.)
              fprintf(out, " -e%g", D_ERATE);
            if (D_OLEN != 1000)
              fprintf(out, " -l%d", D_OLEN);
            if (TRACE != 100)
              fprintf(out, " -s%d", TRACE);
            if (strcmp(TMPDIR, "/tmp"))
              fprintf(out, " -P%s", TMPDIR);
            if (NTHREADS != 4)
              fprintf(out, " -T%d", NTHREADS);
      
            fprintf(out, " -mdust -mtan");    // add dust and tan track
            if (usepath)
              fprintf(out, " %s/%s.%d %s/%s.%d", pwd, root, i, pwd, root, i);
            else
              fprintf(out, " %s.%d %s.%d", root, i, root, i);
            fprintf(out, "\n");
          }
      
        // 07 create overlap based repeat (rep) track - step 2 check .las files
      
        njobs = (lblock - fblock) / BUNIT + 1;
        if (VERBOSE)
          fprintf(stderr, "# 07 Repeat Masking: Check all .las - jobs (%d) \n", njobs);
      
        if (PROJECT != NULL)
          { fclose(out);
            sprintf(name, "%s.07.CHECK.REP", PROJECT);
            out = fopen(name, "w");
          }
      
        for (i = fblock; i <= lblock; i = j)
          { j = i + BUNIT;
            if (j > lblock + 1)
              j = lblock + 1;
      
            fprintf(out, "LAcheck -S");
            if (VERBOSE)
              fprintf(out, " -v");
      
            if (usepath)
              fprintf(out, " %s/%s", pwd, root);
            else
              fprintf(out, " %s", root);
      
            for (k = i; k < j; k++)
              fprintf(out, " %s.%d.%s.%d", root, k, root, k);
            fprintf(out, "\n");
          }
      
        // 08 create overlap based repeat (rep) track  - step 3 TANmask
      
        if (PROJECT != NULL)
          { fclose(out);
            sprintf(name, "%s.08.REPMASK", PROJECT);
            out = fopen(name, "w");
          }
      
        if (VERBOSE)
          fprintf(stderr, "# 08 Repeat masking: REPmask - jobs (%d)\n", njobs);
      
        for (i = fblock; i <= lblock; i = j)
          { j = i + BUNIT;
            if (j > lblock + 1)
              j = lblock + 1;
      
            fprintf(out, "REPmask -c%d", COVERAGE);
            if (VERBOSE)
              fprintf(out, " -v");
      
            if (usepath)
              fprintf(out, " %s/%s", pwd, root);
            else
              fprintf(out, " %s", root);
      
            for (k = i; k < j; k++)
              fprintf(out, " %s.%d.%s.%d", root, k, root, k);
            fprintf(out, "\n");
          }
      
        // 09 create overlap based repeat (rep) track  - step 4 remove .las files
      
        if (VERBOSE)
          fprintf(stderr, "# 09 Repeat masking: Remove .las files\n");
      
        if (PROJECT != NULL)
          { fclose(out);
            sprintf(name, "%s.09.RM.REP", PROJECT);
            out = fopen(name, "w");
          }
      
        for (i = fblock; i <= lblock; i = j)
          { j = i + BUNIT;
            if (j > lblock + 1)
              j = lblock + 1;
      
            fprintf(out, "rm");
            for (k = i; k < j; k++)
              fprintf(out, " %s.%d.%s.%d.las", root, k, root, k);
            fprintf(out, "\n");
          }

        if (FORCE && tblock >= 1)
          { fprintf(out, "rm %s.%d.%d.las\n", root, tblock, BLOCKID);
            tblock = 0;
          }
      }
  
    // DASCOVER PROPER
  
    njobs = (lblock - (tblock+1)) / BUNIT + 1;
    if (VERBOSE)
      fprintf(stderr, "# 10 DAScover: daligner - jobs (%d)\n", njobs);
  
    if (PROJECT != NULL)
      { fclose(out);
        sprintf(name, "%s.10.DALIGNER", PROJECT);
        out = fopen(name, "w");
      }
  
    for (i = tblock+(BLOCKID <= tblock); i <= lblock; i = j)
      { if (i == tblock)
          j = tblock+1;
        else
          { j = i + BUNIT;
            if (j > lblock + 1)
              j = lblock + 1;
            if (j == i+1 && i == BLOCKID)
              continue;
          }
  
        fprintf(out, "daligner -mdust -mtan -mrep");
        if (VERBOSE)
          fprintf(out, " -v");
        if (D_KMER != 14)
          fprintf(out, " -k%d", D_KMER);
        if (D_BAND != 6)
          fprintf(out, " -w%d", D_BAND);
        if (D_HITS != 35)
          fprintf(out, " -h%d", D_HITS);
        if (D_ERATE > 0.)
          fprintf(out, " -e%g", D_ERATE);
        if (D_OLEN != 1000)
          fprintf(out, " -l%d", D_OLEN);
        if (D_FREQ > 0)
          fprintf(out, " -t%d", D_FREQ);
        if (TRACE != 100)
          fprintf(out, " -s%d", TRACE);
        if (strcmp(TMPDIR, "/tmp"))
          fprintf(out, " -P%s", TMPDIR);
        if (NTHREADS != 4)
          fprintf(out, " -T%d", NTHREADS);
        if (i != tblock)
          fprintf(out," -A");
  
        if (usepath)
          fprintf(out, " %s/%s.%d", pwd, root, BLOCKID);
        else
          fprintf(out, " %s.%d", root, BLOCKID);
  
        for (k = i; k < j; k++)
          { if (k == BLOCKID)
              continue;
            if (i == tblock)
              k = BLOCKID;
            if (usepath)
              fprintf(out, " %s/%s.%d", pwd, root, k);
            else
              fprintf(out, " %s.%d", root, k);
          }
        fprintf(out, "\n");
      }
  
    //  Check .las files (optional)
  
    if (VERBOSE)
      fprintf(stderr, "# 11 DAScover: Check new .las files - jobs (%d)\n", njobs);
  
    if (PROJECT != NULL)
      { fclose(out);
        sprintf(name, "%s.11.CHECK.DAL", PROJECT);
        out = fopen(name, "w");
      }
  
    for (i = tblock+1; i <= lblock; i = j)
      { j = i + BUNIT;
        if (j > lblock + 1)
          j = lblock + 1;
  
        fprintf(out, "LAcheck -S");
        if (VERBOSE)
          fprintf(out, " -v");
  
        if (usepath)
          fprintf(out, " %s/%s", pwd, root);
        else
          fprintf(out, " %s", root);
  
        fprintf(out, " %s.%d.%s.@%d-%d", root, BLOCKID, root, i, j - 1);
        fprintf(out, "\n");
      }
  
    // LAmerge
  
    if (VERBOSE)
      fprintf(stderr, "# 12 DAScover: LAmerge final .las\n");
  
    if (PROJECT != NULL)
      { fclose(out);
        sprintf(name, "%s.12.MERGE", PROJECT);
        out = fopen(name, "w");
      }
  
    fprintf(out, "LAmerge");
    if (VERBOSE)
      fprintf(out, " -v");
  
    if (usepath)
      fprintf(out, " %s/%s.%d.%d", pwd, root, lblock, BLOCKID);
    else
      fprintf(out, " %s.%d.%d", root, lblock, BLOCKID);

    if (tblock >= 1)
      { if (usepath)
          fprintf(out, " %s/%s.%d.%d", pwd, root, tblock, BLOCKID);
        else
          fprintf(out, " %s.%d.%d", root, tblock, BLOCKID);
      }

    fprintf(out, " %s.%d.%s.@%d-%d", root, BLOCKID, root, tblock+1, lblock);
    fprintf(out, "\n");
  
    // check final las file
  
    if (VERBOSE)
      fprintf(stderr, "# 13 DAScover: Check final merged .las\n");
  
    if (PROJECT != NULL)
      { fclose(out);
        sprintf(name, "%s.13.CHECK.MRG", PROJECT);
        out = fopen(name, "w");
      }
  
    fprintf(out, "LAcheck -S");
    if (VERBOSE)
      fprintf(out, " -v");
 
    if (usepath)
      fprintf(out, " %s/%s %s/%s.%d.%d", pwd, root, pwd, root, lblock, BLOCKID);
    else
      fprintf(out, " %s %s.%d.%d", root, root, lblock, BLOCKID);
    fprintf(out, "\n");
  
    // 14 DAScover - step 4 remove DB.BLOCKID.DB.k.las files
  
    if (VERBOSE)
      fprintf(stderr, "# 14 DAScover: Remove intermediate .las files\n");
  
    if (PROJECT != NULL)
      { fclose(out);
        sprintf(name, "%s.14.RM.DAL", PROJECT);
        out = fopen(name, "w");
      }
  
    if (tblock >= 1)
      { if (usepath)
          printf("rm %s/%s.%d.%d.las\n",pwd,root,tblock,BLOCKID);
        else
          printf("rm %s.%d.%d.las\n",root,tblock,BLOCKID);
      }

    for (i = fblock; i <= lblock; i = j)
      { j = i + BUNIT;
        if (j > lblock + 1)
          j = lblock + 1;

        fprintf(out, "rm");
        for (k = i; k < j; k++)
          fprintf(out, " %s.%d.%s.%d.las", root, BLOCKID, root, k);
        fprintf(out, "\n");
      }
  
    // DAScover
  
    if (VERBOSE)
      fprintf(stderr, "# 15 DAScover\n");
  
    if (PROJECT != NULL)
      { fclose(out);
        sprintf(name, "%s.15.DASCOVER", PROJECT);
        out = fopen(name, "w");
      }
  
    fprintf(out, "DAScover  -v");
    if (usepath)
      fprintf(out, " %s/%s %s/%s.%d.%d.las", pwd, root, pwd, root, lblock, BLOCKID);
    else
      fprintf(out, " %s %s.%d.%d.las", root, root, lblock, BLOCKID);
   
    fprintf(out, "\n");
  
    if (PROJECT != NULL)
      fclose(out);
  }

  exit (0);
}
