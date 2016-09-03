
# DA MASKER: The Dazzler Repeat Masking Suite

## _Author:  Gene Myers_
## _First:   April 10, 2016_

<p style="text-align: justify;">
Scrubbing is complicated by the presence of repeats.  We currently handle this by soft-masking all tandem and interspersed repeats in the input data when computing overlaps.  This implies that reads that are completely repetitive sequence are not scrubbed.  This is typically a small part of the data set and a portion thereof that is currently not correctly assembled by any assembly system at the current time, and therefore the masking is of minor consequence.  Eventually, these completely masked reads will be analyzed in downstream processes that will attempt to resolve ultra-long (15Kbp or more) repeats.
</p>

<p style="text-align: justify;">
The masking suite therefore consists of several programs that in combination can be used to produce repeat masks for a data set as follows:
</p>

```
1.  REPmask [-v] [-m<track(rep)>] -c<int> <subject:db> <overlaps:las> ...
```

<p style="text-align: justify; margin-left: 40px;">
This command takes as input a database <source> and a sequence of sorted local alignments blocks, <overlaps>, produced by a daligner run for said database.  Note carefully that <source> must always refer to the entire DB, only <overlaps> can involve a block number.
</p>

<p style="text-align: justify; margin-left: 50px;">
REPmask examines each pile for an A-read and determines the intervals that are covered -c or more times by LAs.  This set of intervals is output as a repeat mask for A in an interval track with default name .rep, that can be overridden with the -m option.  If the -v option is set, then the number of intervals and total base pairs in intervals is printed.
</p>

```
2. tander [-v] [-k<int(12)>] [-w<int(4)>] [-h<int(35)>] [-T<int(4)>]
               [-e<double(.70)] [-l<int(1000)] [-s<int(100)]
               <path:db|dam> ...
```

...This program is a variation of daligner tailored to the task of comparing each read against itself (and only those comparisons).   As such each block or DB serves as both the source and target, and the -b, -A, -I, -t, -M, -H, and -m options are irrelevant.  The remaining options are exactly as for daligner (see here).

...For each subject block, say X, this program produces just 4NTHREAD files X.T#.las where -T is the number of threads (by default 4) and where all the alignments do not involve complementing the B-read (which is also the A-read).  These should then be sorted and merged with LAsort and LAmerge as for example performed by the script generator HPC.TANmask described below.

```
3. TANmask [-v] [-l<int(500)>] [-m<track(tan)>] <subject:db> <overlaps:las> ...
```

...This command takes as input a database <source> and a sequence of sorted local alignments blocks, <overlaps>, produced by a datander run for said database.  Note carefully that <source> must always refer to the entire DB, only <overlaps> can involve a block number.

...TANmask examines each pile for an A-read and finds those self-LAs whose two alignment intervals overlap and for which the union of these two intervals is -l bases or longer.  Each of these regions signals a tandem element in A of length -l or greater, and a disjoint list of these is built.  This set of intervals is output as a tandem mask for A in an interval track with default name .tan, that can be overridden with the -m option.  If the -v option is set, then the number of intervals and total base pairs in intervals is printed.

```
4. HPC.REPmask [-vbd]
               [-t<int>] [-w<int(6)>] [-l<int(1000)>] [-s<int(100)>]
               [-M<int>] [-B<int(4)>] [-D<int( 250)>] [-T<int(4)>] [-f<name>] 
               [-k<int(14)>] [-h<int(35)>] [-e<double(.70)>] [-m<track>]+
               -g<int> -c<int> <reads:db|dam> [<first:int>[-<last:int>]]
```

...HPC.REPmask writes a UNIX shell script to the standard output that consists of a sequence of commands that effectively runs daligner, LAsort, and LAmerge to compare all consecutive groups of -g blocks against each other and then applies REPmask to the result in order to generate a repeat mask for the database <path>.   For example if -g is 3, then it will compare blocks 1-3 against each other, blocks 4-6 against each other, and so on with daligner, and then sort and merge so that all the alignments with an A-read in block i are in the file <path>.i.las (e.g. <path>.2.las will contain all alignments where the A-read is in block 2 and the B-read is in blocks 1, 2, or 3).  Thereafter "REPmask -c -mrep<-g> <path> <path>.i.las" is run for every block i, resulting in a .rep<-g> block for each mask that can then be combined with Catrack into a single track for the entire DB.

The data base must have been previously split by DBsplit and all options, except -B, -D, -d, and -f are passed through to the calls to daligner. The defaults for these parameters are as for daligner. The -v flag, for verbose-mode, is also passed to all calls to LAsort and LAmerge.  The -d and -f parameters are explained later.  The -B and -D options control the form of the script generated by HPC.REPmask as follows.  The -B option determines the number of block comparisons per daligner job and the -D option determines the fan-in for the hierarchical LAmerge process exactly as for HPC.daligner.  For a database divided into N sub-blocks, the calls to daligner will produce a total of 2g^2(N/g)T .las files when daligner is run with T threads as specified by the -T option.  These will then be sorted and merged into Ng sorted .las files, g^2 for each of the N/g block groups. These are then merged in ceil(log_D g) phases where the number of files decreases geometrically in -D until there is 1 file per block. So at the end one has N sorted .las files, one per block.

...If the integers <first> and <last> are missing then the script produced is for every block in the database. If <first> is present then HPC.REPmask produces an incremental script that performs as much of the task as possible for blocks <first> through <last> (<last> = <first> if not present) where it is assumed that the script has been called previously up to block <first>-1.  If <last> is not evenly divisible by -g then the script performs the block comparisons and merges for the last incomplete group but does not yet invoke REPmask as all the necessary comparisons have not been made.  Symmetrically, if the initial part of the range completes a previously incomplete block group, then the script generates calls to complete the comparisons for the initial group and produces the repeat-masks for that group.

...The command script output by HPC.REPmask and other HPC.<x> programs consists of command blocks each of which begins with a comment line (begins with #) followed by a potentially long list of lines each containing a shell command.  Command blocks whose comment mentions "jobs" and gives the number of said in parenthesis, we call parallel blocks because each command line in the block can be sent to a node in a cluster for independent execution, i.e. none of the commands in a block depend on another in the block.  The remaining command blocks we call house-keeping blocks because they can be executed by the shell on the launch/server node and the commands are either checking the integrity of .las files with LAcheck, or removing intermediate files with rm. Each block should be performed in the order given and should complete before the next block is performed.

...If the -f option is set, then each command block is written to a file with a name of the form <name>.#.<description> where <name> is specified by the user in the -f option argument, # gives the order in which the command block in the given file is to be performed in relation to other command block files, and <description> is a (very) short symbolic reminder of what the block is doing.  For example, "HPC.REPmask -g3 -c10 -fJOBS DB" would produce the files:

```
  JOBS.01.OVL
  JOBS.02.SORT
  JOBS.03.CHECK.OPT
  JOBS.04.RM
  JOBS.MERGE
  JOBS.06.CHECK.OPT
  JOBS.07.RM
  JOBS.08.MASK
  JOBS.09.RM
```

...The number of command blocks varies as it depends on the number of merging rounds required in the external sort of the .las files.  The files with the suffix .OPT are optional and need not be executed albeit we highly recommend that one run all the CHECK blocks.

...The -d option requests scripts that organize files into a collection of sub-directories so as not to overwhelm the underlying OS for large genomes.  For a DB divided into N blocks and the daligner calls in the script will produce 2gNT .las-files where T is the number of threads specified by the -T option passed to daligner (default is 4).  With the -d option set, N sub-directories (with respect to the directory HPC.daligner is called in) of the form "temp<i>" for i from 1 to N are created in an initial command block, and then all intermediate files are placed in those sub-directories, with a maximum of g(2T+1) files appearing in any sub-directory at any given point in the process.

```
5. HPC.TANmask [-vd] [-k<int(12)>] [-w<int(4)>] [-h<int(35)>] [-T<int(4)>]
                     [-e<double(.70)] [-l<int(1000)] [-s<int(100)] [-f<name>]
                     <reads:db|dam> [<first:int>[-<last:int>]]
```

...HPC.TANmask writes a UNIX shell script to the standard output that runs datander on all relevant blocks of the supplied DB, then sorts and merges the resulting alignments into a single .las for each block, and finally calls TANmask on each LA block to produce a .tan mask for each block that can be merge into a single track for the entire DB with Catrack.

...All option arguments are passed through to datander except for the -d and -f options which serve the same role as for HPC.REPmask above.  The -v option is passed to all programs in the script, and the -l option is also passed to TANmask.  If the integers <first> and <last> are missing then the script produced is for every block in the database <reads>. If <first> is present then HPC.TANmask produces a script that produces .tan tracks for blocks <first> through <last> (<last> = <first> if not present).
