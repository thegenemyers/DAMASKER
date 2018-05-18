
# Damasker: The Dazzler Repeat Masking Suite

## _Author:  Gene Myers_
## _First:   April 10, 2016_

For typeset documentation, examples of use, and design philosophy please go to
my [blog](https://dazzlerblog.wordpress.com/command-guides/damasker-commands).

All programs add suffixes (e.g. .db, .las) as needed.
For the commands that take multiple .las files as arguments, i.e. REPmask
and TANmask, one can place a @-sign in the name, which is then interpreted as the sequence of files
obtained by replacing the @-sign by 1, 2, 3, ... in sequence until a number is reached for which
no file exists with that name.  One can also place a @-sign followed by an integer, say, i, in which
case the sequence starts at i.  Lastly, one can also place @i-j where i and j are integers, in
which case the sequence is from i to j, inclusive.

Scrubbing is complicated by the presence of repeats.  We currently handle this by soft-masking all
tandem and interspersed repeats in the input data when computing overlaps.  This implies that reads
that are completely repetitive sequence are not scrubbed.  This is typically a small part of the
data set and a portion thereof that is currently not correctly assembled by any assembly system at
the current time, and therefore the masking is of minor consequence.  Eventually, these completely
masked reads will be analyzed in downstream processes that will attempt to resolve ultra-long
(15Kbp or more) repeats.

The masking suite therefore consists of several programs that in combination can be used to
produce repeat masks for a data set as follows:

```
1.  REPmask [-v] [-n<track(rep)>] -c<int> <subject:db> <overlaps:las> ...
```

This command takes as input a database \<source\> and a sequence of sorted local alignments blocks, \<overlaps\>, produced by a daligner run for said database.  Note carefully that \<source\> must always refer to the entire DB, only \<overlaps\> can involve a block number.

REPmask examines each pile for an A-read and determines the intervals that are covered -c or more times by LAs.  This set of intervals is output as a repeat mask for A in an interval track with default name .rep, that can be overridden with the -n option.  If the -v option is set, then the number of intervals and total base pairs in intervals is printed.

```
2. datander [-v] [-k<int(12)>] [-w<int(4)>] [-h<int(35)>] [-T<int(4)>]
                 [-e<double(.70)>] [-l<int(1000)>] [-s<int(100)>] [-P<dir(/tmp)>]
                 <path:db|dam> ...
```

This program is a variation of daligner tailored to the task of comparing each read against itself (and only those comparisons).   As such each block or DB serves as both the source and target, and the -b, -A, -I, -t, -M, -H, and -m options are irrelevant.  The remaining options are exactly as for daligner (see here).  For each subject block, say X, this program produces a single file TAN.X.las where all the alignments do not involve complementing the B-read (which is also the A-read).

```
3. TANmask [-v] [-l<int(500)>] [-n<track(tan)>] <subject:db> <overlaps:las> ...
```

This command takes as input a database \<source\> and a sequence of sorted local alignments blocks, \<overlaps\>, produced by a datander run for said database.  Note carefully that \<source\> must always refer to the entire DB, only \<overlaps\> can involve a block number.

TANmask examines each pile for an A-read and finds those self-LAs whose two alignment intervals overlap and for which the union of these two intervals is -l bases or longer.  Each of these regions signals a tandem element in A of length -l or greater, and a disjoint list of these is built.  This set of intervals is output as a tandem mask for A in an interval track with default name .tan, that can be overridden with the -n option.  If the -v option is set, then the number of intervals and total base pairs in intervals is printed.

```
4. HPC.REPmask [-vbd]
               [-t<int>] [-w<int(6)>] [-l<int(1000)>] [-s<int(100)>] [-M<int>]
               [-n<name(rep-g)>] [-P<dir(/tmp)>] [-B<int(4)>] [-T<int(4)>] [-f<name>] 
               [-k<int(14)>] [-h<int(35)>] [-e<double(.70)>] [-m<track>]+
               -g<int> -c<int> <reads:db|dam> [<first:int>[-<last:int>]]
```

HPC.REPmask writes a UNIX shell script to the standard output that consists of a sequence of
commands that effectively runs daligner and LAmerge to compare all consecutive groups
of -g blocks against each other and then applies REPmask to the result in order to generate a
repeat mask for the database \<path\>.   For example if -g is 3, then it will compare blocks
1-3 against each other, blocks 4-6 against each other, and so on with daligner, and then sort
and merge so that all the alignments with an A-read in block i are in the file \<path\>.i.las
(e.g. \<path\>.2.las will contain all alignments where the A-read is in block 2 and the B-read
is in blocks 1, 2, or 3).  Thereafter "REPmask \<-c\> \<-n\> \<path\> \<path\>.i.las" is run
for every block i, resulting in a .\<-n\> block track for each block that can then be combined with
Catrack into a single track for the entire DB.

The data base must have been previously split by DBsplit and all options, except -B, -d, and -f are passed through to the calls to daligner or REPmask as appropriate. The defaults for these parameters are as for daligner and REPmask. The -v flag, for verbose-mode, is passed to all commands.  The -d and -f parameters are explained later.  The -B option controls the form of the script generated by HPC.REPmask as follows.  The -B option determines the number of block comparisons per daligner job.
For a database divided into N sub-blocks, the calls to daligner will produce a total of gN .las files, g<sup>2</sup> for each of the N/g block groups. These are then merged so that there is 1 file per block. So at the end one has N sorted .las files, one per block.

If the integers \<first\> and \<last\> are missing then the script produced is for every block in the database. If \<first\> is present then HPC.REPmask produces an incremental script that performs as much of the task as possible for blocks \<first\> through \<last\> (\<last\> = \<first\> if not present) where it is assumed that the script has been called previously up to block \<first\>-1.  If \<last\> is not evenly divisible by -g then the script performs the block comparisons and merges for the last incomplete group but does not yet invoke REPmask as all the necessary comparisons have not been made.  Symmetrically, if the initial part of the range completes a previously incomplete block group, then the script generates calls to complete the comparisons for the initial group and produces the repeat-masks for that group.

The command script output by HPC.REPmask and other HPC.\<x\> programs consists of command blocks each of which begins with a comment line (begins with #) followed by a potentially long list of lines each containing a shell command.  Command blocks whose comment mentions "jobs" and gives the number of said in parenthesis, we call parallel blocks because each command line in the block can be sent to a node in a cluster for independent execution, i.e. none of the commands in a block depend on another in the block.  The remaining command blocks we call house-keeping blocks because they can be executed by the shell on the launch/server node and the commands are either checking the integrity of .las files with LAcheck, or removing intermediate files with rm. Each block should be performed in the order given and should complete before the next block is performed.

If the -f option is set, then each command block is written to a file with a name of the form \<name\>.#.\<description\> where \<name\> is specified by the user in the -f option argument, # gives the order in which the command block in the given file is to be performed in relation to other command block files, and \<description\> is a (very) short symbolic reminder of what the block is doing.  For example, "HPC.REPmask -g3 -c10 -fJOBS DB" would produce the files:

```
  JOBS.01.OVL
  JOBS.02.CHECK.OPT
  JOBS.03.MERGE
  JOBS.04.RM
  JOBS.05.MASK
  JOBS.06.RM
```

The number of command blocks varies as it depends on the number of merging rounds required in the external sort of the .las files.  The files with the suffix .OPT are optional and need not be executed albeit we highly recommend that one run the CHECK block.

The -d option requests scripts that organize files into a collection of sub-directories so as not to overwhelm the underlying OS for large genomes.  For a DB divided into N blocks and the daligner calls in the script will produce 2gNT .las-files where T is the number of threads specified by the -T option passed to daligner (default is 4).  With the -d option set, N sub-directories (with respect to the directory HPC.daligner is called in) of the form "temp\<i\>" for i from 1 to N are created in an initial command block, and then all intermediate files are placed in those sub-directories, with a maximum of g(2T+1) files appearing in any sub-directory at any given point in the process.

```
5. HPC.TANmask [-v] [-k<int(12)>] [-w<int(4)>] [-h<int(35)>] [-T<int(4)>] [-P<dir(/tmp)>]
                    [-n<name(tan)>] [-e<double(.70)>] [-l<int(1000)>] [-s<int(100)>] [-f<name>]
                    <reads:db|dam> [<first:int>[-<last:int>]]
```

HPC.TANmask writes a UNIX shell script to the standard output that runs datander on all relevant blocks of the supplied DB, then sorts and merges the resulting alignments into a single .las for each block, and finally calls TANmask on each LA block to produce a tandem mask with name \<-n\> for each block that can be merge into a single track for the entire DB with Catrack.

All option arguments are passed through to datander or TANmask except for -l which is passed to both, and except for the -f option which serves the same role as for HPC.REPmask above.  The -v option is passed to all programs in the script.  If the integers \<first\> and \<last\> are missing then the script produced is for every block in the database \<reads\>. If \<first\> is present then HPC.TANmask produces a script that produces .tan tracks for blocks \<first\> through \<last\> (\<last\> = \<first\> if not present).

```
6. HPC.DAScover [-vlF] [-U(w<int(64)> |t<double(.2)> |m<int(10)> |b) ]
                       [-S(k<int(12)> |w<int(4)> |h<int(35)> |e<double(.7)> |l<int(500)>) ]
                       [-L(k<int(14)> |w<int(6)> |h<int(35)> |e<double(.7)> |l<int(1000)> |t<int>) ]
                       [-c<int(10)>] [-s<int(100)] [-M<int>] [-P<dir(/tmp)>] [-T<int(4)>]
                       [-B<int(4)>] [-f<name>] <reads:db|dam> [<target:int(1)>]
```

HPC.DAScover writes a UNIX shell script to the standard output that incrementally updates the .las
file for a given target block against all others in a growing database.   As such there are files
that persist between executions of HPC.DAScover scripts.  After executing the script for target
block, say T, on a data base DB with N blocks, the following files will be present (in the
directory containing the database): (a) dust, tan, and rep block tracks for blocks 1 to N-1, and
(b) DB.<N-1>.T.las that contains all the LAs between reads in block 1 and all reads in blocks
1 to N-1.  The last block, N, is not compared because one cannot be certain it will not change
when new data is added to DB.  The -l option forces that last block to be included, but should
only be invoked when one is certain that no more data will be added to the database.  If necessary,
one can produce a script that starts from the beginning by setting the -F option.

One can produce scripts that operate on *different* target blocks for the same DB.  The desired
target block is the optional second argment to HPC.DAScover.  The production of the block tracks
for the DB are common to all  computations, and not performed again for each target block.  Only
the comparisons of the target block versus all other blocks are required.

The UNIX shell script produced by HPC.DAScover invokes many Dazzler commands including DBdust,
datander, REPmask, daligner, and LAmerge.  Many of the options to HPC.DAScover are directed to
these underlying commands as follows.
All the options beginning with -U are passed to DBdust with the U removed.  Similarly, -S options
are passed to datander, and -L options are directed at all daligner calls.  The -M option sets
the memory limit for all daligner calls.  The -P option sets the scratch directory for daligner,
datander, and LAmerge.  The -T option sets the number of threads and the -s options sets the
trace point spacing for all daligner and datander calls.
The -c option is passed to REPmask calls.

The -B option controls the number of blocks compared or analyzed in any given command line
call (except for LAmerge).  And as for other HPC script generators, the -f\<n\> option directs
the script into a series of files whose names begin with \<n\> that should then be performed
in sequence.
