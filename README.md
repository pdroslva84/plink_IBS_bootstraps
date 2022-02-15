# plink_IBS_bootstraps.py

A script to calculate pairwise distance matrices on a PLINK dataset with bootstrap replicates.


Performs bootstrap replicates of a PLINK dataset (bed, bim, bam) [http://www.cog-genomics.org/plink/1.9/], calculates a pairwise Identity by State (IBS) distance matrix for each replicate using the "--distance square 1-ibs flat-missing" command and combines the output into a single PHYLIP file (multiple matrix file for input in 'neighbor') [http://evolution.genetics.washington.edu/phylip/doc/neighbor.html]

Requires:
- PLINK 1.9 (in system path)
- pyplink (https://github.com/lemieuxl/pyplink)

To set up a conda environment with these requirements:

```
$ conda create --name plink_bootstraps
$ source activate plink_bootstrsaps
$ conda install plink -c bioconda
$ conda install pyplink -c http://statgen.org/wp-content/uploads/Softwares/pyplink
```


Run the script with the command:
```
$ python plink_IBS_bootstraps.py <plink dataset prefix> <number of desired bootstrap replicates>
```

This script was written in Python 2, which is no longer maintained.