#!/usr/bin/env python

###########################
# plink_IBS_bootstraps.py #
###########################

'''
a script to calculate pairwise distance matrices on a PLINK dataset with bootstrap replicates

performs bootstrap replicates of a PLINK dataset (bed, bim, bam)
[http://www.cog-genomics.org/plink/1.9/],
calculates a pairwise Identity by State (IBS) distance matrix for each replicate
using the "--distance square 1-ibs flat-missing" command
and combines the output into a single PHYLIP file (multiple matrix file for input in 'neighbor')
[http://evolution.genetics.washington.edu/phylip/doc/neighbor.html]

requires:
- PLINK 1.9 (in system path)
- pyplink (https://github.com/lemieuxl/pyplink)

to set up a conda environment with these requirements:
$ conda create --name plink_bootstraps
$ source activate plink_bootstrsaps
$ conda install plink -c bioconda
$ conda install pyplink -c http://statgen.org/wp-content/uploads/Softwares/pyplink

run the script with the command:
$ python plink_IBS_bootstraps.py <plink dataset prefix> <number of desired bootstrap replicates>

'''

import argparse
from pyplink import PyPlink
import random
import subprocess
import numpy as np


def get_resampled_loci(number_markers):
    '''
    returns a list of loci positions resampled with replacement
    '''
    bootstrapped_list = [random.randint(0, number_markers-1) for _ in range(number_markers)]
    bootstrapped_list.sort()
    return bootstrapped_list


### parsing command line arguments
parser = argparse.ArgumentParser()
parser.add_argument('bfile', help='plink binary (bed, bim, fam) dataset prefix')
parser.add_argument('nboot', help='number of bootstrap replicates', type=int)
args = parser.parse_args()


### perform plink bootstrap replicates
with PyPlink(args.bfile) as bed:
    bim = bed.get_bim() # returns pandas.DataFrame of bim file

    nmarkers = bed.get_nb_markers()
    nsamples = bed.get_nb_samples()
    print("### Loaded {0} markers and {1} samples...".format(nmarkers, nsamples))
        
    for rep in range(args.nboot):
        print("### Performing bootstrap replicate {0} of {1}...".format(rep+1, args.nboot))
        rep_list = get_resampled_loci(nmarkers) # gets a list of resampled markers
        
        rep_basename = "rep" + str(rep)
        with PyPlink(rep_basename, "w") as outbed, open(rep_basename + ".bim", "w") as outbim:

            for marker_position in rep_list:                
                bed.seek(marker_position)
                marker, genotypes = bed.next()
                # write marker info to outbim
                marker_info = bim.loc[marker]
                marker_line = '\t'.join([_ for _ in map(str, marker_info.tolist())]) + '\n'
                outbim.write(marker_line)
                # write genotypes to outbed
                outbed.write_genotypes(genotypes)

        # create copy of fam file with rep_basename.fam
        cp_cmd_line = "cp {0}.fam {1}.fam".format(args.bfile, rep_basename)
        subprocess.call(cp_cmd_line.split())
        # call plink to compute distance matrix on rep dataset
        plink_cmd_line = "plink --dog --bfile {st} --distance square 1-ibs flat-missing --out {st}".format(st=rep_basename)
        subprocess.call(plink_cmd_line.split(), stdout=open("/dev/null", "w"))
        # clean created files, leaving only matrices (.mdist)
        clean_cmd_line = "rm {st}.bed {st}.bim {st}.fam {st}.nosex {st}.log {st}.mdist.id".format(st=rep_basename)
        subprocess.call(clean_cmd_line.split())
 

### concatenate and format all generated matrices for input to PHYLIP

# get sample names from .fam file
samples = []
with open(args.bfile+".fam", "r") as fam:
    for line in fam:
        sample_id = line.strip().split()[1]
        samples.append(sample_id)

with open("infile", "w") as outf:

    for rep in range(args.nboot):
        # load each replicate matrix as a numpy array
        rep_file = "rep"+str(rep)+".mdist"
        rep_matrix = np.loadtxt(rep_file)
        
        # write matrix in PHYLIP format
        outf.write("    {}\n".format(str(len(samples))))
        for i, sample in enumerate(samples):
            matrix_line = np.array2string(rep_matrix[i,], separator='  ', threshold=5000, formatter={'float': lambda x: "{0:0.6f}".format(x)})[1:-1].replace('\n', '')
            out_line = sample.ljust(12) + matrix_line + "\n"
            outf.write(out_line)
        outf.write("\n")