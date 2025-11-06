# MMAPPR2
[![Build Status](https://travis-ci.org/kjohnsen/MMAPPR2.svg?branch=master)](https://travis-ci.org/kjohnsen/MMAPPR2)

## Mutation Mapping Analysis Pipeline for Pooled RNA-Seq
### Kyle Johnsen, Nathaniel Jenkins, Jonathon Hill

### Introduction
MMAPPR2 maps mutations resulting from pooled RNA-seq data from the F2
cross of forward genetic screens. Its predecessor is described in a paper published
in Genome Research (Hill et al. 2013). MMAPPR2 accepts aligned BAM files as well as
a reference genome as input, identifies loci of high sequence disparity between the
control and mutant RNA sequences, predicts variant effects, 
and outputs a ranked list of candidate mutations.

[See vignette for instructions](vignettes/MMAPPR2.Rmd)

Publication for the [original MMAPPR](http://genome.cshlp.org/content/23/4/687.full.pdf)

## Installation Notes
MMAPPR2 depends on Samtools to function. It must be installed and in the PATH to be found by the appropriate functions.

### Installing Samtools
Instructions to install samtools can be found at https://github.com/samtools/samtools and installation instructions are in the INSTALL file included with samtools.

