TEToolkit
=========

    TEToolkit is composed of two tools, TEpeaks and TEtranscripts, each described in
    its own section below.

    *NOTE* Both programs rely on specially curated GTF files, which are not
    packaged with this software due to their size. Please go to 
    http://hammelllab.labsites.cshl.edu/software#TEToolkit 
    for instructions to download the curated GTF files.

    TEpeak takes ChIP-seq (and similar data) alignment files (BAM or BED),
    identiifes narrow peaks, and is also able to do differential analysis over
    peaks of two sets of libraries. It is an extension of MACS by adding the
    funcionality of taking into account multi-reads, another normalization
    method, bin correlation, and differential analysis. The differential
    analysis is performed using DESeq. 

    TEtranscripts takes RNA-seq (and similar data) and annotates reads to both
    genes & transposable elements. It then performs differential analysis using
    DESeq.


+ [Github Page](https://github.com/mhammell-laboratory/tetoolkit)
+ [Pypi Page](https://pypi.python.org/pypi/TEToolkit)
+ [MHammell Lab](http://hammelllab.labsites.cshl.edu/software)

Created by Ying Jin, Eric Paniagua, Oliver Tam & Molly Hammell, February 2014

Copyright (C) 2014-2015 Ying Jin, Eric Paniagua, Oliver Tam & Molly Hammell
Contact: Ying Jin (yjin@cshl.edu)

Requirements
------------

    Python:     2.6.x or 2.7.x (not tested in Python 3.x)
    pysam    (tested using version 0.8.2.1)
    R:          2.15.x or greater
    DESeq:      1.5.x or greater


Installation
------------

    1) Download compressed tarball.
    2) Unpack tarball.
    3) Navigate into unpacked directory.
    4) Run the following:
       $ python setup.py install

    If you want to install locally (e.g. /local/home/usr),
       run this command instead:
    $ python setup.py install --prefix /local/home/usr

    *NOTE* In the above example, you must add /local/home/usr/bin
    to the PATH variable, and 
    /local/home/usr/lib/python2.X/site-packages 
    to the PYTHONPATH variable, where python2.X refers to the 
    python version (e.g. python2.7 if using python version 2.7.x).


================================
TEpeaks
================================

Usage
-----

    usage: TEpeaks -t treatment sample [treatment sample ...] 
                        -c control sample [control sample ...]
                        --tinput treatment input
                        --cinput control input
                        -s genome  
                        [optional arguments]

    Required arguments:
      -t | --treatment [treatment sample 1 treatment sample 2...]
         _Sample files in group 1 (e.g. treatment/mutant), separated by space_
         _Sample files in group 2 (e.g. control/wildtype), separated by space_
      --tinput    treatment input 
      -s genome  (hg: human19, mm: mouse9, dm: dm3)

    Optional arguments:
      -c | --control [control sample 1 control sample 2 ...]
      --cinput  control input
      --format [input file format]
         Input file format: BAM or BED. DEFAULT: BAM
      --project [name]      Name of this project. DEFAULT: TEpeak_out
      -p | --padj [pvalue]
         FDR cutoff for significance. DEFAULT: 1e-5
      -n | --norm [normalization]
         Normalization method : sd (library size),
                                bc (bin correlation). DEFAULT: sd
      -r | --step           step size. DEFAULT: 100
      -a | --auto           auto detect shiftsize. DEFAULT: False
      -d | --diff           require differential analysis
      -g | --gap            maximum allowed gap. DEFAULT: 1000
      -f | --fragsize       fragment size. DEFAULT: 200
      --lmfold              lower bound of fold change for modeling shipsize.
                            DEFAULT: 10
      --umfold              upper bound of fold change for modeling shiftsize.
                            DEFAULT: 30
      --minread             minimal reads of a peak. DEFAULT: 5
      --mode                TE counting mode. 'uniq' consider uniq-reads only. 'multi' distribute to all alignments. DEFAULT: multi
      --wig                 generate wiggle file for peaks (normalize to
                                10 million reads in total(library size))
      -h | --help           help info


Example Command Lines
---------------------

    TEpeaks --format BAM -t S1.bam --tinput S1input.bam -s mm -n sd --mode multi

    TEpeaks --format BAM -t S1.bam S2.bam -c C1.bam C2.bam  --tinput S1input.bam  --cinput C1input.bam -s mm -n sd --diff --mode multi



================================
TEtranscripts
================================

Usage
-----

    usage: TEtranscripts -t treatment sample [treatment sample ...] 
                         -c control sample [control sample ...]
                         --GTF genic-GTF-file
                         --TE TE-GTF-file 
                         [optional arguments]

    Required arguments:
      -t | --treatment [treatment sample 1 treatment sample 2...]
         Sample files in group 1 (e.g. treatment/mutant), separated by space
      -c | --control [control sample 1 control sample 2 ...]
         Sample files in group 2 (e.g. control/wildtype), separated by space
      --GTF genic-GTF-file  GTF file for gene annotations
      --TE TE-GTF-file      GTF file for transposable element annotations

    Optional arguments:

      *Input/Output options*
      --format [input file format]
         Input file format: BAM or SAM. DEFAULT: BAM
      --stranded [option]   Is this a stranded library? (yes, no, or reverse).
                            DEFAULT: yes.
      --sortByPos           Input file is sorted by chromosome position.
      --project [name]      Prefix used for output files (e.g. project name)
                            DEFAULT: TEtranscript_out

      *Analysis options*
      --mode [TE counting mode]
         How to count TE:
            uniq        (unique mappers only)
            multi       (distribute among all alignments).
         DEFAULT: uniq
      --minread [min_read] read count cutoff. DEFAULT: 1
      -L | --fragmentLength [fragLength]
         Average length of fragment used for paired-end sequencing
         DEFAULT: estimated from paired-end input alignment file
      -n | --norm [normalization]
         Normalization method : DESeq_default (default normalization method of DESeq), TC (total annotated read counts), quant (quantile normalization). 
         DEFAULT: DESeq_default
      -i | --iteration 
         maximum number of iterations used to optimize multi-reads assignment. DEFAULT: 0
      -p | --padj [pvalue]
         FDR cutoff for significance. DEFAULT: 0.05
      -f | --foldchange [foldchange]
         Fold-change ratio (absolute) cutoff for differential expression. 
         DEFAULT: 1

      *Other options*
      -h | --help
         Show help message
      --verbose [number]
         Set verbose level.
           0: only show critical messages
           1: show additional warning messages
           2: show process information
           3: show debug messages
         DEFAULT: 2
      --version
         Show program's version and exit

    *NOTE* BAM files must be either unsorted or sorted by queryname. If the BAM files are sorted by position, please use the '--sortByPos--' option

Example Command Lines
---------------------

    *** If BAM files are unsorted, or sorted by queryname: ***

    TEtranscripts --format BAM --mode multi -t RNAseq1.bam RNAseq2.bam -c CtlRNAseq1.bam CtlRNAseq.bam --project sample_nosort_test

    *** If BAM files are sorted by coordinates/position: ***

    TEtranscripts --sortByPos --format BAM --mode multi -t RNAseq1.bam RNAseq2.bam -c CtlRNAseq1.bam CtlRNAseq.bam --project sample_sorted_test

======================
Copying & distribution
======================

TEtranscripts and TEpeaks are part of TEToolKit.

TEToolKit is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but *WITHOUT ANY WARRANTY*; without even the implied warranty of
*MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE*.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with TEToolKit.  If not, see http://www.gnu.org/licenses/.


