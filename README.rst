TEtranscripts
=============

Version: 2.2.3

*NOTE* TEtranscripts and TEcount rely on specially curated GTF files, which are not
packaged with this software due to their size. Please go to
`our website <http://hammelllab.labsites.cshl.edu/software#TEtranscripts>`_
for instructions to download the curated GTF files.

TEtranscripts and TEcount takes RNA-seq (and similar data) and annotates reads to both
genes & transposable elements. TEtranscripts then performs differential analysis using
DESeq2.


`Github Page <https://github.com/mhammell-laboratory/TEtranscripts>`_

`Pypi Page <https://pypi.python.org/pypi/TEtranscripts>`_

`MHammell Lab <http://hammelllab.labsites.cshl.edu/software>`_

Created by Ying Jin, Eric Paniagua, Oliver Tam & Molly Hammell, February 2014

Copyright (C) 2014-2021 Ying Jin, Eric Paniagua, Talitha Forcier, Oliver Tam & Molly Hammell

Contact: Oliver Tam (tam@cshl.edu) or Talitha Forcier (talitha@cshl.edu)

Requirements
------------

Python:     2.7.x or 3.2.x (tested on Python 2.7.11 and 3.7.7)

pysam:      0.9.x or greater

R:          2.15.x or greater

DESeq2:     1.10.x or greater


Installation
------------

1. Download compressed tarball.
2. Unpack tarball.
3. Navigate into unpacked directory.
4. Run the following::

    $ python setup.py install

If you want to install locally (e.g. /local/home/usr),
run this command instead::

    $ python setup.py install --prefix /local/home/usr

*NOTE* In the above example, you must add

    /local/home/usr/bin

to the PATH variable, and

     /local/home/usr/lib/pythonX.Y/site-packages

to the PYTHONPATH variable, where X refers to the major
python version, and Y refers to the minor python version.
(e.g. `python2.7` if using python version 2.7.x, and
`python3.6` if using python version 3.6.x)


Alternative Singularity Installation for HPC
--------------------------------------------

Many High Performance Compunting clusters (HPCs) have
access to singularity which allows for the download and
exectuation of containers, TEtranscripts also has a
container through docker, it can be downloaded by
singularity thusly:

``singularity pull tetranscripts.sif docker://mhammelllab/tetranscripts:latest``

Execution is then through singularity as well:

``singularity exec tetranscripts.sif TEtranscripts -t <treatment sample> -c <control sample> --GTF <genic-GTF-file> --TE <TE-GTF-file>``

TEtranscripts
=============

Usage
-----

::

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
      --stranded [option]   Is this a stranded library? (no, forward, or reverse).
                 no      -  Library is unstranded
                 forward -  "Second-strand" cDNA library
                            (e.g. QIAseq stranded)
                 reverse -  "First-strand" cDNA library
                            (e.g. Illumina TruSeq stranded)
                            DEFAULT: no.
      --sortByPos           Input file is sorted by chromosome position.
      --project [name]      Prefix used for output files (e.g. project name)
                            DEFAULT: TEtranscript_out
      --outdir [directory]  Directory for output files.
                            DEFAULT: current directory

      *Analysis options*
      --mode [TE counting mode]
         How to count TE:
            uniq        (unique mappers only)
            multi       (distribute among all alignments).
         DEFAULT: multi
      --minread [min_read] read count cutoff. DEFAULT: 1
      -L | --fragmentLength [fragLength]
         Average length of fragment used for single-end sequencing
         DEFAULT: For paired-end, estimated from the input alignment file. For single-end, ignored by default.
      -i | --iteration
         maximum number of iterations used to optimize multi-reads assignment. DEFAULT: 100
      -p | --padj [pvalue]
         FDR cutoff for significance. DEFAULT: 0.05
      -f | --foldchange [foldchange]
         Fold-change ratio (absolute) cutoff for differential expression.
         DEFAULT: 1

      *DESeq1 compatibility options*
      --DESeq
         Use DESeq (instead of DESeq2) for differential analysis.
      -n | --norm [normalization]
         Normalization method : DESeq_default (default normalization method of DESeq), TC (total annotated read counts), quant (quantile normalization). Only applicable if DESeq is used instead of DESeq2.
         DEFAULT: DESeq_default

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

*NOTE* BAM files must be either unsorted or sorted by queryname. If the BAM files are sorted by position, please use the :code:`--sortByPos` option


Example Command Lines
---------------------

If BAM files are unsorted, or sorted by queryname::

    TEtranscripts --format BAM --mode multi -t RNAseq1.bam RNAseq2.bam -c CtlRNAseq1.bam CtlRNAseq.bam --project sample_nosort_test

If BAM files are sorted by coordinates/position::

    TEtranscripts --sortByPos --format BAM --mode multi -t RNAseq1.bam RNAseq2.bam -c CtlRNAseq1.bam CtlRNAseq.bam --project sample_sorted_test

Cluster Usage Recommendation
----------------------------

In our experience, we recommend around 20-30Gb of memory for analyzing human samples (hg19) with around 20-30 million mapped reads when running on a cluster.


TEcount
=======

Usage
-----

::

    usage: TEcount -b RNAseq BAM
                   --GTF genic-GTF-file
                   --TE TE-GTF-file
                   [optional arguments]

    Required arguments:
      -b | --BAM alignment-file  RNAseq alignment file (BAM preferred)
      --GTF genic-GTF-file       GTF file for gene annotations
      --TE TE-GTF-file           GTF file for transposable element annotations

    Optional arguments:

      *Input/Output options*
      --format [input file format]
         Input file format: BAM or SAM. DEFAULT: BAM
      --stranded [option]   Is this a stranded library? (no, forward, or reverse).
                 no      -  Library is unstranded
                 forward -  "Second-strand" cDNA library
                            (e.g. QIAseq stranded)
                 reverse -  "First-strand" cDNA library
                            (e.g. Illumina TruSeq stranded)
                            DEFAULT: no.
      --sortByPos           Input file is sorted by chromosome position.
      --project [name]      Prefix used for output files (e.g. project name)
                            DEFAULT: TEcount_out
      --outdir [directory]  Directory for output files.
                            DEFAULT: current directory

      *Analysis options*
      --mode [TE counting mode]
         How to count TE:
            uniq        (unique mappers only)
            multi       (distribute among all alignments).
         DEFAULT: multi
      -L | --fragmentLength [fragLength]
         Average length of fragment used for single-end sequencing
         DEFAULT: For paired-end, estimated from the input alignment file. For single-end, ignored by default.
      -i | --iteration
         maximum number of iterations used to optimize multi-reads assignment. DEFAULT: 100

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

*NOTE* BAM files must be either unsorted or sorted by queryname. If the BAM files are sorted by position, please use the :code:`--sortByPos` option


Example Command Lines
---------------------

If BAM files are unsorted, or sorted by queryname::

    TEcount --format BAM --mode multi -b RNAseq.bam --project sample_nosort_test

If BAM files are sorted by coordinates/position::

    TEtranscripts --sortByPos --format BAM --mode multi -b RNAseq.bam --project sample_sorted_test

Cluster Usage Recommendations
-----------------------------

TEcount is better suited than TEtranscripts for usage in the cluster environment, as each sample (e.g. replicates of an experiment) can be quantified on separate nodes. The output can then be merged into a single count table for differential analysis.
In our experience, we recommend around 20-30Gb of memory for analyzing human samples (hg19) with around 20-30 million mapped reads when running on a cluster.


Recommendations for TEtranscripts input files
=============================================

TEtranscripts can perform transposable element quantification from alignment results (e.g. BAM files) generated from a variety of programs.
Given the variety of experimental systems, we could not provide an optimal alignment strategy for every approach. Therefore,
we recommend that users identify the optimal parameters for their particular genome and alignment program in order to get the best
results.

When optimizing the alignment parameters, we recommend taking these points into consideration:

*Allowing sufficient number of multi-mappers during alignment*

Most alignment programs provide only 1 alignment per read by default. We recommend reporting multiple alignments per read. We have found
that reporting a maximum of 100 alignments per read provides an optimal compromise between the size of the alignment file and recovery
of multi-mappers in many genome builds. However, we highly suggest that users optimize this parameter for their particular experiment,
as this could significantly improve the quality of transposable element quantification.

*Optimizing alignment parameters for non-reference strains*

It is common that the specific laboratory strains used in an experiment contains genomic variations not present in the reference strain.
While this can be mitigated through allowing mismatches during alignments, certain lab strains (e.g. Drosophila melanogaster) have
diverged significantly from the reference genomes. We highly recommend that users should refine their alignment procedures to better
account for the expected variations between their lab strains and the reference genome, which will accordingly improve their analysis
with TEtranscripts. Users can also align to a custom genome build specific to their organism, though they would need GTF annotations for
genes and transposable elements that are compatible with their custom genome in order to utilize TEtranscripts. Please contact us if you
require advice in generating these annotation files.

*Paired end sequencing input*

For paired-end libraries, it is recommended that only alignments from properly paired reads are present in the input BAM file. I.e., each read 1 alignment should only have a single read 2 alignment. For example, if read 1 matched 3 genomic locations (A, B, C), then if read 2 also match 3 genomic locations (A', B', C'), then all three pairs of alignments could be used (and should be in the BAM file). However, if alignment C of read 1 was matched with more than one alignment of read 2 (e.g. C' and C*), then alignment C should be discarded (as there are unmatched alignments between read 1 and read 2). `STAR <https://github.com/alexdobin/STAR>`_ only outputs properly paired alignments by default, while `Bowtie2 <http://bowtie-bio.sourceforge.net/bowtie2/index.shtml>`_ requires the :code:`--no-mixed` parameter to be used.

*Specific recommendations when using STAR*

`STAR <https://github.com/alexdobin/STAR>`_ utilizes two parameters for optimal identification of multi-mappers `--outFilterMultimapNmax` and `--outAnchorMultimapNmax`.
The author of STAR recommends that `--winAnchorMultimapNmax` should be set at twice the value used in `--outFilterMultimapNmax`,
but no less than 50. In our study, we used the same number for both parameters (100), and found negligible differences in identifying
multi-mappers. Upon further discussion with the author of STAR, we recommend that setting the same value for `--winAnchorMultimapNmax`
and `--outFilterMultimapNmax`, though we highly suggest users test multiple values of `--winAnchorMultimapNmax` to identify the
optimal value for their experiment.


Copying & distribution
======================

TEtranscripts and TEcount are part of `TEToolkit suite <http://hammelllab.labsites.cshl.edu/software/>`_.

TEtranscripts is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but *WITHOUT ANY WARRANTY*; without even the implied warranty of
*MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE*.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with TEtranscripts.  If not, see `this website <http://www.gnu.org/licenses/>`_.

Citation
======================

If using the software in a publication, please cite the
`following <https://pubmed.ncbi.nlm.nih.gov/26206304/>`_:

Jin Y, Tam OH, Paniagua E, Hammell M. (2015) TEtranscripts: a package
for including transposable elements in differential expression
analysis of RNA-seq datasets. Bioinformatics. 31(22):3593-9.
