[![MIT license](https://img.shields.io/badge/license-MIT-green.svg)](https://github.com/compbiocore/BioCookbook.jl/blob/master/LICENSE)

# BioCookbook.jl

A collection of Julia scripts for bioinformatics in use by CBC. We have:

 * `concatenate.jl`: Concatenates multiple sequence alignments (as FASTAs) into a supermatrix. Outputs a FASTA file. We have a Jupyter notebook showing that it is [about 6x faster](https://github.com/compbiocore/BioCookbook.jl/blob/master/benchmarks.ipynb) than the equivalent BioPython function.
 * `fastq_screen_parse.jl`: Go through the results of [FASTQscreen](https://www.bioinformatics.babraham.ac.uk/projects/fastq_screen/) and output either the intersection or the union of the results.
 * `extractReads.jl`: Extract reads from a region in a BAM file.
 * `clipanalysis.jl`: Analysis of soft and hard clipped reads in a BAM file.