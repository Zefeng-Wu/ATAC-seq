#!/usr/bin/env Rscript

# 1. Main Functions and Steps:

# Functions： Read alignment files director contaning chromain modification alignment files (*.bam file), genomic annotation (.gtf) file and gene expression (table) file.Then it converts them into normalized read counts matrix where each
# row stands for a gene while each column is the chromatin moidification intensity. Then, machine learning methods are used to model gene expression based on the levels of those chromatin modifications.

# example
# USAGE:  Rscript PGEP.r -b bam_files_director -o output_file_name

# parameters
# Options
# -i/--input  bam file directers
# -o/--output output name
# -f/--fragment 

options(warn = -1)

# 2. library dependence insepct, if not, install.
package_list <- c("optparse","ggplot2","ATACseqQC","GenomicFeatures","rtracklayer")

for(p in package_list){
  if(!suppressWarnings(suppressMessages(require(p, character.only = TRUE, quietly = TRUE, warn.conflicts = FALSE)))){
    install.packages(p, repos="http://cran.r-project.org")
    suppressWarnings(suppressMessages(library(p, character.only = TRUE, quietly = TRUE, warn.conflicts = FALSE)))
  }
}
# another install source 
if (FALSE){
  # Bioconductor install
  source("https://bioconductor.org/biocLite.R")
  biocLite(c("bamsignals","GenomicFeatures"))
}

# Clean enviroment object
rm(list=ls()) 

# Load essential packages
library(optparse)
library(ATACseqQC)
library(GenomicFeatures)
library(rtracklayer)

# Parse commond lines
if (TRUE){
  option_list <- list(
    make_option(c("-b", "--bams"), type = "character", default = NULL,
                help = "Input bam file", metavar = "character"),
    make_option(c("-o", "--output"), type = "character", default = NULL,
                help = "output directory or prefix"),
    make_option(c("-f", "--fragmentsize"), type = "character", default = NULL,
                help = "output fragment size distribution directory or prefix")
    
    
  )
  opt_parser = OptionParser(option_list=option_list)
  opts = parse_args(opt_parser)
  # Show and verify your input information 
  print(paste("The bam files is", opts$bams,  sep = ""))
  print(paste("The output file prefix is ", opts$output, sep = ""))
  print(paste("The output fragment size pdf is ", opts$fragmentsize, sep = ""))
}

## verify the parameters existence 
if (is.null(opts$b)){
  print_help(opt_parser)
  stop("Please input the bam file director ", call.=FALSE)
}

if (is.null(opts$o)){
  print_help(opt_parser)
  stop("Please input out file prefix ", call.=FALSE)
}

if (is.null(opts$f)){
  print_help(opt_parser)
  stop("Please input fragment size out file prefix ", call.=FALSE)
}

 
#######
tags <- c("AS", "XN", "XM", "XO", "XG", "NM", "MD", "YS", "YT")
## shift the coordinates of 5'ends of alignments in the bam file
gal <- readBamFile(opts$b, tag=tags,asMates=TRUE)
#shift the GAlignmentsLists by 5' ends. All reads aligning to the positive strand will be offset by +4bp, and all reads aligning to the negative strand will be offset -5bp by default.
gal1 <- shiftGAlignmentsList(gal)
export(gal1, opts$o)
#### fragment  size output
pdf(opts$f/basename(opts$b), width =10, height=8) 
fragSizeDist(bamFiles=opts$b, bamFIles.labels=basename(opts$b))
message("ATACseqQC annalysis is finished!")


