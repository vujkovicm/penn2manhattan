#!/bin/bash
#$ -cwd

# transform pennCNV output (to match Illumina's cnvPartition output)
qsub -l mem_free=64G -l h_vmem=64G -cwd -b y 'R CMD BATCH "--args cases.rawcnv" /mydir/penn2temp.R'
qsub -l mem_free=64G -l h_vmem=64G -cwd -b y 'R CMD BATCH "--args controls.rawcnv" /mydir/penn2temp.R'

# perform statistics and gwas-like Manhattan visualisation
qsub -l mem_free=64G -l h_vmem=64G -cwd -b y 'R CMD BATCH "--args 1" /mydir/temp2man.R'
qsub -l mem_free=64G -l h_vmem=64G -cwd -b y 'R CMD BATCH "--args 2" /mydir/temp2man.R'
qsub -l mem_free=64G -l h_vmem=64G -cwd -b y 'R CMD BATCH "--args 3" /mydir/temp2man.R'
qsub -l mem_free=64G -l h_vmem=64G -cwd -b y 'R CMD BATCH "--args 4" /mydir/temp2man.R'
qsub -l mem_free=64G -l h_vmem=64G -cwd -b y 'R CMD BATCH "--args 5" /mydir/temp2man.R'
qsub -l mem_free=64G -l h_vmem=64G -cwd -b y 'R CMD BATCH "--args 6" /mydir/temp2man.R'
qsub -l mem_free=64G -l h_vmem=64G -cwd -b y 'R CMD BATCH "--args 7" /mydir/temp2man.R'
qsub -l mem_free=64G -l h_vmem=64G -cwd -b y 'R CMD BATCH "--args 8" /mydir/temp2man.R'
qsub -l mem_free=64G -l h_vmem=64G -cwd -b y 'R CMD BATCH "--args 9" /mydir/temp2man.R'
qsub -l mem_free=64G -l h_vmem=64G -cwd -b y 'R CMD BATCH "--args 10" /mydir/temp2man.R'
qsub -l mem_free=64G -l h_vmem=64G -cwd -b y 'R CMD BATCH "--args 11" /mydir/temp2man.R'
qsub -l mem_free=64G -l h_vmem=64G -cwd -b y 'R CMD BATCH "--args 12" /mydir/temp2man.R'
qsub -l mem_free=64G -l h_vmem=64G -cwd -b y 'R CMD BATCH "--args 13" /mydir/temp2man.R'
qsub -l mem_free=64G -l h_vmem=64G -cwd -b y 'R CMD BATCH "--args 14" /mydir/temp2man.R'
qsub -l mem_free=64G -l h_vmem=64G -cwd -b y 'R CMD BATCH "--args 15" /mydir/temp2man.R'
qsub -l mem_free=64G -l h_vmem=64G -cwd -b y 'R CMD BATCH "--args 16" /mydir/temp2man.R'
qsub -l mem_free=64G -l h_vmem=64G -cwd -b y 'R CMD BATCH "--args 17" /mydir/temp2man.R'
qsub -l mem_free=64G -l h_vmem=64G -cwd -b y 'R CMD BATCH "--args 18" /mydir/temp2man.R'
qsub -l mem_free=64G -l h_vmem=64G -cwd -b y 'R CMD BATCH "--args 19" /mydir/temp2man.R'
qsub -l mem_free=64G -l h_vmem=64G -cwd -b y 'R CMD BATCH "--args 20" /mydir/temp2man.R'
qsub -l mem_free=64G -l h_vmem=64G -cwd -b y 'R CMD BATCH "--args 21" /mydir/temp2man.R'
qsub -l mem_free=64G -l h_vmem=64G -cwd -b y 'R CMD BATCH "--args 22" /mydir/temp2man.R'
