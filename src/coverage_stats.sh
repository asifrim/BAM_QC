#!/bin/bash

INPUT=$1
OUTPUTDIR=$2
BED_FILE=/nfs/users/nfs_a/as33/Projects/BAM_QC/data/GENCODE17_CDS_V5_captured.bed.gz
DDG2P_FILE=/nfs/users/nfs_a/as33/Projects/BAM_QC/data/hgnc_genes.bed_20_03_2015.txt.gz

TMPDIR=/tmp/$(date +%s | sha256sum | base64 | head -c 32)
echo $TMPDIR
mkdir $TMPDIR

module add hgi/bedtools/latest

HISTFILE=$TMPDIR/$(basename $INPUT).coverage.hist.gz
bedtools coverage -abam $INPUT -b $BED_FILE -hist | gzip > $HISTFILE
PROBESUMMARYFILE=$OUTPUTDIR/$(basename $INPUT).probe_coverage
GENESUMMARYFILE=$OUTPUTDIR/$(basename $INPUT).gene_coverage
ALLSUMMARYFILE=$OUTPUTDIR/$(basename $INPUT).all_coverage   
python /nfs/users/nfs_a/as33/Projects/BAM_QC/src/process_coverage_hist.py $HISTFILE $PROBESUMMARYFILE $GENESUMMARYFILE $ALLSUMMARYFILE $DDG2P_FILE
gzip $PROBESUMMARYFILE
gzip $GENESUMMARYFILE
rm -rf $TMPDIR

