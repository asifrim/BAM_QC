#!/bin/bash

INPUT=$1
OUTPUTDIR=$2
BED_FILE=/nfs/users/nfs_a/as33/Projects/BAM_QC/data/GENCODE17_CDS_V5_captured.bed
TMPDIR=/tmp/$(date +%s | sha256sum | base64 | head -c 32)
echo $TMPDIR
mkdir $TMPDIR

module add hgi/bedtools/latest

HISTFILE=$TMPDIR/$(basename $INPUT).coverage.hist.gz
bedtools coverage -abam $INPUT -b $BED_FILE -hist | gzip > $HISTFILE
zcat $TMPDIR/$(basename $INPUT).coverage.hist | head
SUMMARYFILE=$OUTPUTDIR/$(basename $INPUT).coverage 
python /nfs/users/nfs_a/as33/Projects/BAM_QC/src/process_coverage_hist.py $HISTFILE $SUMMARYFILE
# rm -rf $TMPDIR

