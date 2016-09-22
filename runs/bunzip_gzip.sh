#!usr/bin/env sh

for f in *.bz2
do 
tar -vxjf $f
file="${f%\.tar*}"
gzip $file/*.fastq
done