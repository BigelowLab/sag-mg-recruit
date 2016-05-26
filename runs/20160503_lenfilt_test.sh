#!/usr/bin/env bash

# mp_recruit.py on 

mgs='LVP-2 LVP-4 LVP-5 LVP-8'

log=mp_recruit.log


parallel --retries 1 --bar --load 70% --joblog $log --jobs 2 python mp_recruit.py jrmr \
--outdir /mnt/stepanauskas_nfs/julia/fragdev/20160503/len_filt \
--fq1 /mnt/stepanauskas_nfs/julia/fragdev/20160431/{1}.extendedFrags.fastq.gz \
--refdir /mnt/stepanauskas_nfs/julia/Test_FragRecruitment/Masked_genomes/ --len_filter 150 ::: $mgs


