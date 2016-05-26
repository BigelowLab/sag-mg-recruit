#!/usr/bin/env bash

# mp_recruit.py on 

mgs='LVP-2 LVP-4 LVP-5 LVP-8'

log=mp_recruit.log


parallel --retries 1 --bar --load 70% --joblog $log --jobs 2 python mp_recruit.py jrmr \
--fq2 /mnt/stepanauskas_nfs/mpachiadaki/From_Ramunas/Vents_Sievert_150902/Metagenomes/QC_reads/{1}_trimmed_R2.fastq.gz \
--outdir ./test2 {1} /mnt/stepanauskas_nfs/mpachiadaki/From_Ramunas/Vents_Sievert_150902/Metagenomes/QC_reads/{1}_trimmed_R1.fastq.gz \
/mnt/stepanauskas_nfs/julia/Test_FragRecruitment/Masked_genomes/ ::: $mgs


