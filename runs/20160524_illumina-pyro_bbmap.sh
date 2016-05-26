#!/usr/bin/env bash

module load parallel
module load bbmap
module load samtools

mg1=./Test_FragRecruitment/Sakinaw/Hallam_metagenomes/2079.6.1746.fastq.gz 
mg2=./Test_FragRecruitment/Sakinaw/Hallam_metagenomes/GZFH_454.fastq.gz
gnms=./Test_FragRecruitment/Sakinaw/Hallam_SAGs_AAA255/Masked_genomes_AAA255/testfirst/
outdir=./fragdev/20160525/

python mp_recruit.py jrmr --fq1 $mg1 --outdir $outdir --join --len_filter 150 --refdir $gnms

python mp_recruit.py jrmr --fq1 $mg2 --outdir $outdir --len_filter 150 --refdir $gnms


mg1joined=./fragdev/20160525/2079.6.1746.fastq.extendedFrags_gt150.fastq.gz


python ~/scripts/rcov2.py print_cov --fastq $mg2lenfilt --reference {1} --outdir $outdir  

#python rcov2.py print_cov \
#--fastq /mnt/stepanauskas_nfs/julia/Test_FragRecruitment/Sakinaw/Hallam_metagenomes/2079.6.1746.fastq.gz \
#--contigs /mnt/stepanauskas_nfs/julia/Test_FragRecruitment/Sakinaw/Hallam_SAGs_AAA255/Masked_genomes_AAA255/testfirst/AAA255A6_96912.fasta \
#--outdir /mnt/stepanauskas_nfs/julia/fragdev/20160524

# this run was still running at the end of the day in a screen session
python ~/scripts/mp_recruit.py jrmr --fq1 ../Test_FragRecruitment/Sakinaw/Hallam_metagenomes/2079.6.1746.fastq.gz \
--outdir 20160525 --refdir ../Test_FragRecruitment/Sakinaw/Hallam_SAGs_AAA255/Masked_genomes_AAA255/testfirst/ \
--join --len_filter 150 &> 0525_jrmr2.log

# this run finished:
python ~/scripts/mp_recruit.py jrmr --fq1 ../Test_FragRecruitment/Sakinaw/Hallam_metagenomes/GZFH_454.fastq.gz \
--outdir 20160525 --refdir ../Test_FragRecruitment/Sakinaw/Hallam_SAGs_AAA255/Masked_genomes_AAA255/testfirst/ \
--len_filter 150 &> 0525_jrmr3.log

log=0525_rcov1.log

outdir=./fragdev/20160525/
mg2lenfilt=/mnt/stepanauskas_nfs/julia/Test_FragRecruitment/Sakinaw/Hallam_metagenomes/GZFH_454_gt150.fastq.gz

## running this at the end of the day in a screen session
parallel --retries 1 --bar --load 70% --joblog $log --jobs 2 python ~/scripts/rcov2.py \
print_cov --fastq $mg2lenfilt --reference {1} --outdir $outdir :::: firstlist.txt