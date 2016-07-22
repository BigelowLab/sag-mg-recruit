# sag-mg-recruit
Read recruitment pipeline for comparison of sag abundance across metagenomes

To install, clone repo from github, move into directory and type:

python setup.py install

## Dependencies:
> python 2.7

> Flash

> bwa

> samtools

> bedtools


After installation, for instructions on how to run type:

```
sag-mg-recruit -h
```

which will return:

```
Usage: sag-mg-recruit [OPTIONS] INPUT_MG_TABLE INPUT_SAG_TABLE

Options:
  --outdir TEXT      directory location to place output files
  --cores INTEGER    number of cores to run on, default=8
  --mmd FLOAT        for join step: mismatch density for join step in mg
                     processing, default=0.05
  --mino INTEGER     for join step: minimum overlap for join step,
                     default=35bp
  --maxo INTEGER     for join step: maximum overlap for join step,
                     default=150bp
  --minlen INTEGER   for metagenomes: minimum read length, default=150bp
  --pctid INTEGER    for alignment: minimum percent identity to keep within
                     overlapping region, default=95
  --overlap INTEGER  for alignment: percent read that must overlap with
                     reference sequence to keep, default=95
  --log TEXT         name of log file, else, log sent to standard out
  -h, --help         Show this message and exit.
```

### Input details:

```input_mg_table``` should be a csv file created using the mg_template.xlsx within the program directory.  The columns are:
- name: desired name used to refer to each metagenome; any string without spaces, special characters or "."
- mg_f: path to forward reads for metagenome
- mg_r: path to reverse metagenomic reads (if there are none, put "None" in this column)
- wgs_technology: specify whether the library was sequenced with either illumina (enter "illumina") or 454-pyrosequencing (enter "pyro")
- join: designate whether you want the metagenome to be joined or not; either True or False

See mg_template.xlsx for example table.  Input your own information and save as a .csv file.

```intput_sag_table``` is a csv formatted table with the following columns:
- sag_name: any string with no spaces or '.'
- fasta_file: None if sag should be masked, otherwise path to input fasta to be processed
- gbk_file: path to SAG's annotated gbk file if mask = True
- mask: boolean indicating whether the SAG should have the 16/23S sequences masked (TRUE) or not (FALSE)

See sag_template.xlsx for example table.  Input your own information in excel and save as a .csv file.

### Example input script:
for help:

```sag-mg-recruit -h```

to run with 95% identity alignment, 40 cores, minimum read length of 100:

```
sag-mg-recruit --outdir <path to output dir> --cores 40 --minlen 100  --pctid 95 --log recruitment.log <input mg table> <input sag table>
```



### Output:

The program will create a new output directory with three sub-directories:
```coverage```, ```mgs``` and ```sags``` which will contain various output files created during the recruitment process.  

The final output table will be located in the main output directory, called "summary_table_pctidXX_minlenXX_overlapXX.txt".

That table has one row per mg-sag pair with the following columns:
- sag: SAG name
- metagenome: metagenome name
- Percent_scaffolds_with_any_coverage: percent of SAG scaffolds with any coverage
- Percent_of_reference_bases_covered: percent of SAG bases covered by at least one read
- Average_coverage: mean coverage across the SAG
- total_reads_recruited: total number of metagenomic reads recruited to the SAG
- mg_wgs_technology: metagenome sequencing technology used, designated in the input_mg_table
- mg_read_count: total number of reads used in the processed metagenome
- sag_completeness: % sag completeness as determined by checkm
- sag_total_bp: number of bp in input SAG fasta file
- sag_size_mbp: SAG size in megabasepairs (mbp)
- reads_per_mbp: number of reads recruited per SAG mbp
- prop_mgreads_per_mbp: proprtion metagenomic reads recruited per SAG mbp


#### Changing alignment parameters

After you've run sag-mg-recruit, you can re-run the analysis with changed ```--pctid``` or ```--minlen``` parameters much more quickly than the original run by using the same input files, the same output file but only changing the ```--pctid``` and/or ```--minlen``` parameters. 

For example, say your first run was:

```
sag-mg-recruit --outdir /home/smroutput --cores 40 --minlen 150  --pctid 95 --overlap 95 --log recruitment.log mgtbl.txt sagtbl.txt
```

Re-run with different minlen and pctid parameters:
```
sag-mg-recruit --outdir /home/smroutput --cores 40 --minlen 100  --pctid 90 --overlap 100 --log recruitment.log mgtbl.txt sagtbl.txt
```

**Note** a quirk of one of the dependencies (FLASH) which is used for the join step: if there are metagenomic reads that need to be joined, make sure the output directory is specificied as a relative rather than absolute path.
