# mg-sag-recruitment
Read recruitment pipeline 2.0 for comparison of sag abundance across metagenomes

To install, clone repo from github, move into directory and type:

python setup.py install

After installation, for instructions on how to run type:

```
sag-mg-recruit -h
```

which will return:

```
Usage: sag-mg-recruit [OPTIONS] INPUT_MG_TABLE INPUT_SAG_LIST

Options:
  --outdir TEXT     directory location to place output files
  --cores INTEGER   number of cores to run on
  --mmd FLOAT       for join step: mismatch density for join step in mg
                    processing
  --mino INTEGER    for join step: minimum overlap for join step
  --maxo INTEGER    for join step: maximum overlap for join step
  --minlen INTEGER  for metagenomes: minimum read length
  --pctid INTEGER   for alignment: minimum percent identity to keep
  --log TEXT        name of log file for this process
  -h, --help        Show this message and exit.
```

```input_mg_table``` should be a csv file created using the mg_template.xlsx within the program directory.  The columns are:
- mg_f: path to forward reads for metagenome
- mg_r: path to reverse metagenomic reads (if there are none, put "None" in this column)
- wgs_technology: specify whether the library was sequenced with either illumina (enter "illumina") or 454-pyrosequencing (enter "pyro")
- join: designate whether you want the metagenome to be joined or not; either True or False

```intput_sag_list``` is a text file containing a list of file paths to SAGs in fasta format to be used by this process, entered by line.

The program will create a new output directory with three sub-directories:
```coverage```, ```mgs``` and ```sags``` which will contain various output files created during the recruitment process.  

The final output table will be located in the main output directory, called "summary_table.txt".

That table has one row per mg-sag pair with the following columns:
- sag: SAG name, which is determined as all text within the SAG file name before the first "."
- metagenome: metagenome name, determined in the same way as SAG name
- Percent_scaffolds_with_any_coverage: percent of SAG scaffolds with any coverage
- Percent_of_reference_bases_covered: percent of SAG bases covered by at least one read
- Average_coverage: mean coverage across the SAG
- total_reads_recruited: total number of metagenomic reads recruited to the SAG
- mg_wgs_technology: metagenome sequencing technology used, designated in the input_mg_table
- mg_read_count: total number of reads used in the processed metagenome
- sag_completeness: % sag completeness as determined by checkm
- sag_total_bp: number of bp in input SAG fasta file
- sag_calculated_length: length of SAG genome given the sag completeness and sag_total_bp
- prop_mg_recruited: proportion of metagenomic reads recruited to SAG genome, corrected based on SAG completeness
- prop_mg_adjusted: proportion of metagenomic reads adjusted based on sequencing technology.


