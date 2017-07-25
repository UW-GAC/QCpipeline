QC Pipeline
Step numbers refer to the SOP (version "transition_revised")
Python scripts must be run on cluster head node, as they submit cluster jobs
Run "<script.py> --help" to see options and config parameters for each script
for R function help, use "help(function_name)" within R

2) Create project directory
create_project_dir.sh ProjectName user group
(where "user" is the first analyst working on the project,
and "group" is geneva/olga/cidr/etc.)


4-5) Create scan and snp annotation files
(save as ScanAnnotationDataFrame and SnpAnnotationDataFrame)

5b) Check that the Sample ID - file mapping is correct using the
config file:
> R -q --vanilla --args create_datafiles.config outfile < sample_id_from_files.R


6) Create GDS files
test 5 samples first:
> create_datafiles.py --email user@uw.edu --test create_datafiles.config
(where "user@uw.edu" is your email address)

run all:
> create_datafiles.py --email user@uw.edu create_datafiles.config
Check output to make sure creation was successful and all checks were passed
Pay attention to sections marked "MANUAL REVIEW"

Quality score will only be included in "xy" file if raw_qCol is not NA.

Use option "--checkPlink" to check a CIDR-generated PLINK file against
the newly created GDS file.

To create files in batches and combine:
> create_datafiles.py --email user@uw.edu --batches 5 create_datafiles.config
(creates 5 separate GDS files with subsets of the data)
> create_datafiles.py --email user@uw.edu --batches 5 --combine_batches create_datafiles.config
(combines batches into single file)
> rm batch*.gds
(delete batches after verifying that combine was successful)


9) Gender check (heterozygosity and mean intensity)
> gender_check.py --email user@uw.edu gender.config


10) Missing call rate
"round2" in config file should be FALSE
> missing.py --email user@uw.edu missing.config


NA) Chromosome anomalies (need missing call rate first)
test first:
> chrom_anomalies.py --email user@uw.edu chrom_anom.config 1 10 5 --baf --loh --stats
(first 10 scans in 2 batches of 5 scans each)

> chrom_anomalies.py --email user@uw.edu chrom_anom.config start end by --baf --loh --stats
(where "start end by" indicates how many scans to run at a time,
 e.g. "1 2000 500" means run scans 1-2000 in batches of 500)

Use the --maf option to exclude SNPs below a given MAF (e.g. 0.05):
> chrom_anomalies.py --email user@uw.edu --maf 0.05 chrom_anom.config 1 10 5
This will run the allele frequency calculation first.

Note that if allele frequency and/or BAF variance files specified in
the config file already exist, they will not be recalculated.

Use the flags "--baf", "--loh", and "--stats" to run components of the
pipeline separately.  Note that running loh requires output from baf,
and running stats requires output from both baf and loh.  (Thus it
does not make sense to use options --baf --stats together without
--loh.)


11) Batch quality checks
> batch.py --email user@uw.edu batch.config


14) IBD (SNP selection, run IBD, plots, inbreeding coefficients)
> ibd.py --email user@uw.edu ibd.config
SNPs for IBD are selected with LD pruning.
The LD pruning file will be used if it already exists, otherwise it
will be created.
To calculate individual inbreeding coefficients, use option
"--inbreed",  (Not valid for KING.)
If you would like to use multithreading when available, add the
"--multithread" option when calling ibd.py. "--multithread 1-4" will
use the maximum number of cores available (between 1-4), and 
"--multithread 4" will request 4 cores. Do not request more than 8 cores.

For family studies, provide a family variable in the config file.


NA) Sample quality check
> sample_qualty.py --email user@uw.edu sample_quality.config


16a) Plate layout maps
"annot_scan_plateCol" and "annot_scan_wellCol" should refer to the plates
on which the samples were shipped to the genotyping center.
> plate_layout.py --email user@uw.edu plate_layout.config


19) Recalculate missing call rates
"round2" in config file should be TRUE
> missing.py --email user@uw.edu missing.config


24) Create subject-level GDS genotype file with anomalies filtered
> geno_filt_subset.py --email user@uw.edu subset.config


25) PCA
First round, unduplicated study samples + external hapmaps:
> pca.py --email user@uw.edu pca.config --combined

Make sure both datasets are using the same alleles (e.g.. TOP, PLUS).
Set columns to use for alleles A and B in the config file (defaults to 
"alleleA" and "alleleB" for both datasets).

Second round, unrelated study samples:
> pca.py --email user@uw.edu pca.config

Selecting which samples should be included is left to the user: make
vectors of scanID and save as RData files, then give the path to
these files in the config file.

For subsequent runs with individual ethnic groups, make a new
"study_unrelated.RData" file and update this parameter in the configuration
> pca.py --email user@uw.edu pca.config

The LD pruning file will be used if it already exists, otherwise it
will be created.
If you would like to use multithreading when available, add the
"--multithread" option when calling pca.py. "--multithread 1-4" will
use the maximum number of cores available (between 1-4), and
"--multithread 4" will request 4 cores. Do not request more than 8 cores.


30) HWE
> hwe.py --email user@uw.edu hwe.config
(repeat with different config files for mutiple ethnic groups)
Chomosomes 1-23 will be run in parallel. To specify a chromosome range:
> hwe.py --email user@uw.edu hwe.config start_chrom end_chrom


31) Allele frequency, duplicate discordance and Mendelian errors
> snp_filt.py --email user@uw.edu snp_filt.config
use option --MAconc for minor allele concordance 

Allele frequency, mendelian errors, and dupSNP will be run using the filtered subject-level GDS.

use option --dupSNP for duplicate SNP discordance - snp annotation
needs a column with integer ids for duplicate SNPs (and NA for
singletons)

use option --mendClust to make cluster plots binned by number of
Mendelian errors (for family studies)

set corr.by.snp=TRUE in the config file to calculate correlation by
SNP (this is slow)


32) Allele frequency and heterozygosity by ethnic group and sex. Run using the filtered subject-level GDS.
> snp_filt_ethn.py --email user@uw.edu snp_filt_ethn.config
(repeat with different config files for mutiple ethnic groups)


34)
Association tests:

To do association tests and plotting (including QQ, Manhattan and cluster plots) in one shot:
> assoc.py assoc.config start_chrom end_chrom --assoc --merge --plotQQManh --plotClust --queue gcc.q --email netID@uw.edu

To run step by step:
> assoc.py assoc.config start_chrom end_chrom --assoc (or --merge/--plotQQManh/--plotClust) --email netID@uw.edu

--merge generates a combined file for each model

--plotQQManh and --plotClust plot QQ, Manhattan, and cluster plots of
  SNPs on the chromosomes specified with plot_chroms in the  config file 

--queue gcc.q (or other queues) specifies the queue to use. The
  default is all.q. Check qstat first to decide which queue to use for
  --assoc
        
If there are categorical covariates in the model(s), specify them in the config file so they will be converted into factors automatically.


dbGaP files:
make PLINK files
> plink.py  --email user@uw.edu --filtered plink.config
The unfiltered plink file will be made from the sample-level GDS
using only subj.plink samples.
If "--filtered" argument is given, the script will ALSO create a
filtered plink file from the subject-level GDS.
