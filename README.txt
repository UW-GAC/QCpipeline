QC Pipeline
Step numbers refer to the SOP
Python scripts must be run on pearson, as they submit cluster jobs
Run "<script.py> --help" to see options and config parameters for each script
see QCpipeline-manual.pdf for R function help

1) Create project directory
create_project_dir.sh ProjectName user group
(where "user" is the first analyst working on the project,
and "group" is geneva/olga/cidr/etc.)


2-8) Create scan and snp annotation files
(save as ScanAnnotationDataFrame and SnpAnnotationDataFrame)

3e) Check that the Sample ID - file mapping is correct using the
netcdf config file:
> R -q --vanilla --args ncdf.config outfile < sample_id_from_files.R

9-10) Create NetCDF and GDS files
test 5 samples first:
> netcdf.py --email user@uw.edu --test ncdf.config
(where "user@uw.edu" is your email address)

run all:
> netcdf.py --email user@uw.edu ncdf.config
Check output to make sure creation was successful and all checks were passed
Pay attention to sections marked "MANUAL REVIEW"

Quality score will only be included in "xy" file if raw_qCol is not NA.

Use option "--checkPlink" to check a CIDR-generated PLINK file (with
A/B coding) against the newly created netCDF file.


13-15) Gender check (heterozygosity and mean intensity)
> gender_check.py --email user@uw.edu gender.config


17-19) Missing call rate
"round2" in config file should be FALSE
> missing.py --email user@uw.edu missing.config


20) Chromosome anomalies (need missing call rate first)
test first:
> chrom_anomalies.py --email user@uw.edu chrom_anom.config 1 10 5 \
--baf --loh --stats
(first 10 scans in 2 batches of 5 scans each)

> chrom_anomalies.py --email user@uw.edu chrom_anom.config start end by \
--baf --loh --stats
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


21) Batch quality checks (allele frequency test and plots)
> batch.py --email user@uw.edu batch.config
Default is chisq test (--type chisq).  
For Fisher test:
> batch.py --type fisher --email user@uw.edu batch.config


25-27) IBD (SNP selection, run IBD, plots, inbreeding coefficients)
> ibd.py --email user@uw.edu ibd.config
SNPs for IBD are selected with LD pruning.
The LD pruning file will be used if it already exists, otherwise it
will be created.
If you would like to use multithreading when available, add the
"--multithread" option when calling ibd.py. "--multithread 1-4" will
use the maximum number of cores available (between 1-4), and 
"--multithread 4" will request 4 cores. Do not request more than 8 cores.

27) Plate layout maps
"annot_scan_plateCol" and "annot_scan_wellCol" should refer to the plates
on which the samples were shipped to the genotyping center.
> plate_layout.py --email user@uw.edu plate_layout.config

28) Sample quality check
> sample_qualty.py --email user@uw.edu sample_quality.config


34) Recalculate missing call rates
"round2" in config file should be TRUE
> missing.py --email user@uw.edu missing.config


38) Create subject-level NetCDF or GDS genotype file with anomalies filtered
> geno_filt_subset.py --email user@uw.edu ncdf_subset.config


39) PCA
First round, unduplicated study samples + external hapmaps:
> pca.py --email user@uw.edu pca.config --combined
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


42a)
iii) HWE
> hwe.py --email user@uw.edu hwe.config
(repeat with different config files for mutiple ethnic groups)

To run chromosomes in parallel:
> hwe.py --email user@uw.edu hwe.config start end by
where start, end, by are chromosome numbers.  'by' is optional, if
omitted chromosomes will be run individually.  Values of 'end' > 23 are
ignored (in this case 'end' will be set to 23).


iv-vi) Allele frequency, duplicate discordance and Mendelian errors
> snp_filt.py --email user@uw.edu snp_filt.config
use option --MAconc for minor allele concordance 

Allele frequency, mendelian errors, and dupSNP will be run using the filtered subject-level netCDF.

use option --dupSNP for duplicate SNP discordance - snp annotation
needs a column with integer ids for duplicate SNPs (and NA for
singletons)

use option --mendClust to make cluster plots binned by number of
Mendelian errors (for family studies)

set corr.by.snp=TRUE in the config file to calculate correlation by
SNP (this is slow)


vii-viii) Allele frequency and heterozygosity by ethnic group and sex. Run using the filtered subject-level netcdf.
> snp_filt_ethn.py --email user@uw.edu snp_filt_ethn.config
(repeat with different config files for mutiple ethnic groups)


43d-e)
Association tests:

# To do association tests and plotting (including QQ, Manhattan and cluster plots) in one shot:
/projects/geneva/geneva_sata/GCC_code/QCpipeline/assoc.py \
/projects/geneva/geneva_sata/GCC_code/QCpipeline/config_examples/assoc.config \
start_chrom end_chrom --assoc --merge --plotQQManh --plotClust --queue gcc.q --email netID@uw.edu

# To run step by step:
/projects/geneva/geneva_sata/GCC_code/QCpipeline/assoc.py \
/projects/geneva/geneva_sata/GCC_code/QCpipeline/config_examples/assoc.config
start_chrom end_chrom --assoc (or --merge/--plotQQManh/--plotClust) --email netID@uw.edu

--merge generates a combined file for each model

--plotQQManh and --plotClust plot QQ, Manhattan, and cluster plots of
  SNPs on the chromosomes specified with plot_chroms in the  config file 

--queue gcc.q (or other queues) specifies the queue to use. The
  default is gcc.q. Check qstat first to decide which queue to use for
  --assoc
        
# If there are categorical covariates in the model(s), specify them in the config file so they will be converted into factors automatically.


dbGaP files:
3) make PLINK files
> plink.py  --email user@uw.edu --filtered plink.config
The unfiltered plink file will be made from the sample-level netCDF/GDS
using only subj.plink samples.
If "--filtered" argument is given, the script will ALSO create a
filtered plink file from the subject-level netCDF/GDS.
