QC Pipeline
Step numbers refer to the SOP
Python scripts must be run on pearson, as they submit cluster jobs
see QCpipeline-manual.pdf for R function help

1) Create project directory
create_project_dir.sh ProjectName user
(where "user" is the first analyst working on the project)


2-8) Create scan and snp annotation files
(save as ScanAnnotationDataFrame and SnpAnnotationDataFrame)


9-10) Create NetCDF and GDS files
test 5 samples first:
> python netcdf.py --email user@uw.edu --test ncdf.config
(where "user@uw.edu" is your email address)

run all:
> python netcdf.py --email user@uw.edu ncdf.config
Check output to make sure creation was successful and all checks were passed
Pay attention to sections marked "MANUAL REVIEW"


13-15) Gender check (heterozygosity and mean intensity)
> python gender_check.py --email user@uw.edu gender.config


17-19) Missing call rate
"round2" in config file should be FALSE
> python missing.py --email user@uw.edu missing.config


20) Chromosome anomalies (need missing call rate first)
test first:
> python chrom_anomalies.py --email user@uw.edu chrom_anom.config 1 10 5
(first 10 scans in 2 batches of 5 scans each)

> python chrom_anomalies.py --email user@uw.edu chrom_anom.config start end by
(where "start end by" indicates how many scans to run at a time,
 e.g. "1 2000 500" means run scans 1-2000 in batches of 500)

Use the --maf option to exclude SNPs below a given MAF (e.g. 0.05):
> python chrom_anomalies.py --email user@uw.edu --maf 0.05 chrom_anom.config 1 10 5
This will run the allele frequency calculation first.

Note that if allele frequency and/or BAF variance files specified in
the config file already exist, they will not be recalculated.


21) Batch quality checks (allele frequency test and plots)
> python batch.py --email user@uw.edu batch.config
Default is chisq test (--type chisq).  
For Fisher test:
> python batch.py --type fisher --email user@uw.edu batch.config


25-27) IBD (allele frequency, SNP selection, run IBD, plots,
            inbreeding coefficients)
> python ibd.py --email user@uw.edu ibd.config

The allele frequency file will be used if it already exists, otherwise
it will be created.


28) Sample quality check
> python sample_qualty.py --email user@uw.edu sample_quality.config


33) Recalculate missing call rates
"round2" in config file should be TRUE
> python missing.py --email user@uw.edu missing.config


37) Create subject-level NetCDF and GDS genotype files with anomalies filtered
> python netcdf_subset.py --email user@uw.edu ncdf_subset.config


38) PCA
First round, unduplicated study samples + external hapmaps:
> python pca.py --email user@uw.edu pca.config --combined
Second round, unrelated study samples:
> python pca.py --email user@uw.edu pca.config

Selecting which samples should be included is left to the user: make
vectors of scanID and save as RData files, then give the path to
these files in the config file.

For subsequent runs with individual ethnic groups, make a new
"study_unrelated.RData" file and update this parameter in the configuration
> python pca.py --email user@uw.edu pca.config

The LD pruning file will be used if it already exists, otherwise it
will be created.


41a)
iii) HWE
> python hwe.py --email user@uw.edu hwe.config
(repeat with different config files for mutiple ethnic groups)

iv-vi) Allele frequency, duplicate discordance and Mendelian errors
> python snp_filt.py --email user@uw.edu snp_filt.config

vii-viii) Allele frequency and heterozygosity by ethnic group and sex
> python snp_filt_ethn.py --email user@uw.edu snp_filt_ethn.config
(repeat with different config files for mutiple ethnic groups)


42d-e)
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
