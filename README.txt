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
python netcdf.py --email user@uw.edu --test ncdf.config
(where "user@uw.edu" is your email address)

run all:
python netcdf.py --email user@uw.edu ncdf.config
Check output to make sure creation was successful and all checks were passed
Pay attention to sections marked "MANUAL REVIEW"

13-15) Gender check (heterozygosity and mean intensity)
python gender_check.py --email user@uw.edu gender.config

17-18) Missing call rate
"round2" in config file should be FALSE
python missing.py --email user@uw.edu missing.config

19) Chromosome anomalies (need missing call rate first)
test first:
python chrom_anomalies.py --email user@uw.edu chrom_anom.config 1 10 5
(first 10 scans in 2 batches of 5 scans each)

python chrom_anomalies.py --email user@uw.edu --skipSD chrom_anom.config start end by
(where "start end by" indicates how many scans to run at a time,
 e.g. "1 2000 500" means run scans 1-2000 in batches of 500)
("skipSD" means skip calculating the BAF standard deviation, since
this was already done in the test run)

20) Batch quality checks (allele frequency test and plots)
python batch.py --email user@uw.edu batch.config

24-26) IBD (allele frequency, SNP selection, run IBD, plots,
            inbreeding coefficients)
python ibd.py --email user@uw.edu ibd.config

27) Sample quality check
python sample_qualty.py --email user@uw.edu sample_quality.config

32) Recalculate missing call rates
"round2" in config file should be TRUE
python missing.py --email user@uw.edu missing.config

36) Create subject-level NetCDF genotype file
python netcdf_subset.py --email user@uw.edu ncdf_subset.config

38) PCA
First, run 2 rounds of PCA: 1) with external HapMaps and
2) unrelated study samples
python pca.py --email user@uw.edu pca.config --combined
(where "--combined" option means run PCA with external HapMaps in
additional to unrelated study samples)

For subsequent runs with individual ethnic groups, make a new
"study_unrelated.RData" file and note in the configuration
python pca.py --email user@uw.edu pca.config

42a)
iii) HWE
python hwe.py --email user@uw.edu hwe.config
(repeat with different config files for mutiple ethnic groups)

iv-vi) Allele frequency, duplicate discordance and Mendelian errors
python snp_filt.py --email user@uw.edu snp_filt.config

vii-viii) Allele frequency and heterozygosity by ethnic group and sex
python snp_filt_ethn.py --email user@uw.edu snp_filt_ethn.config
(repeat with different config files for mutiple ethnic groups)
