QC Pipeline
Step numbers refer to the SOP

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
Check output to make sure creation was successful

13-15) Gender check (heterozygosity and mean intensity)
python gender_check.py --email user@uw.edu gender.config

19-20) Missing call rate
python missing.py --email user@uw.edu missing.config

18) Chromosome anomalies (need missing call rate first)
test first:
python chrom_anomalies.py --email user@uw.edu chrom_anom.config 1 10 5
(first 10 scans in 2 batches of 5 scans each)

python chrom_anomalies.py --email user@uw.edu --skipSD chrom_anom.config start end by
(where "start end by" indicates how many scans to run at a time,
 e.g. "1 2000 500" means run scans 1-2000 in batches of 500)
("skipSD" means skip calculating the BAF standard deviation, since
this was already done in the test run)
