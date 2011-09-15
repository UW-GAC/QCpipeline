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
