QC Pipeline

1) Create project directory
create_project_dir.sh ProjectName sdmorris

2) Create scan and snp annotation files

3) Create NetCDF files
test 5 samples first:
python netcdf.py --email sdmorris@uw.edu --test ncdf.config
run all:
python netcdf.py --email sdmorris@uw.edu ncdf.config
Check output to make sure creation was successful

4) Create GDS file
python gds.py --email sdmorris@uw.edu ncdf.config
Check output to make sure creation was successful!

5) Gender check (heterozygosity and mean intensity)
python gender_check.py --email sdmorris@uw.edu gender.config
