#! /usr/local/bin/python2.7

"""Create NetCDF and GDS files"""

import sys
import os
import subprocess
from optparse import OptionParser

usage = """%prog [options] config

Create genotype, XY intensity, and BAF/LRR netCDF files, 
and genotype GDS file.

Required config parameters:
annot_scan_file  scan annotation file
annot_snp_file   snp annotation file
array_build      genotyping build (e.g., hg19)
array_name       name of the genotyping array
raw_path         path to raw data
nc_bl_file       output BAF/LRR netCDF file
nc_xy_file       output XY intensity netCDF file (quality score optional)
nc_geno_file     output genotype netCDF file
gds_geno_file    output genotype GDS file

Required for "checkPlink" option:
plink_prefix     prefix of PLINK file to check

Optional config parameters [default]:
annot_scan_fileCol    [file]                 column of raw data file in scan annotation
annot_scan_nameCol    [Sample.Name]          column of raw data sample name in scan annotation
annot_snp_nameCol     [rsID]                 column of raw data snp name in snp annotation
raw_scanNameInFile    [1]                    1=repeated in column, -1=embedded in column heading, 0=none
raw_sepType           [,]                    column separator in raw data (e.g. ",", "\t")
raw_skipNum           [11]                   number of rows to skip (including header)
raw_colTotal          [19]                   total number of columns in raw data files
raw_snpCol            [1]                    column number with snp name
raw_sampleCol         [2]                    column number with scan name
raw_a1Col             [10]                   column number with allele 1
raw_a2Col             [11]                   column number with allele 2
raw_genoCol			  [NA]					 column number with diploid genotype of allele1 and alelle2
raw_qCol              [NA]                   column number with quality score (NA to omit)
raw_xCol              [14]                   column number with X intensity
raw_yCol              [15]                   column number with Y intensity
raw_bafCol            [18]                   column number with BAF
raw_lrrCol            [19]                   column number with LRR
nc_bl_checkFile       [nc_bl_check.RData]    output file for BAF/LRR netCDF check
nc_bl_diagFile        [nc_bl_diag.RData]     output file for BAF/LRR netCDF creation
nc_geno_checkFile     [nc_geno_check.RData]  output file for genotype netCDF check
nc_geno_diagFile      [nc_geno_diag.RData]   output file for genotype netCDF creation
nc_xy_checkFile       [nc_xy_check.RData]    output file for XY netCDF check
nc_xy_diagFile        [nc_xy_diag.RData]     output file for XY netCDF creation
out_plink_logfile     [plink_check.log]      output file for PLINK check"""
parser = OptionParser(usage=usage)
parser.add_option("-p", "--pipeline", dest="pipeline",
                  default="/projects/geneva/geneva_sata/GCC_code/QCpipeline",
                  help="pipeline source directory")
parser.add_option("-e", "--email", dest="email", default=None,
                  help="email address for job reporting")
parser.add_option("-q", "--queue", dest="qname", default="gcc.q", 
                  help="cluster queue name [default %default]")
parser.add_option("-t", "--test", dest="test",
                  action="store_true", default=False,
                  help="test with first 5 scans only")
parser.add_option("-o", "--overwrite", dest="overwrite",
                  action="store_true", default=False,
                  help="overwrite existing files")
parser.add_option("--checkPlink", dest="plink",
                  action="store_true", default=False,
                  help="check AB-coded PLINK file")
(options, args) = parser.parse_args()

if len(args) != 1:
    parser.error("incorrect number of arguments")

config = args[0]
pipeline = options.pipeline
email = options.email
test = options.test
overwrite = options.overwrite
qname = options.qname
plink = options.plink

sys.path.append(pipeline)
import QCpipeline

if test:
    testStr = "test"
else:
    testStr = ""

configdict = QCpipeline.readConfig(config)

if not overwrite:
    for file in (configdict['nc_geno_file'], configdict['nc_xy_file'],
                 configdict['nc_bl_file'], configdict['gds_geno_file']):
        if os.path.exists(file):
            sys.exit(file + "already exists; use -o flag to overwrite")

driver = os.path.join(pipeline, "runRscript.sh")

jobid = dict()
for job in ["ncdf_geno", "ncdf_xy", "ncdf_bl"]:
    rscript = os.path.join(pipeline, "R", job + ".R")
    jobid[job] = QCpipeline.submitJob(job, driver, [rscript, config, testStr], queue=qname, email=email)
    
if not test:
    job = "gds_geno"
    rscript = os.path.join(pipeline, "R", job + ".R")
    jobid[job] = QCpipeline.submitJob(job, driver, [rscript, config], holdid=[jobid['ncdf_geno']], queue=qname, email=email)
    
    if plink:
        holdid = [jobid['ncdf_geno']]

        # convert bed to ped
        prefix = configdict['plink_prefix']
        if not os.path.exists(prefix + ".ped"):
            if os.path.exists(prefix+".fam"):
                subprocess.call(["mv", prefix+".fam", prefix+".fam.orig"])
            subprocess.call(["cp", prefix+".fam.unr", prefix+".fam"])
            job = "plink_bed2ped"
            arglist = ["--noweb", "--bfile", prefix, "--recode", "--out", prefix]
            jobid[job] = QCpipeline.submitJob(job, "plink", arglist, options="-b y -j y -cwd", queue=qname, email=email)
            holdid.append(jobid['plink_bed2ped'])

        job = "plink_check"
        rscript = os.path.join(pipeline, "R", job + ".R")
        jobid[job] = QCpipeline.submitJob(job, driver, [rscript, config, "ABcoding"], holdid=holdid, queue=qname, email=email)
        
