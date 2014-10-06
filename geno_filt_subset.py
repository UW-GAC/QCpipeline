#! /usr/local/bin/python2.7

"""Filter and subset a NetCDF or GDS file"""

import sys
import os
import subprocess
from optparse import OptionParser

usage = """%prog [options] config

Create subset netCDF or GDS file with chromosome anomalies filtered
and scans excluded.

Required config parameters:
annot_scan_file  scan annotation file
annot_snp_file   snp annotation file
in_file          input genotype netCDF/GDS file
out_file         output genotype netCDF/GDS file

Optional config parameters [default]:
chrom_anom_file       [NA]    data frame of chromosome anomalies, with columns scanID, chromosome, left.base, right.base, whole.chrom, filter
filterYinF            [TRUE]  filter Y chromosome for females?
scan_include_file     [NA]    vector of scanID to include (NA=all)"""
parser = OptionParser(usage=usage)
parser.add_option("-p", "--pipeline", dest="pipeline",
                  default="/projects/geneva/gcc-fs2/GCC_Code/QCpipeline",
                  help="pipeline source directory")
parser.add_option("-e", "--email", dest="email", default=None,
                  help="email address for job reporting")
parser.add_option("-q", "--queue", dest="qname", default="gcc.q", 
                  help="cluster queue name [default %default]")
(options, args) = parser.parse_args()

if len(args) != 1:
    parser.error("incorrect number of arguments")

config = args[0]
pipeline = options.pipeline
email = options.email
qname = options.qname

sys.path.append(pipeline)
import QCpipeline

driver = os.path.join(pipeline, "runRscript.sh")

jobid = dict()
job = "subset"
rscript = os.path.join(pipeline, "R", job + ".R")
jobid[job] = QCpipeline.submitJob(job, driver, [rscript, config], queue=qname, email=email)
