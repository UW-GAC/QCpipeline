#! /usr/local/bin/python2.7

"""Discordance between two data sets"""

import sys
import os
import subprocess
from optparse import OptionParser

usage = """%prog [options] config

Calculate discordance between two data sets.

Required config parameters:
annot_scan_file_1  scan annotation file for dataset 1
annot_scan_file_2  scan annotation file for dataset 2
annot_snp_file_1   snp annotation file for dataset 1
annot_snp_file_2   snp annotation file for dataset 2
geno_file_1        genotype file (netCDF or GDS) for dataset 1
geno_file_2        genotype file (netCDF or GDS) for dataset 2
out_prefix         prefix for output files

Optional config parameters [default]:
annot_scan_subjCol_1  [subjectID]  column of matching subject ID in scan annotation 1
annot_scan_subjCol_2  [subjectID]  column of matching subject ID in scan annotation 2
annot_snp_snpCol_1    [rsID]       column of matching snp ID in snp annotation 1
annot_snp_snpCol_2    [rsID]       column of matching snp ID in snp annotation 2
scan_exclude_file_1   [NA]         vector of scanID to exclude from dataset 1
scan_exclude_file_2   [NA]         vector of scanID to exclude from dataset 2
snp_include_file      [NA]         vector of matching snp ID to include
"""
parser = OptionParser(usage=usage)
parser.add_option("-p", "--pipeline", dest="pipeline",
                  default="/projects/geneva/geneva_sata/GCC_code/QCpipeline",
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

for type in ["all", "minor", "miss2fail"]:
    job = "dup_disc_2sets"
    rscript = os.path.join(pipeline, "R", job + ".R")
    jobid[type] = QCpipeline.submitJob("disc_"+type, driver, [rscript, config, type], queue=qname, email=email)

    job = "disc_snp_bins"
    rscript = os.path.join(pipeline, "R", job + ".R")
    jobid[type+"_summary"] = QCpipeline.submitJob("disc_"+type+"_summary", driver, [rscript, config, type], holdid=jobid[type], queue=qname, email=email)
