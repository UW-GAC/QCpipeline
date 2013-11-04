#! /usr/local/bin/python2.7

"""SNP filters: allele frequency and heterozygosity by ethnic group and sex"""

import sys
import os
import subprocess
from optparse import OptionParser

usage = """%prog [options] config

SNP filters by ethnic group and sex:
1) allele frequency
2) heterozygosity

Required config parameters:
annot_scan_file  scan annotation file
nc_geno_file     genotype netCDF file

Optional config parameters [default]:
scan_exclude_file   [NA]                 vector of scanID to exclude (all but one ethnicity)
out_afreq_file      [allele_freq.RData]  output file for allele frequency
out_het_file        [het_by_snp.RData]   output file for heterozygosity by snp"""
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
job = "allele_freq"
rscript = os.path.join(pipeline, "R", job + ".R")
jobid[job] = QCpipeline.submitJob(job, driver, [rscript, config], queue=qname, email=email)

job = "het_by_snp"
rscript = os.path.join(pipeline, "R", job + ".R")
jobid[job] = QCpipeline.submitJob(job, driver, [rscript, config], queue=qname, email=email)

