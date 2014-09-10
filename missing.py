#! /usr/local/bin/python2.7

"""Missing call rate by SNP and scan"""

import sys
import os
import subprocess
from optparse import OptionParser

usage = """%prog [options] config

Calculate missing call rate by SNP and scan.

Required config parameters:
annot_scan_file  scan annotation file
annot_snp_file   snp annotation file
geno_file        genotype file (netCDF or GDS)

Optional config parameters [default]:
annot_snp_IntensityOnlyCol  [NA]                 column of intensity-only in snp annotation
round2                      [TRUE]               compute e2 and n2 (as well as e1 and n1)?
scan_exclude_file           [NA]                 vector of scanID to exclude
snp_exclude_file            [NA]                 vector of snpID to exclude
out_e1_file                 [missing.e1.RData]   output missing.e1 data file
out_e1_hist                 [missing_e1.pdf]     output missing.e1 histogram
out_e2_file                 [missing.e2.RData]   output missing.e2 data file
out_e2_hist                 [missing_e2.pdf]     output missing.e2 histogram
out_n1_file                 [missing.n1.RData]   output missing.n1 data file
out_n2_file                 [missing.n2.RData]   output missing.n2 data file
out_snp_summary             [snp_summary.RData]  output snp summary table"""
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

job = "missing"
rscript = os.path.join(pipeline, "R", job + ".R")
jobid = QCpipeline.submitJob(job, driver, [rscript, config], queue=qname, email=email)

