#! /usr/local/bin/python2.7

"""Heterozygosity and mean intensity calculations for gender check"""

import sys
import os
import subprocess
from optparse import OptionParser

usage = """%prog [options] config

Calculate heterozygosity and mean intensity by scan and chromosome.
Plot X vs Y intensity, X and Y intensity vs X heterozyogsity, and
X vs autosomal heterozygosity with annotated males and females color-coded.

Required config parameters:
annot_scan_file  scan annotation file
annot_snp_file   snp annotation file
nc_geno_file     genotype netCDF file
nc_xy_file       XY intensity netCDF file

Optional config parameters [default]:
annot_scan_sexCol  [sex]                column of annotated sex (M/F) in scan annotation
out_het_file       [het_by_scan.RData]  output heterozygosity by scan and chromosome
out_inten_file     [mean_inten.RData]   output mean intensity by scan and chromosome
out_pdf            [sex_check.pdf]      output sex check plot"""
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
for job in ["het_by_scan", "mean_inten"]:
    rscript = os.path.join(pipeline, job + ".R")
    jobid[job] = QCpipeline.submitJob(job, driver, [rscript, config], queue=qname, email=email)

job = "gender_plot"
rscript = os.path.join(pipeline, job + ".R")
jobid[job] = QCpipeline.submitJob(job, driver, [rscript, config], holdid=[jobid['het_by_scan'], jobid['mean_inten']], queue=qname, email=email)
