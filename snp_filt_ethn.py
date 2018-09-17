#! /usr/local/bin/python2.7

"""SNP filters: allele frequency and heterozygosity by ethnic group and sex"""

import QCpipeline
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
geno_file        subject-level genotype file (netCDF or GDS)

Optional config parameters [default]:
scan_exclude_file   [NA]                 vector of scanID to exclude (all but one ethnicity)
out_afreq_file      [allele_freq.RData]  output file for allele frequency
out_het_file        [het_by_snp.RData]   output file for heterozygosity by snp"""
parser = OptionParser(usage=usage)
parser.add_option("-e", "--email", dest="email", default=None,
                  help="email address for job reporting")
parser.add_option("-q", "--queue", dest="qname", default="gcc.q", 
                  help="cluster queue name [default %default]")
parser.add_option("-o", "--options", dest="qsubOptions", default="",
                  help="additional options to pass to qsub, excluding -hold_jid, -N, -m e -M, -N, and -q")
(options, args) = parser.parse_args()

if len(args) != 1:
    parser.error("incorrect number of arguments")

config = args[0]
email = options.email
qname = options.qname
qsubOptions = options.qsubOptions

pipeline = os.path.dirname(os.path.abspath(sys.argv[0]))
driver = os.path.join(pipeline, "runRscript.sh")

jobid = dict()
job = "allele_freq"
rscript = os.path.join(pipeline, "R", job + ".R")
jobid[job] = QCpipeline.submitJob(job, driver, [rscript, config], queue=qname, email=email, qsubOptions=qsubOptions)

job = "het_by_snp"
rscript = os.path.join(pipeline, "R", job + ".R")
jobid[job] = QCpipeline.submitJob(job, driver, [rscript, config], queue=qname, email=email, qsubOptions=qsubOptions)

