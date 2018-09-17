#! /usr/local/bin/python2.7

"""Discordance between two data sets"""

import QCpipeline
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
annot_scan_subjCol_1  [subjectID]         column of matching subject ID in scan annotation 1
annot_scan_subjCol_2  [subjectID]         column of matching subject ID in scan annotation 2
annot_snp_snpCol_1    [rsID]              column of matching snp ID in snp annotation 1 (if matching on name)
annot_snp_snpCol_2    [rsID]              column of matching snp ID in snp annotation 2 (if matching on name)
match_snps_on         [position,alleles]  how to match snps (position, alleles, and/or name)
scan_exclude_file_1   [NA]                vector of scanID to exclude from dataset 1
scan_exclude_file_2   [NA]                vector of scanID to exclude from dataset 2
snp_exclude_file_1    [NA]                vector of snpID to exclude from dataset 1
snp_exclude_file_2    [NA]                vector of snpID to exclude from dataset 2
snp_include_file      [NA]                vector of matching snp ID to include
out_summary_prefix    [NA]                prefix for summary output files
summary_include_file  [NA]                vector of snp ID to include in summary
"""
parser = OptionParser(usage=usage)
parser.add_option("-e", "--email", dest="email", default=None,
                  help="email address for job reporting")
parser.add_option("-q", "--queue", dest="qname", default="gcc.q", 
                  help="cluster queue name [default %default]")
parser.add_option("--overall", dest="overall",
                  action="store_true", default=False,
                  help="overall discordance")
parser.add_option("--minor", dest="minor",
                  action="store_true", default=False,
                  help="minor allele discordance")
parser.add_option("--miss1fail", dest="miss1fail",
                  action="store_true", default=False,
                  help="minor allele discordance where missing in dataset 1 counts as discordance")
parser.add_option("--miss2fail", dest="miss2fail",
                  action="store_true", default=False,
                  help="minor allele discordance where missing in dataset 2 counts as discordance")
parser.add_option("--senspec", dest="senspec",
                  action="store_true", default=False,
                  help="minor allele sensitivity and specificity")
parser.add_option("--options", dest="qsubOptions", default="",
                  help="additional options to pass to qsub, excluding -hold_jid, -N, -m e -M, -N, and -q")
(options, args) = parser.parse_args()

if len(args) != 1:
    parser.error("incorrect number of arguments")

config = args[0]
email = options.email
qname = options.qname
qsubOptions = options.qsubOptions

# which tests to run?
tests = []
if options.overall:
    tests.append("all")
if options.minor:
    tests.append("minor")
if options.miss1fail:
    tests.append("miss1fail")
if options.miss2fail:
    tests.append("miss2fail")
if options.senspec:
    tests.append("senspec")

pipeline = os.path.dirname(os.path.abspath(sys.argv[0]))

driver = os.path.join(pipeline, "runRscript.sh")

jobid = dict()

for type in tests:
    job = "dup_disc_2sets"
    rscript = os.path.join(pipeline, "R", job + ".R")
    jobid[type] = QCpipeline.submitJob("disc_"+type, driver, [rscript, config, type], queue=qname, email=email, qsubOptions=qsubOptions)

    # bin by MAF and cluster separation - requires these columns in SNP annotation
    # job = "disc_snp_bins"
    # rscript = os.path.join(pipeline, "R", job + ".R")
    # jobid[type+"_summary"] = QCpipeline.submitJob("disc_"+type+"_summary", driver, [rscript, config, type], holdid=jobid[type], queue=qname, email=email, qsubOptions=qsubOptions)
