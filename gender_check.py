#! /usr/local/bin/python2.7

"""Heterozygosity and mean intensity calculations for gender check"""

import QCpipeline
import sys
import os
import subprocess
from optparse import OptionParser

usage = """%prog [options] config

Calculate heterozygosity and mean intensity by scan and chromosome.
Plot X vs Y intensity, X and Y intensity vs X heterozyogsity, and
X vs autosomal heterozygosity with annotated males and females color-coded.
Plot BAF/LRR for X and XY SNPs for sex chromosome anomalies identified in comments.

Required config parameters:
annot_scan_file  scan annotation file
annot_snp_file   snp annotation file
build            genome build (hg18 or hg19)
geno_file        genotype file (netCDF or GDS)
xy_file          XY intensity file (netCDF or GDS)
bl_file          BAF/LRR file (netCDF or GDS)

Optional config parameters [default]:
annot_scan_commentCol       [Comment]            column of comments in scan annotation incuding sex chrom anomaly codes (XXX, XXY, XYY, XO)
annot_scan_localIDCol       [local.scanID]       column of local scanID in scan annotation (for anom plot titles)
annot_scan_sexCol           [sex]                column of annotated sex (M/F/NA) in scan annotation
annot_snp_IntensityOnlyCol  [NA]                 column of intensity-only in snp annotation
out_het_file                [het_by_scan.RData]  output heterozygosity by scan and chromosome
out_inten_file              [mean_inten.RData]   output mean intensity by scan and chromosome
out_pdf                     [sex_check.pdf]      output sex check plot
out_sexchrom_prefix         [sexchrom_anom]      output prefix for sex chrom anomaly plots
out_autosome_prefix         [autosome]           output prefix for mean autosome intensity plots"""
parser = OptionParser(usage=usage)
parser.add_option("-e", "--email", dest="email", default=None,
                  help="email address for job reporting")
parser.add_option("-q", "--queue", dest="qname", default="gcc.q", 
                  help="cluster queue name [default %default]")
(options, args) = parser.parse_args()

if len(args) != 1:
    parser.error("incorrect number of arguments")

config = args[0]
email = options.email
qname = options.qname

pipeline = os.path.dirname(os.path.abspath(sys.argv[0]))

driver = os.path.join(pipeline, "runRscript.sh")

jobid = dict()
for job in ["het_by_scan", "mean_inten", "sexchrom_anom_plot"]:
    rscript = os.path.join(pipeline, "R", job + ".R")
    jobid[job] = QCpipeline.submitJob(job, driver, [rscript, config], queue=qname, email=email)

job = "gender_plot"
rscript = os.path.join(pipeline, "R", job + ".R")
jobid[job] = QCpipeline.submitJob(job, driver, [rscript, config], holdid=[jobid['het_by_scan'], jobid['mean_inten']], queue=qname, email=email)
