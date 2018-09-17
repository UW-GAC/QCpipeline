#! /usr/local/bin/python2.7

"""Sample quality checks"""

import QCpipeline
import sys
import os
import subprocess
from optparse import OptionParser

usage = """%prog [options] config

Make BAF/LRR plots for 4 types of samples: BAF SD outliers, 
BAF asymmetry > 0.05, heterozygosity by race outliers, and
high IBD connectivity samples.

Required config parameters:
annot_scan_file     scan annotation file
annot_scan_raceCol  column of race in scan annotation
annot_snp_file      snp annotation file
baf_mean_file       mean of BAF file (created by chrom_anomalies.py)
baf_sd_file         SD of BAF file (created by chrom_anomalies.py)
bl_file             BAF/LRR file (netCDF or GDS)

Optional config parameters [default]:
annot_scan_hetACol    [het.A]             column of autosomal heterozygosity in scan annotation
annot_snp_missingCol  [missing.n1]        column of missing call rate in snp annotation
ibd_con_file          [NA]                vector of scanID with high IBD connectivity
plot_all_unknown      [TRUE]              should all samples with unknown race be plotted?
race_unknown          [NA]                value of race corresponding to unknown or missing
range_het             [1.5]               multiple of interquartile range for selecting SD outliers
range_sd              [1.5]               multiple of interquartile range for selecting het outliers
out_flagged_file      [qual_check.RData]  output file with flagged scanIDs  
out_plot_prefix       [qual_check]        output prefix for BAF/LRR plots
out_baf_sd_boxplot    [baf_sd.pdf]        output BAF SD boxplot
out_het_boxplot       [het_outl.pdf]      output heterozygosity by race boxplot"""
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

job = "sample_quality"
rscript = os.path.join(pipeline, "R", job + ".R")
jobid = QCpipeline.submitJob(job, driver, [rscript, config], queue=qname, email=email, qsubOptions=qsubOptions)
