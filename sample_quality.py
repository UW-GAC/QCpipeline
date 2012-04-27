#! /usr/local/bin/python2.7

"""Sample quality checks"""

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
nc_bl_file          BAF/LRR netCDF file

Optional config parameters [default]:
annot_scan_hetACol        [het.A]              column of autosomal heterozygosity in scan annotation
annot_snp_missingCol      [missing.n1]         column of missing call rate in snp annotation
ibd_con_file              [NA]                 vector of scanID with high IBD connectivity
plot_all_unknown          [TRUE]               should all samples with unknown race be plotted?
race_unknown              [NA]                 value of race corresponding to unknown or missing
range_het                 [1.5]                multiple of interquartile range for selecting SD outliers
range_sd                  [1.5]                multiple of interquartile range for selecting het outliers
out_baf_asym              [baf_aysm.RData]     output BAF asymmetry > 0.05 file
out_baf_asym_plot_prefix  [baf_asym]           output prefix for BAF/LRR plots of BAF asymmetry > 0.05
out_baf_sd_boxplot        [baf_sd.pdf]         output BAF SD boxplot
out_baf_sd_outliers       [baf_sd_outl.RData]  output BAF SD outliers file
out_baf_sd_plot_prefix    [baf_sd_outl]        output prefix for BAF/LRR plots of BAF SD outliers
out_het_boxplot           [het_outl.pdf]       output heterozygosity by race boxplot
out_het_outliers          [het_outl.RData]     output heterozygosity by race outliers file
out_het_plot_prefix       [het_outl]           output prefix for BAF/LRR plots of het outliers
out_ibd_plot_prefix       [high_con]           output prefix for BAF/LRR plots of high IBD connectivity"""
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

job = "sample_quality"
rscript = os.path.join(pipeline, job + ".R")
jobid = QCpipeline.submitJob(job, driver, [rscript, config], queue=qname, email=email)
