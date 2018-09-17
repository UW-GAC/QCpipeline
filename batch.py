#! /usr/local/bin/python2.7

"""Batch quality checks"""

import QCpipeline
import sys
import os
import subprocess
from optparse import OptionParser

usage = """%prog [options] config

Make plots to evaluate batch effects.
HapMaps are always excluded.

Required config parameters:
annot_scan_file     scan annotation file
annot_scan_raceCol  column of race in scan annotation
inten_file          mean intensity file (created by gender_check.py)

Optional config parameters [default]:
annot_scan_batchCol         [Sample.Plate]                 column of batch in scan annotation
annot_scan_hapmapCol        [geno.cntl]                    column of hapmap (0/1) in scan annotation
annot_scan_missAutoCol      [miss.e1.auto]                 column of autosomal MCR in scan annotation
annot_scan_redoCol          [Redo.Processing.Plate]        column of redo plate (Y/N) in scan annotation
scan_exclude_file           [NA]                           vector of scanID to exclude
out_hist_plot               [batch_nscan_hist.pdf]         output nscan histogram
out_inten_plot              [batch_chr1inten.pdf]          output chrom 1 intensity boxplot
out_mcr_plot                [batch_mcr.pdf]                output MCR boxplot
out_meanmcr_nscan_plot      [batch_meanmcr_nscan.pdf]      output MCR vs nscan plot"""

## docs longer needed if we are not running allele frequency test
# Run batch test of allele frequency (chisq or fisher) and make plots.
# geno_file           genotype file (netCDF or GDS)
# out_batch_prefix            [batch_test]                   output prefix for test results
# out_lambda_race_plot        [batch_lambda_race.pdf]        output lambda vs race frac plot
# out_meanchisq_nscan_plot    [batch_meanchisq_nscan.pdf]    output chisq vs nscan plot
# out_meanchisq_race_plot     [batch_meanchisq_race.pdf]     output chisq vs race frac plot
# out_meanmcr_meanchisq_plot  [batch_meanmcr_meanchisq.pdf]  output MCR vs chisq plot
# out_meanmcr_meanor_plot     [batch_meanmcr_meanor.pdf]     output MCR vs OR plot
# out_meanor_nscan_plot       [batch_meanor_nscan.pdf]       output OR vs nscan plot
# out_meanor_race_plot        [batch_meanor_race.pdf]        output OR vs race frac plot

parser = OptionParser(usage=usage)
parser.add_option("-e", "--email", dest="email", default=None,
                  help="email address for job reporting")
parser.add_option("-q", "--queue", dest="qname", default="gcc.q", 
                  help="cluster queue name [default %default]")
parser.add_option("-o", "--options", dest="qsubOptions", default="",
                  help="additional options to pass to qsub, excluding -hold_jid, -N, -m e -M, -N, and -q")
# parser.add_option("-t", "--type", dest="type", default="fisher",
#                   help="test type [fisher (default) or chisq]")
(options, args) = parser.parse_args()

if len(args) != 1:
    parser.error("incorrect number of arguments")

config = args[0]
email = options.email
#type = options.type
qname = options.qname
qsubOptions = options.qsubOptions

pipeline = os.path.dirname(os.path.abspath(sys.argv[0]))

driver = os.path.join(pipeline, "runRscript.sh")

jobid = dict()

# job = "batch_test"
# rscript = os.path.join(pipeline, "R", job + ".R")
# jobid[job] = QCpipeline.submitJob(job, driver, [rscript, config, type], queue=qname, email=email, qsubOptions=qsubOptions)

job = "batch_plots"
rscript = os.path.join(pipeline, "R", job + ".R")
# jobid[job] = QCpipeline.submitJob(job, driver, [rscript, config, type], holdid=[jobid['batch_test']], queue=qname, email=email, qsubOptions=qsubOptions)
jobid[job] = QCpipeline.submitJob(job, driver, [rscript, config], queue=qname, email=email, qsubOptions=qsubOptions)
