#! /usr/local/bin/python2.7

"""Hardy-Weinberg Equilibrium test"""

import sys
import os
import subprocess
from optparse import OptionParser

usage = """%prog [options] config chromStart chromEnd

Test for deviations from Hardy-Weinberg Equilbrium and plot results.
chromStart, chromEnd are chromosome numbers to run in parallel (values > 23 are ignored).
Cluster plots are randomly sampled SNPs from bins of pvalue.  There are
3 sets of 9 SNPs from each bin and an additional 9*7 SNPs from the bins 
flanking 1e-4.

Required config parameters:
annot_scan_file    scan annotation file
annot_snp_file     snp annotation file
geno_file          genotype file (netCDF or GDS, filtered subject-level recommended)
samp_geno_file     sample-level genotype file (netCDF or GDS) for plots
samp_xy_file       sample-level XY intensity file (netCDF or GDS) for plots
scan_include_file  vector of scanID to include in HWE

Optional config parameters [default]:
annot_snp_missingCol  [missing.n1]     column of missing call rate in snp annotation
annot_snp_rsIDCol     [rsID]           column of rsID in snp annotation
out_clust_prefix      [hwe_clust]      output prefix for cluster plots
out_hwe_prefix        [hwe]            output prefix for HWE results
out_inbrd_plot        [hwe_inbrd.pdf]  output histogram of inbreeding coefficients
out_maf_plot          [hwe_maf.png]    output plot of MAF vs p value
out_qq_plot           [hwe_qq.png]     output QQ plots (autosomal and X)
out_sim_prefix        [hwe_sim]        prefix for simulated HWE results
nsnp_simulate         [50000]          max number of SNPs to simulate"""
parser = OptionParser(usage=usage)
parser.add_option("-p", "--pipeline", dest="pipeline",
                  default="/projects/geneva/gcc-fs2/GCC_Code/QCpipeline",
                  help="pipeline source directory")
parser.add_option("-e", "--email", dest="email", default=None,
                  help="email address for job reporting")
parser.add_option("-q", "--queue", dest="qname", default="all.q", 
                  help="cluster queue name [default %default]")
(options, args) = parser.parse_args()

if len(args) < 1 or len(args) > 3:
    parser.error("incorrect number of arguments")

config = args[0]
pipeline = options.pipeline
email = options.email
qname = options.qname

if len(args) > 1:
    cStart = args[1] 
    cEnd = args[2] 
    if int(cEnd) > 23:
        cEnd = "23"
else:
    cStart = "1"
    cEnd = "23"
    
if int(cStart) <= int(cEnd):
    arrayRange = ":".join([cStart, cEnd])
else:
    sys.exit("cEnd is smaller than cStart")

sys.path.append(pipeline)
import QCpipeline

driver = os.path.join(pipeline, "runRscript.sh")
driver_array = os.path.join(pipeline, "runRscript_array.sh")

jobid = dict()

job = "hwe"
rscript = os.path.join(pipeline, "R", job + ".R")
jobid[job] = QCpipeline.submitJob(job, driver_array, [rscript, config], queue=qname, email=email, arrayRange=arrayRange)

holdid = [jobid["hwe"].split(".")[0]]
job = "hwe_combine"
rscript = os.path.join(pipeline, "R", job + ".R")
jobid[job] = QCpipeline.submitJob(job, driver, [rscript, config, cStart, cEnd], holdid=holdid, queue=qname, email=email)


holdid = [jobid['hwe_combine']]
job = "hwe_simulate"
rscript = os.path.join(pipeline, "R", job + ".R")
jobid[job] = QCpipeline.submitJob(job, driver, [rscript, config], holdid=holdid, queue=qname, email=email)

job = "hwe_plots"
rscript = os.path.join(pipeline, "R", job + ".R")
holdid = [jobid["hwe_simulate"]]
jobid[job] = QCpipeline.submitJob(job, driver, [rscript, config], holdid=holdid, queue=qname, email=email)

