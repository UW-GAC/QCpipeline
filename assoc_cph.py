#! /usr/local/bin/python2.7

"""Association tests"""

import QCpipeline
import sys
import os
import subprocess
from optparse import OptionParser

usage = """%prog [options] config chromStart chromEnd

Preliminary association tests with the following steps:
1) Option "assoc": run selected chromosomes from chromStart to chromEnd in parallel
2) Option "merge": merge results from chromosomes chromStart to chromEnd
3) Option "plotQQManh": make QQ and manhattan plots
4) Option "plotClust": make genotype cluster plots of top hits

Required config parameters:
annot_scan_file    scan annotation file
annot_snp_file     snp annotation file
covars             covariates for each model, quoted and space delimited (e.g., "sex race")
covars_as_factor   covariates to be cast as factors, quoted and space delimited (e.g., "sex race")
model_type         model type (survival)
event              event (outcome) variable
time_to_event      time to event variable
geno_file          genotype file (netCDF or GDS, filtered subject-level recommended)
samp_geno_file     sample-level genotype file (netCDF or GDS) for plots
samp_xy_file       sample-level XY intensity file (netCDF or GDS) for plots

Optional config parameters [default]:    
annot_snp_filtCol       [quality.filter]  column for quality filter (T/F) in snp annotation
annot_snp_rsIDCol       [rsID]            column for rsID in snp annotation
block_size              [5000]            block size for reading in SNPs
effect_allele           [minor]           effect allele ("minor" or "alleleA")
gene_action             [additive]        gene action
ivar                    [NA]              interaction variable
maf.filter.type         [snp.specific]    type of MAF filter to apply ("absolute" or "snp.specific")
maf.absolute.threshold  [0.02]            absolute MAF filter threshold to apply to plots: MAF > x
maf.survival.threshold  [75]              survival snp.specific filter threshold: 2*MAF*(1-MAF)*N.events > x
plot_chroms             [NA]              chromosomes to plot, quoted and space delimited, ":" for range (e.g., "1:22 23 25") (NA=all)
scan_exclude            [NA]              vector of scanID to exclude
signif_line             [5e-08]           genome-wide significance level for manhattan plot
out_assoc_prefix        [assoc]           output prefix for association results       
out_plot_prefix         [assoc]           output prefix for plots"""
parser = OptionParser(usage=usage)
parser.add_option("-e", "--email", dest="email", default=None,
                  help="email address for job reporting")
parser.add_option("-q", "--queue", dest="qname", default="all.q", 
                  help="cluster queue name [default %default]")
parser.add_option("-a", "--assoc", dest="assoc",
                  action="store_true", default=False,
                  help="run association tests for chroms chromStart - chromEnd")
parser.add_option("-m", "--merge", dest="merge",
                  action="store_true", default=False,
                  help="merge output for chromosomes chromStart - chromEnd")
parser.add_option("-o", "--options", dest="qsubOptions", default="",
                  help="additional options to pass to qsub, excluding -hold_jid, -N, -m e -M, -N, and -q")
parser.add_option("--plotQQManh", dest="plotqq",
                  action="store_true", default=False,
                  help="QQ and Manhattan plots")
parser.add_option("--plotClust", dest="plotcl",
                  action="store_true", default=False,
                  help="cluster plots")
(options, args) = parser.parse_args()

email = options.email
assoc = options.assoc
merge = options.merge
plotq = options.plotqq
plotc = options.plotcl
qname = options.qname
qsubOptions = options.qsubOptions

# 3 arguments (config file, starting chrom, and end chrom) when assoc = True
if len(args) < 1:
    parser.error("incorrect number of arguments")
if (len(args) != 3) and (assoc or merge):
    parser.error("incorrect number of arguments")

config = args[0]
if len(args) > 1:
    cStart = args[1] 
    cEnd = args[2] 

pipeline = os.path.dirname(os.path.abspath(sys.argv[0]))

configdict = QCpipeline.readConfig(config)

driver = os.path.join(pipeline, "runRscript.sh")
driver_array = os.path.join(pipeline, "runRscript_array.sh")

jobid = dict()

if assoc:
    # run by chrom
    job = "assoc_cph"
    rscript = os.path.join(pipeline, "R", job + ".R")
    # range of chroms
    if int(cStart) <= int(cEnd):
        arrayRange = ":".join([cStart, cEnd])
    else:
        sys.exit("cEnd is smaller than cStart")

    jobid[job] = QCpipeline.submitJob(job, driver_array, [rscript, config], queue=qname, email=email, qsubOptions=qsubOptions, arrayRange=arrayRange)

if merge:
    job = "assoc_combine"
    rscript = os.path.join(pipeline, "R", job + ".R")
    # generate holdid list 
    if assoc: # need to wait till association tests finish running
        holdid = [jobid["assoc_cph"].split(".")[0]]
        jobid[job] = QCpipeline.submitJob(job, driver, [rscript, config, cStart, cEnd], holdid=holdid, queue=qname, email=email, qsubOptions=qsubOptions)
    else: # assoc == False means association tests were run in a previous run and will not carry over holdids
        jobid[job] = QCpipeline.submitJob(job, driver, [rscript, config, cStart, cEnd], queue=qname, email=email, qsubOptions=qsubOptions)
        
if plotq:
    job = "plot_qq_manh"
    rscript = os.path.join(pipeline, "R", job + ".R")
    # generate holdid list 
    if merge:
        holdid = [jobid["assoc_combine"]]
        jobid[job] = QCpipeline.submitJob(job, driver, [rscript, config], holdid=holdid, queue=qname, email=email, qsubOptions=qsubOptions)
    else:
        jobid[job] = QCpipeline.submitJob(job, driver, [rscript, config], queue=qname, email=email, qsubOptions=qsubOptions)

if plotc:
    job = "plot_cluster"
    rscript = os.path.join(pipeline, "R", job + ".R")
    # generate holdid list 
    holdid = []
    if merge:
        holdid = [jobid["assoc_combine"]]
        jobid[job] = QCpipeline.submitJob(job, driver, [rscript, config], holdid=holdid, queue=qname, email=email, qsubOptions=qsubOptions)
    else:
        jobid[job] = QCpipeline.submitJob(job, driver, [rscript, config], queue=qname, email=email, qsubOptions=qsubOptions)

