#! /usr/local/bin/python2.7

"""Association tests"""

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
covar.list         file with list of covariates for each model (e.g., list(c("sex"), c("sex", "race")))
covars_as_factor   covariates to be cast as factors, quoted and space delimited (e.g., "sex race")
gene_action        gene actions for each model (usually additive), quoted and space delimited
model_type         model types (logistic or linear), quoted and space delimited
outcome            outcome variable for logistic regression (0=control, 1=case), quoted and space delimited
nc_geno_file       genotype netCDF file (filtered subject-level recommended)
nc_samp_geno_file  sample-level genotype netCDF file for plots
nc_samp_xy_file    sample-level XY intensity netCDF file for plots

Optional config parameters [default]:    
annot_snp_filtCol  [quality.filter]  column for quality filter (T/F) in snp annotation
annot_snp_rsIDCol  [rsID]            column for rsID in snp annotation
maf.filter         [0.02]            MAF filter threshold to apply to plots
plot_chroms        [NA]              integer vector of chromosomes to plot (NA=all)
scan_chrom_filter  [NA]              scan-chromosome filter matrix
scan_exclude       [NA]              vector of scanID to exclude
signif_line        [5e-08]           genome-wide significance level for manhattan plot
out_assoc_prefix   [assoc]           output prefix for association results       
out_plot_prefix    [assoc]           output prefix for plots"""
parser = OptionParser(usage=usage)
parser.add_option("-p", "--pipeline", dest="pipeline",
                  default="/projects/geneva/geneva_sata/GCC_code/QCpipeline",
                  help="pipeline source directory")
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
parser.add_option("--plotQQManh", dest="plotqq",
                  action="store_true", default=False,
                  help="QQ and Manhattan plots")
parser.add_option("--plotClust", dest="plotcl",
                  action="store_true", default=False,
                  help="cluster plots")
(options, args) = parser.parse_args()

config = args[0]
pipeline = options.pipeline
email = options.email
assoc = options.assoc
merge = options.merge
plotq = options.plotqq
plotc = options.plotcl
qname = options.qname

# 3 arguments (config file, starting chrom,and end chrom) when assoc = True
# require 3 arguments all the time (config, start chrom, and end chrom)
if (len(args) == 3):
    cStart = int(args[1]) 
    cEnd = int(args[2]) 

print cStart, cEnd

sys.path.append(pipeline)
import QCpipeline

configdict = QCpipeline.readConfig(config)

driver = os.path.join(pipeline, "runRscript.sh")

jobid = dict()

if assoc:
    # run by chrom
    job = "run.assoc"
    rscript = os.path.join(pipeline, "R", job + ".R")
    # range of chroms
    chroms = []
    if cStart <= cEnd:
        chroms = range(cStart, cEnd+1)
    else:
        sys.exit("cEnd is smaller than cStart")

    for ichrom in chroms:
        jobid[job+"."+str(ichrom)] = QCpipeline.submitJob(job+".chrom"+str(ichrom), driver, [rscript, config, str(ichrom)], queue=qname, email=email)

if merge:
    job = "merge.chroms"
    rscript = os.path.join(pipeline, "R", job + ".R")
    # generate holdid list 
    if assoc: # need to wait till association tests finish running
        holdid = []
        for i in range(cStart, cEnd+1):
            holdid = holdid + [jobid["run.assoc." + str(i)]]
        #print "hold id for merge: "
        #print holdid
        jobid[job] = QCpipeline.submitJob(job, driver, [rscript, config, str(cStart), str(cEnd)], holdid=holdid, queue=qname, email=email)
    else: # assoc == False means association tests were run in a previous run and will not carry over holdids
        jobid[job] = QCpipeline.submitJob(job, driver, [rscript, config, str(cStart), str(cEnd)], queue=qname, email=email)
    #print jobid
        
if plotq:
    job = "plot.qq.manh"
    rscript = os.path.join(pipeline, "R", job + ".R")
    # generate holdid list 
    if merge:
        holdid = [jobid["merge.chroms"]]
        #print "hold id for plotq: "
        #print holdid
        jobid[job] = QCpipeline.submitJob(job, driver, [rscript, config], holdid=holdid, queue=qname, email=email)
    else:
        jobid[job] = QCpipeline.submitJob(job, driver, [rscript, config], queue=qname, email=email)

if plotc:
    job = "plot.cluster"
    rscript = os.path.join(pipeline, "R", job + ".R")
    # generate holdid list 
    holdid = []
    if merge:
        holdid = [jobid["merge.chroms"]]
        #print "hold id for plotc: "
        #print holdid
    if plotq:
        holdid = holdid + [jobid["plot.qq.manh"]]
        #print "hold id for plotc: "
        #print holdid
    if len(holdid)> 0:    
        jobid[job] = QCpipeline.submitJob(job, driver, [rscript, config], holdid=holdid, queue=qname, email=email)
    else:
        jobid[job] = QCpipeline.submitJob(job, driver, [rscript, config], queue=qname, email=email)

