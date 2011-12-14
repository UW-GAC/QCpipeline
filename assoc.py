#! /usr/local/bin/python2.7

"""Association tests"""
# Example:
# python assoc.py /projects/geneva/gcc-fs2/MitchellPak/j7shen/results/analysis/assoc/assoc.config 1 26 --covarsex --sex sex --assoc --merge --plotQQManh --plotClust

import sys
import os
import subprocess
from optparse import OptionParser

usage = """python %prog [options] config chromStart chromEnd"""
parser = OptionParser(usage=usage)
parser.add_option("-p", "--pipeline", dest="pipeline",
                  default="/projects/geneva/geneva_sata/GCC_code/QCpipeline",
                  help="pipeline source directory")
# this is temporary 
# parser.add_option("-l", "--localpipe", dest="localPipe",
#                  default="/projects/geneva/gcc-fs2/MitchellPak/j7shen/R_MitchellPak",
#                  help="local pipeline source directory")
parser.add_option("-e", "--email", dest="email", default=None,
                  help="email address for job reporting")
# to run Y separately without sex in the model
parser.add_option("-c", "--covarsex", dest="covarsex",
                  action="store_true", default=False,
                  help="to run chrom Y without sex in the model if sex is a covariate")
parser.add_option("-x", "--sex", dest="sex", default="sex",
                  help="variable name for sex")
parser.add_option("-a", "--assoc", dest="assoc",
                  action="store_true", default=False,
                  help="run association tests for chroms btw chromStart and chromEnd")
parser.add_option("-m", "--merge", dest="merge",
                  action="store_true", default=False,
                  help="merge output for all chromosomes 1-26")
parser.add_option("--plotQQManh", dest="plotqq",
                  action="store_true", default=False,
                  help="QQ and Manhattan plots")
parser.add_option("--plotClust", dest="plotcl",
                  action="store_true", default=False,
                  help="cluster plots")
parser.add_option("--queue", dest="qname",
                  default="gcc.q", help="cluster plots")
(options, args) = parser.parse_args()

config = args[0]
pipeline = options.pipeline
#localpipe = options.localPipe
email = options.email
covarsex = options.covarsex
sex = options.sex
assoc = options.assoc
merge = options.merge
plotq = options.plotqq
plotc = options.plotcl
qname = options.qname

# 3 arguments (config file, starting chrom,and end chrom) when assoc = True
#if assoc:
#    if (len(args) != 3):
#        parser.error("incorrect number of arguments")
#    else:
#        cStart = max(1,int(args[1])) # lower bound = 1 
#        cEnd = min(26,int(args[2])) # upper bound = 26
# 1 argument when assoc = False
#if ((not assoc) & (merge | plotq | plotc) & (len(args) != 1)):

# require 3 arguments all the time (config, start chrom, and end chrom)
if (len(args) == 3):
    cStart = max(1,int(args[1])) # lower bound = 1 
    cEnd = min(26,int(args[2])) # upper bound = 26

    
sys.path.append(pipeline)
import QCpipeline

configdict = QCpipeline.readConfig(config)

driver = os.path.join(pipeline, "runRscript.sh")

jobid = dict()

if assoc:
    # run by chrom
    job = "run.assoc"
    rscript = os.path.join(pipeline, job + ".R")
    # range of chroms
    chroms = []
    if cStart <= cEnd:
        chroms = range(cStart, cEnd+1)
    else:
        sys.exit("cEnd is smaller than cStart")

    if (cStart == 26 & cEnd == 26):
       jobid[job+".26"] = QCpipeline.submitJob(job+".chrom26", driver, [rscript, config, "26"], queue=qname, email=email)
    elif (cEnd < 25 | (cEnd >= 25 & (not covarsex))): # no need to single out y if cEnd is lt 25 or covarsex = False
        for ichrom in chroms:
            jobid[job+"."+str(ichrom)] = QCpipeline.submitJob(job+".chrom"+str(ichrom), driver, [rscript, config, str(ichrom)], queue=qname, email=email)
    else: # run y without sex 
        jobid[job+".25"] = QCpipeline.submitJob(job+".chrom25", driver, [rscript, config, "25", sex], queue=qname, email=email)
        del chroms[chroms.index(25)]
        for ichrom in chroms:
            jobid[job+"."+str(ichrom)] = QCpipeline.submitJob(job+".chrom"+str(ichrom), driver, [rscript, config, str(ichrom)], queue=qname, email=email)
    #print jobid


if merge:
    job = "merge.chroms"
    rscript = os.path.join(pipeline, job + ".R")
    # generate holdid list 
    if assoc: # need to wait till association tests finish running
        holdid = []
        for i in range(cStart, cEnd+1):
            holdid = holdid + [jobid["run.assoc." + str(i)]]
        #print "hold id for merge: "
        #print holdid
        jobid[job] = QCpipeline.submitJob(job, driver, [rscript, config, str(cStart), str(cEnd)], holdid=holdid, queue="gcc.q", email=email)
    else: # assoc == False means association tests were run in a previous run and will not carry over holdids
        jobid[job] = QCpipeline.submitJob(job, driver, [rscript, config, str(cStart), str(cEnd)], queue="gcc.q", email=email)
    #print jobid
        
if plotq:
    job = "plot.qq.manh"
    rscript = os.path.join(pipeline, job + ".R")
    # generate holdid list 
    if merge:
        holdid = [jobid["merge.chroms"]]
        #print "hold id for plotq: "
        #print holdid
        jobid[job] = QCpipeline.submitJob(job, driver, [rscript, config], holdid=holdid, queue="gcc.q", email=email)
    else:
        jobid[job] = QCpipeline.submitJob(job, driver, [rscript, config], queue="gcc.q", email=email)

if plotc:
    job = "plot.cluster"
    rscript = os.path.join(pipeline, job + ".R")
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
        jobid[job] = QCpipeline.submitJob(job, driver, [rscript, config], holdid=holdid, queue="gcc.q", email=email)
    else:
        jobid[job] = QCpipeline.submitJob(job, driver, [rscript, config], queue="gcc.q", email=email)

