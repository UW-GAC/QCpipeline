#! /usr/local/bin/python2.7

"""Association tests"""
# Example:
# python assoc.py /projects/geneva/gcc-fs2/MitchellPak/j7shen/results/analysis/assoc/assoc.config 1 26 --covarsex --assoc --merge --plotQQManh --plotClust

import sys
import os
import subprocess
from optparse import OptionParser

usage = """python %prog [options] config chromStart chromEnd"""
parser = OptionParser(usage=usage)
parser.add_option("-p", "--pipeline", dest="pipeline",
                  default="/projects/geneva/geneva_sata/GCC_code/QCpipeline",
                  help="pipeline source directory")
parser.add_option("-e", "--email", dest="email", default=None,
                  help="email address for job reporting")
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
parser.add_option("-q", "--queue", dest="qname",
                  default="gcc.q", help="cluster queue name")
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
    rscript = os.path.join(pipeline, job + ".R")
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

