#! /usr/local/bin/python2.7

"""Hardy-Weinberg Equilibrium test"""

import sys
import os
import subprocess
from optparse import OptionParser

usage = """python %prog [options] config start end by"""
parser = OptionParser(usage=usage)
parser.add_option("-p", "--pipeline", dest="pipeline",
                  default="/projects/geneva/geneva_sata/GCC_code/QCpipeline",
                  help="pipeline source directory")
parser.add_option("-e", "--email", dest="email", default=None,
                  help="email address for job reporting")
parser.add_option("-q", "--queue", dest="qname",
                  default="gcc.q", help="cluster queue name")
(options, args) = parser.parse_args()

if len(args) < 1 or len(args) > 4:
    parser.error("incorrect number of arguments")

config = args[0]
pipeline = options.pipeline
email = options.email
qname = options.qname

# if start and end arguments provided, split by chromosome
split = False
if len(args) > 1:
    split = True
    start = int(args[1])
    end = int(args[2])
    if end > 23:
        end = 23

    if len(args) > 3:
        by = int(args[3])
    else:
        by = 1

    cStart = range(start, end+1, by)
    cEnd = [x+by-1 for x in cStart]
    if cEnd[-1] > 23:
        cEnd[-1] = 23

sys.path.append(pipeline)
import QCpipeline

driver = os.path.join(pipeline, "runRscript.sh")

jobid = dict()
job = "hwe"
rscript = os.path.join(pipeline, job + ".R")
if split:
    jobid[job] = []
    for c in range(0,len(cStart)):
        jobid[job].append(QCpipeline.submitJob(job+str(cStart[c]), driver, [rscript, config, str(cStart[c]), str(cEnd[c])], queue=qname, email=email))
else:
    jobid[job] = QCpipeline.submitJob(job, driver, [rscript, config], queue=qname, email=email)

if split:
    job = "hwe_merge"
    rscript = os.path.join(pipeline, job + ".R")
    jobid[job] = QCpipeline.submitJob(job, driver, [rscript, config, str(start), str(end), str(by)], holdid=jobid['hwe'], queue=qname, email=email)

job = "hwe_plots"
rscript = os.path.join(pipeline, job + ".R")
if split:
    holdid = [jobid['hwe_merge']]
else:
    holdid=jobid['hwe']
jobid[job] = QCpipeline.submitJob(job, driver, [rscript, config], holdid=holdid, queue=qname, email=email)

