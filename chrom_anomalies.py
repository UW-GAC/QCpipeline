#! /usr/local/bin/python2.7

"""Chromosome anomalies"""

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
parser.add_option("-m", "--maf", dest="maf", default=None,
                  help="do not use SNPs with MAF below this threshold")
parser.add_option("-q", "--queue", dest="qname",
                  default="gcc.q", help="cluster queue name")
(options, args) = parser.parse_args()

if len(args) != 4:
    parser.error("incorrect number of arguments")

config = args[0]
start = int(args[1])
end = int(args[2])
by = int(args[3])
pipeline = options.pipeline
email = options.email
maf = options.maf
qname = options.qname

sys.path.append(pipeline)
import QCpipeline

configdict = QCpipeline.readConfig(config)

driver = os.path.join(pipeline, "runRscript.sh")

jobid = dict()

# skip BAF if file already exists
waitSD = False
if os.path.exists(configdict['out_baf_med_file']):
    print "using BAF SD file " + configdict['out_baf_med_file']
else:
    job = "baf_variance"
    rscript = os.path.join(pipeline, job + ".R")
    jobid[job] = QCpipeline.submitJob(job, driver, [rscript, config], queue=qname, email=email)
    waitSD = True

# if MAF threshold is needed, calculate allele freq
useMAF = False
waitAfreq = False
if maf is not None:
    useMAF = True
    if os.path.exists(configdict['out_afreq_file']):
        print "using allele freq file " + configdict['out_afreq_file']
    else:
        job = "allele_freq"
        rscript = os.path.join(pipeline, job + ".R")
        jobid[job] = QCpipeline.submitJob(job, driver, [rscript, config], queue=qname, email=email)
        waitAfreq = True

job = "anom_baf"
rscript = os.path.join(pipeline, job + ".R")
jobid[job] = []
istart = start
iend = start + by - 1
while end >= istart:
    if useMAF:
        rargs = [rscript, config, str(istart), str(iend), maf]
    else:
        rargs = [rscript, config, str(istart), str(iend)]

    if waitSD and waitAfreq:
        holdid = [jobid["baf_variance"], jobid["allele_freq"]]
    elif waitSD:
        holdid = [jobid["baf_variance"]]
    elif waitAfreq:
        holdid = [jobid["allele_freq"]]
    else:
        holdid = None

    jobid[job].append(QCpipeline.submitJob(job+str(istart), driver, rargs, holdid=holdid, queue=qname, email=email))

    istart = istart + by
    iend = istart + by - 1
    if (iend > end):
        iend = end

job = "anom_combine"
rscript = os.path.join(pipeline, job + ".R")
jobid[job + "_baf"] = QCpipeline.submitJob(job+"_baf", driver, [rscript, config, "BAF", str(end), str(by)], holdid=jobid['anom_baf'], queue=qname, email=email)


job = "anom_loh"
rscript = os.path.join(pipeline, job + ".R")
jobid[job] = []
istart = start
iend = start + by - 1
while end >= istart:
    if useMAF:
        rargs = [rscript, config, str(istart), str(iend), maf]
    else:
        rargs = [rscript, config, str(istart), str(iend)]

    jobid[job].append(QCpipeline.submitJob(job+str(istart), driver, rargs, holdid=[jobid["anom_combine_baf"]], queue=qname, email=email))

    istart = istart + by
    iend = istart + by - 1
    if (iend > end):
        iend = end

job = "anom_combine"
rscript = os.path.join(pipeline, job + ".R")
jobid[job + "_loh"] = QCpipeline.submitJob(job+"_loh", driver, [rscript, config, "LOH", str(end), str(by)], holdid=jobid['anom_loh'], queue=qname, email=email)


job = "anom_stats"
rscript = os.path.join(pipeline, job + ".R")
if useMAF:
    rargs = [rscript, config, maf]
else:
    rargs = [rscript, config]
jobid[job] = QCpipeline.submitJob(job, driver, rargs, holdid=[jobid["anom_combine_loh"]], queue=qname, email=email)
