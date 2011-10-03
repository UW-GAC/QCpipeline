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
parser.add_option("-s", "--skipSD", dest="skipSD", 
                  action="store_true", default=False,
                  help="skip BAF variance calculations (use only if median BAF SD output file already exists)")
(options, args) = parser.parse_args()

if len(args) != 4:
    parser.error("incorrect number of arguments")

config = args[0]
start = int(args[1])
end = int(args[2])
by = int(args[3])
pipeline = options.pipeline
email = options.email
skipSD = options.skipSD

sys.path.append(pipeline)
import QCpipeline

driver = os.path.join(pipeline, "runRscript.sh")

jobid = dict()

# make this optional in case we run script repeatedly for testing
if not skipSD:
    job = "baf_variance"
    rscript = os.path.join(pipeline, job + ".R")
    jobid[job] = QCpipeline.submitJob(job, driver, [rscript, config], email=email)


job = "anom_baf"
rscript = os.path.join(pipeline, job + ".R")
jobid[job] = []
istart = start
iend = start + by - 1
while end >= istart:
    if not skipSD:
        holdid = [jobid["baf_variance"]]
    else:
        holdid = None
    jobid[job].append(QCpipeline.submitJob(job+str(istart), driver, [rscript, config, str(istart), str(iend)], holdid=holdid, email=email))

    istart = istart + by
    iend = istart + by - 1
    if (iend > end):
        iend = end

job = "anom_combine"
rscript = os.path.join(pipeline, job + ".R")
jobid[job + "_baf"] = QCpipeline.submitJob(job+"_baf", driver, [rscript, config, "BAF", str(end), str(by)], holdid=jobid['anom_baf'], email=email)


job = "anom_loh"
rscript = os.path.join(pipeline, job + ".R")
jobid[job] = []
istart = start
iend = start + by - 1
while end >= istart:
    jobid[job].append(QCpipeline.submitJob(job+str(istart), driver, [rscript, config, str(istart), str(iend)], holdid=[jobid["anom_combine_baf"]], email=email))

    istart = istart + by
    iend = istart + by - 1
    if (iend > end):
        iend = end

job = "anom_combine"
rscript = os.path.join(pipeline, job + ".R")
jobid[job + "_loh"] = QCpipeline.submitJob(job+"_loh", driver, [rscript, config, "LOH", str(end), str(by)], holdid=jobid['anom_loh'], email=email)


job = "anom_stats"
rscript = os.path.join(pipeline, job + ".R")
jobid[job] = QCpipeline.submitJob(job, driver, [rscript, config], holdid=[jobid["anom_combine_loh"]], email=email)
