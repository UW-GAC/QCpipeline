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

if email is not None:
    emailStr = "-m e -M " + email
else:
    emailStr = ""

driver = os.path.join(pipeline, "runRscript.sh")

jobid = dict()

# make this optional in case we run script repeatedly for testing
if not skipSD:
    job = "baf_variance"
    rscript = os.path.join(pipeline, job + ".R")
    qsub = "qsub %s -N %s %s %s %s" % (emailStr, job, driver, rscript, config)
    process = subprocess.Popen(qsub, shell=True, stdout=subprocess.PIPE)
    pipe = process.stdout
    qsubout = pipe.readline()
    jobid[job] = qsubout.split()[2]
    print qsubout


job = "anom_baf"
rscript = os.path.join(pipeline, job + ".R")
jobid[job] = []
istart = start
iend = start + by - 1
while end >= istart:
    if not skipSD:
        qsub = "qsub -hold_jid %s %s -N %s_%i %s %s %s %s %s" % (jobid["baf_variance"], emailStr, job, istart, driver, rscript, config, istart, iend)
    else:
        qsub = "qsub %s -N %s_%i %s %s %s %s %s" % (emailStr, job, istart, driver, rscript, config, istart, iend)
    process = subprocess.Popen(qsub, shell=True, stdout=subprocess.PIPE)
    pipe = process.stdout
    qsubout = pipe.readline()
    jobid[job].append(qsubout.split()[2])
    print qsubout

    istart = istart + by
    iend = istart + by - 1
    if (iend > end):
        iend = end

job = "anom_combine"
rscript = os.path.join(pipeline, job + ".R")
qsub = "qsub -hold_jid %s %s -N %s_baf %s %s %s BAF %s %s" % (','.join(jobid["anom_baf"]), emailStr, job, driver, rscript, config, end, by)
process = subprocess.Popen(qsub, shell=True, stdout=subprocess.PIPE)
pipe = process.stdout
qsubout = pipe.readline()
jobid[job + "_baf"] = qsubout.split()[2]
print qsubout


job = "anom_loh"
rscript = os.path.join(pipeline, job + ".R")
jobid[job] = []
istart = start
iend = start + by - 1
while end >= istart:
    qsub = "qsub -hold_jid %s %s -N %s_%i %s %s %s %s %s" % (jobid["anom_combine_baf"], emailStr, job, istart, driver, rscript, config, istart, iend)
    process = subprocess.Popen(qsub, shell=True, stdout=subprocess.PIPE)
    pipe = process.stdout
    qsubout = pipe.readline()
    jobid[job].append(qsubout.split()[2])
    print qsubout

    istart = istart + by
    iend = istart + by - 1
    if (iend > end):
        iend = end

job = "anom_combine"
rscript = os.path.join(pipeline, job + ".R")
qsub = "qsub -hold_jid %s %s -N %s_loh %s %s %s LOH %s %s" % (','.join(jobid["anom_loh"]), emailStr, job, driver, rscript, config, end, by)
process = subprocess.Popen(qsub, shell=True, stdout=subprocess.PIPE)
pipe = process.stdout
qsubout = pipe.readline()
jobid[job + "_loh"] = qsubout.split()[2]
print qsubout


job = "anom_stats"
rscript = os.path.join(pipeline, job + ".R")
qsub = "qsub -hold_jid %s %s -N %s %s %s %s" % (jobid["anom_combine_loh"], emailStr, job, driver, rscript, config)
process = subprocess.Popen(qsub, shell=True, stdout=subprocess.PIPE)
pipe = process.stdout
qsubout = pipe.readline()
jobid[job] = qsubout.split()[2]
print qsubout
