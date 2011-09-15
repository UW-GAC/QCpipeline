"""Heterozygosity and mean intensity calculations for gender check"""

import sys
import os
import subprocess
from optparse import OptionParser

usage = """python %prog [options] config"""
parser = OptionParser(usage=usage)
parser.add_option("-p", "--pipeline", dest="pipeline",
                  default="/projects/geneva/geneva_sata/GCC_code/QCpipeline",
                  help="pipeline source directory")
parser.add_option("-e", "--email", dest="email", default=None,
                  help="email address for job reporting")
(options, args) = parser.parse_args()

if len(args) != 1:
    parser.error("incorrect number of arguments")

config = args[0]
pipeline = options.pipeline
email = options.email

if email is not None:
    emailStr = "-m e -M " + email
else:
    emailStr = ""

driver = os.path.join(pipeline, "runRscript.sh")

jobid = dict()
for job in ["het_by_scan", "mean_inten"]:
    rscript = os.path.join(pipeline, job + ".R")
    qsub = "qsub %s -N %s %s %s %s" % (emailStr, job, driver, rscript, config)
    process = subprocess.Popen(qsub, shell=True, stdout=subprocess.PIPE)
    pipe = process.stdout
    qsubout = pipe.readline()
    jobid[job] = qsubout.split()[2]
    print qsubout

job = "gender_plot"
rscript = os.path.join(pipeline, job + ".R")
qsub = "qsub -hold_jid %s,%s %s -N %s %s %s %s" % (jobid['het_by_scan'], jobid['mean_inten'], emailStr, job, driver, rscript, config)
retcode = subprocess.call(qsub, shell=True)
if retcode != 0:
    sys.exit("job submission failed for command\n" + qsub)
    
