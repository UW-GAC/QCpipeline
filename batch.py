#! /usr/local/bin/python2.7

"""Batch quality checks"""

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
parser.add_option("-t", "--type", dest="type", default="chisq",
                  help="test type (chisq or fisher)")
parser.add_option("-q", "--queue", dest="qname",
                  default="gcc.q", help="cluster queue name")
(options, args) = parser.parse_args()

if len(args) != 1:
    parser.error("incorrect number of arguments")

config = args[0]
pipeline = options.pipeline
email = options.email
type = options.type
qname = options.qname

sys.path.append(pipeline)
import QCpipeline

driver = os.path.join(pipeline, "runRscript.sh")

jobid = dict()
if (type == "chisq"):
    job = "batch_chisq"
elif (type == "fisher"):
    job = "batch_fisher"
else:
    sys.exit("test type must be chisq or fisher")

rscript = os.path.join(pipeline, job + ".R")
jobid["batch_test"] = QCpipeline.submitJob(job, driver, [rscript, config], queue=qname, email=email)

job = "batch_plots"
rscript = os.path.join(pipeline, job + ".R")
jobid[job] = QCpipeline.submitJob(job, driver, [rscript, config, type], holdid=[jobid['batch_test']], queue=qname, email=email)
