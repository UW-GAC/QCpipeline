#! /usr/local/bin/python2.7

"""Sample quality checks"""

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

sys.path.append(pipeline)
import QCpipeline

driver = os.path.join(pipeline, "runRscript.sh")

job = "sample_quality"
rscript = os.path.join(pipeline, job + ".R")
jobid = QCpipeline.submitJob(job, driver, [rscript, config], email=email)
