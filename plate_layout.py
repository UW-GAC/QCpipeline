#! /usr/local/bin/python2.7

"""Plate layout checks"""

import sys
import os
import subprocess
from optparse import OptionParser

usage = """%prog [options] config

Generate plate map layout plots to assess sample identity issues.

Required config parameters:
XXX

Optional config parameters [default]:
XXX
"""
parser = OptionParser(usage=usage)
parser.add_option("-p", "--pipeline", dest="pipeline",
                  default="/projects/geneva/geneva_sata/GCC_code/QCpipeline",
                  help="pipeline source directory")
parser.add_option("-e", "--email", dest="email", default=None,
                  help="email address for job reporting")
parser.add_option("-q", "--queue", dest="qname", default="gcc.q", 
                  help="cluster queue name [default %default]")
(options, args) = parser.parse_args()

if (len(args) != 1):
    parser.error("incorrect number of arguments")

config = args[0]
pipeline = options.pipeline
email = options.email
qname = options.qname

sys.path.append(pipeline)
import QCpipeline

driver = os.path.join(pipeline, "runRscript.sh")


jobid = dict()

job = "plate_layout"
rscript = os.path.join(pipeline, "R", job + ".R")
jobid["plate_layout"] = QCpipeline.submitJob(job, driver, [rscript, config], queue=qname, email=email)