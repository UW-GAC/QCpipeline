#! /usr/local/bin/python2.7

"""Create filtered and unfiltered PLINK files"""

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
parser.add_option("-q", "--queue", dest="qname",
                  default="gcc.q", help="cluster queue name")
parser.add_option("-f", "--filtered", dest="filt",
                  action="store_true", default=False,
                  help="create filtered PLINK file")
(options, args) = parser.parse_args()

if len(args) != 1:
    parser.error("incorrect number of arguments")

config = args[0]
pipeline = options.pipeline
email = options.email
qname = options.qname
filt = options.filt

sys.path.append(pipeline)
import QCpipeline

configdict = QCpipeline.readConfig(config)
plinkfile = configdict["out_plink_prefix"]

driver = os.path.join(pipeline, "runRscript.sh")

jobid = dict()

job = "plink_unfiltered"
rscript = os.path.join(pipeline, job + ".R")
jobid[job] = QCpipeline.submitJob(job, driver, [rscript, config], queue=qname, email=email)

job = "plink_unfiltered_bed"
arglist = ["--noweb", "--make-bed", "--file", plinkfile, "--out", plinkfile]
jobid[job] = QCpipeline.submitJob(job, "plink", arglist, options="-b y -j y -cwd",
                                  holdid=[jobid["plink_unfiltered"]], queue=qname, email=email)

job = "plink_unfiltered_tar"
arglist = ["-cvzf", plinkfile+".tar.gz", plinkfile+".bed", plinkfile+".bim", plinkfile+".fam"]
jobid[job] = QCpipeline.submitJob(job, "tar", arglist, options="-b y -j y -cwd",
                                  holdid=[jobid["plink_unfiltered_bed"]], queue=qname, email=email)

if filt:
    plinkfile = configdict["out_plink_prefix"] + "_filtered"
    job = "plink_filtered"
    rscript = os.path.join(pipeline, job + ".R")
    jobid[job] = QCpipeline.submitJob(job, driver, [rscript, config], queue=qname, email=email)

    job = "plink_filtered_bed"
    arglist = ["--noweb", "--make-bed", "--file", plinkfile, "--out", plinkfile]
    jobid[job] = QCpipeline.submitJob(job, "plink", arglist, options="-b y -j y -cwd",
                                      holdid=[jobid["plink_filtered"]], queue=qname, email=email)

    job = "plink_filtered_tar"
    arglist = ["-cvzf", plinkfile+".tar.gz", plinkfile+".bed", plinkfile+".bim", plinkfile+".fam"]
    jobid[job] = QCpipeline.submitJob(job, "tar", arglist, options="-b y -j y -cwd",
                                      holdid=[jobid["plink_filtered_bed"]], queue=qname, email=email)
