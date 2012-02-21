#! /usr/local/bin/python2.7

"""SNP filters: allele frequency, duplicate discordance, Mendelian errors"""

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
parser.add_option("--MAconc", dest="mac",
                  action="store_true", default=False,
                  help="minor allele concordance")
parser.add_option("--dupSNP", dest="dupsnp",
                  action="store_true", default=False,
                  help="duplicate SNP discordance")
(options, args) = parser.parse_args()

if len(args) != 1:
    parser.error("incorrect number of arguments")

config = args[0]
pipeline = options.pipeline
email = options.email
qname = options.qname
mac = options.mac
dupsnp = options.dupsnp

sys.path.append(pipeline)
import QCpipeline

driver = os.path.join(pipeline, "runRscript.sh")

jobid = dict()
job = "allele_freq"
rscript = os.path.join(pipeline, job + ".R")
jobid[job] = QCpipeline.submitJob(job, driver, [rscript, config], queue=qname, email=email)

job = "dup_disc"
rscript = os.path.join(pipeline, job + ".R")
jobid[job] = QCpipeline.submitJob(job, driver, [rscript, config], holdid=[jobid["allele_freq"]], queue=qname, email=email)

if mac:
    job = "dup_disc_maf"
    rscript = os.path.join(pipeline, job + ".R")
    jobid[job] = QCpipeline.submitJob(job, driver, [rscript, config], holdid=[jobid["allele_freq"]], queue=qname, email=email)

job = "snp_filt_plots"
rscript = os.path.join(pipeline, job + ".R")
holdid = [jobid["dup_disc"]]
if mac:
    holdid.append(jobid["dup_disc_maf"])
jobid[job] = QCpipeline.submitJob(job, driver, [rscript, config], holdid=holdid, queue=qname, email=email)

job = "mendel_err"
rscript = os.path.join(pipeline, job + ".R")
jobid[job] = QCpipeline.submitJob(job, driver, [rscript, config], queue=qname, email=email)

if dupsnp:
    job = "dup_snps"
    rscript = os.path.join(pipeline, job + ".R")
    jobid[job] = QCpipeline.submitJob(job, driver, [rscript, config], queue=qname, email=email)
