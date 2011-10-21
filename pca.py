#! /usr/local/bin/python2.7

"""Principal Component Analysis"""

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
parser.add_option("-c", "--combined", dest="combined",
                  action="store_true", default=False,
                  help="make combined dataset")
parser.add_option("-s", "--skipLD", dest="skipLD",
                  action="store_true", default=False,
                  help="skip LD pruning (use only if LD pruning output file already exists)")
(options, args) = parser.parse_args()

if len(args) != 1:
    parser.error("incorrect number of arguments")

config = args[0]
pipeline = options.pipeline
email = options.email
combined = options.combined
skipLD = options.skipLD

sys.path.append(pipeline)
import QCpipeline

driver = os.path.join(pipeline, "runRscript.sh")

jobid = dict()

if combined:
    job = "dup_disc_ext"
    rscript = os.path.join(pipeline, job + ".R")
    jobid[job] = QCpipeline.submitJob(job, driver, [rscript, config], email=email)

    job = "combine_gds"
    rscript = os.path.join(pipeline, job + ".R")
    jobid[job] = QCpipeline.submitJob(job, driver, [rscript, config], email=email)

if not skipLD:
    job = "ld_pruning"
    rscript = os.path.join(pipeline, job + ".R")
    jobid[job] = QCpipeline.submitJob(job, driver, [rscript, config], email=email)

if combined:
    job = "pca_combined"
    rscript = os.path.join(pipeline, job + ".R")
    holdid = [jobid['dup_disc_ext'], jobid['combine_gds']]
    if not skipLD:
        holdid.append(jobid['ld_pruning'])
    jobid[job] = QCpipeline.submitJob(job, driver, [rscript, config], holdid=holdid, email=email)

    job = "pca_combined_plots"
    rscript = os.path.join(pipeline, job + ".R")
    jobid[job] = QCpipeline.submitJob(job, driver, [rscript, config], holdid=[jobid['pca_combined']], email=email)


job = "pca_study"
rscript = os.path.join(pipeline, job + ".R")
if not skipLD:
    holdid = [jobid['ld_pruning']]
else:
    holdid = None
jobid[job] = QCpipeline.submitJob(job, driver, [rscript, config], holdid=holdid, email=email)

job = "pca_study_plots"
rscript = os.path.join(pipeline, job + ".R")
jobid[job] = QCpipeline.submitJob(job, driver, [rscript, config], holdid=[jobid['pca_study']], email=email)
