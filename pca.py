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
(options, args) = parser.parse_args()

if len(args) != 1:
    parser.error("incorrect number of arguments")

config = args[0]
pipeline = options.pipeline
email = options.email
combined = options.combined

sys.path.append(pipeline)
import QCpipeline

configdict = QCpipeline.readConfig(config)

driver = os.path.join(pipeline, "runRscript.sh")

jobid = dict()

if combined:
    job = "dup_disc_ext"
    rscript = os.path.join(pipeline, job + ".R")
    jobid[job] = QCpipeline.submitJob(job, driver, [rscript, config], email=email)

    job = "combine_gds"
    rscript = os.path.join(pipeline, job + ".R")
    jobid[job] = QCpipeline.submitJob(job, driver, [rscript, config], email=email)

# skip LD if file already exists
waitLD = False
if os.path.exists(configdict['out_pruned_file']):
    print "using LD pruned file " + configdict['out_pruned_file']
else:
    job = "ld_pruning"
    if combined:
        holdid = [jobid['dup_disc_ext']]
    else:
        holdid = None
    rscript = os.path.join(pipeline, job + ".R")
    jobid[job] = QCpipeline.submitJob(job, driver, [rscript, config], holdid=holdid, email=email)
    waitLD = True

if combined:
    job = "pca_combined"
    rscript = os.path.join(pipeline, job + ".R")
    holdid = [jobid['dup_disc_ext'], jobid['combine_gds']]
    if waitLD:
        holdid.append(jobid['ld_pruning'])
    jobid[job] = QCpipeline.submitJob(job, driver, [rscript, config], holdid=holdid, email=email)

    job = "pca_plots"
    rscript = os.path.join(pipeline, job + ".R")
    jobid[job] = QCpipeline.submitJob(job+"_combined", driver, [rscript, config, "combined"], holdid=[jobid['pca_combined']], email=email)

else:
    job = "pca_study"
    rscript = os.path.join(pipeline, job + ".R")
    if waitLD:
        holdid = [jobid['ld_pruning']]
    else:
        holdid = None
    jobid[job] = QCpipeline.submitJob(job, driver, [rscript, config], holdid=holdid, email=email)

    job = "pca_plots"
    rscript = os.path.join(pipeline, job + ".R")
    jobid[job] = QCpipeline.submitJob(job+"_study", driver, [rscript, config, "study"], holdid=[jobid['pca_study']], email=email)
