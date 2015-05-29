#! /usr/local/bin/python2.7

"""Create NetCDF and GDS files"""

import sys
import os
import subprocess
from optparse import OptionParser

usage = """%prog [options] config

Creates GDS files from IMPUTE2 files.
Required config parameters:

Optional config parameters [default]:
"""
parser = OptionParser(usage=usage)
parser.add_option("-e", "--email", dest="email", default=None,
                  help="email address for job reporting")
parser.add_option("-q", "--queue", dest="qname", default="all.q", 
                  help="cluster queue name [default %default]")
parser.add_option("-p", "--pipeline", dest="pipeline",
                  default=None,
                  help="test pipeline path")
parser.add_option("-o", "--overwrite", dest="overwrite",
                  action="store_true", default=False,
                  help="overwrite existing files")
#parser.add_option("-v", "--verbose", dest="verbose",
#                  action="store_true", default=TRUE,
#                  help="verbose mode (print all qsub messages)")
(options, args) = parser.parse_args()

if len(args) != 1:
    parser.error("incorrect number of arguments")

config = args[0]
#pipeline = options.pipeline
email = options.email
pipeline = options.pipeline
overwrite = options.overwrite
qname = options.qname
verbose = True

sys.path.append(pipeline)
import QCpipeline as qcp


# can get chromosomes here, later.
configdict = qcp.readConfig(config)

jobid = dict()


chromosomes = range(1,24)

# check if directory exists, and if not, created it
if not os.path.exists(configdict["out_gds_dir"]):
    os.makedirs(configdict["out_gds_dir"])

# check if directory is empty -- not the best way but.
if (os.listdir(configdict["out_gds_dir"]) != []) and not overwrite:
    sys.exit("Directory " + configdict["out_gds_dir"] + " is not empty.")


# check if any files exist - chromosomes 1-23
for chromosome in chromosomes:
    
    filename = configdict["out_gds_dir"] + "/" + configdict["out_gds_prefix"] + "_chr-%d.gds" % chromosome

    if (not overwrite) & os.path.exists(filename):
        sys.exit(filename + " already exists; use -o flag to overwrite all files")


# check if any files exist - other and failed files
for chromosome in ["other", "failed"]:
    
    filename = configdict["out_gds_dir"] + "/" + configdict["out_gds_prefix"] + "_chr-" + chromosome + ".gds"

    if (not overwrite) & os.path.exists(filename):
        sys.exit(filename + " already exists; use -o flag to overwrite all files")



driver = os.path.join(pipeline, "runRscript_array.sh")
    
args = ""
job = "gds_create_imputed"
rscript = os.path.join(pipeline, "R", "gds_create_imputed.R")

# array job for imputed chromosomes
jobid[job] = qcp.submitJob(job, driver, [rscript, config, args], queue=qname, email=email,
                                    verbose=verbose, arrayRange="1:23")

## OTHER AND FAILED commented out for now
# print "Submitting job for other and failed SNPs"
# 
# # submit "other" snps (non-imputed)
# holdid = [jobid["gds_create_imputed"].split(".")[0]] # because it's from an array job
# args = " ".join(["other", testStr])
# job = "gds_create_other_chr-other"
# rscript = os.path.join(pipeline, "R", "gds_create_other.R")
# jobid[job] = qcp.submitJob(job, driver, [rscript, config, args], queue=qname, email=email,
#                                     holdid=holdid, verbose=verbose)
# 
# # submit "failed" snps
# holdid = [jobid["gds_create_other_chr-other"]] # hold on "other" snps since we need the snp annotation
# args = " ".join(["failed", testStr])
# job = "gds_create_other_chr-failed"
# rscript = os.path.join(pipeline, "R", "gds_create_other.R")
# jobid[job] = qcp.submitJob(job, driver, [rscript, config, args], queue=qname, email=email,
#                                     holdid=holdid, verbose=verbose)

holdid = [jobid["gds_create_imputed"].split(".")[0]]
job = "gds_imputed_cleanup"
rscript = os.path.join(pipeline, "R", "gds_imputed_cleanup.R")

jobid[job] = qcp.submitJob(job, driver, [rscript, config, args], queue=qname, email=email,
                                    holdid=holdid, verbose=verbose)
