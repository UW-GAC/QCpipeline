#! /usr/local/bin/python2.7

"""Create NetCDF and GDS files"""

import QCpipeline as qcp
import sys
import os
import subprocess
from optparse import OptionParser

usage = """%prog [options] config chromosomes

Creates GDS files from IMPUTE2 files, ie:
%prog my.config 1:23

Required config parameters:

Optional config parameters [default]:
"""
parser = OptionParser(usage=usage)
parser.add_option("-e", "--email", dest="email", default=None,
                  help="email address for job reporting")
parser.add_option("-q", "--queue", dest="qname", default="all.q", 
                  help="cluster queue name [default %default]")
parser.add_option("--overwrite", dest="overwrite",
                  action="store_true", default=False,
                  help="overwrite existing files")
parser.add_option("--options", dest="qsubOptions", default="",
                  help="additional options to pass to qsub, excluding -hold_jid, -N, -m e -M, -N, and -q")
# parser.add_option("-t", "--test", dest="test",
                  # action="store_true", default=False,
                  # help="test with chromosomes 21 and 22 only")
#parser.add_option("-v", "--verbose", dest="verbose",
#                  action="store_true", default=TRUE,
#                  help="verbose mode (print all qsub messages)")
(options, args) = parser.parse_args()

if len(args) != 2:
    parser.error("incorrect number of arguments")

config = args[0]
email = options.email
overwrite = options.overwrite
qname = options.qname
verbose = True
qsubOptions = options.qsubOptions
# test = options.test

chromosomeString = args[1]

pipeline = os.path.dirname(os.path.abspath(sys.argv[0]))

# can get chromosomes here, later.
configdict = qcp.readConfig(config)

jobid = dict()

# check if directory exists, and if not, created it
if not os.path.exists(configdict["out_gds_dir"]):
    os.makedirs(configdict["out_gds_dir"])

# check if directory is empty -- not the best way but.
if (os.listdir(configdict["out_gds_dir"]) != []) and not overwrite:
    sys.exit("Directory " + configdict["out_gds_dir"] + " is not empty.")


# check if any files exist - chromosomes 1-23
for chromosome in qcp.parseChromosomes([chromosomeString]):
    
    filename = configdict["out_gds_dir"] + "/" + configdict["out_gds_prefix"] + "_chr-%s.gds" % chromosome

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
# borrowed from below - not sure it's right (SN)

# array job for imputed chromosomes
jobid[job] = qcp.submitJob(job, driver, [rscript, config, args, chromosomeString], queue=qname, email=email, qsubOptions=qsubOptions,
                                    verbose=verbose, arrayRange=chromosomeString)

## OTHER AND FAILED commented out for now
# print "Submitting job for other and failed SNPs"
# 
# # submit "other" snps (non-imputed)
# holdid = [jobid["gds_create_imputed"].split(".")[0]] # because it's from an array job
# args = " ".join(["other", testStr])
# job = "gds_create_other_chr-other"
# rscript = os.path.join(pipeline, "R", "gds_create_other.R")
# jobid[job] = qcp.submitJob(job, driver, [rscript, config, args], queue=qname, email=email, qsubOptions=qsubOptions,
#                                     holdid=holdid, verbose=verbose)
# 
# # submit "failed" snps
# holdid = [jobid["gds_create_other_chr-other"]] # hold on "other" snps since we need the snp annotation
# args = " ".join(["failed", testStr])
# job = "gds_create_other_chr-failed"
# rscript = os.path.join(pipeline, "R", "gds_create_other.R")
# jobid[job] = qcp.submitJob(job, driver, [rscript, config, args], queue=qname, email=email, qsubOptions=qsubOptions,
#                                     holdid=holdid, verbose=verbose)

holdid = [jobid["gds_create_imputed"].split(".")[0]]
job = "gds_imputed_cleanup"
rscript = os.path.join(pipeline, "R", "gds_imputed_cleanup.R")

jobid[job] = qcp.submitJob(job, driver, [rscript, config, args], queue=qname, email=email, qsubOptions=qsubOptions,
                                    holdid=holdid, verbose=verbose)
