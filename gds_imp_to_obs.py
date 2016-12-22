#! /usr/local/bin/python2.7

"""Create GDS files"""

import QCpipeline as qcp
import sys
import os
import subprocess
from optparse import OptionParser


usage = """%prog [options] config

Creates GDS files of observed SNPs from a set of imputed files.
Required config parameters:

Optional config parameters [default]:
"""
parser = OptionParser(usage=usage)
parser.add_option("-e", "--email", dest="email", default=None,
                  help="email address for job reporting")
parser.add_option("-q", "--queue", dest="qname", default="all.q", 
                  help="cluster queue name [default %default]")
parser.add_option("-o", "--overwrite", dest="overwrite",
                  action="store_true", default=False,
                  help="overwrite existing files")
(options, args) = parser.parse_args()

if len(args) != 1:
    parser.error("incorrect number of arguments")

config = args[0]
email = options.email
overwrite = options.overwrite
qname = options.qname
verbose = True


pipeline = os.path.dirname(os.path.abspath(sys.argv[0]))


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





driver = os.path.join(pipeline, "runRscript_array.sh")
    
args = ""
job = "gds_imp_to_obs"
rscript = os.path.join(pipeline, "R", "gds_imp_to_obs.R")

# array job for imputed chromosomes
jobid[job] = qcp.submitJob(job, driver, [rscript, config, args], queue=qname, email=email,
                                    verbose=verbose, arrayRange="1:23")

