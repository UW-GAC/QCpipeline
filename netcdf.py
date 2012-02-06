#! /usr/local/bin/python2.7

"""Create NetCDF and GDS files"""

import sys
import os
from optparse import OptionParser

usage = """python %prog [options] config"""
parser = OptionParser(usage=usage)
parser.add_option("-p", "--pipeline", dest="pipeline",
                  default="/projects/geneva/geneva_sata/GCC_code/QCpipeline",
                  help="pipeline source directory")
parser.add_option("-e", "--email", dest="email", default=None,
                  help="email address for job reporting")
parser.add_option("-t", "--test", dest="test",
                  action="store_true", default=False,
                  help="test with first 5 scans only")
parser.add_option("-o", "--overwrite", dest="overwrite",
                  action="store_true", default=False,
                  help="overwrite existing files")
parser.add_option("-q", "--queue", dest="qname",
                  default="gcc.q", help="cluster queue name")
(options, args) = parser.parse_args()

if len(args) != 1:
    parser.error("incorrect number of arguments")

config = args[0]
pipeline = options.pipeline
email = options.email
test = options.test
overwrite = options.overwrite
qname = options.qname

sys.path.append(pipeline)
import QCpipeline

if test:
    testStr = "test"
else:
    testStr = ""

if not overwrite:
    configdict = QCpipeline.readConfig(config)
    for file in (configdict['nc_geno_file'], configdict['nc_xy_file'],
                 configdict['nc_bl_file'], configdict['gds_geno_file']):
        if os.path.exists(file):
            sys.exit(file + "already exists; use -o flag to overwrite")

driver = os.path.join(pipeline, "runRscript.sh")

jobid = dict()
for job in ["ncdf_geno", "ncdf_xy", "ncdf_bl"]:
    rscript = os.path.join(pipeline, job + ".R")
    jobid[job] = QCpipeline.submitJob(job, driver, [rscript, config, testStr], queue=qname, email=email)
    
if not test:
    job = "gds_geno"
    rscript = os.path.join(pipeline, job + ".R")
    jobid[job] = QCpipeline.submitJob(job, driver, [rscript, config], holdid=[jobid['ncdf_geno']], queue=qname, email=email)
    
