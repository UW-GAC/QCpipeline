"""Create GDS file from NetCDF file"""

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

if email is not None:
    emailStr = "-m e -M " + email
else:
    emailStr = ""

nameStr = "gds_geno"
driver = os.path.join(pipeline, "runRscript.sh")
rscript = os.path.join(pipeline, nameStr + ".R")
qsub = "qsub %s -N %s %s %s %s" % (emailStr, nameStr, driver, rscript, config)
retcode = subprocess.call(qsub, shell=True)
if retcode != 0:
    sys.exit("job submission failed for command\n" + qsub)
