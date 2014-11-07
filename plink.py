#! /usr/local/bin/python2.7

"""Create PLINK file from netCDF or GDS"""

import sys
import os
import subprocess
from optparse import OptionParser

usage = """%prog [options] config

Create subject-level PLINK files.  Any samples with subj.plink=FALSE are excluded.
1) Make ped and map files
2) Convert to bed/bim/fam
3) Create tar.gz file with bed/bim/fam

Required config parameters:
annot_scan_file       scan annotation file (with columns subj.plink, family, father, mother)
annot_snp_file        snp annotation file
geno_file             genotype file (netCDF or GDS)
out_log_prefix        output prefix for log files
out_plink_prefix      output prefix for plink files

Optional config parameters [default]:
annot_scan_subjectCol  [subjectID]  column of subjectID in scan annotation
annot_snp_alleleACol   [alleleA]    column of allele A in snp annotation
annot_snp_alleleBCol   [alleleB]    column of allele B in snp annotation
annot_snp_rsIDCol      [rsID]       column of rsID in snp annotation"""
parser = OptionParser(usage=usage)
parser.add_option("-p", "--pipeline", dest="pipeline",
                  default="/projects/geneva/gcc-fs2/GCC_Code/QCpipeline",
                  help="pipeline source directory")
parser.add_option("-e", "--email", dest="email", default=None,
                  help="email address for job reporting")
parser.add_option("-q", "--queue", dest="qname", default="gcc.q", 
                  help="cluster queue name [default %default]")
(options, args) = parser.parse_args()

if len(args) != 1:
    parser.error("incorrect number of arguments")

config = args[0]
pipeline = options.pipeline
email = options.email
qname = options.qname

sys.path.append(pipeline)
import QCpipeline

configdict = QCpipeline.readConfig(config)

driver = os.path.join(pipeline, "runRscript.sh")
plinkfile = configdict["out_plink_prefix"]

jobid = dict()

job = "plink_from_geno"
rscript = os.path.join(pipeline, "R", job + ".R")
jobid[job] = QCpipeline.submitJob(job, driver, [rscript, config], queue=qname, email=email)

job = "plink_bed"
arglist = ["--noweb", "--make-bed", "--file", plinkfile, "--out", plinkfile]
jobid[job] = QCpipeline.submitJob(job, "plink", arglist, options="-b y -j y -cwd",
                                  holdid=[jobid["plink_from_geno"]], queue=qname, email=email)

job = "plink_tar"
(path, file) = os.path.split(plinkfile)
if path == "":
    path = "."
arglist = ["-cvz", "--directory", path, "-f", plinkfile+".tar.gz", file+".bed", file+".bim", file+".fam"]
jobid[job] = QCpipeline.submitJob(job, "tar", arglist, options="-b y -j y -cwd",
                                  holdid=[jobid["plink_bed"]], queue=qname, email=email)
