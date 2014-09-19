#! /usr/local/bin/python2.7

"""Create filtered and unfiltered PLINK files"""

import sys
import os
import subprocess
from optparse import OptionParser

usage = """%prog [options] config

Create subject-level PLINK files.
1) Make ped and map files
2) Convert to bed/bim/fam
3) Create tar.gz file with bed/bim/fam
An unfiltered file is always created from the sample-level netCDF 
using subj.plink samples.
If the "filtered" option is given, a filtered file is also created
from the subject-level netCDF.

Required config parameters:
annot_scan_file       scan annotation file (with columns subj.plink, family, father, mother)
annot_snp_file        snp annotation file
samp_geno_file        sample-level genotype file (netCDF or GDS)
out_log_prefix        output prefix for log files
out_plink_prefix      output prefix for plink files

Required for "filtered" option:
subj_geno_file        subject-level genotype file (netCDF or GDS)

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
rscript = os.path.join(pipeline, "R", job + ".R")
jobid[job] = QCpipeline.submitJob(job, driver, [rscript, config], queue=qname, email=email)

job = "plink_unfiltered_bed"
arglist = ["--noweb", "--make-bed", "--file", plinkfile, "--out", plinkfile]
jobid[job] = QCpipeline.submitJob(job, "plink", arglist, options="-b y -j y -cwd",
                                  holdid=[jobid["plink_unfiltered"]], queue=qname, email=email)

job = "plink_unfiltered_tar"
(path, file) = os.path.split(plinkfile)
if path == "":
    path = "."
arglist = ["-cvz", "--directory", path, "-f", plinkfile+".tar.gz", file+".bed", file+".bim", file+".fam"]
jobid[job] = QCpipeline.submitJob(job, "tar", arglist, options="-b y -j y -cwd",
                                  holdid=[jobid["plink_unfiltered_bed"]], queue=qname, email=email)

if filt:
    plinkfile = configdict["out_plink_prefix"] + "_filtered"
    job = "plink_filtered"
    rscript = os.path.join(pipeline, "R", job + ".R")
    jobid[job] = QCpipeline.submitJob(job, driver, [rscript, config], queue=qname, email=email)

    job = "plink_filtered_bed"
    arglist = ["--noweb", "--make-bed", "--file", plinkfile, "--out", plinkfile]
    jobid[job] = QCpipeline.submitJob(job, "plink", arglist, options="-b y -j y -cwd",
                                      holdid=[jobid["plink_filtered"]], queue=qname, email=email)

    job = "plink_filtered_tar"
    (path, file) = os.path.split(plinkfile)
    if path == "":
        path = "."
    arglist = ["-cvz", "--directory", path, "-f", plinkfile+".tar.gz", file+".bed", file+".bim", file+".fam"]
    jobid[job] = QCpipeline.submitJob(job, "tar", arglist, options="-b y -j y -cwd",
                                      holdid=[jobid["plink_filtered_bed"]], queue=qname, email=email)
