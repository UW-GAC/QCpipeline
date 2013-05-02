#! /usr/local/bin/python2.7

"""Identity By Descent"""

import sys
import os
import subprocess
from optparse import OptionParser

usage = """%prog [options] config

Identity by Descent with the following steps:
1) Select SNPs with LD pruning  (unless file already exists)
   (from autosomal, non-monomorphic, missing<0.05)
2) IBD calculations
3) Assign observed relationships and plot results
4) Calculate individual inbreeding coefficients

Required config parameters:
annot_scan_file    scan annotation file
annot_snp_file     snp annotation file
gds_geno_file      genotype GDS file
out_ibd_kc32_file  output file with data frame of pairs with KC > 1/32

Optional config parameters [default]:
annot_scan_subjectCol      [subjectID]               column of subjectID in scan annotation
annot_snp_missingCol       [missing.n1]              column of missing call rate in snp annotation
exp_rel_file               [NA]                      file with data frame of expected relationships
ibd_method                 [MoM]                     IBD method (MoM, MLE, or KING)
ld_r_threshold             [0.32]                    r threshold for LD pruning (0.32 = sqrt(0.1))
ld_win_size                [10]                      size of sliding window for LD pruning (in Mb)
maf_threshold              [0.05]                    minimum MAF for non-monomorphic SNPs to consider in LD pruning
scan_ibd_include_file      [NA]                      vector of scanID to include in IBD (NA=all)
scan_pruning_include_file  [NA]                      vector of scanID to include in LD pruning (NA=all)
out_ibd_con_file           [ibd_connectivity.RData]  output connectivity data file
out_ibd_con_plot           [ibd_connectivity.pdf]    output connectivity plot
out_ibd_exp_plot           [ibd_expected.pdf]        output IBD plot color-coded by expected relationships
out_ibd_file               [ibd.RData]               output file of full IBD results
out_ibd_obs_plot           [ibd_observed.pdf]        output IBD plot color-coded by observed relationships
out_ibd_rel_file           [ibd_obsrel.RData]        output file of observed relationships
out_ibd_unexp_plot         [ibd_unexpected.pdf]      output IBD plot color-coded by expected relationships, and different symbols for unexpected with KC > 0.1
out_ibd_unobs_dup_file     [ibd_unobs_dup.RData]     output file of expected but unobserved duplicates
out_ibd_unobs_rel_file     [ibd_unobs_rel.RData]     output file of expected but unobserved relationships
out_inbrd_file             [inbreed_coeff.RData]     output file of inbreeding coefficients
out_snp_file               [ibd_snp_sel.RData]       output file of pruned snps  (if file exists, pruning step is skipped)"""
parser = OptionParser(usage=usage)
parser.add_option("-p", "--pipeline", dest="pipeline",
                  default="/projects/geneva/geneva_sata/GCC_code/QCpipeline",
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

jobid = dict()

# skip LD if file already exists
waitLD = False
if os.path.exists(configdict['out_snp_file']):
    print "using SNPs in " + configdict['out_snp_file']
else:
    job = "ibd_snp_sel"
    rscript = os.path.join(pipeline, "R", job + ".R")
    jobid[job] = QCpipeline.submitJob(job, driver, [rscript, config], queue=qname, email=email)
    waitLD = True

job = "ibd"
if waitLD:
    holdid = [jobid['ibd_snp_sel']]
else:
    holdid = None
rscript = os.path.join(pipeline, "R", job + ".R")
jobid[job] = QCpipeline.submitJob(job, driver, [rscript, config], holdid=holdid, queue=qname, email=email)

job = "ibd_plots"
rscript = os.path.join(pipeline, "R", job + ".R")
jobid[job] = QCpipeline.submitJob(job, driver, [rscript, config], holdid=[jobid['ibd']], queue=qname, email=email)

job = "inbreed_coeff"
if waitLD:
    holdid = [jobid['ibd_snp_sel']]
else:
    holdid = None
rscript = os.path.join(pipeline, "R", job + ".R")
jobid[job] = QCpipeline.submitJob(job, driver, [rscript, config], holdid=holdid, queue=qname, email=email)
