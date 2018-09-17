#! /usr/local/bin/python2.7

"""Identity By Descent"""

import QCpipeline
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
4) Calculate individual inbreeding coefficients (if method is MoM)

Required config parameters:
annot_scan_file    scan annotation file
gds_geno_file      genotype GDS file
out_snp_file       output file of pruned snps  (if file exists, pruning step is skipped)

Optional config parameters [default]:
annot_scan_familyCol       [NA]                      column of familyID in scan annotation; if not NA, used when ibd_method="KING" to calculate kinship coefficients within families differently
annot_scan_subjectCol      [subjectID]               column of subjectID in scan annotation
exp_rel_file               [NA]                      file with data frame of expected relationships
ibd_method                 [KING]                    IBD method (MoM, MLE, or KING)
ld_r_threshold             [0.32]                    r threshold for LD pruning (0.32 = sqrt(0.1))
ld_win_size                [10]                      size of sliding window for LD pruning (in Mb)
maf_threshold              [0.05]                    minimum MAF for non-monomorphic SNPs to consider in LD pruning
scan_ibd_include_file      [NA]                      vector of scanID to include in IBD (NA=all)
scan_pruning_include_file  [NA]                      vector of scanID to include in LD pruning (NA=all)
snp_pruning_include_file   [NA]                      vector of snpID to include in LD pruning (NA=all)
unexpected_threshold       [deg2]                    threshold for defining relationships as unexpected (deg2 or deg3)
out_ibd_con_file           [ibd_connectivity.RData]  output connectivity data file
out_ibd_con_plot           [ibd_connectivity.pdf]    output connectivity plot
out_ibd_exp_plot           [ibd_expected.pdf]        output IBD plot color-coded by expected relationships
out_ibd_file               [ibd.RData]               output file of full IBD results
out_ibd_kc32_file          [ibd_kc32.RData]          output file with data frame of pairs with KC > 1/32
out_ibd_obs_plot           [ibd_observed.pdf]        output IBD plot color-coded by observed relationships
out_ibd_rel_file           [ibd_obsrel.RData]        output file of observed relationships
out_ibd_unexp_plot         [ibd_unexpected.pdf]      output IBD plot color-coded by expected relationships, and different symbols for unexpected with KC > 0.1
out_ibd_unobs_dup_file     [ibd_unobs_dup.RData]     output file of expected but unobserved duplicates
out_ibd_unobs_rel_file     [ibd_unobs_rel.RData]     output file of expected but unobserved relationships
out_inbrd_file             [inbreed_coeff.RData]     output file of inbreeding coefficients"""
parser = OptionParser(usage=usage)
parser.add_option("-e", "--email", dest="email", default=None,
                  help="email address for job reporting")
parser.add_option("-q", "--queue", dest="qname", default="gcc.q", 
                  help="cluster queue name [default %default]")
parser.add_option("-m", "--multithread", dest="multithread", default=None,
                  help="number of cores to use; either a number (e.g, 1) or a range of numbers (e.g., 1-4) [default 1 core]")
parser.add_option("-i", "--inbreed", dest="inbreed",
                  action="store_true", default=False,
                  help="calculate individual inbreeding coefficients (not for use with KING)")
parser.add_option("-o", "--options", dest="qsubOptions", default="",
                  help="additional options to pass to qsub, excluding -hold_jid, -N, -m e -M, -N, and -q")
(options, args) = parser.parse_args()

if len(args) != 1:
    parser.error("incorrect number of arguments")

config = args[0]
email = options.email
qname = options.qname
multithread = options.multithread
inbreed = options.inbreed
qsubOptions = options.qsubOptions

pipeline = os.path.dirname(os.path.abspath(sys.argv[0]))

configdict = QCpipeline.readConfig(config)

driver = os.path.join(pipeline, "runRscript.sh")

jobid = dict()


if multithread is not None:
    optionsMulti = qsubOptions + " -pe local " + multithread
else:
    optionsMulti = qsubOptions


# skip LD if file already exists
# multithreading not implemented in SNPRelate code for snpgdsLDpruning, so don't pass optionsMulti
waitLD = False
if os.path.exists(configdict['out_snp_file']):
    print "using SNPs in " + configdict['out_snp_file']
else:
    job = "ibd_snp_sel"
    rscript = os.path.join(pipeline, "R", job + ".R")
    jobid[job] = QCpipeline.submitJob(job, driver, [rscript, config], queue=qname, email=email, qsubOptions=qsubOptions)
    waitLD = True

job = "ibd"
if waitLD:
    holdid = [jobid['ibd_snp_sel']]
else:
    holdid = None
rscript = os.path.join(pipeline, "R", job + ".R")
jobid[job] = QCpipeline.submitJob(job, driver, [rscript, config], holdid=holdid, queue=qname, email=email, qsubOptions=optionsMulti)

job = "ibd_plots"
rscript = os.path.join(pipeline, "R", job + ".R")
jobid[job] = QCpipeline.submitJob(job, driver, [rscript, config], holdid=[jobid['ibd']], queue=qname, email=email, qsubOptions=qsubOptions)

# don't run inbreed_coeff for KING.
if inbreed:
    job = "inbreed_coeff"
    if waitLD:
        holdid = [jobid['ibd_snp_sel']]
    else:
        holdid = None
    rscript = os.path.join(pipeline, "R", job + ".R")
    jobid[job] = QCpipeline.submitJob(job, driver, [rscript, config], holdid=holdid, queue=qname, email=email, qsubOptions=qsubOptions)
