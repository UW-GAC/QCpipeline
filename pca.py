#! /usr/local/bin/python2.7

"""Principal Component Analysis"""

import sys
import os
import subprocess
from optparse import OptionParser

usage = """%prog [options] config

Principal Component Analysis with the following steps:
1) If using "combined" option:
  a) create combined GDS file from two netCDF files
  b) check duplicate discordance between datasets
2) LD pruning to select SNPs (unless file already exists)
  - uses scans from study_unrelated_file
3) PCA calculations
  a) if combined, scans will be from study_unduplicated_file
     plus unrelated hapmaps from study controls and external
  b) if not combined, scans will be from study_unrelated_file only
4) Plot results

Required config parameters:
annot_scan_file       scan annotation file
annot_scan_raceCol    column of race in scan annotation
annot_snp_file        snp annotation file
build                 genome build (hg18 or hg19)
gds_geno_file         genotype GDS file (filtered subject-level recommended)
study_unrelated_file  vector of scanID for PCA and LD pruning (no hapmaps)
out_corr_file         output file of PC-SNP correlations
out_pca_file          output file of PCA results
out_pruned_file       output file of pruned snps (if file exists, pruning step is skipped)

Required for "combined" option:
ext_annot_scan_file       external dataset scan annotation file
ext_annot_snp_file        external dataset snp annotation file
ext_nc_geno_file          external dataset genotype netCDF file
nc_geno_file              genotype netCDF file (filtered subject-level recommended)
study_unduplicated_file   vector of scanID from study for combined PCA
out_comb_scan_annot_file  output combined scan annotation
out_comb_snp_annot_file   output combined snp annotation
out_comb_gds_geno_file    output combined GDS file
out_disc_file             output duplicate discordance file

Optional config parameters [default]:
annot_scan_ethnCol           [NA]                     column of ethnicity in scan annotation
annot_scan_hapmapCol         [geno.cntl]              column of hapmap (0/1) in scan annotation
annot_scan_subjectCol        [subjectID]              column of subjectID in scan annotation
annot_scan_unrelCol          [unrelated]              column of unrelated (T/F) in scan annotation
annot_snp_rsIDCol            [rsID]                   column of rsID in snp annotation
comb_scan_exclude_file       [NA]                     vector of scanID (study and/or external) to exclude from combined GDS
comb_snp_exclude_file        [NA]                     vector of snpID (in study annotation) to exclude from combined GDS
ext_annot_scan_raceCol       [pop.group]              column of race in external scan annotation
ext_annot_scan_subjectCol    [subjectID]              column of subjectID in external scan annoatation
ext_annot_scan_unrelCol      [unrelated]              column of unrelated (T/F) in external scan annoatation
ext_annot_snp_rsIDCol        [rsID]                   column of rsID in external snp annotation
ld_r_threshold               [0.32]                   r threshold for LD pruning (0.32 = sqrt(0.1))
ld_win_size                  [10]                     size of sliding window for LD pruning (in Mb)
num_evs_to_plot              [12]                     number of eigenvectors for correlation and scree plots
snp_pruning_include_file     [NA]                     vector of snpID to include in LD pruning (NA=all)
out_corr_plot_prefix         [pca_corr]               output prefix for correlation plots (all SNPs)
out_corr_pruned_plot_prefix  [NA]                     output prefix for correlation plots (pruned SNPs only)
out_dens_plot                [pca_dens.pdf]           output plot of EV2 vs EV1 with density sidebars
out_disc_plot                [dup_disc_ext.pdf]       output duplicate discordance plot
out_ev12_plot                [pca_ev12.pdf]           output plot of EV2 vs EV1
out_pairs_plot               [pca_pairs.png]          output pairs plot of EV 1-4
out_scree_plot               [pca_scree.pdf]          output scree plot

Additional parameters:
each value for race should be a parameter with an associated color
optionally, each value for ethnicity can be a parameter with an associated symbol
e.g.,
CEU           blue
YRI           red
Not_Hispanic  1
Hispanic      4"""
parser = OptionParser(usage=usage)
parser.add_option("-p", "--pipeline", dest="pipeline",
                  default="/projects/geneva/geneva_sata/GCC_code/QCpipeline",
                  help="pipeline source directory")
parser.add_option("-e", "--email", dest="email", default=None,
                  help="email address for job reporting")
parser.add_option("-q", "--queue", dest="qname", default="gcc.q", 
                  help="cluster queue name [default %default]")
parser.add_option("-c", "--combined", dest="combined",
                  action="store_true", default=False,
                  help="make combined dataset")
parser.add_option("-m", "--multithread", dest="multithread", default=None,
                  help="number of cores to use; either a number (e.g, 1) or a range of numbers (e.g., 1-4) [default 1 core]")
(options, args) = parser.parse_args()

if len(args) != 1:
    parser.error("incorrect number of arguments")

config = args[0]
pipeline = options.pipeline
email = options.email
combined = options.combined
qname = options.qname
multithread = options.multithread


sys.path.append(pipeline)
import QCpipeline

configdict = QCpipeline.readConfig(config)

driver = os.path.join(pipeline, "runRscript.sh")

jobid = dict()


if multithread is not None:
    optionsMulti = "-pe local " + multithread
else:
    optionsMulti = ""


if combined:
    job = "combine_gds"
    rscript = os.path.join(pipeline, "R", job + ".R")
    jobid[job] = QCpipeline.submitJob(job, driver, [rscript, config], queue=qname, email=email)

    job = "dup_disc_ext"
    rscript = os.path.join(pipeline, "R", job + ".R")
    jobid[job] = QCpipeline.submitJob(job, driver, [rscript, config], holdid=[jobid['combine_gds']], queue=qname, email=email)

# skip LD if file already exists
# multithreading not implemented in SNPRelate code for snpgdsLDpruning, so don't pass optionsMulti
waitLD = False
if os.path.exists(configdict['out_pruned_file']):
    print "using LD pruned file " + configdict['out_pruned_file']
else:
    job = "ld_pruning"
    if combined:
        holdid = [jobid['dup_disc_ext']]
    else:
        holdid = None
    rscript = os.path.join(pipeline, "R", job + ".R")
    jobid[job] = QCpipeline.submitJob(job, driver, [rscript, config], holdid=holdid, queue=qname, email=email)
    waitLD = True

if combined:
    job = "pca_combined"
    rscript = os.path.join(pipeline, "R", job + ".R")
    holdid = [jobid['dup_disc_ext'], jobid['combine_gds']]
    if waitLD:
        holdid.append(jobid['ld_pruning'])
    jobid[job] = QCpipeline.submitJob(job, driver, [rscript, config], holdid=holdid, queue=qname, email=email, options=optionsMulti)

    job = "pca_plots"
    rscript = os.path.join(pipeline, "R", job + ".R")
    jobid[job] = QCpipeline.submitJob(job+"_combined", driver, [rscript, config, "combined"], holdid=[jobid['pca_combined']], queue=qname, email=email)

else:
    job = "pca_study"
    rscript = os.path.join(pipeline, "R", job + ".R")
    if waitLD:
        holdid = [jobid['ld_pruning']]
    else:
        holdid = None
    jobid[job] = QCpipeline.submitJob(job, driver, [rscript, config], holdid=holdid, queue=qname, email=email, options=optionsMulti)

    job = "pca_plots"
    rscript = os.path.join(pipeline, "R", job + ".R")
    jobid[job] = QCpipeline.submitJob(job+"_study", driver, [rscript, config, "study"], holdid=[jobid['pca_study']], queue=qname, email=email)
