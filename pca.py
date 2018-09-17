#! /usr/local/bin/python2.7

"""Principal Component Analysis"""

import QCpipeline
import sys
import os
import subprocess
from optparse import OptionParser

usage = """%prog [options] config

Principal Component Analysis with the following steps:
1) If using "combined" option:
  a) create combined GDS file from two netCDF/GDS files
    (default is to remove discordant SNPs)
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
out_pruned_file       output file of pruned snps (if file exists, pruning step is skipped)

Required for "combined" option:
ext_annot_scan_file       external dataset scan annotation file
ext_annot_snp_file        external dataset snp annotation file
ext_geno_file             external dataset genotype netCDF or GDS file
study_unduplicated_file   vector of scanID from study for combined PCA
out_comb_prefix           prefix for output combined GDS file and annotation

Optional config parameters [default]:
annot_scan_ethnCol           [NA]                   column of ethnicity in scan annotation
annot_scan_hapmapCol         [geno.cntl]            column of hapmap (0/1) in scan annotation
annot_scan_subjectCol        [subjectID]            column of subjectID in scan annotation
annot_scan_unrelCol          [unrelated]            column of unrelated (T/F) in scan annotation
annot_snp_alleleACol         [alleleA]              column of allele A in snp annotation
annot_snp_alleleBCol         [alleleB]              column of allele B in snp annotation
annot_snp_rsIDCol            [rsID]                 column of rsID in snp annotation
comb_scan_exclude_file       [NA]                   vector of scanID (study and/or external) to exclude from combined GDS
comb_snp_exclude_file        [NA]                   vector of snpID (in study annotation) to exclude from combined GDS
ext_annot_scan_raceCol       [pop.group]            column of race in external scan annotation
ext_annot_scan_subjectCol    [subjectID]            column of subjectID in external scan annotation
ext_annot_scan_unrelCol      [unrelated]            column of unrelated (T/F) in external scan annotation
ext_annot_snp_alleleACol     [alleleA]              column of allele A in external snp annotation
ext_annot_snp_alleleBCol     [alleleB]              column of allele B in external snp annotation
ext_annot_snp_rsIDCol        [rsID]                 column of rsID in external snp annotation
include_study_hapmaps        [TRUE]                 logical for whether to include study HapMaps in combined PCA
ld_r_threshold               [0.32]                 r threshold for LD pruning (0.32 = sqrt(0.1))
ld_win_size                  [10]                   size of sliding window for LD pruning (in Mb)
num_evs_to_plot              [12]                   number of eigenvectors for correlation and scree plots
remove_discordant            [TRUE]                 logical for whether to remove discordant SNPs from combined file
snp_pruning_include_file     [NA]                   vector of snpID to include in LD pruning (NA=all)
out_pca_file                 [pca.RData]            output file of PCA results
out_corr_file                [pca_corr.RData]       output file of PC-SNP correlations
out_corr_plot_prefix         [pca_corr]             output prefix for correlation plots (all SNPs)
out_corr_pruned_plot_prefix  [NA]                   output prefix for correlation plots (pruned SNPs only)
out_dens_plot                [pca_dens.pdf]         output plot of EV2 vs EV1 with density sidebars
out_ev12_plot                [pca_ev12.pdf]         output plot of EV2 vs EV1
out_ev12_plot_hapmap         [pca_ev12_hapmap.pdf]  output plot of EV2 vs EV1 for hapmaps only (combined only)
out_ev12_plot_study          [pca_ev12_study.pdf]   output plot of EV2 vs EV1 for study subjects only (combined only)
out_pairs_plot               [pca_pairs.png]        output pairs plot of EV 1-4
out_scree_plot               [pca_scree.pdf]        output scree plot
out_parcoord_plot            [pca_parcoord.png]     output PNG parallel coordinates plot of EV 1-8
parcoord_vars                [NA]                   scan annotation variables for parallel coordinates plots (study run only)
out_parcoord_prefix          [pca_parcoord]         prefix for PNG parcoord plots of parcoord_vars (study run only)

Additional parameters:
each value for race should be a parameter with an associated color
optionally, each value for ethnicity can be a parameter with an associated symbol
e.g.,
CEU           blue
YRI           red
Not_Hispanic  1
Hispanic      4"""
parser = OptionParser(usage=usage)
parser.add_option("-e", "--email", dest="email", default=None,
                  help="email address for job reporting")
parser.add_option("-q", "--queue", dest="qname", default="gcc.q", 
                  help="cluster queue name [default %default]")
parser.add_option("-c", "--combined", dest="combined",
                  action="store_true", default=False,
                  help="make combined dataset")
parser.add_option("-m", "--multithread", dest="multithread", default=None,
                  help="number of cores to use; either a number (e.g, 1) or a range of numbers (e.g., 1-4) [default 1 core]")
parser.add_option("-o", "--options", dest="qsubOptions", default="",
                  help="additional options to pass to qsub, excluding -hold_jid, -N, -m e -M, -N, and -q")
(options, args) = parser.parse_args()

if len(args) != 1:
    parser.error("incorrect number of arguments")

config = args[0]
email = options.email
combined = options.combined
qname = options.qname
multithread = options.multithread
qsubOptions = options.qsubOptions


pipeline = os.path.dirname(os.path.abspath(sys.argv[0]))

configdict = QCpipeline.readConfig(config)

driver = os.path.join(pipeline, "runRscript.sh")

jobid = dict()


if multithread is not None:
    optionsMulti = qsubOptions + " -pe local " + multithread
else:
    optionsMulti = qsubOptions

if combined:
    type = "combined"
else:
    type = "study"

if combined:
    job = "combine_gds"
    rscript = os.path.join(pipeline, "R", job + ".R")
    jobid[job] = QCpipeline.submitJob(job, driver, [rscript, config], queue=qname, email=email, qsubOptions=qsubOptions)

# skip LD if file already exists
# multithreading not implemented in SNPRelate code for snpgdsLDpruning, so don't pass optionsMulti
waitLD = False
if os.path.exists(configdict['out_pruned_file']):
    print "using LD pruned file " + configdict['out_pruned_file']
else:
    job = "ld_pruning"
    if combined:
        holdid = [jobid['combine_gds']]
    else:
        holdid = None
    rscript = os.path.join(pipeline, "R", job + ".R")
    jobid[job] = QCpipeline.submitJob(job, driver, [rscript, config, type], holdid=holdid, queue=qname, email=email, qsubOptions=qsubOptions)
    waitLD = True

if combined:
    job = "pca_combined"
    rscript = os.path.join(pipeline, "R", job + ".R")
    holdid = [jobid['combine_gds']]
    if waitLD:
        holdid.append(jobid['ld_pruning'])
    jobid['pca'] = QCpipeline.submitJob(job, driver, [rscript, config], holdid=holdid, queue=qname, email=email, qsubOptions=optionsMulti)

else:
    job = "pca_study"
    rscript = os.path.join(pipeline, "R", job + ".R")
    if waitLD:
        holdid = [jobid['ld_pruning']]
    else:
        holdid = None
    jobid['pca'] = QCpipeline.submitJob(job, driver, [rscript, config], holdid=holdid, queue=qname, email=email, qsubOptions=optionsMulti)

job = "pca_plots"
rscript = os.path.join(pipeline, "R", job + ".R")
jobid[job] = QCpipeline.submitJob(job+"_"+type, driver, [rscript, config, type], holdid=[jobid['pca']], queue=qname, email=email, qsubOptions=qsubOptions)
