#! /usr/local/bin/python2.7

"""SNP filters: allele frequency, duplicate discordance, Mendelian errors"""

import sys
import os
import subprocess
from optparse import OptionParser

usage = """%prog [options] config

SNP filters:
1) Allele frequency
2) Duplicate sample discordance
  a) discordance over all alleles
  b) option "MAconc": discordance excluding major homozygotes
3) Mendelian errors
  a) option "mendClust" makes cluster plots
4) Option "dupSNP": Duplicate SNP discordance
5) Plot summary and results

Required config parameters:
annot_scan_file  scan annotation file
annot_snp_file   snp annotation file
geno_file        sample-level genotype file (netCDF or GDS)
geno_subj_file   subject-level genotype file (netCDF or GDS)

Required for "mendClust" option:
xy_file         XY intensity file (netCDF or GDS)
mend.bin.start  start for error bins, e.g., "0 1 6"
mend.bin.end    end for error bins, e.g., "1 5 10"

Optional config parameters [default]:
annot_scan_hapmapCol      [geno.cntl]           column of hapmap (0/1) in scan annotation
annot_scan_subjectCol     [subjectID]           column of subjectID in scan annotation
annot_snp_dupSnpCol       [dup.pos.id]          column of dup snp id in snp annotation
annot_snp_missingCol      [missing.n1]          column of missing call rate in snp annotation
annot_snp_rsIDCol         [rsID]                column of rsID in snp annotation
corr.by.snp               [FALSE]               compute correlation by SNP? (slow)
maf.bin                   [0.01]                bin for plotting results by minor allele frequency
disc_scan_exclude_file    [NA]                  vector of scanID to exclude from dup sample discord
dupsnp_scan_exclude_file  [NA]                  vector of scanID to exclude from dup snp discord
mend_scan_exclude_file    [NA]                  vector of scanID to exclude from mendelian errors
scan_exclude_file         [NA]                  vector of scanID to exclude from allele frequency
out_afreq_file            [allele_freq.RData]   output file for allele frequency
out_disc_file             [dup_disc.RData]      output file for duplicate sample discordance
out_disc_maf_file         [dup_disc_maf.RData]  output file for minor allele discordance
out_disc_plot             [dup_disc.pdf]        output plot of duplicate sample discordance
out_dupsnp_file           [dup_snps.RData]      output file for duplicate snp discordance
out_ma_conc_plot          [snp_ma_conc.pdf]     output plot of minor allele concordance
out_maf_autosomes_plot    [maf_aut_hist.pdf]    output histogram of MAF for autosomes
out_maf_plot              [maf_hist.pdf]        output hisotogram of MAF for all snps
out_maf_xchrom_plot       [maf_x_hist.pdf]      output histogram of MAF for X chrom snps
out_mend_clust_prefix     [mendel_clust]        output prefix for mendelian error cluster plots
out_mend_file             [mendel_err.RData]    output file for mendelian errors
out_snp_conc_plot         [snp_conc.pdf]        output plot of snp concordance
out_snp_corr_plot         [snp_corr.pdf]        output plot of snp correlation"""
parser = OptionParser(usage=usage)
parser.add_option("-p", "--pipeline", dest="pipeline",
                  default="/projects/geneva/gcc-fs2/GCC_Code/QCpipeline",
                  help="pipeline source directory")
parser.add_option("-e", "--email", dest="email", default=None,
                  help="email address for job reporting")
parser.add_option("-q", "--queue", dest="qname", default="gcc.q", 
                  help="cluster queue name [default %default]")
parser.add_option("--MAconc", dest="mac",
                  action="store_true", default=False,
                  help="minor allele concordance")
parser.add_option("--dupSNP", dest="dupsnp",
                  action="store_true", default=False,
                  help="duplicate SNP discordance")
parser.add_option("--mendClust", dest="mendclust",
                  action="store_true", default=False,
                  help="make cluster plots for Mendelian errors")
(options, args) = parser.parse_args()

if len(args) != 1:
    parser.error("incorrect number of arguments")

config = args[0]
pipeline = options.pipeline
email = options.email
qname = options.qname
mac = options.mac
dupsnp = options.dupsnp
mendclust = options.mendclust

sys.path.append(pipeline)
import QCpipeline

driver = os.path.join(pipeline, "runRscript.sh")

jobid = dict()
job = "allele_freq"
rscript = os.path.join(pipeline, "R", job + ".R")
jobid[job] = QCpipeline.submitJob(job, driver, [rscript, config], queue=qname, email=email)

job = "dup_disc"
rscript = os.path.join(pipeline, "R", job + ".R")
jobid[job] = QCpipeline.submitJob(job, driver, [rscript, config], holdid=[jobid["allele_freq"]], queue=qname, email=email)

if mac:
    job = "dup_disc_maf"
    rscript = os.path.join(pipeline, "R", job + ".R")
    jobid[job] = QCpipeline.submitJob(job, driver, [rscript, config], holdid=[jobid["allele_freq"]], queue=qname, email=email)

job = "snp_filt_plots"
rscript = os.path.join(pipeline, "R", job + ".R")
holdid = [jobid["dup_disc"]]
if mac:
    holdid.append(jobid["dup_disc_maf"])
jobid[job] = QCpipeline.submitJob(job, driver, [rscript, config], holdid=holdid, queue=qname, email=email)

job = "mendel_err"
rscript = os.path.join(pipeline, "R", job + ".R")
jobid[job] = QCpipeline.submitJob(job, driver, [rscript, config], queue=qname, email=email)

if mendclust:
    job = "mendel_plots"
    rscript = os.path.join(pipeline, "R", job + ".R")
    jobid[job] = QCpipeline.submitJob(job, driver, [rscript, config], holdid=[jobid["mendel_err"]], queue=qname, email=email)

if dupsnp:
    job = "dup_snps"
    rscript = os.path.join(pipeline, "R", job + ".R")
    jobid[job] = QCpipeline.submitJob(job, driver, [rscript, config], queue=qname, email=email)
