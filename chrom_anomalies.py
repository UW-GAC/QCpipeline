#! /usr/local/bin/python2.7

"""Chromosome anomalies"""

import sys
import os
import subprocess
from optparse import OptionParser

usage = """%prog [options] config start end by

Detect chromosome anomalies with the following steps:
1) Find BAF variance (median BAF SD)
2) If using MAF threshold, calculate allele frequency
3) Detect BAF anomalies, running samples in parallel with start,end,by
4) Combine BAF anomalies
5) Detect LOH anomalies, running samples in parallel with start,end,by
6) Combine LOH anomalies
7) Calculate anomaly stats and make plots of anomalies > thresh.indiv 
   (or sum on a chromosome > thresh.sum)

Required config parameters:
annot_scan_file  scan annotation file
annot_snp_file   snp annotation file
build            genome build (hg18 or hg19)
nc_bl_file       BAF/LRR netCDF file
nc_geno_file     genotype netCDF file
project          project name (prepend to output files)
out_anom_dir     output directory for anomaly data
out_plot_dir     output directory for anomaly plots

Optional config parameters [default]:
annot_snp_IntensityOnlyCol  [NA]                            column of intensity-only in snp annotation
annot_snp_missingCol        [missing.n1]                    column of missing call rate in snp annotation
chromXY                     [FALSE]                         find anomalies in pseudoautosomal region?  
plot.win                    [1]                             size of plot window (multiple of anomaly length)
scan_exclude_file           [NA]                            vector of scanID to exclude
thresh.indiv                [10]                            threshold for sum of all anomalies on a chromosome (Mb)
thresh.sum                  [5]                             threshold for plotting individual anomalies (Mb)  
out_afreq_file              [allele_freq.RData]             output (or existing) allele frequency file
out_baf_mean_file           [baf_mean_by_scan_chrom.RData]  output BAF mean file
out_baf_med_file            [median_baf_sd_by_scan.RData]   output (or existing) median BAF SD file
out_baf_sd_file             [baf_sd_by_scan_chrom.RData]    output BAF SD file
out_eligible_snps           [snps_eligible.RData]           output file of eligible SNPs
out_plot_prefix             [long_plot]                     output prefix for anomaly plots
"""
parser = OptionParser(usage=usage)
parser.add_option("-p", "--pipeline", dest="pipeline",
                  default="/projects/geneva/geneva_sata/GCC_code/QCpipeline",
                  help="pipeline source directory")
parser.add_option("-e", "--email", dest="email", default=None,
                  help="email address for job reporting")
parser.add_option("-q", "--queue", dest="qname", default="gcc.q", 
                  help="cluster queue name [default %default]")
parser.add_option("-m", "--maf", dest="maf", default=None,
                  help="do not use SNPs with MAF below this threshold")
parser.add_option("--baf", dest="baf",
                  action="store_true", default=False,
                  help="Run BAF detection")
parser.add_option("--loh", dest="loh",
                  action="store_true", default=False,
                  help="Run LOH detection")
parser.add_option("--stats", dest="stats",
                  action="store_true", default=False,
                  help="Run anom stats")
(options, args) = parser.parse_args()

if len(args) != 4:
    parser.error("incorrect number of arguments")

config = args[0]
start = int(args[1])
end = int(args[2])
by = int(args[3])
pipeline = options.pipeline
email = options.email
maf = options.maf
qname = options.qname
runbaf = options.baf
runloh = options.loh
runstats = options.stats

sys.path.append(pipeline)
import QCpipeline

configdict = QCpipeline.readConfig(config)

driver = os.path.join(pipeline, "runRscript.sh")

jobid = dict()

# skip SD if file already exists
waitSD = False
if os.path.exists(configdict['out_baf_med_file']):
    print "using BAF SD file " + configdict['out_baf_med_file']
else:
    job = "baf_variance"
    rscript = os.path.join(pipeline, "R", job + ".R")
    jobid[job] = QCpipeline.submitJob(job, driver, [rscript, config], queue=qname, email=email)
    waitSD = True

# if MAF threshold is needed, calculate allele freq
useMAF = False
waitAfreq = False
if maf is not None:
    useMAF = True
    if os.path.exists(configdict['out_afreq_file']):
        print "using allele freq file " + configdict['out_afreq_file']
    else:
        job = "allele_freq"
        rscript = os.path.join(pipeline, "R", job + ".R")
        jobid[job] = QCpipeline.submitJob(job, driver, [rscript, config], queue=qname, email=email)
        waitAfreq = True

if runbaf:
    job = "anom_baf"
    rscript = os.path.join(pipeline, "R", job + ".R")
    jobid[job] = []

    if waitSD and waitAfreq:
        holdid = [jobid["baf_variance"], jobid["allele_freq"]]
    elif waitSD:
        holdid = [jobid["baf_variance"]]
    elif waitAfreq:
        holdid = [jobid["allele_freq"]]
    else:
        holdid = None

    istart = start
    iend = start + by - 1
    while end >= istart:
        if useMAF:
            rargs = [rscript, config, str(istart), str(iend), maf]
        else:
            rargs = [rscript, config, str(istart), str(iend)]

        jobid[job].append(QCpipeline.submitJob(job+str(istart), driver, rargs, holdid=holdid, queue=qname, email=email))

        istart = istart + by
        iend = istart + by - 1
        if (iend > end):
            iend = end

    job = "anom_combine"
    rscript = os.path.join(pipeline, "R", job + ".R")
    jobid[job + "_baf"] = QCpipeline.submitJob(job+"_baf", driver, [rscript, config, "BAF", str(end), str(by)], holdid=jobid['anom_baf'], queue=qname, email=email)


if runloh:
    job = "anom_loh"
    rscript = os.path.join(pipeline, "R", job + ".R")
    jobid[job] = []
    
    if runbaf:
        holdid = [jobid["anom_combine_baf"]]
    else:
        holdid = None

    istart = start
    iend = start + by - 1
    while end >= istart:
        if useMAF:
            rargs = [rscript, config, str(istart), str(iend), maf]
        else:
            rargs = [rscript, config, str(istart), str(iend)]

        jobid[job].append(QCpipeline.submitJob(job+str(istart), driver, rargs, holdid=holdid, queue=qname, email=email))

        istart = istart + by
        iend = istart + by - 1
        if (iend > end):
            iend = end

    job = "anom_combine"
    rscript = os.path.join(pipeline, "R", job + ".R")
    jobid[job + "_loh"] = QCpipeline.submitJob(job+"_loh", driver, [rscript, config, "LOH", str(end), str(by)], holdid=jobid['anom_loh'], queue=qname, email=email)


if runstats:
    job = "anom_stats"
    rscript = os.path.join(pipeline, "R", job + ".R")

    if runloh:
        holdid = [jobid["anom_combine_loh"]]
    else:
        holdid = None

    if useMAF:
        rargs = [rscript, config, maf]
    else:
        rargs = [rscript, config]
    jobid[job] = QCpipeline.submitJob(job, driver, rargs, holdid=holdid, queue=qname, email=email)
