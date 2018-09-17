#! /usr/local/bin/python2.7

"""Chromosome anomalies"""

import QCpipeline
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
annot_scan_file   scan annotation file
annot_snp_file    snp annotation file
build             genome build (hg18 or hg19)
bl_file           BAF/LRR file (netCDF or GDS)
geno_file         genotype file (netCDF or GDS)
project           project name (prepend to output files)
out_anom_dir      output directory for anomaly data
out_plot_dir      output directory for anomaly plots
out_baf_med_file  output (or existing) median BAF SD file

Required for "--maf" option:
out_afreq_file    output (or existing) allele frequency file

Optional config parameters [default]:
annot_snp_IntensityOnlyCol  [NA]                            column of intensity-only in snp annotation
annot_snp_missingCol        [missing.n1]                    column of missing call rate in snp annotation
chromXY                     [FALSE]                         find anomalies in pseudoautosomal region?  
plot.win                    [1]                             size of plot window (multiple of anomaly length)
scan_exclude_file           [NA]                            vector of scanID to exclude
snp_exclude_file            [NA]                            vector of snpID to exclude
thresh.indiv                [5]                             threshold for plotting individual anomalies (Mb)
thresh.sum                  [10]                            threshold for sum of all anomalies on a chromosome (Mb)  
out_baf_mean_file           [baf_mean_by_scan_chrom.RData]  output BAF mean file
out_baf_sd_file             [baf_sd_by_scan_chrom.RData]    output BAF SD file
out_eligible_snps           [snps_eligible.RData]           output file of eligible SNPs
out_plot_prefix             [long_plot]                     output prefix for anomaly plots
"""
parser = OptionParser(usage=usage)
parser.add_option("-e", "--email", dest="email", default=None,
                  help="email address for job reporting")
parser.add_option("-q", "--queue", dest="qname", default="all.q", 
                  help="cluster queue name [default %default]")
parser.add_option("-m", "--maf", dest="maf", default=None,
                  help="do not use SNPs with MAF below this threshold [default %default]")
parser.add_option("--baf", dest="baf",
                  action="store_true", default=False,
                  help="Run BAF detection")
parser.add_option("--loh", dest="loh",
                  action="store_true", default=False,
                  help="Run LOH detection")
parser.add_option("--stats", dest="stats",
                  action="store_true", default=False,
                  help="Run anom stats")
parser.add_option("-o", "--options", dest="qsubOptions", default="",
                  help="additional options to pass to qsub, excluding -hold_jid, -N, -m e -M, -N, and -q")
(options, args) = parser.parse_args()

email = options.email
maf = options.maf
qname = options.qname
runbaf = options.baf
runloh = options.loh
runstats = options.stats
qsubOptions = options.qsubOptions

if len(args) < 1:
    parser.error("incorrect number of arguments")
elif len(args) != 4 and (runbaf or runloh):
    parser.error("incorrect number of arguments")

config = args[0]
if len(args) > 1:
    start = int(args[1])
    end = int(args[2])
    by = int(args[3])

pipeline = os.path.dirname(os.path.abspath(sys.argv[0]))

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
    jobid[job] = QCpipeline.submitJob(job, driver, [rscript, config], queue=qname, email=email, qsubOptions=qsubOptions)
    waitSD = True

# if MAF threshold is needed, calculate allele freq
useMAF = False
waitAfreq = False
if maf is not None and float(maf) > 0:
    useMAF = True
    if os.path.exists(configdict['out_afreq_file']):
        print "using allele freq file " + configdict['out_afreq_file']
    else:
        job = "allele_freq"
        rscript = os.path.join(pipeline, "R", job + ".R")
        jobid[job] = QCpipeline.submitJob(job, driver, [rscript, config], queue=qname, email=email, qsubOptions=qsubOptions)
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

        jobid[job].append(QCpipeline.submitJob(job+str(istart), driver, rargs, holdid=holdid, queue=qname, email=email, qsubOptions=qsubOptions))

        istart = istart + by
        iend = istart + by - 1
        if (iend > end):
            iend = end

    job = "anom_combine"
    rscript = os.path.join(pipeline, "R", job + ".R")
    rargs = [rscript, config, "BAF", str(end), str(by)]
    jobid[job + "_baf"] = QCpipeline.submitJob(job+"_baf", driver, rargs, holdid=jobid['anom_baf'], queue=qname, email=email, qsubOptions=qsubOptions)


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
        rargs = [rscript, config, str(istart), str(iend)]
        jobid[job].append(QCpipeline.submitJob(job+str(istart), driver, rargs, holdid=holdid, queue=qname, email=email, qsubOptions=qsubOptions))

        istart = istart + by
        iend = istart + by - 1
        if (iend > end):
            iend = end

    job = "anom_combine"
    rscript = os.path.join(pipeline, "R", job + ".R")
    rargs = [rscript, config, "LOH", str(end), str(by)]
    jobid[job + "_loh"] = QCpipeline.submitJob(job+"_loh", driver, rargs, holdid=jobid['anom_loh'], queue=qname, email=email, qsubOptions=qsubOptions)


if runstats:
    job = "anom_stats"
    rscript = os.path.join(pipeline, "R", job + ".R")

    if runloh:
        holdid = [jobid["anom_combine_loh"]]
    else:
        holdid = None

    rargs = [rscript, config]
    jobid[job] = QCpipeline.submitJob(job, driver, rargs, holdid=holdid, queue=qname, email=email, qsubOptions=qsubOptions)
