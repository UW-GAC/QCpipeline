#! /usr/local/bin/python2.7

"""Mixed model association tests"""

import QCpipeline as qcp
import sys
import os
import subprocess
import argparse


description = "Run association tests using the mixed model."
parser = argparse.ArgumentParser(description=description, formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument("config", help="config file", nargs=1)
parser.add_argument("chromosomes", help="chromosomes to run (e.g. 1:23 or 1 3 5)", nargs="*") # nargs="*" means match 0 to any number
parser.add_argument("-e", "--email", dest="email", default=None,
                    help="email address for job reporting")
parser.add_argument("-q", "--queue", dest="queue", default="all.q", 
                    help="cluster queue name")
parser.add_argument("-v", "--varcomp", dest="varcomp", 
                    action="store_true", default=False,
                    help="estimate variance components (null model)")
parser.add_argument("-a", "--assoc", dest="assoc",
                    action="store_true", default=False,
                    help="run association tests for chroms chromStart - chromEnd")
parser.add_argument("-o", "--options", dest="qsubOptions", default="",
                    help="additional options to pass to qsub, excluding -hold_jid, -N, -m e -M, -N, and -q")
parser.add_argument("--plotQQManh", dest="plotqq",
                    action="store_true", default=False,
                    help="QQ and Manhattan plots")

arguments = parser.parse_args()

# check assoc/chromosomes
if arguments.assoc and len(arguments.chromosomes) < 1:
    print "chromosomes must be given if --assoc is specified."
    parser.print_usage()

pipeline = os.path.dirname(os.path.abspath(sys.argv[0]))

# parse chromosomes
chromosomes = qcp.parseChromosomes(arguments.chromosomes)

config_file = arguments.config[0]
configdict = qcp.readConfig(config_file)


# other checks:
# snp segment file
if arguments.assoc:
    
    try:
        segmap = qcp.getSegmentMapping(configdict["snp_segment_file"])
    except KeyError:
        print "'snp_segment_file' required in config file"
        sys.exit(0)
    except IOError:
        print "snp_segment_file " + configdict["snp_segment_file"] + " not found"
        sys.exit(0)
    except ValueError:
        print "'snp_segment_file format is incorrect"
        sys.exit(0)
    except:
        raise


# set up bash scripts to run R scripts on the cluster
driver = os.path.join(pipeline, "runRscript.sh")
driver_array = os.path.join(pipeline, "runRscript_array.sh")

jobid = dict()

qsub_kwargs = {"email": arguments.email,
              "queue": arguments.queue}
qsub_kwargs_no_email = dict((k, qsub_kwargs[k]) for k in qsub_kwargs.keys() if k not in ("email"))

holdid = None

if arguments.varcomp:
    job = "assoc_mixed_varcomp"
    rscript = os.path.join(pipeline, "R", job + ".R")
    
    args = config_file
    
    jobid[job] = qcp.submitJob(job, driver, [rscript, args], **qsub_kwargs)
    
    holdid = jobid[job]
    
    
if arguments.assoc:

    assoc_script = "assoc_mixed_segment"
    assoc_rscript = os.path.join(pipeline, "R", assoc_script + ".R")
    
    comb_script = "assoc_mixed_combine"
    comb_rscript = os.path.join(pipeline, "R", comb_script + ".R")
    
    print "submitting jobs for chromosomes:", ", ".join(chromosomes) + "\n"
    
    # run by chromosome
    # run on all by default? snp segment mapping?
    for chromosome in chromosomes:
        seg_range = qcp.getChromSegments(configdict["snp_segment_file"], chromosome)
        array_range = ":".join(str(x) for x in seg_range)
        
        job = assoc_script + "_chr" + str(chromosome)
        
        qsub_id = qcp.submitJob(job, driver_array, [assoc_rscript, config_file, str(chromosome)], holdid=holdid, arrayRange=array_range, verbose=False, **qsub_kwargs_no_email)
        jobid[job] = qsub_id.split(".")[0]
        
        holdid_comb = jobid[job]
        
        # combine segment files
        job = comb_script + "_chr" + str(chromosome)
        qsub_id = qcp.submitJob(job, driver, [comb_rscript, config_file, str(chromosome)], holdid=holdid_comb, verbose=False, **qsub_kwargs)
        jobid[job] = qsub_id
        
    holdid = [value for (key, value) in jobid.iteritems() if "combine" in key]

    
if arguments.plotqq:
    script = "plot_qq_manh"
    rscript = os.path.join(pipeline, "R", script + ".R")
    
    qsub_id = qcp.submitJob(script, driver, [rscript, config_file], holdid=holdid, verbose=True, **qsub_kwargs)
    
    
    
