"""Utility functions for QC pipeline"""

import sys
import subprocess
import itertools
import numpy as np

def readConfig(file):
    """Read a pipeline config file.

    Usage: 
    config = readconfig(file)

    Arguments: 
    file - name of config file to read

    Returns: 
    dictionary with config values

    """

    config = dict()
    f = open(file, 'r')
    for line in f:
        if line.startswith('#'):
            continue

        tmp = line.split()
        if len(tmp) == 2:
            (key, value) = tmp
            config[key] = value
            
    f.close()
    return config


def submitJob(job, cmd, args, queue="gcc.q", holdid=None, email=None, qsubOptions="",
              arrayRange=None, verbose=True):
    """Sumbit a pipeline job.

    Usage: 
    jobid = submitJob(job, cmd, args, queue="gcc.q", holdid=None, email=None, qsubOptions="", arrayRange=None, verbose=True)

    Arguments:
    job - name of job
    cmd - command to execute
    args - list of arguments to cmd
    queue - compute cluster queue to submit job into
    holdid - list of job ids that must be complete before this job is run
    email - email address to notify when job is complete
    qsubOptions - additional options to pass to qsub
    arrayRange - specified, range of array jobs to pass to qsub (ie 1-23)
    verbose - Print out stdout from qsub?

    Returns:
    id of the submitted job

    """

    nameStr = "-N " + job
    
    queueStr = "-q " + queue
    
    if holdid is not None and holdid != []:
        if isinstance(holdid, str):
            holdid = [holdid]
        holdStr = "-hold_jid " + ",".join(holdid)
    else:
        holdStr = ""

    if email is not None:
        emailStr = "-m e -M " + email
    else:
        emailStr = ""

    if arrayRange is not None:
        arrayStr = "-t " + arrayRange
    else:
        arrayStr = ""
        
    argStr = " ".join(args)

    qsub = "qsub -S /bin/bash %s %s %s %s %s %s %s %s" % (qsubOptions, nameStr, arrayStr, holdStr, queueStr, emailStr, cmd, argStr)
    
    #print qsub # this is helpful to uncomment for occasional debugging
    
    process = subprocess.Popen(qsub, shell=True, stdout=subprocess.PIPE)
    pipe = process.stdout
    qsubout = pipe.readline()
    jobid = qsubout.split()[2]
    
    if verbose:
        print qsubout
    
    return jobid




# key for sorting chromosomes with integers first and other/failed last.
chromosomeSortKey = lambda item: (int(item.partition(" ")[0]) if item[0].isdigit() else float('inf'), item)

# function to be used as a decorator to sort chromosomes!
# put @sortChromosomes on the line above any function definition that returns chromosomes
def sortChromosomes(func):

    def sorter(*args, **kwargs):
        #return sorted(*args, key=chromosomeSortKey)
        x = func(*args, **kwargs)
        return sorted(x, key=chromosomeSortKey)

    return sorter



def getChromSegments(map_file, chromosome):
    
    x = getSegmentMapping(map_file)
    
    segments = x["segment"][x["chrom"] == int(chromosome)]
    
    return (min(segments), max(segments))



def getSegmentMapping(map_file):
    """ """
    x = np.genfromtxt(map_file, delimiter=",", names=True)
    
    # change data type of chrom and segment to integers
    dt = x.dtype
    names = dt.names
    # figure out how to get i that corresponds to "chrom"
    dt = dt.descr
    dt[x.dtype.names.index("chrom")] = ("chrom", "int")
    dt[x.dtype.names.index("segment")] = ("segment", "int")
    data = np.array(x, dtype=dt)

    return data

# chromosomeList is a list of strings, ie ["1", "2", "5:10"]
# does not deal with "other" or "failed yet"
@sortChromosomes
def parseChromosomes(chromosomeList):
    """Parses list of chromosomes/chromosome ranges to run on (e.g., ['1', '2', '5:7'] returns ['1', '2', '5', '6', '7']). No error checking for chromosomes in the correct range (1-23) or non-integer chromosomes in a:b statement."""
    # an example of what the list comprehension and itertools calls do:
    # given a list:
    #  ["1", "2", "6:8"]
    # returns ["1", "2", "6", "7", "8"]
    tmp = [[x[0]] if len(x) == 1 else range(int(x[0]), int(x[1])+1) for x in [r.split(":") for r in chromosomeList]]
    chromosomes = list(itertools.chain.from_iterable(tmp))
    
    chromosomes = [str(x) for x in chromosomes]
    return chromosomes


