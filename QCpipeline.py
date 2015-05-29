"""Utility functions for QC pipeline"""

import sys
import subprocess

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
