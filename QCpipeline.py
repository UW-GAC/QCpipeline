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


def submitJob(job, cmd, args, holdid=None, email=None):
    """Sumbit a pipeline job.

    Usage: 
    jobid = submitJob(job, cmd, args, holdid=None, email=None)

    Arguments:
    job - name of job
    cmd - command to execute
    args - list of arguments to cmd
    holdid - list of job ids that must be complete before this job is run
    email - email address to notify when job is complete

    Retuns:
    id of the submitted job

    """

    nameStr = "-N " + job

    if holdid is not None:
        if isinstance(holdid, str):
            holdid = [holdid]
        holdStr = "-hold_jid " + ",".join(holdid)
    else:
        holdStr = ""

    if email is not None:
        emailStr = "-m e -M " + email
    else:
        emailStr = ""

    argStr = " ".join(args)

    qsub = "qsub %s %s %s %s %s" % (nameStr, holdStr, emailStr, cmd, argStr)
    process = subprocess.Popen(qsub, shell=True, stdout=subprocess.PIPE)
    pipe = process.stdout
    qsubout = pipe.readline()
    jobid = qsubout.split()[2]
    print qsubout
    return jobid
