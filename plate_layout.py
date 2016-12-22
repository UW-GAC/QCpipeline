#! /usr/local/bin/python2.7

"""Plate layout checks"""

import QCpipeline
import sys
import os
import subprocess
from optparse import OptionParser

usage = """%prog [options] config

Generate plate map layout plots to assess sample identity issues.

Required config parameters:
annot_scan_file     scan annotation file
ibd_obsrel_file     output from ibd of observed relatives

Optional config parameters [default]:
annot_scan_subjectCol       [subjectID]             name of subjectID column in scan annotation
annot_scan_sexCol           [sex]                   name of genetic sex column in scan annotation
annot_scan_annotSexCol      [annot.sex]             name of annotated sex column in scan annotation
annot_scan_annotSexMale     [M]                     coding for annotated sex = Male in scan annotation
annot_scan_annotSexFemale   [F]                     coding for annotated sex = Female in scan annotation
annot_scan_annotSexUnknown  [NA]                    coding for annotated sex = Unknown in scan annotation
annot_scan_plateCol         [Sample.Plate]          name of plate column in scan annotation
annot_scan_wellCol          [Sample.Well]           name of well column in scan annotation
ibd_unobsdup_file           [NA]                    output file from IBD with unobserved duplicates
min_num_problems            [1]                     minimum number of problems per plate required to display in pdf
out_annot_file              [out_scan_annot.RData]  output data frame with id problems listed
out_plate_plot              [plate_layout.pdf]      output pdf of plates with id problems shown
out_plate_file              [plate_layout.RData]    output file showing which plate goes with which page in pdf
scan_hilite_file            [NA]                    vector of user-specified scanIDs to higlight in pdf
scan_contaminated_file      [NA]                    vector of contaminated scanIDs to show in pdf
"""
parser = OptionParser(usage=usage)
parser.add_option("-e", "--email", dest="email", default=None,
                  help="email address for job reporting")
parser.add_option("-q", "--queue", dest="qname", default="gcc.q", 
                  help="cluster queue name [default %default]")
(options, args) = parser.parse_args()

if (len(args) != 1):
    parser.error("incorrect number of arguments")

config = args[0]
email = options.email
qname = options.qname

pipeline = os.path.dirname(os.path.abspath(sys.argv[0]))

driver = os.path.join(pipeline, "runRscript.sh")


jobid = dict()

job = "plate_layout"
rscript = os.path.join(pipeline, "R", job + ".R")
jobid["plate_layout"] = QCpipeline.submitJob(job, driver, [rscript, config], queue=qname, email=email)
