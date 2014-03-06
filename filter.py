#!/usr/bin/python
import sys

def selectRuns(uidlst, run_dat, filteredF):
    uids = dict([(l.strip().split()[0], None,) for l in open(uidlst, "rU")])
    rundat = open(run_dat, "rU")
    fout = open(filteredF, "w")
    for l in rundat:
	ident = l.strip().split()[0]
	if ident not in uids:
	    continue
	fout.write(l)
    fout.close()
    
if len(sys.argv) != 4:
    print "USAGE: %s <UID LIST> <RUN DAT> <OUTPUT>"%(sys.argv[0])
    sys.exit(1)
    
selectRuns(sys.argv[1], sys.argv[2], sys.argv[3])
