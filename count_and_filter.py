#!/usr/bin/python

import sys


"""
#1. inputs will include a list of UIDS or hyperstats file, and the iterations data.
#2. the record counts in the iteration rows will be computed by:
a. split row on \t
b. remove UID column
c. divide remaining columns by the total record count found in the UID
"""

if len(sys.argv) != 5:
    print "USAGE: %s <UID list/hyperstats file> <Iteration data> <filtered output> <counts output>"%(sys.argv[0])
    sys.exit(0)

UIDs = dict([ (i.strip().split("\t")[0], None ) for i in open(sys.argv[1])])

counts = {}
fout = open(sys.argv[3], "w")

for l in open(sys.argv[2]):
    l = l.strip()
    ls = l.split("\t")
    if ls[0] not in UIDs:
        continue
    print >> fout, l
    splitID = ls[0].split("_")
    congruent = int(splitID[-6])
    total = int(splitID[-7])
    reclen = (len(ls) - 1) / total
    
    for x in ls[1: (congruent * reclen + 1) :reclen]:
        if x not in counts:
            counts[x] = 0
        counts[x] += 1
        print "x"
    print "-"*80
fout.close()
fout = open(sys.argv[4], "w")
print >> fout, "\n".join([ "%s\t%s"%(j,k) for j, k in counts.iteritems()])
fout.close()
