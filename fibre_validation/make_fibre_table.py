#!/usr/bin/python
import sys
import ROOT
from rat import ratdb
fname=sys.argv[1]
db = ROOT.RAT.DB.Get()
db.LoadDefaults()
output = open("TELLIE_FIBRES.txt",'w')
output.write("#Fibre x y z u v w\n")
output.write("#-----------------\n")
with open(fname) as f:
    for line in f:
        fibre=line.rstrip('\n')
        print fibre
        entry = db.GetLink("FIBRE",fibre)
        x = entry.GetD("x")
        y = entry.GetD("y")
        z = entry.GetD("z")
        u = entry.GetD("u")
        v = entry.GetD("v")
        w = entry.GetD("w")
        output.write("%6s %lf %lf %lf %lf %lf %lf\n" % (fibre, x, y, z, u, v, w))
    f.close()
output.close()
