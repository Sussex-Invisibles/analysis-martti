#!/usr/bin/python
import sys
import rat
from ROOT import RAT
from rat import ratdb
fname=sys.argv[1]
db = RAT.DB.Get() #ratdb.RATDBConnector('postgres://snoplus@pgsql.snopl.us:5400/ratdb')
run = RAT.DS.Run()
output = open("TELLIE_PCA.txt",'w')
output.write("#Fibre Channel Run IPW Photons\n")
output.write("#-----------------------------\n")
with open(fname) as f:
    for line in f:
        num=int(line.rstrip('\n'))
        run.SetRunID(num)
        db.BeginOfRun(run)
        info = db.GetLink("TELLIE_RUN") #db.fetch(obj_type="TELLIE_RUN", run)
        subRunInfo = info.GetJSON("sub_run_info")
        print num
        channel = subRunInfo[0]['channel'].getInteger()
        # Catch exceptions where DB has no 'fibre' entry
        if   num == 100873: fibre = "FT019A"
        elif num == 101929: fibre = "FT093A"
        elif num == 101933: fibre = "FT094A"
        elif num == 101944: fibre = "FT047A"
        elif num == 102157: fibre = "FT101A"
        else: fibre = subRunInfo[0]['fibre'].getString()
        ipw     = subRunInfo[0]['pulse_width'].getInteger()
        # Catch exceptions where DB has non-integral 'photons' entry
        if num == 101705 or num == 101905: photons = subRunInfo[0]['photons'].getReal()
        else: photons = subRunInfo[0]['photons'].getInteger()
        output.write("%6s %2d %6d %5d %6d\n" % (fibre, channel, num, ipw, int(photons)))
    f.close()
output.close()

