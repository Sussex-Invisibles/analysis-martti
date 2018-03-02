#!/usr/bin/python
import sys
import rat
from ROOT import RAT
#from rat import ratdb
fname=sys.argv[1]
db = RAT.DB.Get() #ratdb.RATDBConnector('postgres://snoplus@pgsql.snopl.us:5400/ratdb')
run = RAT.DS.Run()
output = open("TELLIE_delays_auto.txt",'w')
output.write("#Fibre Run trigger_delay fibre_delay\n")
output.write("#-------------------------------------\n")
with open(fname) as f:
    for line in f:
        num=int(line.rstrip('\n'))
        run.SetRunID(num)
        db.BeginOfRun(run)
        info = db.GetLink("TELLIE_RUN") #db.fetch(obj_type="TELLIE_RUN", run)
        subRunInfo = info.GetJSON("sub_run_info")
        print num
        # Catch exceptions where DB has no 'fibre' entry
        if   num == 100873: fibre = "FT019A"
        elif num == 101929: fibre = "FT093A"
        elif num == 101933: fibre = "FT094A"
        elif num == 101944: fibre = "FT047A"
        elif num == 102157: fibre = "FT101A"
        else: fibre = subRunInfo[0]['fibre'].getString()
        trigger_delay = subRunInfo[0]['trigger_delay'].getInteger()
        try: fibre_delay   = subRunInfo[0]['fibre_delay'].getReal()
        except: fibre_delay   = subRunInfo[0]['fibre_delay'].getInteger()
        output.write("%6s %6d %4d %5.2f\n" % (fibre, num, trigger_delay, fibre_delay))
    f.close()
output.close()

