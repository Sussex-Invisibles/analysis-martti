#!/usr/bin/python
import sys
import rat
from ROOT import RAT
fname=sys.argv[1]
db = RAT.DB.Get() #ratdb.RATDBConnector('postgres://snoplus@pgsql.snopl.us:5400/ratdb')
#db = ratdb.RATDBConnector('postgres://snoplus@pgsql.snopl.us:5400/ratdb')
run = RAT.DS.Run()
output = open("TELLIE_delays.txt",'w')
output.write("#Run Fibre_delay Trig_delay PCA_offset\n")
output.write("#-------------------------------------\n")
with open(fname) as f:
    l=0
    for line in f:
        num=int(line.rstrip('\n'))
        l = l+1
        print "\n#%d - Getting info for run %d..." % (l, num)
        run.SetRunID(num)
        db.BeginOfRun(run)

        # Get information from TELLIE run table
        info = db.GetLink("TELLIE_RUN") #db.fetch(obj_type="TELLIE_RUN", run)
        subRunInfo = info.GetJSON("sub_run_info")
        try:
            fibre_delay = subRunInfo[0]['fibre_delay'].getReal()
        except Exception:
            fibre_delay = subRunInfo[0]['fibre_delay'].getInteger()
        trig_delay = subRunInfo[0]['trigger_delay'].getInteger()
        if   (num==101929): fibre_name = "FT093A"
        elif (num==101933): fibre_name = "FT094A"
        elif (num==101944): fibre_name = "FT047A"
        elif (num==102157): fibre_name = "FT101A"
        else: fibre_name = subRunInfo[0]['fibre'].getString()

        # Get information from Mark's fibre offsets table
        run.SetRunID(103000)  # crazy DB hack cause of multiple versions
        db.BeginOfRun(run)
        offsets = db.GetLink("TELLIE_PCA_FIBRE_OFFSETS",fibre_name[2:-1])
        this_run = offsets.GetI("run_number_used")
        pca_offset = offsets.GetD("fibre_pca_offset")

        output.write("%6d %4.1f %3d %8.3f\n" % (num, fibre_delay, trig_delay, pca_offset))
    f.close()
output.close()

