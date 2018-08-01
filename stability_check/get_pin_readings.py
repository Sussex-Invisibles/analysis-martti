#!/usr/bin/python
import sys
import rat
from ROOT import RAT
fname=sys.argv[1]
db = RAT.DB.Get() #ratdb.RATDBConnector('postgres://snoplus@pgsql.snopl.us:5400/ratdb')
run = RAT.DS.Run()
output = open("PIN_readings_auto.txt",'w')
output.write("#Fibre Run Subrun IPW PIN RMS\n")
output.write("#----------------------------\n")
with open(fname) as f:
    l=0
    for line in f:
        num=int(line.rstrip('\n'))
        l = l+1
        print "\n#%d - Getting info for run %d..." % (l, num)
        run.SetRunID(num)
        db.BeginOfRun(run)
        info = db.GetLink("TELLIE_RUN") #db.fetch(obj_type="TELLIE_RUN", run)
        subRunInfo = info.GetJSON("sub_run_info")
        # Find fibre name
        fibre = subRunInfo[0]['fibre'].getString()
        # Find number of subruns
        nsubruns = 0
        while 1:
            try:
                channel = subRunInfo[nsubruns]['channel'].getInteger()
                nsubruns = nsubruns+1
            except:
                channel = 0
                print 'Found %d subruns.\n' % nsubruns
                break
        # Get PIN readings for each subrun
        for sr in range(0,nsubruns):
            ipw = subRunInfo[sr]['pulse_width'].getInteger()
            pinval = subRunInfo[sr]['pin_value'].getInteger()
            try:
                pinrms = subRunInfo[sr]['pin_rms'].getReal()
            except:
                pinrms = subRunInfo[sr]['pin_rms'].getInteger()
            output.write("%6s %6d %3d %5d %5d %6.2f\n" % (fibre, num, sr, ipw, pinval, pinrms))
    f.close()
output.close()

