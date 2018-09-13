#!/usr/bin/python
import sys
import rat
from ROOT import RAT
#from rat import ratdb
fname=sys.argv[1]
db = RAT.DB.Get() #ratdb.RATDBConnector('postgres://snoplus@pgsql.snopl.us:5400/ratdb')
run = RAT.DS.Run()
output = open("TELLIE_PCA_auto.txt",'w')
output.write("#Fibre Channel Run IPW Photons PIN RMS\n")
output.write("#-------------------------------------\n")
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
        #if num == 101705 or num == 101905: photons = subRunInfo[0]['photons'].getReal()
        #else: photons = subRunInfo[0]['photons'].getInteger()
        try: photons = subRunInfo[0]['photons'].getReal()
        except: photons = subRunInfo[0]['photons'].getInteger()
        pinval=0.
        pinrms=0.
        for sr in range(0,39):
            #if num == 101852: continue
            try: pinval += subRunInfo[sr]['pin_value'].getInteger()
            except: pinval += subRunInfo[sr]['pin_value'].getReal()
            try: pinrms += subRunInfo[sr]['pin_rms'].getReal()
            except: pinrms += subRunInfo[sr]['pin_rms'].getInteger()
        pinval /= 40.
        pinrms /= 40.
        output.write("%6s %2d %6d %5d %6d %5d %5d\n" % (fibre, channel, num, ipw, int(photons), pinval, int(pinrms)))
    f.close()
output.close()

