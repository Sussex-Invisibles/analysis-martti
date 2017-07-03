import ROOT
from rat import ratdb
db = ratdb.RATDBConnector('postgres://snoplus@pgsql.snopl.us:5400/ratdb')
result = db.fetch(obj_type="TELLIE_RUN", run=15414)
intensity_subrun0 = result[0]['data']['sub_run_info'][0]['photons']
print "Intensity is %d" % intensity_subrun0

ROOT.RAT.DB.Get().LoadDefaults()
fibre = "FT079A"
entry = ROOT.RAT.DB.Get().GetLink("FIBRE", fibre)
vector = ROOT.TVector3(entry.GetD("x"), entry.GetD("y"), entry.GetD("z"))
print "Position is ( %.1f | %.1f | %.1f )" % (vector.X(), vector.Y(), vector.Z())

