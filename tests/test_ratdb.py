import ROOT
from rat import ratdb
db = ratdb.RATDBConnector('postgres://snoplus@pgsql.snopl.us:5400/ratdb')
result = db.fetch(obj_type="TELLIE_RUN", run=102315)
trigger_delay = result[0]['data']['sub_run_info'][0]['trigger_delay']
fibre_delay = result[0]['data']['sub_run_info'][0]['fibre_delay']
print "Fibre delay %.1fns, trigger delay %dns" % (fibre_delay, trigger_delay)

ROOT.RAT.DB.Get().LoadDefaults()
fibre = "FT019A"
entry = ROOT.RAT.DB.Get().GetLink("FIBRE", fibre)
vector = ROOT.TVector3(entry.GetD("x"), entry.GetD("y"), entry.GetD("z"))
print "Position is ( %.1f | %.1f | %.1f )mm" % (vector.X(), vector.Y(), vector.Z())

