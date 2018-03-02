// ---------------------------------------------------------
// Goal:            Plot TELLIE PIN readings vs subruns
// Author:          Martti Nirkko, 07/08/2017
// Compile & run:   clear && g++ -g -o pin_readings.exe pin_readings.C `root-config --cflags --libs` -I$RATROOT/include/libpq -I$RATROOT/include -L$RATROOT/lib -lRATEvent_Linux && ./pin_readings.exe
// ---------------------------------------------------------
#include "../HelperFunc.C"

// Global constants
const int RUN_CLUSTER = 1;  // whether running on cluster (0=local)
const int VERBOSE = 0;      // verbosity flag
const int IS_MC = 0;        // Monte-Carlo flag 
  const int NSUBRUNS = 200; // Number of subruns

// Initialise functions
void get_subrun_nhit(int, int*, int*, bool, bool);

// Main program
using namespace std;
int main(int argc, char** argv) {
  
  // Test flag (0 = process all runs, else run number)
  const int TEST = (RUN_CLUSTER) ? 0 : 1;

  // Loop over all fibres in list
  string input = "PIN_readings_long.txt";
  ifstream in(input.c_str());
  if (!in) { cerr<<"Failed to open "<<input<<endl; exit(1); }
  //TH2F *nhits = new TH2F("nhits","NHit/event vs PIN reading",100,0,1e4,100,0,100);
  string line;
  float subrunpin[NSUBRUNS]={0}, subrunrms[NSUBRUNS]={0};
  float avgnhit[NSUBRUNS]={0};
  int run, subrun, pin;
  float rms;
  int nruns=0, nsubruns=0;
  for (int hdr=0; hdr<2; hdr++) {
    getline(in,line);      // header
  }
  int lastrun = 0;
  TGraphErrors *g[2];
  g[0] = new TGraphErrors();
  g[1] = new TGraphErrors();
  while (true) {
    in >> run >> subrun >> pin >> rms;
    if (!in.good()) break;
    if (TEST && lastrun!=0 && run!=lastrun) continue; // only want one run
    if (run != lastrun) nruns++;
    subrunpin[subrun] = pin;
    subrunrms[subrun] = rms;
    cout<<run<<","<<subrun<<","<<pin<<","<<rms<<endl;
    if (nruns==1) g[nruns-1]->SetPoint(subrun,subrun,pin);
    else g[nruns-1]->SetPoint(subrun,subrun+0.5,pin);
    g[nruns-1]->SetPointError(subrun,0,rms);
    if (run == lastrun) continue; // only do this once per run
    lastrun = run;
    cout<<"nruns="<<nruns<<", nsubruns="<<nsubruns<<endl;
    if (VERBOSE) printf("%6d %2d %5d %.2f\n", run, subrun, pin, rms);
    int subrunevs[NSUBRUNS], subrunhits[NSUBRUNS];
    //get_subrun_nhit(run, subrunevs, subrunhits, (bool)IS_MC, (bool)TEST);
    for (int sr=0; sr<NSUBRUNS; sr++) {
      if(subrunevs[sr]==0) continue;
      avgnhit[sr] = (float)subrunhits[sr]/subrunevs[sr];
      if (VERBOSE) printf("Run %6d subrun %2d has %4d events with an average nhit of %5.2f.\n", run, sr, subrunevs[sr], avgnhit[sr]);
      nsubruns++;
    }
  }
  printf("Ran over %d subruns.\n",nsubruns);
  if (nsubruns==0) {
    cerr<<"*** ERROR *** No input files found."<<endl;
    return 1; 
  }

  // Get zero points
  TGraph *zero = new TGraph();
  int nilpoints = 0;
  for (int i=0; i<NSUBRUNS; i++) {
    if (subrunpin[i]>0) continue;
    zero->SetPoint(nilpoints,i,0);
    nilpoints++;
  }
  
  // Plot graph
  TCanvas *c = new TCanvas("c","",800,600);
  c->SetGrid();
  string title = "TELLIE internal readout stability (FT010A);Subrun number;PIN reading";
  c->DrawFrame(-0.025*NSUBRUNS,0,1.025*NSUBRUNS,2000,title.c_str())->GetYaxis()->SetTitleOffset(1.4);
  //g->SetTitle(title.c_str());
  g[0]->SetMarkerStyle(7);
  g[0]->SetMarkerColor(2);
  g[0]->SetLineColor(2);
  g[0]->Draw("P same");
  g[1]->SetMarkerStyle(7);
  g[1]->SetMarkerColor(4);
  g[1]->SetLineColor(4);
  g[1]->Draw("P same");
  zero->SetMarkerStyle(8);
  zero->SetMarkerColor(4);
  zero->Draw("P same");
  TLegend *leg = new TLegend(0.74,0.77,0.89,0.89);
  leg->AddEntry(g[0],"Master","LP");
  leg->AddEntry(g[1],"Slave","LP");
  leg->Draw();
  c->Print("pin_vs_subrun.pdf");
  c->Print("pin_vs_subrun.png");
  c->Close();

  // End of program
  if (g) delete g[0];
  if (g) delete g[1];
  if (c) delete c;
  return 0;
}

// Returns the fitted light position for a given fibre/run
void get_subrun_nhit(int run, int *srevts, int *srhits, bool isMC=false, bool TEST=false) {

  // ********************************************************************
  // Initialisation
  // ********************************************************************
  
  // Check files for given run
  if(!TEST) printf("*****\n");
  //printf("Checking files for run %d... ", run);
  //string fpath = (RUN_CLUSTER) ? "/lustre/scratch/epp/neutrino/snoplus/TELLIE_PCA_RUNS_PROCESSED" : "/home/nirkko/Software/SNOP/data";
  string fpath = (RUN_CLUSTER) ? "/lustre/scratch/epp/neutrino/snoplus/TELLIE_TEST_RUNS" : "/home/nirkko/Software/SNOP/data";
  string fname = "";
  ifstream f;
  for (int pass=3;pass>=0;pass--) {
    //fname = Form("%s/Analysis_r0000%d_s000_p00%d.root",fpath.c_str(),run,pass);
    fname = Form("%s/SNOP_0000%d_00%d.l2.zdab",fpath.c_str(),run,pass);
    f.open(fname.c_str());
    if (f.good()) break;
  }
  string out = Form("./output/PCA_%d.out",run);
  string img = Form("./images/PCA_%d.pdf",run);
  //ifstream f(fname.c_str());
  ifstream g(out.c_str());
  ifstream h(img.c_str());
  int scanned_file = 0;
  if (!TEST && h.good()) {  // file downloaded and processed
    printf("already processed! Skipping run %d.\n",run);
    return;
  } else if(g.good()) {     // file extracted, but not processed
    printf("not processed! Generating plots for run %d.\n",run);
    scanned_file = 1;
  } else if(!f.good()) {    // file not downloaded
    printf("not downloaded! Skipping run %d.\n",run);
    return;
  } else {                  // file downloaded, but not processed
    printf("OK. Extracting data for run %d.\n",run);
  }

  // Initialise RAT
  RAT::DU::DSReader dsreader(fname);
  const RAT::DU::ChanHWStatus& chs = RAT::DU::Utility::Get()->GetChanHWStatus();
  const RAT::DU::PMTInfo& pmtinfo = RAT::DU::Utility::Get()->GetPMTInfo();
  const int NPMTS = pmtinfo.GetCount();

  // ********************************************************************
  // Sum PMT hit counts for entire run
  // ********************************************************************
    int subrunevtcount[NSUBRUNS] = {0};
    int subrunhitcount[NSUBRUNS] = {0};
    for(int iEntry=0; iEntry<dsreader.GetEntryCount(); iEntry++) {
      const RAT::DS::Entry& ds = dsreader.GetEntry(iEntry);
      int subrunid = ds.GetSubRunID();             // sub run this event belongs to
      
      // Loop over triggered events in each entry
      for(int iEv=0; iEv<ds.GetEVCount(); iEv++) { // mostly 1 event per entry
        const RAT::DS::EV& ev = ds.GetEV(iEv);
        
        // Trigger type
        int trig = ev.GetTrigType();
        if (!(trig & 0x8000)) continue;            // EXT trigger only
        
        // Event observables
        int nhitscleaned = ev.GetNhitsCleaned();   // calibrated PMT hits with removed crosstalk
        subrunhitcount[subrunid-1] += nhitscleaned;
        subrunevtcount[subrunid-1]++;
        
      } // event loop
    } // entry loop
    
    for (int sr=0; sr<NSUBRUNS; sr++) {
      srevts[sr] = subrunevtcount[sr];
      srhits[sr] = subrunhitcount[sr];
    }
    return;

}
