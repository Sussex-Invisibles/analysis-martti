// ---------------------------------------------------------
// Goal:            Check TELLIE fibre stability over time
// Author:          Martti Nirkko, 13/06/2018
// Compile & run:   clear && g++ -g -o stability_check.exe stability_check.C `root-config --cflags --libs` -I$RATROOT/include/libpq -I$RATROOT/include -L$RATROOT/lib -lRATEvent_Linux && ./stability_check.exe
// ---------------------------------------------------------

// Includes and helper functions
#include "../HelperFunc.C"

// Global constants
const int RUN_CLUSTER = 1;  // whether running on cluster (0=local)
const int VERBOSE = 0;      // verbosity flag
const int IS_MC = 0;        // MC flag
const int MAXSUBRUNS = 20;

// Initialise functions
void get_subrun_nhit(int, int*, int*, bool, bool);

// Main program
using namespace std;
int main(int argc, char** argv) {
  
  // Test flag (0 = process all runs, else run number)
  const int TEST = (RUN_CLUSTER) ? 0 : 1;
  
  const int NFIBRES = 5;
  const string FIBRES[NFIBRES] = {"FT010A","FT033A","FT048A","FT079A","FT090A"};
  
  const int NSETS = 4;
  TGraph *g[NSETS][NFIBRES] = {{NULL}};

  string output = "output.root";
  ifstream out(output.c_str());
  TFile *outfile = NULL;
  bool extracted = false;
  if (out) extracted = true;

  if (extracted) {
    cout << "Extracting from file " << output << "..." << endl;  
    outfile = new TFile(output.c_str(),"READ");
    for (int iset=0; iset<NSETS; iset++) {
      for (int ifib=0; ifib<NFIBRES; ifib++) {
        string gname = Form("g_%d_%d",iset,ifib);
        g[iset][ifib] = (TGraph*)outfile->Get(gname.c_str());
        if (g[iset][ifib]) cout << "Loaded graph " << gname << " with " << g[iset][ifib]->GetN() << " points." << endl;
        else cout << "ERROR - could not find graph " << gname << endl;
      }
    }
  
  } else {

    outfile = new TFile(output.c_str(),"RECREATE");
    for (int iset=0; iset<NSETS; iset++) {
      for (int ifib=0; ifib<NFIBRES; ifib++) {
        g[iset][ifib] = new TGraph();
      }
    }

    // Loop over all fibres in list
    string input1 = "PIN_readings_auto.txt";
    ifstream in1(input1.c_str());
    if (!in1) { cerr<<"Failed to open "<<input1<<endl; exit(1); }
    //TH2F *nhits = new TH2F("nhits","NHit/event vs PIN reading",100,0,1e4,100,0,100);
    string line, fibre;
    int set, fib;
    int run, subrun, pin;
    float rms;
    float subrunpin[MAXSUBRUNS]={0}, subrunrms[MAXSUBRUNS]={0}, avgnhit[MAXSUBRUNS]={0};
    int subrunevs[MAXSUBRUNS]={0}, subrunhits[MAXSUBRUNS]={0};
    int nruns=0, nsubruns[NSETS][NFIBRES] = {{0}};
    for (int hdr=0; hdr<2; hdr++) {
      getline(in1,line);      // header
    }
    int lastrun = 0;
    while (true) {
      in1 >> fibre >> run >> subrun >> pin >> rms;
      if (!in1.good()) break;
      if (TEST && lastrun!=0 && run!=lastrun) continue; // only want one run
      subrunpin[subrun] = pin;
      subrunrms[subrun] = rms;

      if (run != lastrun) { // only do this once per run
        lastrun = run;
        
        // Determine dataset index
        set = 0;
        if (run > 111000) set = 1;
        if (run > 111800) set = 2;
        if (run > 115000) set = 3;
        
        // Determine fibre index
        fib = 0;
        for (int ifib=0; ifib<NFIBRES; ifib++) {
          if (fibre != FIBRES[ifib]) continue;
          fib = ifib;
          break;
        }
        
        for (int i=0; i<MAXSUBRUNS; i++) { subrunevs[i]=0; subrunhits[i]=0; }
        get_subrun_nhit(run, subrunevs, subrunhits, (bool)IS_MC, (bool)TEST);
    
        nruns++;
      }

      if(subrunevs[subrun]==0) avgnhit[subrun] = 0;
      else avgnhit[subrun] = (float)subrunhits[subrun]/subrunevs[subrun];
      //printf("Set %d, fibre %d (%s), subrun %d, PIN %.2f, %d events with an average %.2f nhits\n",set+1,fib+1,FIBRES[fib].c_str(),subrun,subrunpin[subrun],subrunevs[subrun],avgnhit[subrun]);
      g[set][fib]->SetPoint(nsubruns[set][fib],subrunpin[subrun],avgnhit[subrun]);
      nsubruns[set][fib]++;
    }
    printf("Ran over %d runs.\n",nruns);
    if (nruns==0) {
      cerr<<"*** ERROR *** No input files found."<<endl;
      return 1; 
    }
  
    for (int iset=0; iset<NSETS; iset++) {
      for (int ifib=0; ifib<NFIBRES; ifib++) {
        if (g[iset][ifib]->GetN()==0) continue;
        g[iset][ifib]->Write(Form("g_%d_%d",iset,ifib));
      }
    }

  }
  
  // Plot options
  gErrorIgnoreLevel = kWarning;
  gStyle->SetTitleOffset(1.3,"x");
  gStyle->SetTitleOffset(1.3,"y");
  
  // Plot graphs
  int xmin[NFIBRES] = {600,1100,500,800,700};
  int col[NSETS] = {2,3,4,92};
  TCanvas *c = new TCanvas("c","",1500,1000);
  c->Divide(3,2);
  TVirtualPad *p;
  for (int ifib=0; ifib<NFIBRES; ifib++) {
    p = c->cd(ifib+1);
    p->SetGrid();
    string title = Form("TELLIE fibre %s;PIN reading;Average NHit",FIBRES[ifib].c_str());
    p->DrawFrame(xmin[ifib],0,xmin[ifib]+700,200)->SetTitle(title.c_str());
    //p->SetLogx();
    //p->SetLogy();
    for (int iset=0; iset<NSETS; iset++) {
      if (!g[iset][ifib]) continue;
      if (g[iset][ifib]->GetN()==0) continue;
      g[iset][ifib]->SetMarkerStyle(7);
      g[iset][ifib]->SetMarkerColor(col[iset]);
      g[iset][ifib]->Draw("P same");
    }
  }
  p = c->cd(6);
  TLegend l(0.1,0.1,0.9,0.9);
  if(g[0][0]) l.AddEntry(g[0][0],"23 March 2018","p");
  if(g[1][0]) l.AddEntry(g[1][0],"25 March 2018","p");
  if(g[2][0]) l.AddEntry(g[2][0],"17 April 2018","p");
  if(g[3][0]) l.AddEntry(g[3][0],"10 July 2018","p");
  l.Draw();
  c->Print("stability_check.pdf");
  c->Print("stability_check.png");
  c->Close();

  // End of program
  outfile->Close();
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
  printf("Checking files for run %d... ", run);
  string fpath = (RUN_CLUSTER) ? "/lustre/scratch/epp/neutrino/snoplus/TELLIE/TELLIE_STABILITY_RUNS_PROCESSED" : "/home/nirkko/Software/SNOP/data";
  string fname = "";
  ifstream f;
  for (int pass=3;pass>=0;pass--) {
    fname = Form("%s/Analysis_r0000%d_s000_p00%d.root",fpath.c_str(),run,pass);
    f.open(fname.c_str());
    if (f.good()) break;
    fname = Form("%s/Calibration_r0000%d_s000_p00%d.root",fpath.c_str(),run,pass);
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
    int subrunevtcount[MAXSUBRUNS] = {0};
    int subrunhitcount[MAXSUBRUNS] = {0};
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
        subrunhitcount[subrunid] += nhitscleaned;
        subrunevtcount[subrunid]++;
        
      } // event loop
    } // entry loop
    cout << endl;
    for (int sr=0; sr<MAXSUBRUNS; sr++) {
      srevts[sr] = subrunevtcount[sr];
      srhits[sr] = subrunhitcount[sr];
      //printf("Subrun %d has %d events with %d hits.\n",sr,subrunevtcount[sr],subrunhitcount[sr]);
    }
    return;

}
