// ---------------------------------------------------------
// Goal:            Plot angular response of TELLIE fibres
// Author:          Martti Nirkko, 26/04/2017
// Compile & run:   clear && g++ -o angular.exe angular.C `root-config --cflags --libs` -I$RATROOT/include/libpq -I$RATROOT/include -L$RATROOT/lib -lRATEvent_Linux && ./angular.exe
// ---------------------------------------------------------

// C++ stuff
#include <fstream>
#include <iostream>
#include <stdio.h>
#include <stdlib.h>

// ROOT stuff
#include "TROOT.h"
#include "TSystem.h"
#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TCanvas.h"
#include "TMath.h"
#include "TPad.h"
#include "TPaveStats.h"
#include "TStyle.h"
#include "TLegend.h"

// RAT stuff
#include "RAT/DU/DSReader.hh"
#include "RAT/DU/PMTInfo.hh"
#include "RAT/DU/Utility.hh"
#include "RAT/DS/Run.hh"
#include "RAT/DS/Entry.hh"
#include "RAT/DS/MC.hh"

// Helper functions
#include "../fibre_validation/HelperFunc.C"
#include "../fibre_validation/Xianguo.C"

// Global constants
const int RUN_CLUSTER = 0;  // whether running on cluster (0=local)
const int VERBOSE = 1;      // verbosity flag
const int IS_MC = 0;        // Monte-Carlo flag 

// Initialise functions
float angular(string, int, bool, bool);

// Main program
using namespace std;
int main(int argc, char** argv) {
  
  // Test flag (0 = process all runs, else run number)
  const int TEST = (RUN_CLUSTER) ? 0 : 102315;

  // Loop over all fibres in list
  string input = "../pca_runs/TELLIE_PCA.txt";
  ifstream in(input.c_str());
  if (!in) { cerr<<"Failed to open "<<input<<endl; exit(1); }
  string line, fibre;
  int node, channel, run, ipw, photons, pin, rms;
  float nhit;
  for (int hdr=0; hdr<2; hdr++) {
    getline(in,line);      // header
  }
  // Initialise RAT
  string example = (RUN_CLUSTER) ? "/lustre/scratch/epp/neutrino/snoplus/TELLIE_PCA_RUNS_PROCESSED/Analysis_r0000102315_s000_p000.root" : "/home/m/mn/mn372/Desktop/Analysis_r0000102315_s000_p000.root";
  //RAT::DU::DSReader dsreader(example);
  //const RAT::DU::PMTInfo& pmtinfo = RAT::DU::Utility::Get()->GetPMTInfo();
  //const int NPMTS = pmtinfo.GetCount();
  
  int nfiles=0; 
  while (true) {
    in >> node >> fibre >> channel >> run >> ipw >> photons >> pin >> rms >> nhit;
    if (!in.good()) break;
    if (TEST && TEST!=run) continue; // only want specified run
    //string fname = Form("TELLIE_angular_%s.root", fibre.c_str());
    //if (TEST) fname = "angular.root";
    float avgnhit = angular(fibre, run, IS_MC, TEST);
    printf("%3d %6s %4.1f = %4.1f\n",channel,fibre.c_str(),avgnhit,nhit);
    nfiles++;
  }
  printf("Ran over %d files.\n",nfiles);
  if (nfiles==0) { 
    cerr<<"*** ERROR *** No input files found."<<endl;
    return 1; 
  }
  return 0;
}

// Define macro
float angular(string fibre, int run, bool isMC=false, bool TEST=false) {

  // ********************************************************************
  // Initialisation
  // ********************************************************************
  
  // Check files for given run
  if(!TEST) printf("*****\n");
  printf("Checking files for run %d... ", run);
  string fpath = (RUN_CLUSTER) ? "/lustre/scratch/epp/neutrino/snoplus/TELLIE_PCA_RUNS_PROCESSED" : "/home/m/mn/mn372/Desktop";
  string fname = "";
  ifstream f;
  for (int pass=3;pass>=0;pass--) {
    fname = Form("%s/Analysis_r0000%d_s000_p00%d.root",fpath.c_str(),run,pass);
    f.open(fname.c_str());
    if (f.good()) break;
  }
  //string out = Form("./output/Angular_%s.out",fibre.c_str());
  string img = Form("./images/Angular_%s.pdf",fibre.c_str());
  //ifstream g(out.c_str());
  ifstream h(img.c_str());
  //int scanned_file = 0;
  if (!TEST && h.good()) {  // file downloaded and processed
    printf("already processed! Skipping fibre %s.\n",fibre.c_str());
    return -999.;
  //} else if(g.good()) {     // file extracted, but not processed
  //  printf("not processed! Generating plots for fibre %s.\n",fibre.c_str());
  //  scanned_file = 1;
  } else if(!f.good()) {    // file not downloaded
    printf("not downloaded! Skipping fibre %s.\n",fibre.c_str());
    return -999.;
  } else {                   // file downloaded, but not processed
    printf("OK. Extracting data for fibre %s.\n",fibre.c_str());
  }

  // Initialise RAT
  RAT::DU::DSReader dsreader(fname);
  const RAT::DU::ChanHWStatus& chs = RAT::DU::Utility::Get()->GetChanHWStatus();
  const RAT::DU::PMTInfo& pmtinfo = RAT::DU::Utility::Get()->GetPMTInfo();
  const int NPMTS = pmtinfo.GetCount();
  
  // Get fibre info (from RATDB) 
  RAT::DB *db = RAT::DB::Get();
  //db->LoadDefaults();	  // Already done when calling DU::Utility::Get()
  RAT::DBLinkPtr entry = db->GetLink("FIBRE",fibre);
  TVector3 fibrepos(entry->GetD("x"), entry->GetD("y"), entry->GetD("z")); // position
  TVector3 fibredir(entry->GetD("u"), entry->GetD("v"), entry->GetD("w")); // direction
  TVector3 lightpos = fibrepos + 2*fibrepos.Mag()*fibredir; // projected light spot centre
  if (VERBOSE) cout << "RATDB: fibre " << fibre << ", pos " << printVector(fibrepos) << ", dir " << printVector(fibredir) << endl;
  
  // Initialise histograms
  int NBINS = 50;
  TH2D *h2 = new TH2D("h2",fibre.c_str(),NBINS,0,0.25,NBINS,240,290);
  TH1D *h1 = new TH1D("h1",fibre.c_str(),NBINS,0,0.25);
  
  // Loop over entries
  int totalnhit=0, count=0;
  for(int iEntry=0; iEntry<dsreader.GetEntryCount(); iEntry++) {
    const RAT::DS::Entry& ds = dsreader.GetEntry(iEntry);
    RAT::DS::MC mc;
    if (isMC) mc = ds.GetMC();  // don't initialise this for real data (crashes)
    
    // Loop over triggered events in each entry
    for(int iEv=0; iEv<ds.GetEVCount(); iEv++) {        // mostly 1 event per entry
      const RAT::DS::EV& ev = ds.GetEV(iEv);
      
      // Trigger type
      int trig = ev.GetTrigType();
      if (!(trig & 0x8000)) continue;                   // EXT trigger only
        
      // Event observables
      int nhitscleaned = ev.GetNhitsCleaned();   // calibrated PMT hits with removed crosstalk
      totalnhit += nhitscleaned;
      count++;
      
      // PMT information
      const RAT::DS::CalPMTs& pmts = ev.GetCalPMTs();
      for(int iPMT=0; iPMT<pmts.GetNormalCount(); iPMT++) {
        RAT::DS::PMTCal pmt = pmts.GetNormalPMT(iPMT);
        int pmtID = pmt.GetID();
        if (!chs.IsTubeOnline(pmtID)) continue;             // test CHS
        if (pmt.GetCrossTalkFlag()) continue;               // remove crosstalk
        if (pmt.GetStatus().GetULong64_t(0) != 0) continue; // test PCA
        TVector3 pmtpos = pmtinfo.GetPosition(pmtID);
        TVector3 track = pmtpos-fibrepos;
        double theta = track.Angle(fibredir);               // angle w.r.t. fibre [rad]
        double pmttime = pmt.GetTime();                     // hit time [ns]
        // Fill histogram
        h2->Fill(theta/pi, pmttime);
      } // pmt loop
	  
    } // event loop
  } // entry loop
  
  // Fit angular slices with a gaussian distribution in time
  h2->FitSlicesY();
  TH1D *h2_0 = (TH1D*)gDirectory->Get("h2_0");  // Constant
  TH1D *h2_1 = (TH1D*)gDirectory->Get("h2_1");  // Mean
  TH1D *h2_2 = (TH1D*)gDirectory->Get("h2_2");  // StdDev
  
  // Put mean and RMS into 1D histogram
  for (int b=0; b<NBINS+2; b++) {
    h1->SetBinContent(b, h2_1->GetBinContent(b));
    h1->SetBinError(b, h2_2->GetBinContent(b));
    //if (VERBOSE) printf("#%2d %.4f %6.2f %6.2f\n",b,h1->GetBinCenter(b),h1->GetBinContent(b),h1->GetBinError(b));
  }
  
  TCanvas *c0 = new TCanvas("","",1200,600);
  c0->Divide(2,1);
  gStyle->SetOptStat(1111);
  gStyle->SetStatX(0.89);
  gStyle->SetStatY(0.33);
  gStyle->SetStatW(0.22);
  gStyle->SetStatH(0.14);
  
  c0->cd(1)->SetGrid();
  h2->SetXTitle("Angle [#pi]");
  h2->SetYTitle("Time [ns]");
  h2->Draw("colz");
  h2->GetXaxis()->SetTitleOffset(1.3);
  h2->GetYaxis()->SetTitleOffset(1.5);
  
  c0->cd(2)->SetGrid();
  h1->SetXTitle("Angle [#pi]");
  h1->SetYTitle("Time [ns]");
  h1->SetLineWidth(2);
  h1->SetStats(0);
  h1->Draw("e");
  h1->GetYaxis()->SetRangeUser(240,290);
  h1->GetXaxis()->SetTitleOffset(1.3);
  h1->GetYaxis()->SetTitleOffset(1.5);
  
  // Save canvas and close
  string outfile = "angular";
  if (!TEST) outfile = Form("images/angular_%s",fibre.c_str());
  c0->Print(Form("%s.png",outfile.c_str()));
  c0->Print(Form("%s.pdf",outfile.c_str()));
  c0->Close();
  
  if (c0) delete c0;
  if (h2) delete h2;
   
  return (float)totalnhit/count;
}

