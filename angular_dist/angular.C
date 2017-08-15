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

// Global constants
const double pi = TMath::Pi();
const int RUN_CLUSTER = 0;  // whether running on cluster (0=local)
const int VERBOSE = 1;      // verbosity flag
const int IS_MC = 0;        // Monte-Carlo flag 

// Initialise functions
float angular(string, string, bool, bool);

// Run program
using namespace std;
int main(int argc, char** argv) {
  
  // Test flag (0 = process all runs, else run number)
  const int TEST = (RUN_CLUSTER) ? 0 : 101834;

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
  RAT::DU::DSReader dsreader(example);
  const RAT::DU::PMTInfo& pmtinfo = RAT::DU::Utility::Get()->GetPMTInfo();
  const int NPMTS = pmtinfo.GetCount();
  printf("Initialised DSReader & PMTInfo (%d PMTs).\n",NPMTS);
  
  int nfiles=0; 
  while (true) {
    in >> node >> fibre >> channel >> run >> ipw >> photons >> pin >> rms >> nhit;
    if (!in.good()) break;
    if (TEST && TEST!=run) continue; // only want specified run
    //string fname = Form("TELLIE_angular_%s.root", fibre.c_str());
    //if (TEST) fname = "angular.root";
    float avgnhit = angular(example.c_str(), fibre, IS_MC, TEST);
    printf("%3d %6s %4.1f\n",channel,fibre.c_str(),avgnhit);
    nfiles++;
  }
  printf("Ran over %d files.\n",nfiles);
  if (nfiles==0) { cerr<<"*** ERROR *** No input files found."<<endl; return 1; }
  
  return 0;
}

// Define macro
float angular(string fname, string fibre, bool isMC, bool TEST) {

  // Initialise RAT
  RAT::DU::DSReader dsreader(fname);
  const RAT::DU::PMTInfo& pmtinfo = RAT::DU::Utility::Get()->GetPMTInfo();
  
  // Initialise histogram
  TH2D *hist = new TH2D("hist",fibre.c_str(),50,0,0.5,40,1550,1750);
  hist->SetXTitle("Angle [#pi]");
  hist->SetYTitle("TAC [ ]");
  
  // Loop over entries
  int count=0, hitpmts=0, badpmts=0, init=0;
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
      
      // Skip first entry in real data (TELLIE init command)
      if (!isMC && !init) {
        init=1;
        continue;
      }
      
      // Event observables
      int nhits = ev.GetNhits();			// normal/inward looking PMT hits
      count++;
      
      // Get fibre position and angle
      TVector3 fibrepos, fibredir, lightpos;
      if (isMC) {
        fibrepos = mc.GetMCParticle(0).GetPosition();
        fibredir = mc.GetMCParticle(0).GetMomentum().Unit();
        //double rmc = mcpos.Mag()/1e3;       // radius [m]
        //double amc = mcpos.Angle(mcang);    // angle w.r.t. centre [rad]
      } else {
        // TODO - this is hardcoded for a test file!
        fibrepos.SetXYZ(3650.54,4359.93,-6173.26);
        fibredir.SetXYZ(-0.401593,-0.552760,0.730191);
        lightpos = fibrepos + 2*fibrepos.Mag()*fibredir; // projected light spot centre
      }
      
      // PMT information
      const RAT::DS::UncalPMTs& pmts = ev.GetUncalPMTs();
      for(int iPMT=0; iPMT<pmts.GetCount(); iPMT++) {
        int pmtID = pmts.GetPMT(iPMT).GetID();
        int pmttac = pmts.GetPMT(iPMT).GetTime();
        TVector3 pmtpos = pmtinfo.GetPosition(pmtID);
        TVector3 track = pmtpos-fibrepos;
        double theta = track.Angle(lightpos-fibrepos);
        hist->Fill(theta/pi,pmttac);
        hitpmts++;
      } // pmt loop
	  
    } // event loop
  } // entry loop
  
  gStyle->SetOptStat(1111);
  gStyle->SetStatX(0.89);
  gStyle->SetStatY(0.33);
  gStyle->SetStatW(0.22);
  gStyle->SetStatH(0.14);
  gStyle->SetPadGridX(true);
  gStyle->SetPadGridY(true);
  gStyle->SetTitleXOffset(1.3);
  gStyle->SetTitleYOffset(1.5);
  
  TCanvas *c0 = new TCanvas("","",1200,900);
  hist->Draw("colz");
  hist->GetXaxis()->SetTitleOffset(1.3);
  hist->GetYaxis()->SetTitleOffset(1.5);
  string imgname = Form("images/tellie_%s.png",fibre.c_str());
  if (TEST) imgname = "angular.png";    // overwrite
  c0->Print(imgname.c_str());
  c0->Close();
  
  if (c0) delete c0;
  if (hist) delete hist;
   
  return (float)hitpmts/count;
}

