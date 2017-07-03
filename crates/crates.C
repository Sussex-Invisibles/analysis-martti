// ---------------------------------------------------------
// Goal:            Figure out which TELLIE fibres hit a given crate
// Author:          Martti Nirkko, 26/04/2017
// Compile & run:   clear && g++ -o crates.exe crates.C `root-config --cflags --libs` -I$RATROOT/include/libpq -I$RATROOT/include -L$RATROOT/lib -lRATEvent_Linux && ./crates.exe
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

using namespace std;

// Global constants
const double pi = TMath::Pi();

// Initialise functions
float crates(string, int, bool);

// Run program
int main(int argc, char** argv) {

  // Crate causing problems
  const int BAD_CRATE = 12;
  
  // Monte-Carlo simulation
  const int IS_MC = 1;
  
  // Loop over all fibres in list 
  ifstream in;
  in.open("tellie_fibres.list");
  string fibre, fname;
  float badpmtfrac;
  int channel=0;
  printf("Estimating hit percentages on crate #%d...\n\n", BAD_CRATE);
  printf("%-3s %-6s %-5s\n","Ch.","Fibre","Frac.");
  printf("-----------------\n");
  while (true) {
    in >> fibre;
    if (!in.good()) break;
    string fname = Form("output/TELLIE_test_%s.root", fibre.c_str());  
    channel++;
    badpmtfrac = crates(fname.c_str(), BAD_CRATE, IS_MC);
    printf("%3d %6s %4.1f%%\n",channel,fibre.c_str(),100.*badpmtfrac);
  }

  return 0;
}

// Define macro
float crates(string fname, int bad_crate, bool isMC) {

  // Initialise RAT
  RAT::DU::DSReader dsreader(fname);
  const RAT::DU::PMTInfo& pmtinfo = RAT::DU::Utility::Get()->GetPMTInfo();
  
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
        //cout << "*** INFO *** Skipping entry " << iEntry << " (_init_ trigger)!" << endl;
        continue;
      }
      
      // Event observables
      int nhits = ev.GetNhits();			// normal/inward looking PMT hits
      count++;
      
      // Get some MC variables, if available
      if (isMC) {
        TVector3 mcpos = mc.GetMCParticle(0).GetPosition();
        TVector3 mcang = mc.GetMCParticle(0).GetMomentum().Unit();
        double rmc = mcpos.Mag()/1e3;       // radius [m]
        double amc = mcpos.Angle(mcang);    // angle w.r.t. centre [rad]
      }
      
      // PMT information
      const RAT::DS::UncalPMTs& pmts = ev.GetUncalPMTs();
      for(int iPMT=0; iPMT<pmts.GetCount(); iPMT++) {
        int pmtID = pmts.GetPMT(iPMT).GetID();
        int crate = pmts.GetPMT(iPMT).GetCrate();
        hitpmts++;
        //cout << "PMT #" << pmtID << " is in crate " << crate << endl;
        if (crate==bad_crate) badpmts++;
      } // pmt loop
	  
    } // event loop
  } // entry loop
  
  //cout << endl << "Found " << count << " events, hitting " << Form("%.1f",(float)badpmts/count) << " PMTs in crate " << bad_crate << " on average." << endl << endl;
  
  /*
  gStyle->SetOptStat(110);
  gStyle->SetStatX(0.88);
  gStyle->SetStatY(0.88);
  gStyle->SetStatW(0.28);
  gStyle->SetStatH(0.18);
  gStyle->SetPadGridX(true);
  gStyle->SetPadGridY(true);
  
  TCanvas *c0 = new TCanvas("","",1200,900);
  c0->Print("test.png");
  c0->Close();
  
  if (c0) delete c0;
  */
  return (float)badpmts/hitpmts;
}

