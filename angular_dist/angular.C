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

using namespace std;

// Global constants
const double pi = TMath::Pi();

// Initialise functions
float angular(string, string, bool);

// Run program
int main(int argc, char** argv) {

  // Test on only one file
  const int TEST = 0;
  
  // Monte-Carlo simulation
  const int IS_MC = 1;
  
  // Loop over all fibres in list 
  string fibre, fname;
  if (TEST) {
    fibre = "FT077A";
    fname = "TELLIE_test/output/TELLIE_test_FT077A.root";
    cout << "Average nhits = " << angular(fname, fibre, IS_MC) << endl;
  } else {
    ifstream in;
    in.open("TELLIE_test/tellie_fibres.list");
    float avgnhit;
    int channel=0;
    printf("%-3s %-6s %-4s\n","Ch.","Fibre","NHit");
    printf("----------------\n");
    while (true) {
      in >> fibre;
      if (!in.good()) break;
      string fname = Form("TELLIE_test/output/TELLIE_test_%s.root", fibre.c_str());  
      channel++;
      avgnhit = angular(fname.c_str(), fibre, IS_MC);
      printf("%3d %6s %4.1f\n",channel,fibre.c_str(),avgnhit);
    }
  }
  
  return 0;
}

// Define macro
float angular(string fname, string fibre, bool isMC) {

  // Initialise RAT
  RAT::DU::DSReader dsreader(fname);
  const RAT::DU::PMTInfo& pmtinfo = RAT::DU::Utility::Get()->GetPMTInfo();
  
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
      
      // Get some MC variables, if available
      TVector3 mcpos, mcang;
      if (isMC) {
        mcpos = mc.GetMCParticle(0).GetPosition();
        mcang = mc.GetMCParticle(0).GetMomentum().Unit();
        double rmc = mcpos.Mag()/1e3;       // radius [m]
        double amc = mcpos.Angle(mcang);    // angle w.r.t. centre [rad]
      }
      
      // PMT information
      const RAT::DS::UncalPMTs& pmts = ev.GetUncalPMTs();
      for(int iPMT=0; iPMT<pmts.GetCount(); iPMT++) {
        int pmtID = pmts.GetPMT(iPMT).GetID();
        int pmttac = pmts.GetPMT(iPMT).GetTime();
        TVector3 pmtpos = pmtinfo.GetPosition(pmtID);
        TVector3 track = pmtpos-mcpos;
        double theta = track.Angle(-mcpos);
        hist->Fill(theta,pmttac);
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
  c0->Print(imgname.c_str());
  c0->Close();
  
  if (c0) delete c0;
  if (hist) delete hist;
   
  return (float)hitpmts/count;
}

