// ---------------------------------------------------------
// Goal:            Figure out which TELLIE fibres hit a given crate
// Author:          Martti Nirkko, 26/04/2017
// Compile & run:   clear && g++ -g -o occupancy.exe occupancy.C `root-config --cflags --libs` -I$RATROOT/include/libpq -I$RATROOT/include -L$RATROOT/lib -lRATEvent_Linux && ./occupancy.exe
// ---------------------------------------------------------

// C++ stuff
#include <fstream>
#include <iostream>
#include <stdio.h>
#include <stdlib.h>

// ROOT stuff
#include "TROOT.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TF2.h"
#include "TH2.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TMarker.h"
#include "TMath.h"
#include "TPad.h"
#include "TPaveStats.h"
#include "TStyle.h"
#include "TSystem.h"

// RAT stuff
#include "RAT/DU/DSReader.hh"
#include "RAT/DU/PMTInfo.hh"
#include "RAT/DU/Utility.hh"
#include "RAT/DS/Run.hh"
#include "RAT/DS/Entry.hh"
#include "RAT/DS/MC.hh"

// Helper functions
#include "HelperFunc.C"
#include "Xianguo.C"

// Global constants
const int RUN_CLUSTER=0;    // whether running on cluster (0=local)
const int IS_MC = 0;        // Monte-Carlo flag 

// Initialise functions
void occupancy(string, int, bool, bool);

// Main program
using namespace std;
int main(int argc, char** argv) {
  
  // Test flag (0 = process all runs, else run number)
  const int TEST = (RUN_CLUSTER) ? 0 : 101849;
  
  // Loop over all fibres in list
  string input = (RUN_CLUSTER) ? "../pca_runs/TELLIE_PCA.txt" : "TELLIE_PCA.txt";
  ifstream in(input.c_str());
  if (!in) { cerr<<"Failed to open TELLIE_PCA.txt"<<endl; exit(1); }
  string line, fibre;
  int node, channel, run, ipw, photons;
  float nhit;
  for (int hdr=0; hdr<2; hdr++) {
    getline(in,line);      // header
  }
  while (true) {
    in >> node >> fibre >> channel >> run >> ipw >> photons >> nhit;
    if (!in.good()) break;
    if (TEST && TEST!=run) continue; // only want specified run
    //printf("%6s %2d %6d %5d %6d\n", fibre.c_str(), channel, run, ipw, photons);
    occupancy(fibre, run, IS_MC, TEST);
  }

  return 0;
}

// Returns the fitted light position for a given fibre/run
void occupancy(string fibre, int run, bool isMC=false, bool TEST=false) {

  // ********************************************************************
  // Initialisation
  // ********************************************************************
  
  // Check files for given run
  if(!TEST) printf("*****\n");
  printf("Checking files for run %d... ", run);
  string fpath = (RUN_CLUSTER) ? "/lustre/scratch/epp/neutrino/snoplus/TELLIE_PCA_RUNS_PROCESSED" : "/home/nirkko/Desktop/fibre_validation";
  string fname = Form("%s/Analysis_r0000%d_s000_p000.root",fpath.c_str(),run);
  string out   = Form("./output/PCA_%s.pdf",fibre.c_str());
  ifstream f(fname.c_str());
  ifstream g(out.c_str());
  if (!TEST && g.good()) {   // file downloaded and processed
    printf("already processed! Skipping fibre %s.\n",fibre.c_str());
    return;
  } else  if(!f.good()) {    // file not downloaded
    printf("not downloaded! Skipping fibre %s.\n",fibre.c_str());
    return;
  } else {                   // file downloaded, but not processed
    printf("OK. Processing fibre %s.\n",fibre.c_str());
  }

  // Initialise RAT
  RAT::DU::DSReader dsreader(fname);
  const RAT::DU::PMTInfo& pmtinfo = RAT::DU::Utility::Get()->GetPMTInfo();
  const int NPMTS = pmtinfo.GetCount();
  printf("Initialised DSReader & PMTInfo.\n");
  
  // Initialise histograms
  TH2D *hicos = new TH2D("hicos","PMT positions",200,0,1,200,0,1); // icosahedral
  TH1I *hhits = new TH1I("hhits","PMT hit count",50,0,3e5);
  BinLog(hhits->GetXaxis(),0.3);
  
  // ********************************************************************
  // Sum PMT hit counts for entire run
  // ********************************************************************
  int count=0, avgnhit=0, hitpmts=0, init=0;
  int pmthitcount[NPMTS];
  memset( pmthitcount, 0, NPMTS*sizeof(int) ); // NPMTS only known at runtime
  printf("Looping over events...\n");
  for(int iEntry=0; iEntry<dsreader.GetEntryCount(); iEntry++) {
    if (iEntry % 10000 == 0) printf("%d / %d\n",iEntry,(int)dsreader.GetEntryCount());
    const RAT::DS::Entry& ds = dsreader.GetEntry(iEntry);
    
    // Loop over triggered events in each entry
    for(int iEv=0; iEv<ds.GetEVCount(); iEv++) {        // mostly 1 event per entry
      const RAT::DS::EV& ev = ds.GetEV(iEv);
      // Trigger type
      int trig = ev.GetTrigType();
      if (!(trig & 0x8000)) continue;                   // EXTA trigger only
      // Event observables
      int nhits = ev.GetNhits();			            // normal/inward looking PMT hits
      avgnhit += nhits;
      count++;
      
      // PMT information
      const RAT::DS::UncalPMTs& pmts = ev.GetUncalPMTs();
      for(int iPMT=0; iPMT<pmts.GetCount(); iPMT++) {
        int pmtID = pmts.GetPMT(iPMT).GetID();
        pmthitcount[pmtID]++;
      } // pmt loop
    } // event loop
  } // entry loop
  printf("Finished loop.\n");
  
  // Find threshold for "screamers"
  int HOTLIMIT;
  GetHotLimit(pmthitcount, NPMTS, HOTLIMIT);
  for(int id=0; id<NPMTS; id++) {
    if(pmthitcount[id]==0) hhits->Fill(0.5); // visible on log scale
    else hhits->Fill(pmthitcount[id]);
  }
  printf("Found screamer threshold: %d\n", HOTLIMIT);
  
  // ********************************************************************
  // Get icosahedral projection of detector
  // ********************************************************************
  // Use icosahedral projection functions from HelperFunc.C
  using namespace func;
  TGraph2D *occupy = new TGraph2D();
  TVector2 icospos;
  TVector3 pmtpos;
  int goodpmts=0;
  printf("Generating icosahedral projection... ");
  for(int id=0; id<NPMTS; id++) {
    pmtpos = pmtinfo.GetPosition(id);
    if (pmtpos.Mag()==0) continue; // not a valid PMT
    if (pmtinfo.GetType(id) != 1) continue; // not a normal PMT
    int face; // side of icosahedron containing PMT
    icospos = func::IcosProject(pmtpos,face);
    //for(int j=0; j<pmthitcount[id]; j++) hicos->Fill(icospos.X(),icospos.Y());
    int step = (int)TMath::Floor(pmthitcount[id]/1.e3);
    if (pmthitcount[id] > HOTLIMIT) step=0; // hot PMT
    occupy->SetPoint(goodpmts,icospos.X(),icospos.Y(),step);
    goodpmts++;
  }
  printf("done.\n");
  
  // ********************************************************************
  // Plotting section
  // ********************************************************************
  TCanvas *c0 = new TCanvas("","",1600,800);
  //gStyle->SetOptStat(111111);
  //gStyle->SetTitleOffset(1.4,"z");
  gStyle->SetPadRightMargin(0.15); // for TH2 color scale
  gStyle->SetPadLeftMargin(0.05);
  gStyle->SetPadTopMargin(0.05);
  gStyle->SetPadBottomMargin(0.05);
  
  /*
  TPad *pad0 = new TPad("pad0","Text",0.01,0.75,0.50,0.99);
  TPad *pad1 = new TPad("pad1","Hist",0.50,0.75,1.00,0.99);
  TPad *pad2 = new TPad("pad2","Icos",0.05,0.41,0.95,0.74);
  TPad *pad3 = new TPad("pad3","Hist",0.01,0.01,0.50,0.40);
  TPad *pad4 = new TPad("pad4","Hist",0.50,0.01,1.00,0.40);
  pad0->Draw();
  pad1->Draw();
  pad2->Draw();
  pad3->Draw();
  pad4->Draw();
  */
  
  // PMT hit count histogram
  c0->cd()->SetGrid();
  occupy->SetMarkerStyle(7);
  //hicos->Draw(); // can't draw points over same plot?
  occupy->Draw("pcolz,ah,fb,bb");
  gPad->SetLineColor(0);
  gPad->SetTheta(90);
  gPad->SetPhi(0);
  gPad->Modified();
  gPad->Update();
 
/*
  occupy->GetXaxis()->SetTickLength(0);
  occupy->GetXaxis()->SetLabelOffset(999);
  occupy->GetYaxis()->SetTickLength(0);
  occupy->GetYaxis()->SetLabelOffset(999);
 */
  // Save canvas and close
  string outfile;
  if (!TEST) outfile = Form("output/PCA_%s",fibre.c_str());
  else       outfile = "occupancy";
  c0->Print(Form("%s.png",outfile.c_str()));
  c0->Print(Form("%s.pdf",outfile.c_str()));
  c0->Close();
  
  // Delete pointers created with 'new'
  if(occupy) delete occupy;
  if(hicos) delete hicos;
  if(hhits) delete hhits;
  
}
