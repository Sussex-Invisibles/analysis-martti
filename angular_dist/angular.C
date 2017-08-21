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
  const int TEST = (RUN_CLUSTER) ? 0 : 102315;// 101410;

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
  //string example = (RUN_CLUSTER) ? "/lustre/scratch/epp/neutrino/snoplus/TELLIE_PCA_RUNS_PROCESSED/Analysis_r0000102315_s000_p000.root" : "/home/m/mn/mn372/Desktop/Analysis_r0000102315_s000_p000.root";
  //RAT::DU::DSReader dsreader(example);
  //const RAT::DU::PMTInfo& pmtinfo = RAT::DU::Utility::Get()->GetPMTInfo();
  //const int NPMTS = pmtinfo.GetCount();
  
  int nfiles=0; 
  while (true) {
    in >> node >> fibre >> channel >> run >> ipw >> photons >> pin >> rms >> nhit;
    if (!in.good()) break;
    if (TEST && TEST!=run) continue; // only want specified run
    float avgnhit = angular(fibre, run, IS_MC, TEST);
    if (avgnhit==-2.) { printf("Could not find data for fibre %s!\n",fibre.c_str()); continue; }
    if (avgnhit==-1.) { printf("Nothing to do!\n"); continue; }
    printf("%3d %8s %8.1f nhit (extracted) %8.1f nhit (table)\n",channel,fibre.c_str(),avgnhit,nhit);
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
  //string fpath = (RUN_CLUSTER) ? "/lustre/scratch/epp/neutrino/snoplus/TELLIE_PCA_RUNS_PROCESSED" : "/home/m/mn/mn372/Desktop";
  //string fpath = (RUN_CLUSTER) ? "/home/m/mn/mn372/Software/SNOP/work/data" : "/home/nirkko/Software/SNOP/work/data";
  string fpath = Form("%s/Software/SNOP/work/data",getenv("HOME"));
  string fname = "";
  ifstream f;
  for (int pass=3;pass>=0;pass--) {
    fname = Form("%s/Analysis_r0000%d_s000_p00%d.root",fpath.c_str(),run,pass);
    f.open(fname.c_str());
    if (f.good()) break;
  }
  string out = Form("./output/angular_%s.out",fibre.c_str());
  string img = Form("./images/angular_%s.pdf",fibre.c_str());
  ifstream g(out.c_str());
  ifstream h(img.c_str());
  int scanned_file = 0;
  if (!TEST && h.good()) {  // file downloaded and processed
    printf("already processed! Skipping fibre %s.\n",fibre.c_str());
    return -1.;
  } else if(g.good()) {     // file extracted, but not processed
    printf("not processed! Generating plots for fibre %s.\n",fibre.c_str());
    scanned_file = 1;
  } else if(!f.good()) {    // file not downloaded
    printf("not downloaded! Skipping fibre %s.\n",fibre.c_str());
    return -2.;
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
  const int NBINS = 30;
  const int MAXANG = 30;
  TH1D *hcoarse = new TH1D("hcoarse",fibre.c_str(),500,0,500);
  TH1D *hpmtseg = new TH1D("hpmtseg",fibre.c_str(),NBINS,0,MAXANG);
  
  // First iteration: Get average hit time and number of PMTs in each segment
  int hitpmts=0, allpmts[NPMTS];
  memset( allpmts, 0, NPMTS*sizeof(int) ); // NPMTS only known at runtime
  for(int iEntry=0; iEntry<dsreader.GetEntryCount(); iEntry++) {
    const RAT::DS::Entry& ds = dsreader.GetEntry(iEntry);
    for(int iEv=0; iEv<ds.GetEVCount(); iEv++) {        // mostly 1 event per entry
      const RAT::DS::EV& ev = ds.GetEV(iEv);
      int trig = ev.GetTrigType();
      if (!(trig & 0x8000)) continue;                   // EXT trigger only
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
        double corr = track.Mag()/c_water;                  // light travel time [ns]
        hcoarse->Fill(pmttime-corr);
        if (allpmts[pmtID]==0) { hpmtseg->Fill(theta*180./pi); allpmts[pmtID]++; }
        hitpmts++;
      } // pmt loop
    } // event loop
  } // entry loop
  
  int commontime = hcoarse->GetXaxis()->GetBinLowEdge(hcoarse->GetMaximumBin());
  commontime += 5/2;              // intermediate step
  commontime -= commontime % 5;   // rounded to the nearest multiple of 5
  if (VERBOSE) printf("Most common hit time: %d ns\n",commontime);
  
  // More histograms
  TH1D *h1 = new TH1D("h1",fibre.c_str(),NBINS,0,MAXANG);
  TH2D *h2 = new TH2D("h2",fibre.c_str(),NBINS,0,MAXANG,NBINS,commontime-15,commontime+15);
  TH1D *herr = new TH1D("herr",fibre.c_str(),NBINS,0,MAXANG);
  
  // Second iteration: Fill histograms
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
        double corr = track.Mag()/c_water;                  // light travel time [ns]
        h2->Fill(theta*180./pi, pmttime-corr);              // angle [deg], time [ns]
      } // pmt loop
    } // event loop
  } // entry loop
  
  // Fit angular slices with a gaussian distribution in time
  h2->FitSlicesY();
  TH1D *h2_0 = (TH1D*)gDirectory->Get("h2_0");  // Constant
  TH1D *h2_1 = (TH1D*)gDirectory->Get("h2_1");  // Mean
  TH1D *h2_2 = (TH1D*)gDirectory->Get("h2_2");  // StdDev
  TH1D *h2_3 = (TH1D*)gDirectory->Get("h2_chi2");  // chi^2/ndof  
  TH1D *h2_py = NULL;

  // Put mean and RMS into 1D histogram
  int binscone=0;
  float avgmean=0., avgdev=0.;
  for (int b=0; b<NBINS+2; b++) {
    h1->SetBinContent(b, h2_1->GetBinContent(b));
    h1->SetBinError(b, h2_2->GetBinContent(b));
    h2->ProjectionY("_py",b,b);
    h2_py = (TH1D*)gDirectory->Get("h2_py");
    h2_0->SetBinContent(b, h2_py->GetSum());  // overwrite Gaussian constant with intensity
    h2_1->SetBinContent(b, h2_py->GetMean()); // overwrite Gaussian mean with mean
    h2_2->SetBinContent(b, h2_py->GetRMS());  // overwrite Gaussian sigma with RMS
    herr->SetBinContent(b, h2_py->GetMeanError()); // not *quite* the same as error on Gaussian mean
    //if (VERBOSE) printf("#%2d %10.6f %10.6f\n",b,h2_1->GetBinError(b),h2_py->GetMeanError());
    //if (VERBOSE) printf("#%2d %.4f %6.2f %6.2f\n",b,h1->GetBinCenter(b),h1->GetBinContent(b),h1->GetBinError(b));
    if (h1->GetBinCenter(b)>0 && h1->GetBinCenter(b)<30) {  // within 30 deg cone
      avgmean += h2_1->GetBinContent(b);
      avgdev += h2_2->GetBinContent(b);
      binscone++;
    }
  }
  avgmean/=binscone;
  avgdev/=binscone;
  if (VERBOSE) printf("Average hit time mean (< 30 deg): %6.2f\n", avgmean);
  if (VERBOSE) printf("Average hit time RMS (< 30 deg):  %6.2f\n", avgdev);
  
  // Normalise angular slices to number of PMTs in slice
  double tmp, pmts;
  for (int binx=0; binx<NBINS+2; binx++) {
    pmts = hpmtseg->GetBinContent(binx);
    tmp = h2_0->GetBinContent(binx);
    if (pmts<=0) h2_0->SetBinContent(binx,0);
    else h2_0->SetBinContent(binx,tmp/pmts);
    //if (VERBOSE) printf("Segment %2d contains %4d PMTs.\n",binx,(int)pmts);
    for (int biny=0; biny<NBINS+2; biny++) {
      tmp = h2->GetBinContent(binx,biny);
      if (pmts<=0) h2->SetBinContent(binx,biny,0);
      else h2->SetBinContent(binx,biny,tmp/pmts);
    }
  }
  
  // Plotting options
  gStyle->SetOptStat(0);
  gStyle->SetPadLeftMargin(0.1);
  gStyle->SetPadRightMargin(0.1);
  gStyle->SetPadTopMargin(0.1);
  gStyle->SetPadBottomMargin(0.1);
  gStyle->SetPadBorderSize(0);
   
  // Define canvas and pads
  TCanvas *c0 = new TCanvas("","",1200,900);
  TPad *pad0 = new TPad("pad0","Hist",0.01,0.35,0.48,0.99);
  TPad *pad1 = new TPad("pad1","Hist",0.52,0.35,1.00,0.99);
  TPad *pad2 = new TPad("pad2","Hist",0.01,0.01,0.25,0.33);
  TPad *pad3 = new TPad("pad3","Hist",0.26,0.01,0.50,0.33);
  TPad *pad4 = new TPad("pad4","Hist",0.51,0.01,0.75,0.33);
  TPad *pad5 = new TPad("pad5","Hist",0.76,0.01,1.00,0.33);
  pad0->Draw();
  pad1->Draw();
  pad2->Draw();
  pad3->Draw();
  pad4->Draw();
  pad5->Draw();
  
  // Time vs. angle (2D)
  pad0->cd()->SetGrid();
  pad0->SetRightMargin(0.12);   // for TH2D color scale
  h2->SetTitle(Form("Intensity profile %s (norm.);Angle [deg];Time [ns]",fibre.c_str()));
  h2->Draw("colz");
  h2->GetXaxis()->SetTitleOffset(1.2);
  h2->GetYaxis()->SetTitleOffset(1.5);
  
  // Mean
  pad1->cd()->SetGrid();
  pad1->SetLeftMargin(0.12);
  h2_1->SetTitle("Mean hit time;Angle [deg];Time [ns]");
  h2_1->SetLineWidth(2);
  h2_1->Draw();
  h2_1->GetXaxis()->SetTitleOffset(1.2);
  h2_1->GetYaxis()->SetTitleOffset(1.7);
  h2_1->GetYaxis()->SetRangeUser(0.995*h2_1->GetMinimum(), 1.005*h2_1->GetMaximum());
  
  // Sum
  pad2->cd()->SetGrid();
  pad2->SetLeftMargin(0.2);   // for axis label
  h2_0->SetTitle("Intensity (norm.);Angle [deg];Avg. NHit/PMT [ ]");
  h2_0->SetLineWidth(2);
  h2_0->Draw();
  h2_0->GetXaxis()->SetTitleOffset(1.2);
  h2_0->GetYaxis()->SetTitleOffset(1.7);
  h2_0->GetYaxis()->SetRangeUser(0, 1.1*h2_0->GetMaximum());
  hpmtseg->Scale(10);
  hpmtseg->SetLineWidth(2);
  hpmtseg->SetLineColor(2);
  hpmtseg->SetTitle("NPMTs (#times10)");
  hpmtseg->Draw("same");
  pad2->BuildLegend();
  
  // RMS
  pad3->cd()->SetGrid();
  h2_2->SetTitle("RMS hit time;Angle [deg];Time [ns]");
  h2_2->SetLineWidth(2);
  h2_2->Draw();
  h2_2->GetXaxis()->SetTitleOffset(1.2);
  h2_2->GetYaxis()->SetTitleOffset(1.4);
  h2_2->GetYaxis()->SetRangeUser(0.95*h2_2->GetMinimum(), 1.05*h2_2->GetMaximum());
 
  // Error on mean
  pad4->cd()->SetGrid();
  pad1->SetLeftMargin(0.2);   // for axis label
  herr->SetTitle("Error on mean hit time;Angle [deg];#Delta t_{mean} [ns]");
  herr->SetLineWidth(2);
  herr->Draw("colz");
  herr->GetXaxis()->SetTitleOffset(1.2);
  herr->GetYaxis()->SetTitleOffset(1.7);
  herr->GetYaxis()->SetRangeUser(0,0.1);

  // Chisquare/NDOF
  pad5->cd()->SetGrid();
  h2_3->SetTitle("Gaussian fit #chi^{2}/ndof;Angle [deg];");
  h2_3->SetLineWidth(2);
  h2_3->Draw();
  h2_3->GetXaxis()->SetTitleOffset(1.2);
  h2_3->GetYaxis()->SetRangeUser(0, 1.1*h2_3->GetMaximum());

  // Save canvas and close
  string outfile = "angular";
  if (!TEST) outfile = Form("images/angular_%s",fibre.c_str());
  c0->Print(Form("%s.png",outfile.c_str()));
  c0->Print(Form("%s.pdf",outfile.c_str()));
  c0->Close();
  
  if (c0) delete c0;
  if (h1) delete h1;
  if (h2) delete h2;
  if (herr) delete herr;
  if (hcoarse) delete hcoarse;
  if (hpmtseg) delete hpmtseg;
  
  return (float)totalnhit/count;
}

