// ---------------------------------------------------------
// Goal:            Plot angular response of TELLIE fibres
// Author:          Martti Nirkko, 26/04/2017
// Compile & run:   clear && g++ -o summary.exe summary.C `root-config --cflags --libs` -I$RATROOT/include/libpq -I$RATROOT/include -L$RATROOT/lib -lRATEvent_Linux && ./summary.exe
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
const int RUN_CLUSTER = 1;  // whether running on cluster (0=local)
const int VERBOSE = 0;      // verbosity flag
const int IS_MC = 0;        // Monte-Carlo flag 

// Initialise functions
TH1D* get_histo(string, int, bool, bool);

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
  
  TH1D* temp = NULL;
  TH1D* hist[95] = {NULL};
  int nfiles=0;
  while (true) {
    in >> node >> fibre >> channel >> run >> ipw >> photons >> pin >> rms >> nhit;
    if (!in.good()) break;
    if (TEST && TEST!=run) continue; // only want specified run
    temp = get_histo(fibre, run, IS_MC, TEST);
    if (!temp) continue;
    /* No longer need comparison to mean value
    float meanval=0.;
    int nonzerobins=0;
    for (int bin=1; bin<=temp->GetNbinsX()+1; bin++) { // no underflow/overflow bins
      if (temp->GetBinContent(bin)==0) continue;
      meanval += temp->GetBinContent(bin);
      nonzerobins++;
    }
    meanval /= nonzerobins;
    for (int bin=0; bin<temp->GetNbinsX()+2; bin++) {
      temp->SetBinContent(bin, temp->GetBinContent(bin)-meanval);
    }
    */
    hist[nfiles] = temp;
    nfiles++;
  }
  printf("Ran over %d files.\n",nfiles);
  if (nfiles==0) { 
    cerr<<"*** ERROR *** No input files found."<<endl;
    return 1; 
  }

  // Plotting options
  gStyle->SetOptStat(0);
  gStyle->SetPadLeftMargin(0.1);
  gStyle->SetPadRightMargin(0.1);
  gStyle->SetPadTopMargin(0.1);
  gStyle->SetPadBottomMargin(0.1);
  gStyle->SetPadBorderSize(0);
   
  // Draw summary plot (full range)
  TCanvas *c0 = new TCanvas("","",1200,900);
  c0->SetGrid();
  c0->DrawFrame(0,-5,30,35,"TELLIE angular systematic;Angle w.r.t. fit position [deg];Mean time for each fibre [ns]");
  for (int fib=0; fib<95; fib++) {
    if (!hist[fib]) continue;
    hist[fib]->SetLineWidth(1);
    hist[fib]->SetLineColor(100-fib);
    hist[fib]->Draw("hist same L");
  }
  string imgfile = "summary";
  c0->Print(Form("%s.png",imgfile.c_str()));
  c0->Print(Form("%s.pdf",imgfile.c_str()));
  c0->Close();
 
  if (c0) delete c0;
  if (hist) for (int fib=0; fib<95; fib++) delete hist[fib];
  return 0;
   
  // Draw summary plot (zoomed in)
  c0 = new TCanvas("","",1200,900);
  c0->SetGrid();
  c0->DrawFrame(0,0,30,10,"TELLIE angular systematic;Angle w.r.t. fit position [deg];Mean time offset for each fibre [ns]");
  for (int fib=0; fib<95; fib++) {
    if (!hist[fib]) continue;
    hist[fib]->SetLineWidth(1);
    hist[fib]->SetLineColor(100-fib);
    hist[fib]->Draw("hist same L");
  }
  imgfile = "summary_zoom";
  c0->Print(Form("%s.png",imgfile.c_str()));
  c0->Print(Form("%s.pdf",imgfile.c_str()));
  c0->Close();
  
  // Free memory
  if (c0) delete c0;
  if (hist) for (int fib=0; fib<95; fib++) delete hist[fib];
  return 0;
}

// Get histogram from file created by "angular.C"
TH1D* get_histo(string fibre, int run, bool MC=false, bool TEST=false) {
  printf("Checking files for run %d... ", run);
  string in = Form("./output/angular_%s.root",fibre.c_str());
  ifstream f(in.c_str());
  if(!f.good()) {    // file not processed
    printf("not processed! Skipping fibre %s.\n",fibre.c_str());
    return NULL;
  } else {           // file processed
    printf("OK. Extracting data for fibre %s.\n",fibre.c_str());
  }
  TFile* infile = new TFile(in.c_str(),"READ");
  TH1D* hist = (TH1D*)infile->Get("h2_1");  // Mean time
  return hist;
}
