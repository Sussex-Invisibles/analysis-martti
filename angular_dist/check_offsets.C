// ---------------------------------------------------------
// Goal:            Check the unusual PMT timing offsets
// Author:          Martti Nirkko, University of Sussex (Oct 2017)
// Compile & run:   clear && g++ -o check_offsets.exe check_offsets.C `root-config --cflags --libs` -I$RATROOT/include/libpq -I$RATROOT/include -L$RATROOT/lib -lRATEvent_Linux && ./check_offsets.exe
// ---------------------------------------------------------

// C++ stuff
#include <fstream>
#include <iostream>
#include <stdio.h>
#include <stdlib.h>

// ROOT stuff
#include <TCanvas.h>
#include <TH1.h>
#include <TH2D.h>
#include <TMath.h>
#include <TROOT.h>
#include <TStyle.h>
#include <TVector3.h>

// RAT stuff
#include <RAT/BitManip.hh>
#include "RAT/DU/DSReader.hh"
#include <RAT/DU/PMTInfo.hh>
#include <RAT/DU/Utility.hh>

using namespace std;

// Run time parameters
const int RUN_CLUSTER = 1;  // whether running on cluster (0=local)
const int VERBOSE = 1;      // verbosity flag
const int IS_MC = 0;        // Monte-Carlo flag

// Global constants
const double pi = TMath::Pi();

// Initialise functions
int check_offsets(string, int, TH2D**, const RAT::DU::PMTInfo&, bool, bool);

// Main program
int main(int argc, char** argv) {
  
  // Test flag (0 = process all runs, else run number)
  const int TEST = (RUN_CLUSTER) ? 0 : 102315;// 101410;
  
  // Initialisation
  string example = Form("%s/Software/SNOP/work/data/Analysis_r0000101834_s000_p000.root",getenv("HOME"));
  RAT::DU::DSReader dsreader(example);
  const RAT::DU::PMTInfo& pmtinfo = RAT::DU::Utility::Get()->GetPMTInfo();
  const int NPMTS = pmtinfo.GetCount();
  int NBINS     = 80;
  int MAXOFFSET = 20; // ns
  int MAXHEIGHT = 10; // m
  const int NHIST = 4;
  TH2D *histos[NHIST];
  histos[0] = new TH2D("pmtids","",NBINS,0,NPMTS+1,NBINS,-MAXOFFSET,MAXOFFSET);
  histos[1] = new TH2D("crates","",20,0,20,NBINS+1,-MAXOFFSET,MAXOFFSET);
  histos[2] = new TH2D("height","",NBINS,-MAXHEIGHT,MAXHEIGHT,NBINS,-MAXOFFSET,MAXOFFSET);
  histos[3] = new TH2D("angle","",NBINS,-1,1,NBINS,-MAXOFFSET,MAXOFFSET);

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
  int nfiles=0;
  while (true) {
    in >> node >> fibre >> channel >> run >> ipw >> photons >> pin >> rms >> nhit;
    if (!in.good()) break;
    if (TEST && TEST!=run) continue; // only want specified run
    int result = check_offsets(fibre, run, histos, pmtinfo, IS_MC, TEST);
    if (result==-1) { printf("Could not find data for fibre %s!\n",fibre.c_str()); continue; }
    printf("%6s run %6d has %2d PMTs with unusual hit time offsets (>3*RMS)\n",fibre.c_str(),run,result);
    nfiles++;
  }
  
  // Check output
  printf("Ran over %d files.\n",nfiles);
  if (nfiles==0) { 
    cerr<<"*** ERROR *** No input files found."<<endl;
    return 1; 
  }
  
  // Histogram titles
  histos[0]->SetTitle(Form("Unusually high PMT timing offsets (%d PCA runs);PMT ID;Offset [ns]",nfiles));
  histos[1]->SetTitle("Unusual offsets VS crates;Crate;Offset [ns]");
  histos[2]->SetTitle("Unusual offsets VS position;Z [m];Offset [ns]");
  histos[3]->SetTitle("Unusual offsets VS angle;#phi [#pi];Offset [ns]");
  
  // Plot output
  gStyle->SetOptStat(0);
  TCanvas *c = new TCanvas("c","",1200,900);
  c->Divide(2,2);
  for (int p=0; p<NHIST; p++) {
    c->cd(p+1)->SetGrid();
    c->cd(p+1)->SetLogz();
    histos[p]->GetZaxis()->SetRangeUser(0.5, histos[p]->GetMaximum());
    histos[p]->Draw("colz");
  }
  c->Print("check_offsets.png");
  c->Print("check_offsets.pdf");
  c->Close();
  
  // Free memory
  if (c) delete c;
  for (int p=0; p<NHIST; p++) { if (histos[p]) delete histos[p]; }
  return 0;
}

// Define macro
int check_offsets(string fibre, int run, TH2D **output, const RAT::DU::PMTInfo& pmtinfo, bool isMC=false, bool TEST=false) {
  
  // Open file
  string input = Form("logs/unusual_timing_%d.log",run);
  ifstream in(input.c_str());
  if (!in) return -1;
  string line;
  for (int hdr=0; hdr<2; hdr++) {
    getline(in,line);      // header
  }
  
  // Read file
  int count=0;
  int pmtid, crate;
  float offset;
  TVector3 pmtpos;
  while (true) {
    in >> pmtid >> offset;
    if (!in.good()) break;
    crate = RAT::BitManip::GetCrate(pmtid);
    pmtpos = pmtinfo.GetPosition(pmtid);
    output[0]->Fill(pmtid, offset);
    output[1]->Fill(crate, offset);
    output[2]->Fill(pmtpos.Z()/1e3, offset);
    output[3]->Fill(pmtpos.Phi()/pi, offset);
    count++;
  }
  return count;
  
}

