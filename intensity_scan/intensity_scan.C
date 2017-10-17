// ---------------------------------------------------------
// Goal:            Check TELLIE PIN readings over subruns
// Author:          Martti Nirkko, 07/08/2017
// Compile & run:   clear && g++ -g -o intensity_scan.exe intensity_scan.C `root-config --cflags --libs` -I$RATROOT/include/libpq -I$RATROOT/include -L$RATROOT/lib -lRATEvent_Linux && ./intensity_scan.exe
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

// Global constants
const int RUN_CLUSTER = 1;  // whether running on cluster (0=local)
const int VERBOSE = 1;      // verbosity flag
const int IS_MC = 0;        // Monte-Carlo flag 
const int NRUNS = 16;       // Runs in total
const int NSUBRUNS = 10;    // Subruns per run

// Initialise functions
void get_pin_values(int, float*, float*, bool, bool);
void get_subrun_nhit(int, int*, int*, int*, bool, bool);

// ----------------------------------------------------------------------------

// Main program
using namespace std;
int main(int argc, char** argv) {
  
  // Test flag (0 = process all runs, else run number)
  const int TEST = (RUN_CLUSTER) ? 0 : 104157;

  // Loop over all fibres in list
  string input = "TELLIE_intensity_runs.list";
  ifstream in(input.c_str());
  if (!in) { cerr<<"Failed to open "<<input<<endl; exit(1); }
  string output = "triggers_missing.txt";
  FILE* outfile = fopen(output.c_str(), "w");  // (over)write
  fclose(outfile); // emptied file contents
  //TH2F *nhits = new TH2F("nhits","NHit/event vs PIN reading",100,0,2000,100,0,200);
  TGraph *g[NRUNS] = {NULL};
  int run, nruns=0, nsubruns=0;
  int col[NRUNS]={0};
  float trig[3][NRUNS]={0};
  float meannhit[NRUNS]={0};
  string fibre[NRUNS] = {""};
  while (true) {
    in >> run;
    if (!in.good()) break;
    if (TEST && TEST!=run) continue;
    if (run>=104144 && run<104158) { col[nruns]=3; fibre[nruns]="FT019A"; }
    if (run>=104158 && run<104170) { col[nruns]=2; fibre[nruns]="FT034A"; }
    if (run>=104170 && run<104177) { col[nruns]=4; fibre[nruns]="FT010A"; }
    int subrunevs[NSUBRUNS]={0}, subrunhits[NSUBRUNS]={0}, triggers[3]={0};
    float subrunpins[NSUBRUNS]={0}, subrunrms[NSUBRUNS]={0}, avgnhit[NSUBRUNS]={0};
    get_pin_values(run, subrunpins, subrunrms, (bool)IS_MC, (bool)TEST);
    get_subrun_nhit(run, subrunevs, subrunhits, triggers, (bool)IS_MC, (bool)TEST);
    for (int sr=0; sr<NSUBRUNS; sr++) {
      if(subrunevs[sr]==0) continue;
      avgnhit[sr] = (float)subrunhits[sr]/subrunevs[sr];
      meannhit[nruns] += avgnhit[sr]/NSUBRUNS;
      if (VERBOSE) printf("Run %6d subrun %2d has %4d events, PIN=%4d, RMS=%6.2f, Nhit=%5.2f.\n", run, sr, subrunevs[sr],(int)subrunpins[sr],subrunrms[sr],avgnhit[sr]);
      nsubruns++;
    }
    if (VERBOSE) printf("Run %6d has %5d triggers (%5d MISS, %5d EXTA).\n",run,triggers[0],triggers[1],triggers[2]);
    for (int i=0; i<3; i++) trig[i][nruns] = (float)triggers[i]/1e3;
    g[nruns] = new TGraph(NSUBRUNS, subrunpins, avgnhit);
    nruns++;
  }
  TGraph *gall = new TGraph(nruns, meannhit, trig[0]);
  TGraph *gmiss = new TGraph(nruns, meannhit, trig[1]);
  TGraph *gexta = new TGraph(nruns, meannhit, trig[2]);
  printf("Ran over %d runs and %d subruns.\n",nruns,nsubruns);
  if (nruns==0 || nsubruns==0) {
    cerr<<"*** ERROR *** No input files found."<<endl;
    return 1; 
  }

  // Plot graph
  TCanvas *c = new TCanvas("c","",800,600);
  c->SetGrid();
  c->DrawFrame(400,0,1200,150,"TELLIE intensity sweep;PIN reading;Average NHit");
  int lastcol=0;
  TLegend *leg = new TLegend(0.71,0.71,0.89,0.89);
  for (int runs=0; runs<nruns; runs++) {
    if (!g[runs]) continue;
    g[runs]->SetMarkerStyle(7);
    g[runs]->SetMarkerColor(col[runs]);
    g[runs]->Draw("P same");
    if(lastcol!=col[runs]) { leg->AddEntry(g[runs],fibre[runs].c_str(),"p"); lastcol=col[runs]; }
  }
  leg->Draw("same");
  c->Print("intensity_scan.pdf");
  c->Print("intensity_scan.png");
  c->Close();

  c = new TCanvas("c","",800,600);
  c->SetGrid();
  c->DrawFrame(0,0,150,60,"TELLIE triggers;Average NHit;Trigger count (10^{3})");
  // Sort graphs along x-axis
  gall->Sort();
  gmiss->Sort();
  gexta->Sort();
  // Events with any trigger bit
  gall->SetMarkerStyle(8);
  gall->SetMarkerColor(4);
  gall->SetLineWidth(2);
  gall->SetLineColor(4);
  gall->Draw("LP same");
  // Events with MISS trigger bit
  gmiss->SetMarkerStyle(8);
  gmiss->SetMarkerColor(2);
  gmiss->SetLineWidth(2);
  gmiss->SetLineColor(2);
  gmiss->Draw("LP same");
  // Events with EXTA trigger bit
  gexta->SetMarkerStyle(8);
  gexta->SetMarkerColor(3);
  gexta->SetLineWidth(2);
  gexta->SetLineColor(3);
  gexta->Draw("LP same");
  // Legend
  leg = new TLegend(0.13,0.13,0.29,0.29);
  leg->AddEntry(gall, "ALL", "p");
  leg->AddEntry(gmiss, "MISS", "p");
  leg->AddEntry(gexta, "EXTA", "p");
  leg->Draw("same");
  c->Print("triggers.pdf");
  c->Print("triggers.png");
  c->Close();
  
  // End of program
  for (int runs=0; runs<NRUNS; runs++) { if(g[runs]) delete g[runs]; }
  if (c) delete c;
  return 0;
}

// ----------------------------------------------------------------------------

// Returns the PIN readings for a given run
void get_pin_values(int run, float *srpin, float *srrms, bool isMC=false, bool TEST=false) {
  string input = "PIN_readings.txt";
  ifstream in(input.c_str());
  if (!in) { cerr<<"Failed to open "<<input<<endl; exit(1); }
  string line;
  int check, subrun, pin;
  float rms;
  for (int hdr=0; hdr<2; hdr++) {
    getline(in,line);      // header
  }
  while (true) {
    in >> check >> subrun >> pin >> rms;
    if (!in.good()) break;
    if (check!=run) continue; // only want specified run
    srpin[subrun] = (float)pin;
    srrms[subrun] = rms;
  }
}

// ----------------------------------------------------------------------------

// Returns the fitted light position for a given fibre/run
void get_subrun_nhit(int run, int *srevts, int *srhits, int *triggers, bool isMC=false, bool TEST=false) {

  // Check files for given run
  if(!TEST) printf("*****\n");
  printf("Checking files for run %d... ", run);
  string fpath = Form("%s/Software/SNOP/work/data",getenv("HOME"));
  string fname = "";
  ifstream f;
  for (int pass=3;pass>=0;pass--) {
    fname = Form("%s/Analysis_r0000%d_s000_p00%d.root",fpath.c_str(),run,pass);
    f.open(fname.c_str());
    if (f.good()) break;
  }
  if(!f.good()) {    // file not downloaded
    printf("not downloaded! Skipping run %d.\n",run);
    return;
  } else {           // file downloaded, but not processed
    printf("OK. Extracting data for run %d.\n",run);
  }

  string output = "triggers_missing.txt";
  FILE* outfile = fopen(output.c_str(), "a");  // append
  
  // Initialise RAT
  RAT::DU::DSReader dsreader(fname);
  //const RAT::DU::ChanHWStatus& chs = RAT::DU::Utility::Get()->GetChanHWStatus();
  //const RAT::DU::PMTInfo& pmtinfo = RAT::DU::Utility::Get()->GetPMTInfo();
  //const int NPMTS = pmtinfo.GetCount();
  
  // Loop over entries
  int all=0, miss=0, exta=0;
  int subrunevtcount[NSUBRUNS] = {0};
  int subrunhitcount[NSUBRUNS] = {0};
  for(int iEntry=0; iEntry<dsreader.GetEntryCount(); iEntry++) {
    const RAT::DS::Entry& ds = dsreader.GetEntry(iEntry);
    int subrunid = ds.GetSubRunID();             // sub run this event belongs to

    // Loop over triggered events in each entry
    for(int iEv=0; iEv<ds.GetEVCount(); iEv++) { // mostly 1 event per entry
      const RAT::DS::EV& ev = ds.GetEV(iEv);
      int gtid = ev.GetGTID();
      
      // Trigger type
      all++;
      int trig = ev.GetTrigType();
      if (trig & 0x04000000) {       // MISS trigger bit set
        miss++;
      }
      if (!(trig & 0x8000)) {        // EXTA trigger bit not set
        if (trig != 64 && trig != 1024)
          fprintf(outfile,"Run %d, subrun %d, GTID %d has trigger mask %d\n",run,subrunid,gtid,trig);
        continue;
      }
      exta++; 

      // Event observables
      int nhitscleaned = ev.GetNhitsCleaned();   // calibrated PMT hits with removed crosstalk
      subrunhitcount[subrunid-1] += nhitscleaned;
      subrunevtcount[subrunid-1]++;
      
    } // event loop
  } // entry loop
  
  fclose(outfile);  

  for (int sr=0; sr<NSUBRUNS; sr++) {
    srevts[sr] = subrunevtcount[sr];
    srhits[sr] = subrunhitcount[sr];
  }
  triggers[0] = all;
  triggers[1] = miss;
  triggers[2] = exta;
  return;
}

// ----------------------------------------------------------------------------
