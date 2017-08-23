// ---------------------------------------------------------
// Goal:            Check TELLIE PIN readings over subruns
// Author:          Martti Nirkko, 07/08/2017
// Compile & run:   clear && g++ -g -o pin_readings.exe pin_readings.C `root-config --cflags --libs` -I$RATROOT/include/libpq -I$RATROOT/include -L$RATROOT/lib -lRATEvent_Linux && ./pin_readings.exe
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
//#include "HelperFunc.C"
//#include "Xianguo.C"

// Global constants
const int RUN_CLUSTER = 1;  // whether running on cluster (0=local)
const int VERBOSE = 0;      // verbosity flag
const int IS_MC = 0;        // Monte-Carlo flag 
const int NCOL = 20;        // number of colors (max. 50)
const int NDOTS = 360;      // number of points in circle
const double DIR_CONE = 48; // opening angle to search for direct light (using aperture: 24 deg)
const double REF_CONE = 20; // opening angle to search for reflected light (using aperture: 9.874 deg)

// Initialise functions
void get_subrun_nhit(int, int*, int*, bool, bool);

// Main program
using namespace std;
int main(int argc, char** argv) {
  
  // Test flag (0 = process all runs, else run number)
  const int TEST = (RUN_CLUSTER) ? 0 : 100801;

  // Loop over all fibres in list
  string input1 = "TELLIE_PCA_set1.txt";
  ifstream in1(input1.c_str());
  if (!in1) { cerr<<"Failed to open TELLIE_PCA_set1.txt"<<endl; exit(1); }
  //TH2F *nhits = new TH2F("nhits","NHit/event vs PIN reading",100,0,1e4,100,0,100);
  string line;
  float subrunpin[40], subrunrms[40];
  float avgnhit[40];
  int run, subrun, pin;
  float rms;
  int nsubruns=0;
  for (int hdr=0; hdr<2; hdr++) {
    getline(in1,line);      // header
  }
  int lastrun = 0;
  while (true) {
    in1 >> run >> subrun >> pin >> rms;
    if (!in1.good()) break;
    if (TEST && TEST!=run) continue; // only want specified run
    subrunpin[subrun] = pin;
    subrunrms[subrun] = rms;
    if (run == lastrun) continue; // only do this once per run
    lastrun = run;
    //if (VERBOSE) printf("%6d %2d %5d %.2f\n", run, subrun, pin, rms);
    int subrunevs[40], subrunhits[40];
    get_subrun_nhit(run, subrunevs, subrunhits, (bool)IS_MC, (bool)TEST);
    for (int sr=0; sr<40; sr++) {
      if(subrunevs[sr]==0) continue;
      avgnhit[sr] = (float)subrunhits[sr]/subrunevs[sr];
      //if (VERBOSE) printf("Run %6d subrun %2d has %4d events with an average nhit of %5.2f.\n", run, sr, subrunevs[sr], avgnhit[sr]);
      nsubruns++;
    }
  }
  printf("Ran over %d subruns.\n",nsubruns);
  if (nsubruns==0) {
    cerr<<"*** ERROR *** No input files found."<<endl;
    return 1; 
  }

  // Plot graph
  TGraph *g = new TGraph(20,subrunpin,avgnhit);
  TCanvas *c = new TCanvas("c","",800,600);
  c->SetGrid();
  g->SetTitle("TELLIE run 100801;PIN reading;Average NHit");
  g->SetMarkerStyle(7);
  g->SetMarkerColor(2);
  g->Draw("AP");
  c->Print("pin_readings.pdf");
  c->Print("pin_readings.png");
  c->Close();

  // End of program
  if (g) delete g;
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
  string fpath = (RUN_CLUSTER) ? "/lustre/scratch/epp/neutrino/snoplus/TELLIE_PCA_RUNS_PROCESSED" : "/home/nirkko/Software/SNOP/data";
  string fname = "";
  ifstream f;
  for (int pass=3;pass>=0;pass--) {
    fname = Form("%s/Analysis_r0000%d_s000_p00%d.root",fpath.c_str(),run,pass);
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

  /*
  // Get fibre information (without using RATDB) 
  string fibre_table = "Fibre_Positions_DocDB1730.csv";
  ifstream tab(fibre_table.c_str());
  if (!tab) { cerr<<"Failed to open "<<fibre_table<<endl; exit(1); }
  string line, ff, fn;
  double fx,fy,fz,fu,fv,fw;
  getline(tab,line);      // header
  while (tab.good()) {
    tab >> ff >> fx >> fy >> fz >> fu >> fv >> fw >> fn;
    if (!tab.good()) break;
    if (ff==fibre) break;
  }
  if(ff!=fibre) { cerr << "Failed to find information for fibre " << fibre << endl; exit(1); }
  //cout << Form("%s %f %f %f %f %f %f %s\n",ff.c_str(),fx,fy,fz,fu,fv,fw,fn.c_str()) << endl;
  TVector3 fibrepos(10*fx,10*fy,10*fz); // position of fibre [mm]
  TVector3 fibredir(fu,fv,fw); // direction of fibre
  TVector3 lightpos = fibrepos + 2*fibrepos.Mag()*fibredir; // projected light spot centre 
  //cout << "DocDB: fibre " << fibre << ", pos " << printVector(fibrepos) << ", dir " << printVector(fibredir) << endl;
  */

  // ********************************************************************
  // Sum PMT hit counts for entire run
  // ********************************************************************
    int subrunevtcount[40] = {0};
    int subrunhitcount[40] = {0};
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
        subrunhitcount[subrunid-1] += nhitscleaned;
        subrunevtcount[subrunid-1]++;
        
      } // event loop
    } // entry loop
    
    for (int sr=0; sr<40; sr++) {
      srevts[sr] = subrunevtcount[sr];
      srhits[sr] = subrunhitcount[sr];
    }
    return;

  /*
    // Write PMT hit count to output file (saves time when rerunning)
    FILE *outFile = fopen(out.c_str(),"w");
    TVector3 pmtpos;
    fprintf(outFile, "Run: %d\n", run);             // run number
    fprintf(outFile, "Events: %d\n", count);        // events with EXTA trigger
    fprintf(outFile, "TotalNHit: %d\n", totalnhit); // calibrated PMT hits with removed crosstalk
    for(int id=0; id<NPMTS; id++) {
      pmtpos = pmtinfo.GetPosition(id);
      int in_spot = 0;
      if (pmtpos.Angle(lightpos) < DIR_CONE/180*pi) in_spot = 1;
      else if (pmtpos.Angle(fibrepos) < REF_CONE/180*pi) in_spot = 2;
      fprintf(outFile, "%d %d %d\n", id, pmthitcount[id], in_spot);
    }
    fclose(outFile);
    scanned_file = 1;
  */
   
  // ********************************************************************
  // Plotting section
  // ********************************************************************
  /*
  TCanvas *c0 = new TCanvas("","",1200,1500);
  gStyle->SetTitleOffset(1.2,"xyz");
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
  
  // Save canvas and close
  string outfile;
  if (!TEST) outfile = Form("images/PCA_%s",fibre.c_str());
  else       outfile = "pin_readings";
  c0->Print(Form("%s.png",outfile.c_str()));
  c0->Print(Form("%s.pdf",outfile.c_str()));
  c0->Close();
  */ 
  // Delete pointers created with 'new'
  //if(hhits) delete hhits;
  //if(c0) delete c0;
  
}
