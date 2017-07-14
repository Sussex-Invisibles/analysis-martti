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
#include "TH3.h"
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
const int RUN_CLUSTER=1;    // whether running on cluster (0=local)
const int IS_MC = 0;        // Monte-Carlo flag 

// Initialise functions
void occupancy(string, int, int*, int, bool, bool);
void GetProjection(TGraph2D*, TGraph2D*, int*, const int, const RAT::DU::PMTInfo&);
using namespace std;


// Main program
int main(int argc, char** argv) {
  
  // Test flag (0 = process all runs, else run number)
  const int TEST = (RUN_CLUSTER) ? 0 : 101849;
  //const int TEST = 102278;  

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
  
  string example = "/lustre/scratch/epp/neutrino/snoplus/TELLIE_PCA_RUNS_PROCESSED/Analysis_r0000102315_s000_p000.root";
  // Initialise RAT
  RAT::DU::DSReader dsreader(example);
  const RAT::DU::PMTInfo& pmtinfo = RAT::DU::Utility::Get()->GetPMTInfo();
  const int NPMTS = pmtinfo.GetCount();
  printf("Initialised DSReader & PMTInfo.\n");
  
  int cvrg[NPMTS];
  memset( cvrg, 0, NPMTS*sizeof(int) ); // NPMTS only known at runtime
  for (int id=0; id<NPMTS; id++) cvrg[id]=-99999;
  
  while (true) {
    in >> node >> fibre >> channel >> run >> ipw >> photons >> nhit;
    if (!in.good()) break;
    if (TEST && TEST!=run) continue; // only want specified run
    //printf("%6s %2d %6d %5d %6d\n", fibre.c_str(), channel, run, ipw, photons);
    occupancy(fibre, run, cvrg, NPMTS, IS_MC, TEST);
  }
  TGraph2D *coverage = new TGraph2D();
  TGraph2D *offpmts = new TGraph2D();
  GetProjection(coverage, offpmts, cvrg, NPMTS, pmtinfo);
  printf("Graph has %d good PMTs and %d bad PMTs.\n", coverage->GetN(), offpmts->GetN());
  
  // ********************************************************************
  // Plotting section
  // ********************************************************************
  TCanvas *c0 = new TCanvas("","",1600,800);
  gStyle->SetOptStat(0);
  gStyle->SetLabelOffset(999,"XY");
  gStyle->SetCanvasBorderSize(0);
  gStyle->SetFrameBorderSize(0);
  
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
  
  // Configure pad
  c0->SetGrid();
  c0->SetLeftMargin(0.005);
  c0->SetRightMargin(0.075); // for TH2 color scale
  c0->SetTopMargin(0.02);
  c0->SetBottomMargin(0.02);
  c0->SetBorderSize(0);
  
  // Draw empty 3D histogram as boundary for 2D graphs
  TH3F *hicos = new TH3F("hicos","",10,0,1,10,0,1,10,0,10); // flat map
  hicos->Draw("a,fb,bb");     // suppress axis, front box, back box
  c0->SetTheta(90-0.001);     // view from above
  c0->SetPhi(0+0.001);        // no x-y rotation
  
  // Draw active PMT coverage
  coverage->SetMarkerStyle(7);
  //coverage->SetMinimum(0);                      // TODO - takes >10min (!?)
  //coverage->GetZaxis()->SetRangeUser(0,10);     // TODO - doesn't work
  coverage->Draw("pcolz,a,fb,bb,same");
  
  // Draw inactive PMTs (grey)
  offpmts->SetMarkerStyle(7);
  offpmts->SetMarkerColor(16);
  offpmts->Draw("p,a,fb,bb,same");
  
  // Save canvas and close
  string outfile = "coverage";
  if (TEST) outfile = Form("coverage_%d",TEST);
  c0->Print(Form("%s.png",outfile.c_str()));
  c0->Print(Form("%s.pdf",outfile.c_str()));
  c0->Close();
  
  // Delete pointers
  if(coverage) delete coverage;
  if(hicos) delete hicos;
  if(offpmts) delete offpmts;
  return 0;
}

// Returns the fitted light position for a given fibre/run
void occupancy(string fibre, int run, int* cvrg, int NPMTS, bool isMC=false, bool TEST=false) {

  // ********************************************************************
  // Initialisation
  // ********************************************************************
  
  // Check files for given run
  if(!TEST) printf("-----\n");
  printf("Checking files for run %d... ", run);
  string out = Form("./output/PCA_%s.out",fibre.c_str());
  ifstream g(out.c_str());
  if(g.good()) {            // file extracted
    printf("OK! Adding direct light for fibre %s.\n",fibre.c_str());
  } else {    // file not available
    printf("not available! Skipping fibre %s.\n",fibre.c_str());
    return;
  }
  
  // ********************************************************************
  // Sum PMT hit counts for entire run
  // ********************************************************************
  int pmthitcount[NPMTS], pmtlightcone[NPMTS];
  memset( pmthitcount, 0, NPMTS*sizeof(int) ); // NPMTS only known at runtime
  memset( pmtlightcone, 0, NPMTS*sizeof(int) ); // NPMTS only known at runtime
  int pmtid, pmthits, onspot;
  int checkrun, count, totalnhit;
  string dummy;
  g >> dummy >> checkrun;
  if(!g.good()) printf("*** ERROR *** Bad input - str=%s int=%d\n",dummy.c_str(),checkrun);
  if(checkrun != run) { printf("*** ERROR *** Bad run number %d\n",checkrun); return; }
  g >> dummy >> count;
  if(!g.good()) { printf("*** ERROR *** Bad input - str=%s int=%d\n",dummy.c_str(),count); return; }
  g >> dummy >> totalnhit;
  if(!g.good()) { printf("*** ERROR *** Bad input - str=%s int=%d\n",dummy.c_str(),totalnhit); return; }
  // Print run info here
  printf("*** INFO *** Run %d has %d EXTA events with %d total NHits.\n",checkrun,count,totalnhit);
  float avgnhit = (float)totalnhit/count;  
    
  while (g.good()) {
    g >> pmtid >> pmthits >> onspot;
    if(!g.good()) break;
    pmthitcount[pmtid]=pmthits;
    pmtlightcone[pmtid]=onspot;
  }
  
  // Find threshold for "screamers"
  int HOTLIMIT;
  GetHotLimit(pmthitcount, NPMTS, HOTLIMIT);
  printf("Found screamer threshold: %d\n", HOTLIMIT);
  
  // ********************************************************************
  // Increase counter for PMTs with more than 5000 nhit for this run
  // ********************************************************************
  for(int id=0; id<NPMTS; id++) {
    if (cvrg[id]<0) cvrg[id]=0;                     // valid PMT ID, reset counter
    if (pmtlightcone[id] != 1) continue;            // not in direct light cone
    if (pmthitcount[id] == 0) continue;             // off PMT
    if (pmthitcount[id] > HOTLIMIT) continue;       // hot PMT
    if (pmthitcount[id] >= count/100.) cvrg[id]++;  // good coverage for this fibre
  }
 
}

void GetProjection(TGraph2D *coverage, TGraph2D* offpmts, int* cvrg, const int NPMTS, const RAT::DU::PMTInfo& pmtinfo) {
  
  // Use icosahedral projection function from HelperFunc.C
  using namespace func;
  TVector2 icospos;
  TVector3 pmtpos;
  printf("Generating icosahedral projection... ");
  int goodpmts=0, badpmts=0;
  double offx[NPMTS], offy[NPMTS];
  for(int id=0; id<NPMTS; id++) {
    pmtpos = pmtinfo.GetPosition(id);
    if (pmtpos.Mag()==0) continue;              // not a valid PMT
    if (pmtinfo.GetType(id) != 1) continue;     // not a normal PMT
    int face;                                   // side of PSUP icosahedron
    icospos = func::IcosProject(pmtpos,face);   // PMT position on flatmap
    if(cvrg[id] <= 0) {                         // inactive PMT
      offpmts->SetPoint(badpmts,icospos.X(),icospos.Y(),0.001); // non-zero Z
      offx[badpmts] = icospos.X();
      offy[badpmts] = icospos.Y();
      //printf("Bad PMT #%d at ( %.3f | %.3f )\t",id,icospos.X(),icospos.Y());
      //printf(" position ( %.3f | %.3f | %.3f )\n",id,pmtpos.X(),pmtpos.Y(),pmtpos.Z());
      badpmts++;
    } else {
      if(cvrg[id] > 10) {                       // unrecognised hot PMT
        printf("*** WARNING *** PMT #%d has coverage %d - setting to 10.\n",id,cvrg[id]);
        cvrg[id]=10;
      }
      coverage->SetPoint(goodpmts,icospos.X(),icospos.Y(),cvrg[id]);
      goodpmts++;
    }
  }
  //offpmts = new TGraph(badpmts,offx,offy);
  printf("done.\n");
  
}
