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
#include "../HelperFunc.C"
#include "../Xianguo.C"

// Global constants
const int RUN_CLUSTER = 1;  // whether running on cluster (0=local)
const int VERBOSE = 1;      // verbosity flag
const int IS_MC = 0;        // Monte-Carlo flag 

// Initialise functions
TH1D* get_histo(string, int, bool, bool);

// Main program
using namespace std;
int main(int argc, char** argv) {
  
  // Test flag (0 = process all runs, else run number)
  const int TEST = (RUN_CLUSTER) ? 0 : 102315;// 101410;

  // Initialise input file
  string input = "ANGULAR_FITRESULTS.txt";
  //string input = "../pca_runs/TELLIE_PCA.txt";
  ifstream in(input.c_str());
  if (!in) { cerr<<"Failed to open "<<input<<endl; exit(1); }
  string line;
  for (int hdr=0; hdr<2; hdr++) {
    getline(in,line);      // header
  }

  // Initialise graph
  TGraph* graph[95] = {NULL};
  //TH1D* temp = NULL;
  //TH1D* hist[95] = {NULL};
  TF1* func[95] = {NULL}; 
 
  // Loop over all fibres in list
  string fibre;
  int run, chisq, ndof;
  float a, sa, b, sb;
  //int node, channel, run, ipw, photons, pin, rms;
  //float nhit;
  int nfiles=0;
  while (true) {
    in >> fibre >> run >> a >> sa >> b >> sb >> chisq >> ndof;
    //in >> node >> fibre >> channel >> run >> ipw >> photons >> pin >> rms >> nhit;
    if (!in.good()) break;
    if (TEST && TEST!=run) continue; // only want specified run
    graph[nfiles] = new TGraph();
    graph[nfiles]->SetTitle(fibre.c_str());
    graph[nfiles]->SetPoint(0,a,b);
    //graph->SetPointError(nfiles,sa,sb);
    //temp = get_histo(fibre, run, IS_MC, TEST);
    //if (!temp) continue;
    //hist[nfiles] = temp;
    func[nfiles] = new TF1(fibre.c_str(), "[0] - [1] + [1]/cos(x/180.*pi)", 0, 24);
    func[nfiles]->SetParameters(a,b);
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
  gStyle->SetTitleOffset(1.3,"xy");  
 
  TCanvas *c0 = new TCanvas("","",1200,900);
  c0->Divide(2,2);
  TLatex *tname = new TLatex();
  tname->SetTextSize(0.03);

  // Fit lines (full range)
  c0->cd(1)->SetGrid();
  c0->cd(1)->DrawFrame(0,0,24,40,"TELLIE angular systematic fits;Angle of PMT w.r.t. fitted fibre direction [deg];Mean hit time offset [ns]");
  for (int fib=0; fib<95; fib++) {
    //if (!hist[fib]) continue;
    //hist[fib]->SetLineWidth(1);
    //hist[fib]->SetLineColor(100-fib);
    //hist[fib]->Draw("hist same L");
    if (!func[fib]) continue;
    func[fib]->SetLineWidth(1);
    func[fib]->SetLineColor(100-fib);
    func[fib]->Draw("L same");
   
    // Highlight unusual fibres
    float yval = func[fib]->GetParameter(0);
    if(yval > 10) {
      tname->DrawLatex(10,yval-1.5,func[fib]->GetName());
    }
  }
  // Box for text
  TBox *tbox = new TBox(0.6,33.5,14.8,39);
  tbox->SetLineColor(2);
  tbox->SetFillColor(kYellow-9);
  tbox->Draw("L same");
  // Fit function as text
  TLatex *text = new TLatex(1,38,"Fit function: y = a #plus b#left(#frac{1}{cos(x)} #minus 1#right)");
  text->SetTextAlign(13);
  text->SetTextFont(62);
  text->SetTextSize(0.04);
  text->Draw("same");

  // Fit lines (zoomed in)
  c0->cd(2)->SetGrid();
  c0->cd(2)->DrawFrame(0,0,16,5,"TELLIE angular systematic fits (zoom);Angle of PMT w.r.t. fitted fibre direction [deg];Mean hit time offset [ns]");
  for (int fib=0; fib<95; fib++) {
    //if (!hist[fib]) continue;
    //hist[fib]->SetLineWidth(1);
    //hist[fib]->SetLineColor(100-fib);
    //hist[fib]->Draw("hist same L");
    if (!func[fib]) continue;
    func[fib]->SetLineWidth(1);
    func[fib]->SetLineColor(100-fib);
    func[fib]->Draw("L same");
    
    // Highlight unusual fibres
    float yval = func[fib]->GetParameter(0);
    if(yval > 2.9) {
      tname->DrawLatex(2,yval+0.1,func[fib]->GetName());
    }
  }
  
  // Fit parameters (full range)
  c0->cd(3)->SetGrid();
  c0->cd(3)->DrawFrame(0,-10,30,70,"Fit parameters;Fit parameter a [ns];Fit parameter b [ns]");
  for (int fib=0; fib<95; fib++) {
    graph[fib]->SetMarkerStyle(7);
    graph[fib]->SetMarkerColor(100-fib);
    graph[fib]->Draw("P same");
    
    // Highlight unusual fibres
    double xval,yval;
    graph[fib]->GetPoint(0,xval,yval);
    if(xval > 10 || yval < 0) {
      tname->SetTextAlign(23);
      tname->DrawLatex(xval-0.5,yval-2,graph[fib]->GetTitle());
    }
  }
  
  // Fit parameters (zoomed in)
  c0->cd(4)->SetGrid();
  c0->cd(4)->DrawFrame(1,10,3,70,"Fit parameters (zoom);Fit parameter a [ns];Fit parameter b [ns]");
  for (int fib=0; fib<95; fib++) {
    graph[fib]->SetMarkerStyle(7);
    graph[fib]->SetMarkerColor(100-fib);
    graph[fib]->Draw("P same");
    
    // Highlight unusual fibres
    double xval,yval;
    graph[fib]->GetPoint(0,xval,yval);
    if(xval > 2.9 || yval < 15 || (xval>1.6 && xval<2.1)) {
      tname->SetTextAlign(32);
      tname->DrawLatex(xval-0.02,yval,graph[fib]->GetTitle());
    }
  }
  
  // Save canvas and close
  string imgfile = "summary";
  c0->Print(Form("%s.png",imgfile.c_str()));
  c0->Print(Form("%s.pdf",imgfile.c_str()));
  c0->Close();
  
  // Free memory
  if (c0) delete c0;
  if (func) for (int fib=0; fib<95; fib++) delete func[fib];
  if (graph) for (int fib=0; fib<95; fib++) delete graph[fib];
  //if (hist) for (int fib=0; fib<95; fib++) delete hist[fib];
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
  TH1D* hist = (TH1D*)infile->Get("hmean");  // Mean time
  return hist;
}
