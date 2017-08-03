// ---------------------------------------------------------
// Goal:            Print PMT positions in 3D and 2D (flatmap)
// Author:          Martti Nirkko, 26/04/2017
// Compile & run:   clear && g++ -g -o positions.exe positions.C `root-config --cflags --libs` -I$RATROOT/include/libpq -I$RATROOT/include -L$RATROOT/lib -lRATEvent_Linux && ./positions.exe
// ---------------------------------------------------------

// C++ stuff
#include <fstream>
#include <iostream>
#include <stdio.h>
#include <stdlib.h>

// ROOT stuff
#include "TROOT.h"
#include "TFile.h"

// RAT stuff
#include "RAT/DU/DSReader.hh"
#include "RAT/DU/PMTInfo.hh"
#include "RAT/DU/Utility.hh"
#include "RAT/DS/Run.hh"
#include "RAT/DS/Entry.hh"
#include "RAT/DS/MC.hh"

// Helper functions
#include "HelperFunc.C"

// Global constants
const int RUN_CLUSTER=0;    // whether running on cluster (0=local)
const int IS_MC = 0;        // Monte-Carlo flag 

// Initialise functions
void positions(string, int, bool, bool);

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
    positions(fibre, run, IS_MC, TEST);
  }

  return 0;
}

// Returns the fitted light position for a given fibre/run
void positions(string fibre, int run, bool isMC=false, bool TEST=false) {

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
  
  int type;
  TVector3 pmtpos;
  TVector2 icospos;
  for(int id=0; id<NPMTS; id++) {
    pmtpos = pmtinfo.GetPosition(id);
    type = pmtinfo.GetType(id);
    //if (pmtpos.Mag()==0) continue; // not a valid PMT
    //if (type!=1) continue; // not a normal PMT
    int face; // side of icosahedron containing PMT
    icospos = func::IcosProject(pmtpos,face);
    if (pmtpos.X()==-99999.) icospos.Set(-99999.,-99999.);
    if (type!=1) printf("%4d %d %8lf %8lf %8lf %8lf %8lf\n",id,type,pmtpos.X(),pmtpos.Y(),pmtpos.Z(),icospos.X(),icospos.Y());
  }
  
}
