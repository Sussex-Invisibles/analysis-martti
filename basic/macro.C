// ---------------------------------------------------------
// Goal:            Basic macro to get event information
// Author:          Martti Nirkko, 21/09/2017
// Compile & run:   clear && g++ -o macro.exe macro.C `root-config --cflags --libs` -I$RATROOT/include/libpq -I$RATROOT/include -L$RATROOT/lib -lRATEvent_Linux && ./macro.exe
// ---------------------------------------------------------

// C++ stuff
#include <fstream>

// ROOT stuff
#include "TH2.h"
#include "TCanvas.h"
#include "TMath.h"
#include "TStyle.h"

// RAT stuff
#include "RAT/DU/DSReader.hh"
#include "RAT/DU/PMTInfo.hh"
#include "RAT/DU/Utility.hh"
#include "RAT/DS/Entry.hh"

// Physical constants
const double pi = TMath::Pi();
const double c_vacuum = 299.792458;        // mm/ns (detector units)
const double n_water = 1.33772;            // at 500 nm -> http://www.philiplaven.com/p20.html
const double c_water = c_vacuum/n_water;   // mm/ns

// Global parameters
const int VERBOSE = 1;      // verbosity flag
const int IS_MC = 0;        // Monte-Carlo flag 

// Initialise functions
string printVector(const TVector3&);
int analyseData(string, string, bool);

// *****************************************************************************
// Main program
using namespace std;
int main(int argc, char** argv) {
  
  // Choose run number, and corresponding TELLIE fibre
  int run = 102315;
  string fibre = "FT019A";
  
  // Set file location
  string path = Form("%s/Software/SNOP/work/data", getenv("HOME"));
  string file = Form("%s/Analysis_r0000%d_s000_p000.root", path.c_str(), run);
  
  // Run macro
  int error = analyseData(file, fibre, IS_MC);
  
  // End of main program
  if (error) { 
    cerr << "*** ERROR *** Macro failed!" << endl;
    return 1;
  }
  return 0;
}

// *****************************************************************************
// Displays a vector as string "( x | y | z )"
string printVector(const TVector3& v) {
  string out = Form("( %.3f | %.3f | %.3f )", v.X(),  v.Y(), v.Z());
  return out.c_str();
}

// *****************************************************************************
// Macro which does all the cool stuff
int analyseData(string file, string fibre, bool isMC=false) {

  // Check files for given run
  ifstream f(file.c_str());
  cout << "Opening file: " << file << endl;
  if (!f.good()) { 
    cerr << "Error - cannot open file!" << endl;
    return 1;
  }
  
  // Initialise RAT
  RAT::DU::DSReader dsreader(file);
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
  if (VERBOSE) cout << "RATDB info: fibre " << fibre << ", position " << printVector(fibrepos) << ", direction " << printVector(fibredir) << endl;
  
  // Initialise histograms
  const int NBINS = 30;
  const int MAXANG = 30;    // deg
  const int MINTIME = 180;  // ns
  const int MAXTIME = 210;  // ns
  TH2D *hist = new TH2D("hist",fibre.c_str(),NBINS,0,MAXANG,NBINS,MINTIME,MAXTIME);
  
  // Initialise variables
  int totalnhit = 0;
  int nevents = 0;
  
  // Loop over all entries in file
  for(int iEntry=0; iEntry<dsreader.GetEntryCount(); iEntry++) {
    const RAT::DS::Entry& ds = dsreader.GetEntry(iEntry);
    RAT::DS::MC mc;
    if (isMC) mc = ds.GetMC();  // don't initialise this for real data (crashes)
    
    // Loop over triggered events in each entry
    for(int iEv=0; iEv<ds.GetEVCount(); iEv++) {            // mostly 1 event per entry
      const RAT::DS::EV& ev = ds.GetEV(iEv);
      
      // Trigger type
      int trig = ev.GetTrigType();
      if (!(trig & 0x8000)) continue;                       // EXT trigger only
      
      // Event observables
      int nhitscleaned = ev.GetNhitsCleaned();              // calibrated PMT hits with removed crosstalk
      totalnhit += nhitscleaned;
      nevents++;
      
      // PMT information
      const RAT::DS::CalPMTs& pmts = ev.GetCalPMTs();
      for(int iPMT=0; iPMT<pmts.GetNormalCount(); iPMT++) {
        RAT::DS::PMTCal pmt = pmts.GetNormalPMT(iPMT);
        int pmtID = pmt.GetID();
        
        // Cut out undesirable PMT conditions
        if (!chs.IsTubeOnline(pmtID)) continue;             // test CHS
        if (pmt.GetCrossTalkFlag()) continue;               // remove crosstalk
        if (pmt.GetStatus().GetULong64_t(0) != 0) continue; // test PCA
        
        // Get useful information
        TVector3 pmtpos = pmtinfo.GetPosition(pmtID);       // PMT position
        TVector3 track = pmtpos-fibrepos;                   // light direction
        double theta = track.Angle(fibredir);               // angle w.r.t. fibre [rad]
        double pmttime = pmt.GetTime();                     // hit time [ns]
        double corr = track.Mag()/c_water;                  // light travel time [ns]
        
        // Fill histograms
        hist->Fill(theta*180./pi, pmttime-corr);            // angle [deg], time [ns]
        
      } // pmt loop
    } // event loop
  } // entry loop
  
  float avg_nhit = (float)totalnhit/nevents;
  cout << "Found " << nevents << " events with an average " << avg_nhit << " nhits/event." << endl;
  
  // Plotting options
  gStyle->SetOptStat(0);    // turn off statistics box
  
  // Define canvas and pads
  TCanvas *c0 = new TCanvas("c0","Awesome plot",800,600);
  c0->SetGrid();
  c0->SetRightMargin(0.15);   // for TH2D color scale
  hist->SetTitle(Form("TELLIE fibre profile (%s);Angle [deg];Time [ns]",fibre.c_str()));
  hist->Draw("colz");

  // Save canvas and close
  string imgfile = "macro";
  c0->Print("macro.png");
  c0->Close();
  
  // Free memory created with "new"
  if (c0) delete c0;
  if (hist) delete hist;
  
  // Return value
  return 0;
}

