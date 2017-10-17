// ---------------------------------------------------------
// Goal:            Figure out which TELLIE fibres hit a given crate
// Author:          Martti Nirkko, 26/04/2017
// Compile & run:   clear && g++ -g -o focal_point.exe focal_point.C `root-config --cflags --libs` -I$RATROOT/include/libpq -I$RATROOT/include -L$RATROOT/lib -lRATEvent_Linux && ./focal_point.exe
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
#include "../HelperFunc.C"
#include "../Xianguo.C"

// Global constants
const int RUN_CLUSTER = 1;  // whether running on cluster (0=local)
const int USE_RATDB = 1;    // whether to use RATDB to get fibre positions (1=yes)
const int VERBOSE = 0;      // verbosity flag
const int IS_MC = 0;        // Monte-Carlo flag 
const int NCOL = 20;        // number of colors (max. 50)
const int NDOTS = 360;      // number of points in circle
const double DIR_CONE = 48; // opening angle to search for direct light (using aperture: 24 deg)
const double REF_CONE = 20; // opening angle to search for reflected light (using aperture: 9.874 deg)

// Initialise functions
void focal_point(string, int, int, int, int, float, TVector3*, TVector3*, bool, bool);

// Main program
using namespace std;
int main(int argc, char** argv) {
  
  // Test flag (0 = process all runs, else run number)
  const int TEST = (RUN_CLUSTER) ? 0 : 102157;
  
  // Run over Monte Carlo file
  if (IS_MC) {
    string fibre = "FT040A";	// FT040A, FT079A
    TVector3 dirfit, reffit;
    focal_point(fibre, -999, -999, -999, -999, -999., &dirfit, &reffit, (bool)IS_MC, false);
    cout << "Results for fibre " << fibre << ":" << endl;
    cout << "- Direct light fit (x,y,z) [mm] = " << printVector(dirfit) << endl;
    cout << "- Reflected light fit (x,y,z) [mm] = " << printVector(reffit) << endl;
    return 0;
  }

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
  TVector3 dirfit, reffit;
  FILE *fitresult = fopen("TELLIE_FITRESULTS.txt","w");
  fprintf(fitresult,"#Run Fibre Direct_light(xyz) Reflected_light(xyz)\n");
  fprintf(fitresult,"#------------------------------------------------\n");
  while (true) {
    in >> node >> fibre >> channel >> run >> ipw >> photons >> pin >> rms >> nhit;
    if (!in.good()) break;
    if (TEST && TEST!=run) continue; // only want specified run
    if (VERBOSE) printf("%6s %2d %6d %5d %6d %5d %5d %.2f\n", fibre.c_str(), channel, run, ipw, photons, pin, rms, nhit);
    dirfit.SetXYZ(0,0,0);
    reffit.SetXYZ(0,0,0);
    focal_point(fibre, channel, run, ipw, photons, nhit, &dirfit, &reffit, (bool)IS_MC, (bool)TEST);
    if (dirfit.Mag()==0) continue;
    fprintf(fitresult, "%6d %6s %.3f %.3f %.3f %.3f %.3f %.3f\n", run, fibre.c_str(), dirfit.X(), dirfit.Y(), dirfit.Z(), reffit.X(), reffit.Y(), reffit.Z());
    cout << "Direct light fit " << printVector(dirfit) << endl;
    cout << "Reflected light fit " << printVector(reffit) << endl;
    nfiles++;
  }
  fclose(fitresult);
  printf("Ran over %d files.\n",nfiles);
  if (nfiles==0) { 
    cerr<<"*** ERROR *** No input files found, or nothing to do!"<<endl;
    return 1;
  }
  return 0;
}

// Calculates the fitted light and fibre positions for a given fibre/run
void focal_point(string fibre, int channel, int run, int ipw, int photons, float nhit, TVector3 *dirfit, TVector3 *reffit, bool isMC=false, bool TEST=false) {

  // ********************************************************************
  // Initialisation
  // ********************************************************************
  
  // Check files for given run
  if(!TEST && VERBOSE) printf("*****\n");
  printf("Checking files for run %d... ", run);
  //string fpath = (RUN_CLUSTER) ? "/lustre/scratch/epp/neutrino/snoplus/TELLIE_PCA_RUNS_PROCESSED" : "/home/nirkko/Desktop/fibre_validation";
  string fpath = Form("%s/Software/SNOP/work/data",getenv("HOME"));
  string fname, out, img;
  string mcopt = "";
  ifstream f;
  if (isMC) {
    mcopt = "centre"; // centre, Z+10cm (shifted AV position)
    fpath = Form("%s/Software/SNOP/work/analysis/tellie_jobs/markstringer",getenv("HOME"));
    fname = Form("%s/MC_%s_AV_%s_processed.root",fpath.c_str(),fibre.c_str(),mcopt.c_str());
    f.open(fname.c_str());
    run = 0; // TODO: figure out how to get this from file
    out = Form("./mc/MC_%s_%s.out",fibre.c_str(),mcopt.c_str());
    img = Form("./mc/MC_%s.pdf",fibre.c_str());
  } else {
    for (int pass=3;pass>=0;pass--) {
      fname = Form("%s/Analysis_r0000%d_s000_p00%d.root",fpath.c_str(),run,pass);
      f.open(fname.c_str());
      if (f.good()) break;
    }
    out = Form("./output/PCA_%s.out",fibre.c_str());
    img = Form("./images/PCA_%s.pdf",fibre.c_str());
  }
  ifstream g(out.c_str());
  ifstream h(img.c_str());
  int scanned_file = 0;
  if (!TEST && h.good()) {  // file downloaded and processed
    printf("already processed! Skipping fibre %s.\n",fibre.c_str());
    return;
  } else if (g.good()) {     // file extracted, but not processed
    printf("not processed! Generating plots for fibre %s.\n",fibre.c_str());
    scanned_file = 1;
  } else if (!f.good()) {    // file not downloaded
    printf("not downloaded! Skipping fibre %s.\n",fibre.c_str());
    return;
  } else {                   // file downloaded, but not processed
    printf("OK. Extracting data for fibre %s.\n",fibre.c_str());
  }

  // Initialise RAT
  RAT::DU::DSReader dsreader(fname);
  const RAT::DU::ChanHWStatus& chs = RAT::DU::Utility::Get()->GetChanHWStatus();
  const RAT::DU::PMTInfo& pmtinfo = RAT::DU::Utility::Get()->GetPMTInfo();
  const int NPMTS = pmtinfo.GetCount();

  if (USE_RATDB) {
    // Get fibre info (from RATDB) 
    RAT::DB *db = RAT::DB::Get();
    //db->LoadDefaults();	  // Already done when calling DU::Utility::Get()
    RAT::DBLinkPtr entry = db->GetLink("FIBRE",fibre);
    TVector3 fibrepos(entry->GetD("x"), entry->GetD("y"), entry->GetD("z")); // position
    TVector3 fibredir(entry->GetD("u"), entry->GetD("v"), entry->GetD("w")); // direction
    TVector3 lightpos = fibrepos + 2*fibrepos.Mag()*fibredir; // projected light spot centre
    cout << "RATDB: fibre " << fibre << ", pos " << printVector(fibrepos) << ", dir " << printVector(fibredir) << endl;
  } else {
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
    cout << "DocDB: fibre " << fibre << ", pos " << printVector(fibrepos) << ", dir " << printVector(fibredir) << endl;
  }

  // Initialise histograms
  TH2D *hicos = new TH2D("hicos","PMT positions",200,0,1,200,0,1); // icosahedral
  TH2D *hcoarse = new TH2D("hcoarse","PMT positions",20,-1,1.001,10,-1.001,0); // units of pi
  TH2D *hfineD = new TH2D("hfineD","PMT positions",1000,-10,10,1000,-10,10); // fine grained
  TH2D *hfineR = new TH2D("hfineR","PMT positions",1000,-10,10,1000,-10,10); // fine grained
  TH1I *hhits = new TH1I("hhits","PMT hit count",50,0,3e5);
  BinLog(hhits->GetXaxis(),0.3);
  
  // ********************************************************************
  // Sum PMT hit counts for entire run
  // ********************************************************************
  int pmthitcount[NPMTS], pmtlightcone[NPMTS];
  memset( pmthitcount, 0, NPMTS*sizeof(int) ); // NPMTS only known at runtime
  memset( pmtlightcone, 0, NPMTS*sizeof(int) ); // NPMTS only known at runtime
  int totalnhit=0, directnhit=0, reflectednhit=0, count=0;
  if (scanned_file) {
    int checkrun;
    string dummy;
    g >> dummy >> checkrun;
    if(!g.good()) printf("*** ERROR *** Bad input - str=%s int=%d\n",dummy.c_str(),checkrun);
    if(checkrun != run && !isMC) { printf("*** ERROR *** Bad run number %d\n",checkrun); return; }
    g >> dummy >> count;
    if(!g.good()) { printf("*** ERROR *** Bad input - str=%s int=%d\n",dummy.c_str(),count); return; }
    g >> dummy >> totalnhit;
    if(!g.good()) { printf("*** ERROR *** Bad input - str=%s int=%d\n",dummy.c_str(),totalnhit); return; }
    // Print run info here
    printf("*** INFO *** Run %d has %d EXTA events with %d total NHits.\n",checkrun,count,totalnhit);
    int pmtid, pmthits, onspot;
    while (g.good()) {
      g >> pmtid >> pmthits >> onspot;
      pmthitcount[pmtid]=pmthits;
      pmtlightcone[pmtid]=onspot;
      if(onspot==1) directnhit += pmthits;
      else if(onspot==2) reflectednhit += pmthits;
    }
  } else {
    for(int iEntry=0; iEntry<dsreader.GetEntryCount(); iEntry++) {
      const RAT::DS::Entry& ds = dsreader.GetEntry(iEntry);
      
      // Loop over triggered events in each entry
      for(int iEv=0; iEv<ds.GetEVCount(); iEv++) { // mostly 1 event per entry
        const RAT::DS::EV& ev = ds.GetEV(iEv);
        
        // Trigger type
        int trig = ev.GetTrigType();
        if (!(trig & 0x8000)) continue;            // EXT trigger only
        
        // Event observables
        int nhitscleaned = ev.GetNhitsCleaned();   // calibrated PMT hits with removed crosstalk
        totalnhit += nhitscleaned;
        count++;
        
        // PMT information
        const RAT::DS::CalPMTs& pmts = ev.GetCalPMTs();
        //int eventnhit=0;
        for(int iPMT=0; iPMT<pmts.GetNormalCount(); iPMT++) {
          RAT::DS::PMTCal pmt = pmts.GetNormalPMT(iPMT);
          int pmtID = pmt.GetID();
          if (!chs.IsTubeOnline(pmtID)) continue;                   // test CHS
          if (pmt.GetCrossTalkFlag()) continue;                     // remove crosstalk
          if (pmt.GetStatus().GetULong64_t(0) != 0) continue;       // test PCA
          pmthitcount[pmtID]++;
          //eventnhit++;
        } // pmt loop
        //if(nhitscleaned!=eventnhit) printf("*** WARNING *** - cleaned nhits %d != %d CalPMT nhits!\n",nhitscleaned,eventnhit);
      } // event loop
    } // entry loop

    // Write PMT hit count to output file (saves time when rerunning)
    FILE *outFile = fopen(out.c_str(),"w");
    TVector3 pmtpos;
    fprintf(outFile, "Run: %d\n", run);             // run number
    fprintf(outFile, "Events: %d\n", count);        // events with EXTA trigger
    fprintf(outFile, "TotalNHit: %d\n", totalnhit); // calibrated PMT hits with removed crosstalk
    for(int id=0; id<NPMTS; id++) {
      pmtpos = pmtinfo.GetPosition(id);
      int in_spot = 0;
      if (pmtpos.Angle(lightpos) < DIR_CONE/180.*pi) in_spot = 1;
      else if (pmtpos.Angle(fibrepos) < REF_CONE/180.*pi) in_spot = 2;
      fprintf(outFile, "%d %d %d\n", id, pmthitcount[id], in_spot);
    }
    fclose(outFile);
    scanned_file = 1;
  
  } // else condition
  
  // Find threshold for "screamers"
  int HOTLIMIT;
  GetHotLimit(pmthitcount, NPMTS, HOTLIMIT);
  for(int id=0; id<NPMTS; id++) {
    if(pmthitcount[id]==0) hhits->Fill(0.5); // visible on log scale
    else hhits->Fill(pmthitcount[id]);
  }
 
  // ********************************************************************
  // Get icosahedral projection of detector
  // ********************************************************************
  // Use icosahedral projection functions from HelperFunc.C
  using namespace func;
  int icosN[NCOL+2]={0};
  int pmtface[NPMTS];
  double icosX[NCOL+2][NPMTS], icosY[NCOL+2][NPMTS];
  memset( icosX, 0, (NCOL+2)*NPMTS*sizeof(double) ); // NPMTS only known at runtime
  memset( icosY, 0, (NCOL+2)*NPMTS*sizeof(double) ); // NPMTS only known at runtime
  memset( pmtface, 0, NPMTS*sizeof(int) ); // NPMTS only known at runtime
  
  TVector2 icospos;
  TVector3 pmtpos;
  int nfacepmts[20] = {0}, nfacegoodpmts[20] = {0};
  TVector3 facecenter[20], faceweight[20];
  for(int id=0; id<NPMTS; id++) {
    pmtpos = pmtinfo.GetPosition(id);
    if (pmtpos.Mag()==0) continue; // not a valid PMT
    if (pmtinfo.GetType(id) != 1) continue; // not a normal PMT
    int face; // side of icosahedron containing PMT
    icospos = func::IcosProject(pmtpos,face);
    pmtface[id] = face; // face of PMT (1-20)
    facecenter[face-1] += pmtpos;
    nfacepmts[face-1]++;
    //for(int j=0; j<pmthitcount[id]; j++) hicos->Fill(icospos.X(),icospos.Y());
    int step = (int)TMath::Ceil(pmthitcount[id]/(1.*HOTLIMIT/NCOL))+1;
    if (pmthitcount[id] > HOTLIMIT) step=0; // hot PMT 
    icosX[step][icosN[step]] = icospos.X();
    icosY[step][icosN[step]] = icospos.Y();
    icosN[step]++;
    if (step<2) continue; // hot or off PMT
    faceweight[face-1] += pmtpos*pmthitcount[id];
    nfacegoodpmts[face-1]++;
  }
  int maxface=-1;
  double faceheat[20]={0}, maxfaceheat=-1;
  for(int fc=0; fc<20; fc++) {
    if (nfacepmts[fc]==0 || nfacegoodpmts[fc]==0) continue;
    facecenter[fc] *= 1./nfacepmts[fc];     // all PMTs weighted equally
    faceweight[fc] *= 1./nfacegoodpmts[fc]; // good PMTS weighted by intensity
    faceheat[fc] = faceweight[fc].Mag();
    if (VERBOSE) printf("Face #%2d has intensity %6.2lfM\n",fc+1,faceheat[fc]/1e6);
    if (faceheat[fc]>maxfaceheat) { maxfaceheat=faceheat[fc]; maxface=fc+1; }
  }
  if (VERBOSE) printf("Hot faces: ");
  TVector3 bestguess(0,0,0);
  int hotfaces=0;
  for(int fc=0; fc<20; fc++) {
    if(faceheat[fc]<maxfaceheat/5.) continue; // face is "hot" if intensity is >20% of hottest face
    if (VERBOSE) printf("%d ",fc+1);
    bestguess += faceweight[fc];
    hotfaces++;
  }
  if (hotfaces>0) { bestguess *= 1./hotfaces; }
  else { cerr<<"Something went wrong! No hot faces found."<<endl; return; }
  if (VERBOSE) printf(" -> best guess direction: %s\n",printVector(bestguess.Unit()).c_str());
  
  // Make graphs for icosahedral projection
  TGraph *icos[NCOL+2];
  for (int s=0; s<NCOL+2; s++) {
    int col  = (int)(50.+s*(50./NCOL));
    if (icosN[s]==0) {icos[s]=NULL; continue;}
    icos[s] = new TGraph(icosN[s],icosX[s],icosY[s]);
    icos[s]->SetMarkerStyle(7);
    icos[s]->SetMarkerColor(col);
  }
  if (icos[0]) icos[0]->SetMarkerColor(1);  // hot PMTs
  if (icos[1]) icos[1]->SetMarkerColor(16); // off PMTs
  
  
  // ********************************************************************
  // First iteration: Fill coarse histogram with PMT hit positions
  // ********************************************************************
  int hitpmts=0, init=0;
  double pmtrad = 0;
  TVector3 pmtsum(0,0,0);
  for(int id=0; id<NPMTS; id++) {
    pmtpos = pmtinfo.GetPosition(id);
    if (pmtpos.Mag()==0) continue;            // not a valid PMT
    if (pmtinfo.GetType(id) != 1) continue;   // not a normal PMT
    if (pmthitcount[id] == 0) continue;       // off PMTs
    if (pmthitcount[id] > HOTLIMIT) continue; // hot PMTs
    for (int j=0; j<pmthitcount[id]; j++) {
      hcoarse->Fill(pmtpos.Phi()/pi, -pmtpos.Theta()/pi);  // negative theta (neck on top)
    }
    pmtrad += pmtpos.Mag()*pmthitcount[id];
    pmtsum += pmtpos*pmthitcount[id];
    hitpmts += pmthitcount[id];
  }
  pmtrad /= hitpmts;
  if(fabs(lightpos.Mag()/pmtrad-1.) > 0.01) {  // tolerate 1%
    cout << "*** WARNING *** Projected light position deviates from PSUP sphere!" << endl;
    return;
  }
  if(VERBOSE) printf("PMT radius = %.3lf mm\n",pmtrad);

  // Get the maximum bin as guesstimate for light spot
  bestguess.SetMag(pmtrad);
  TVector3 guess_dir =  bestguess;
  TVector3 guess_ref = -bestguess;  // TODO - improve guesstimate for reflected light?
  
  
  // ********************************************************************
  // Second iteration: take weighted average around estimated light spots
  // ********************************************************************
  TVector3 direct(0.,0.,0.);
  TVector3 reflected(0.,0.,0.);
  int weight, nearmax, farmax, empty1, empty2;
  GetMaxColVal(guess_dir, pmthitcount, NPMTS, nearmax, empty1, pmtinfo);
  GetMaxColVal(guess_ref, pmthitcount, NPMTS, farmax, empty2, pmtinfo);
  
  if(nearmax==0) cout << "*** WARNING *** No good PMTs in direct light cone!" << endl;
  if(farmax==0) cout << "*** WARNING *** No good PMTs in reflected light cone!" << endl;
  
  for(int id=0; id<NPMTS; id++) {
    pmtpos = pmtinfo.GetPosition(id);
    if (pmtpos.Mag()==0) continue;                  // not a valid PMT position
    if (pmtinfo.GetType(id) != 1) continue;         // not a normal PMT
    weight = pmthitcount[id];
    if (weight == 0) continue;                      // off PMTs
    if (weight > HOTLIMIT) continue;                // hot PMTs
    if (pmtpos.Angle(guess_dir) < pi*DIR_CONE/180.) // within cone of direct light spot
      direct += pmtpos*(1.*weight/hitpmts);
    if (pmtpos.Angle(guess_ref) < pi*REF_CONE/180.) // within cone of reflected light spot
      reflected += pmtpos*(1.*weight/hitpmts);
  }
  if (direct.Mag() == 0) {
    cout << "*** WARNING *** No good PMTs in direct light cone! Setting to z-axis." << endl;
    direct.SetXYZ(0,0,1);
  }
  if (reflected.Mag() == 0) { 
    cout << "*** WARNING *** No good PMTs in reflected light cone! Setting to z-axis." << endl;
    reflected.SetXYZ(0,0,1);
  }
  direct.SetMag(pmtrad);
  reflected.SetMag(pmtrad);
  
  // Fill graphs with new view, centered on weighted light spot
  TGraph *gDir[NCOL+2], *gRef[NCOL+2];  // array of graphs (for 2D view)
  TGraph2D *gDir2D = new TGraph2D();    // 2D graph (for 3D view and fit)
  TGraph2D *gRef2D = new TGraph2D();
  FillHemisphere(direct, pmthitcount, NPMTS, gDir, gDir2D, NCOL, nearmax, pmtinfo);
  FillHemisphere(reflected, pmthitcount, NPMTS, gRef, gRef2D, NCOL, farmax, pmtinfo);
  
  // Get rotation angles for rotating view over weighted light spots
  double rot_Z1, rot_X1, rot_Z2, rot_X2;
  GetRotationAngles(direct, rot_Z1, rot_X1);
  GetRotationAngles(reflected, rot_Z2, rot_X2);
  
  
  // ********************************************************************
  // Third iteration: perform 2D Gaussian fit around weighted light spots
  // ********************************************************************
  // Fit 2D graphs with Gaussian surface
  const int NPARS = 4;
  double parDir[NPARS], parRef[NPARS];
  FitLightSpot(gDir2D,pmtrad/1e3,DIR_CONE,parDir); // units: dist [m], ang [deg]
  FitLightSpot(gRef2D,pmtrad/1e3,REF_CONE,parRef);
  cout << "FIT RESULTS:" << endl;
  cout << "- Direct light";
  for (int p=0; p<NPARS; p++) printf(" %.3lf", parDir[p]);
  cout << endl;
  cout << "- Reflected light";
  for (int p=0; p<NPARS; p++) printf(" %.3lf", parRef[p]);
  cout << endl;
  // Translate to point on PSUP sphere
  TVector3 fitDir(1e3*parDir[1],1e3*parDir[2],sqrt(pmtrad*pmtrad-1e6*parDir[1]*parDir[1]-1e6*parDir[2]*parDir[2]));
  TVector3 fitRef(1e3*parRef[1],1e3*parRef[2],sqrt(pmtrad*pmtrad-1e6*parRef[1]*parRef[1]-1e6*parRef[2]*parRef[2]));
  // Standard deviation (Gaussian sigma) translated to angle [deg]
  double sigmaAngDir = asin(1e3*parDir[3]/pmtrad)*180./pi;
  double sigmaAngRef = asin(1e3*parRef[3]/pmtrad)*180./pi;
  
  // Output values (new)
  dirfit->SetXYZ(fitDir.X(), fitDir.Y(), fitDir.Z());
  dirfit->RotateY(rot_Z1);
  dirfit->RotateZ(rot_X1);
  reffit->SetXYZ(fitRef.X(), fitRef.Y(), fitRef.Z());
  reffit->RotateY(rot_Z2);
  reffit->RotateZ(rot_X2);
  double DIRANG = lightpos.Angle(*dirfit);   // angle between expected and fitted light spot
  double REFANG = fibrepos.Angle(*reffit);   // angle between fibre position and fitted light spot
  /*
  // Output values (old)
  double DIRANG = lightpos.Angle(direct);    // angle between expected and weighted light spot
  double REFANG = fibrepos.Angle(reflected); // angle between fibre position and weighted light spot
  dirfit->SetXYZ(direct.X(), direct.Y(), direct.Z());
  reffit->SetXYZ(reflected.X(), reflected.Y(), reflected.Z());
  */
  
  
  // ********************************************************************
  // Create histograms and graphs for output
  // ********************************************************************
  // Get contours around fitted light spots in (phi,theta)
  TVector3 *dotsD[NDOTS], *dotsR[NDOTS];
  DrawCircle(*dirfit, sigmaAngDir, dotsD, NDOTS);
  DrawCircle(*reffit, sigmaAngRef, dotsR, NDOTS);
  double pcontDx[NDOTS], pcontDy[NDOTS];
  double pcontRx[NDOTS], pcontRy[NDOTS];
  for (int d=0; d<NDOTS; d++) {
    int dummy;
    icospos = func::IcosProject(*dotsD[d],dummy);
    pcontDx[d]=icospos.X();
    pcontDy[d]=icospos.Y();
    icospos = func::IcosProject(*dotsR[d],dummy);
    pcontRx[d]=icospos.X();
    pcontRy[d]=icospos.Y();
  }
  TGraph *pcontD = new TGraph(NDOTS,pcontDx,pcontDy);
  TGraph *pcontR = new TGraph(NDOTS,pcontRx,pcontRy);
  pcontD->SetMarkerStyle(6);
  pcontR->SetMarkerStyle(6);

  // Marker styles in ROOT 5.34:
  /* enum EMarkerStyle {kDot=1, kPlus, kStar, kCircle=4, kMultiply=5,
                        kFullDotSmall=6, kFullDotMedium=7, kFullDotLarge=8,
                        kFullCircle=20, kFullSquare=21, kFullTriangleUp=22,
                        kFullTriangleDown=23, kOpenCircle=24, kOpenSquare=25,
                        kOpenTriangleUp=26, kOpenDiamond=27, kOpenCross=28,
                        kFullStar=29, kOpenStar=30, kOpenTriangleDown=32,
                        kFullDiamond=33, kFullCross=34}; */
  int guessMarker = 1;
  int weightMarker = 2;
  int fitMarker = 34;
  int trueMarker = 24;
  
  // Create markers for relevant points (icosahedral view)
  int nope;
  double ptX, ptY;
  icospos = func::IcosProject(guess_dir,nope);
  ptX = icospos.X();
  ptY = icospos.Y();
  TGraph *pGssDir = new TGraph(1,&ptX,&ptY);        // direct light (guess)
  pGssDir->SetMarkerStyle(guessMarker);
  pGssDir->SetMarkerColor(1);
  pGssDir->SetMarkerSize(2);
  
  icospos = func::IcosProject(guess_ref,nope);
  ptX = icospos.X();
  ptY = icospos.Y();
  TGraph *pGssRef = new TGraph(1,&ptX,&ptY);        // reflected light (guess)
  pGssRef->SetMarkerStyle(guessMarker);
  pGssRef->SetMarkerColor(1);
  pGssRef->SetMarkerSize(1.5);
  
  icospos = func::IcosProject(*dirfit,nope);
  ptX = icospos.X();
  ptY = icospos.Y();
  TGraph *pFitDir = new TGraph(1,&ptX,&ptY);          // direct light (fitted)
  pFitDir->SetMarkerStyle(fitMarker);
  pFitDir->SetMarkerColor(1);
  pFitDir->SetMarkerSize(2);
  
  icospos = func::IcosProject(*reffit,nope);
  ptX = icospos.X();
  ptY = icospos.Y();
  TGraph *pFitRef = new TGraph(1,&ptX,&ptY);          // reflected light (fitted)
  pFitRef->SetMarkerStyle(fitMarker);
  pFitRef->SetMarkerColor(1);
  pFitRef->SetMarkerSize(1.5);
  
  // Create markers for relevant points (rotated view)
  lightpos.RotateZ(-rot_X1);
  lightpos.RotateY(-rot_Z1);
  ptX = lightpos.X()/1e3;
  ptY = lightpos.Y()/1e3;
  TGraph *pFibDir = new TGraph(1,&ptX,&ptY);      // expected light position
  pFibDir->SetMarkerStyle(trueMarker);
  pFibDir->SetMarkerColor(1);
  pFibDir->SetMarkerSize(2);
  
  fibrepos.RotateZ(-rot_X2);
  fibrepos.RotateY(-rot_Z2);
  ptX = fibrepos.X()/1e3;
  ptY = fibrepos.Y()/1e3;
  TGraph *pFibPos = new TGraph(1,&ptX,&ptY);      // expected fibre position
  pFibPos->SetMarkerStyle(trueMarker);
  pFibPos->SetMarkerColor(1);
  pFibPos->SetMarkerSize(1.5);
  
  guess_dir.RotateZ(-rot_X1);
  guess_dir.RotateY(-rot_Z1);
  double guessD_rotX = guess_dir.X()/1e3;
  double guessD_rotY = guess_dir.Y()/1e3;
  
  guess_ref.RotateZ(-rot_X2);
  guess_ref.RotateY(-rot_Z2);
  double guessR_rotX = guess_ref.X()/1e3;
  double guessR_rotY = guess_ref.Y()/1e3;
  
  // View is centered over weighted light spot
  const int nil = 0;
  TGraph *pWgtDir = new TGraph(1,&nil,&nil);         // direct light (weighted)
  pWgtDir->SetMarkerStyle(weightMarker);
  pWgtDir->SetMarkerColor(1);
  pWgtDir->SetMarkerSize(2);
  
  TGraph *pWgtRef = new TGraph(1,&nil,&nil);         // reflected light (weighted)
  pWgtRef->SetMarkerStyle(weightMarker);
  pWgtRef->SetMarkerColor(1);
  pWgtRef->SetMarkerSize(1.5);
  
  // No rotation required, fit was performed in rotated view
  double fitD_rotX = fitDir.X()/1e3;
  double fitD_rotY = fitDir.Y()/1e3;
  double fitR_rotX = fitRef.X()/1e3;
  double fitR_rotY = fitRef.Y()/1e3;

  double pcircDx[2][NDOTS], pcircDy[2][NDOTS], pcircRx[2][NDOTS], pcircRy[2][NDOTS];
  int nfrontD=0, nbackD=0, nfrontR=0, nbackR=0;
  for (int d=0; d<NDOTS; d++) {
    // Rotate circle around direct spot
    dotsD[d]->RotateZ(-rot_X1);
    dotsD[d]->RotateY(-rot_Z1);
    if(dotsD[d]->Z() > 0) {
      pcircDx[0][nfrontD]=dotsD[d]->X()/1e3;
      pcircDy[0][nfrontD]=dotsD[d]->Y()/1e3;
      nfrontD++;
    } else {
      pcircDx[1][nbackD]=dotsD[d]->X()/1e3;
      pcircDy[1][nbackD]=dotsD[d]->Y()/1e3;
      nbackD++;
    }
    // Rotate circle around reflected spot
    dotsR[d]->RotateZ(-rot_X2);
    dotsR[d]->RotateY(-rot_Z2);
    if(dotsR[d]->Z() > 0) {
      pcircRx[0][nfrontR]=dotsR[d]->X()/1e3;
      pcircRy[0][nfrontR]=dotsR[d]->Y()/1e3;
      nfrontR++;
    } else {
      pcircRx[1][nbackR]=dotsR[d]->X()/1e3;
      pcircRy[1][nbackR]=dotsR[d]->Y()/1e3;
      nbackR++;
    }
  }
  TGraph *pcircD1=NULL, *pcircD2=NULL, *pcircR1=NULL, *pcircR2=NULL;
  if(nfrontD>0) { pcircD1 = new TGraph(nfrontD,pcircDx[0],pcircDy[0]);
                  pcircD1->SetMarkerStyle(6); }
  if(nbackD>0)  { pcircD2 = new TGraph(nbackD,pcircDx[1],pcircDy[1]);
                  pcircD2->SetMarkerStyle(1); }
  if(nfrontR>0) { pcircR1 = new TGraph(nfrontR,pcircRx[0],pcircRy[0]);
                  pcircR1->SetMarkerStyle(6); }
  if(nbackR>0)  { pcircR2 = new TGraph(nbackR,pcircRx[1],pcircRy[1]);
                  pcircR2->SetMarkerStyle(1); }
  
  
  // ********************************************************************
  // Plotting section
  // ********************************************************************
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
  
  // ----------
  // Run summary
  pad0->cd();
  TLatex *title, *t[6], *v[6];
  float avgnhit = (float)totalnhit/count;
  char* text[6] = { Form(" - Run number:"),
                    Form(" - Fibre name:"),
                    Form(" - Number of events:"),
                    Form(" - Average NHit:"),
                    Form(" - Fit deviation (dir.):"),
                    Form(" - Fit deviation (refl.):")
                  };
  char* vals[6] = { Form("%d",run),
                    Form("%s",fibre.c_str()),
                    Form("%d",count),
                    Form("%.2f",avgnhit),
                    Form("%.2f#circ",DIRANG/pi*180),
                    Form("%.2f#circ",REFANG/pi*180)
                  };
  
  // Highlight unsatisfactory values in bold
  int nexpected = (isMC ? 1e5 : 2e5);
  if (count < 0.99*nexpected) vals[2] = Form("#bf{%s}",vals[2]);     // lost more than 1% of triggers
  if (avgnhit<25 || avgnhit>45) vals[3] = Form("#bf{%s}",vals[3]);   // nhit far outside optimal range (30-40)
  if (DIRANG/pi*180 >= 10.) vals[4] = Form("#bf{%s}",vals[4]); // bad direct light fit
  if (REFANG/pi*180 >= 10.) vals[5] = Form("#bf{%s}",vals[5]); // bad reflected light fit
  
  // Place text in pad
  if (isMC) title = new TLatex(0.05,0.9,"SNO+ TELLIE MC simulation");
  else      title = new TLatex(0.05,0.9,"SNO+ TELLIE PCA data");
  title->SetTextAlign(12);
  title->SetTextFont(62);
  title->SetTextSize(0.12);
  title->Draw();
  for (int l=0; l<6; l++) {
    t[l] = new TLatex(0.05,0.7-0.12*l,text[l]);
    t[l]->SetTextAlign(11);
    t[l]->SetTextFont(82);
    t[l]->SetTextSize(0.08);
    t[l]->Draw();
    v[l] = new TLatex(0.75,0.7-0.12*l,vals[l]);
    v[l]->SetTextAlign(11);
    v[l]->SetTextFont(82);
    v[l]->SetTextSize(0.08);
    v[l]->Draw();
  }
  
  // Indicate colour scale
  TLatex *txtD = new TLatex(9.5,9.5,Form("NHit/PMT #leq %d",nearmax));
  txtD->SetTextAlign(33);
  txtD->SetTextFont(82);
  txtD->SetTextSize(0.04);
  TLatex *txtR = new TLatex(9.5,9.5,Form("NHit/PMT #leq %d",farmax));
  txtR->SetTextAlign(33);
  txtR->SetTextFont(82);
  txtR->SetTextSize(0.04);
  
  // ----------
  // PMT hit count histogram
  pad1->cd()->SetGrid();
  pad1->SetLogx();
  pad1->SetLogy();
  hhits->SetMinimum(0.5);
  hhits->SetMaximum(5e3);
  hhits->SetTitle("PMT hit count;NHit/PMT;NPMTs");
  hhits->SetLineWidth(2);
  hhits->GetXaxis()->SetTitleOffset(1.3);
  hhits->Draw();
  
  // Draw off and hot PMT bins in their respective colour
  TH1I *hhitslo = (TH1I*)hhits->Clone("hhitslo");
  hhitslo->SetFillColor(16);
  hhitslo->GetXaxis()->SetRangeUser(0.4,0.8);
  hhitslo->Draw("same");
  TH1I *hhitshi = (TH1I*)hhits->Clone("hhitshi");
  hhitshi->SetFillColor(1);
  double hotlimedge;
  for (int b=0; b<hhitshi->GetNbinsX(); b++) {
    if(hhitshi->GetBinCenter(b)<HOTLIMIT) continue;
    hotlimedge = hhitshi->GetBinCenter(b);
    break;
  }
  hhitshi->GetXaxis()->SetRangeUser(hotlimedge,2e5);
  hhitshi->Draw("same");
  
  // Draw line indicating "hot PMT" limit
  TLine *lhot = new TLine(HOTLIMIT,0,HOTLIMIT,5e3);
  lhot->SetLineWidth(2);
  lhot->SetLineColor(2);
  lhot->Draw("same");
  
  // ----------
  // Icosahedral projection of detector display
  pad2->cd()->SetGrid();
  hicos->SetTitle("Detector display (PMT hit sum)");
  //hicos->Draw();
  for(int s=0;s<NCOL+2;s++) { 
    if(!icos[s]) continue;
    if(icos[s]->GetN()==0) continue;
    icos[s]->Draw("P same");
  }
  //phot->Draw("P same");
  pGssDir->Draw("P same");
  pGssRef->Draw("P same");
  pFitDir->Draw("P same");
  pFitRef->Draw("P same");
  pcontD->Draw("P same");
  pcontR->Draw("P same");
  
  // ----------
  // View from direct light spot (fitted)
  pad3->cd()->SetGrid();
  hfineD->SetTitle("Direct light (PMT hit sum);X' [m];Y' [m]");
  hfineD->GetXaxis()->SetTitleOffset(1.3);
  hfineD->GetYaxis()->SetTitleOffset(1.4);
  hfineD->Draw("scat");
  hfineD->SetStats(0);
  for(int s=0;s<NCOL+2;s++) { if(!gDir[s]) continue; if(gDir[s]->GetN()==0) continue; gDir[s]->Draw("P same"); }
  pFibDir->Draw("P same");
  pGssDir->DrawGraph(1,&guessD_rotX,&guessD_rotY,"P same");
  pWgtDir->Draw("P same");
  pFitDir->DrawGraph(1,&fitD_rotX,&fitD_rotY,"P same");
  if(pcircD1) pcircD1->Draw("P same");
  if(pcircD2) pcircD2->Draw("P same");
  txtD->Draw();  
  
  // ----------
  // View from reflected light spot (fitted)
  pad4->cd()->SetGrid();
  hfineR->SetTitle("Reflected light (PMT hit sum);X'' [m];Y'' [m]");
  hfineR->GetXaxis()->SetTitleOffset(1.3);
  hfineR->GetYaxis()->SetTitleOffset(1.4);
  hfineR->Draw("scat");
  hfineR->SetStats(0);
  for(int s=0;s<NCOL+2;s++) { if(!gRef[s]) continue; if(gRef[s]->GetN()==0) continue; gRef[s]->Draw("P same"); }
  pFibPos->Draw("P same");
  pGssRef->DrawGraph(1,&guessR_rotX,&guessR_rotY,"P same");
  pWgtRef->Draw("P same");
  pFitRef->DrawGraph(1,&fitR_rotX,&fitR_rotY,"P same");
  if(pcircR1) pcircR1->Draw("P same");
  if(pcircR2) pcircR2->Draw("P same");
  txtR->Draw();  
  
  // ----------
  // Save canvas and close
  string outfile = "focal_point";
  if (!TEST) outfile = Form("images/PCA_%s",fibre.c_str());
  if (isMC) outfile = Form("mc/MC_%s_%s",fibre.c_str(),mcopt.c_str());
  c0->Print(Form("%s.png",outfile.c_str()));
  c0->Print(Form("%s.pdf",outfile.c_str()));
  c0->Close();
  
  // Delete pointers created with 'new'
  if(hicos) delete hicos;
  if(hcoarse) delete hcoarse;
  if(hfineD) delete hfineD;
  if(hfineR) delete hfineR;
  if(hhits) delete hhits;
  if(hhitslo) delete hhitslo;
  if(hhitshi) delete hhitshi;
  //if(gDir2D) delete gDir2D;
  //if(gRef2D) delete gRef2D;
  if (icos) { for(int s=0; s<NCOL+2; s++) delete icos[s]; }
  if (gDir) { for(int s=0; s<NCOL+2; s++) delete gDir[s]; }
  if (gRef) { for(int s=0; s<NCOL+2; s++) delete gRef[s]; }
  if(pcontD) delete pcontD;
  if(pcontR) delete pcontR;
  if(pGssDir) delete pGssDir;
  if(pGssRef) delete pGssRef;
  if(pFibPos) delete pFibPos;
  if(pFibDir) delete pFibDir;
  if(pWgtDir) delete pWgtDir;
  if(pWgtRef) delete pWgtRef;
  if(pcircD1) delete pcircD1;
  if(pcircD2) delete pcircD2;
  if(pcircR1) delete pcircR1;
  if(pcircR2) delete pcircR2;
  if(c0) delete c0;
  if(title) delete title;
  if(t) for(int l=0;l<6;l++) delete t[l];
  if(v) for(int l=0;l<6;l++) delete v[l];
  if(txtD) delete txtD;
  if(txtR) delete txtR;
  if(lhot) delete lhot;
  
}
