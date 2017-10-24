// ---------------------------------------------------------
// Goal:            Plot angular response of TELLIE fibres
// Author:          Martti Nirkko, University of Sussex (Oct 2017)
// Compile & run:   clear && g++ -o angular.exe angular.C `root-config --cflags --libs` -I$RATROOT/include/libpq -I$RATROOT/include -L$RATROOT/lib -lRATEvent_Linux && ./angular.exe
// ---------------------------------------------------------

// C++ stuff
#include <fstream>
#include <iostream>
#include <stdio.h>
#include <stdlib.h>

// ROOT stuff
#include <TFile.h>
#include <TLegend.h>
#include <TPaveStats.h>
#include <TStyle.h>
#include <TSystem.h>

// RAT stuff
#include <RAT/DS/Entry.hh>
#include <RAT/DS/MC.hh>
#include <RAT/DS/Run.hh>
#include <RAT/DU/DSReader.hh>
#include <RAT/DU/PMTInfo.hh>
#include <RAT/DU/Utility.hh>

// Helper functions
#include "../HelperFunc.C"
#include "../Xianguo.C"

// Run time parameters
const int RUN_CLUSTER = 0;  // whether running on cluster (0=local)
const int USE_RATDB = 1;    // whether to use RATDB to get fibre positions (1=yes)
const int VERBOSE = 1;      // verbosity flag
const int IS_MC = 0;        // Monte-Carlo flag 

// Initialise functions
float angular(string, int, bool, bool);

// Main program
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
  
  int nfiles=0; 
  while (true) {
    in >> node >> fibre >> channel >> run >> ipw >> photons >> pin >> rms >> nhit;
    if (!in.good()) break;
    if (TEST && TEST!=run) continue; // only want specified run
    float avgnhit = angular(fibre, run, IS_MC, TEST);
    if (avgnhit==-2.) { printf("Could not find data for fibre %s!\n",fibre.c_str()); continue; }
    if (avgnhit==-1.) { if (VERBOSE) printf("Nothing to do!\n"); continue; }
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
  if(!TEST && VERBOSE) printf("*****\n");
  printf("Checking files for run %d... ", run);
  string fpath = Form("%s/Software/SNOP/work/data",getenv("HOME"));
  string fname = "";
  ifstream f;
  for (int pass=3;pass>=0;pass--) {
    fname = Form("%s/Analysis_r0000%d_s000_p00%d.root",fpath.c_str(),run,pass);
    f.open(fname.c_str());
    if (f.good()) break;
  }
  string out = Form("./output/angular_%s.root",fibre.c_str());
  string img = Form("./images/angular_%s.pdf",fibre.c_str());
  ifstream g(out.c_str());
  ifstream h(img.c_str());
  int scanned_file = 0;
  if (!TEST && h.good()) {        // file downloaded and processed
    printf("already processed! Skipping fibre %s.\n",fibre.c_str());
    return -1.;
  } else if(g.good()) {  // file extracted, but not processed
    printf("not processed! Generating plots for fibre %s.\n",fibre.c_str());
    scanned_file = 1;
  } else if(!f.good()) {          // file not downloaded
    printf("not downloaded! Skipping fibre %s.\n",fibre.c_str());
    return -2.;
  } else {                        // file downloaded, but not processed
    printf("OK. Extracting data for fibre %s.\n",fibre.c_str());
  }

  string outpath = Form("output/angular_%s.root",fibre.c_str());
  TFile *outfile = NULL;
  int totalnhit, count;
  TH1D *hmeantime=NULL, *hrmstime=NULL;
  TH1D *h1=NULL, *herr=NULL, *hcoarse=NULL, *hpmtseg=NULL;
  TH1D *h2_0=NULL, *h2_1=NULL, *h2_2=NULL, *h2_3=NULL, *h2_py=NULL;
  TH2D *h2=NULL;
  TGraphErrors *pmtfits=NULL;
  if (scanned_file) {
    outfile = new TFile(outpath.c_str(),"READ");
    h1 = (TH1D*)outfile->Get("h1"); // Time vs angle
    h2 = (TH2D*)outfile->Get("h2"); // Time vs angle
    pmtfits = (TGraphErrors*)outfile->Get("pmtfits"); // Fitted prompt peaks
    hmeantime = (TH1D*)outfile->Get("hmeantime"); // Binned PMT mean times
    hrmstime = (TH1D*)outfile->Get("hrmstime"); // Binned PMT errors on mean
    herr = (TH1D*)outfile->Get("herr"); // Error on mean time
    hpmtseg = (TH1D*)outfile->Get("hpmtseg"); // Number of PMTs in angular bin
    h2_0 = (TH1D*)outfile->Get("h2_0");  // Intensity
    h2_1 = (TH1D*)outfile->Get("h2_1");  // Mean
    h2_2 = (TH1D*)outfile->Get("h2_2");  // RMS
    h2_3 = (TH1D*)outfile->Get("h2_chi2");  // Gaussian GOF
    totalnhit = h2_0->Integral(); // TODO - check if correct for all cases!
    count = h2_0->GetEntries()*h2_0->GetNbinsX(); // histogram normalised
  } else {
    outfile = new TFile(outpath.c_str(),"RECREATE");

    // Initialise RAT
    RAT::DU::DSReader dsreader(fname);
    const RAT::DU::ChanHWStatus& chs = RAT::DU::Utility::Get()->GetChanHWStatus();
    const RAT::DU::PMTInfo& pmtinfo = RAT::DU::Utility::Get()->GetPMTInfo();
    const int NPMTS = pmtinfo.GetCount();
    
    TVector3 fibrepos(0,0,0);
    TVector3 fibredir(0,0,0);
    TVector3 lightpos(0,0,0);
    if (USE_RATDB) {
      // Get fibre info (from RATDB) 
      RAT::DB *db = RAT::DB::Get();
      //db->LoadDefaults();	  // Already done when calling DU::Utility::Get()
      RAT::DBLinkPtr entry = db->GetLink("FIBRE",fibre);
      fibrepos.SetXYZ(entry->GetD("x"), entry->GetD("y"), entry->GetD("z")); // position of fibre [mm]
      fibredir.SetXYZ(entry->GetD("u"), entry->GetD("v"), entry->GetD("w")); // direction of fibre
      lightpos = fibrepos + 2*fibrepos.Mag()*fibredir;                       // projected light spot centre
      if (VERBOSE)  cout << "RATDB: fibre " << fibre << ", pos " << printVector(fibrepos) << ", dir " << printVector(fibredir) << endl;
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
      fibrepos.SetXYZ(10*fx,10*fy,10*fz);                // position of fibre [mm]
      fibredir.SetXYZ(fu,fv,fw);                         // direction of fibre
      lightpos = fibrepos + 2*fibrepos.Mag()*fibredir;   // projected light spot centre 
      if (VERBOSE) cout << "DocDB: fibre " << fibre << ", pos " << printVector(fibrepos) << ", dir " << printVector(fibredir) << endl;
    }

    // TELLIE fibre/trigger delays & Mark's PCA offsets
    ifstream del("TELLIE_delays.txt");
    if (!del) { cerr<<"ERROR - could not open TELLIE_delays.txt!"<<endl; return -999.; }
    int num, trig_delay;
    float fibre_delay, pca_offset;
    string line;
    for (int hdr=0; hdr<2; hdr++) {
      getline(del,line);      // header
    }
    while (true) {
      del >> num >> fibre_delay >> trig_delay >> pca_offset;
      if (!del.good()) break;
      if (num == run) break;
    }
    if (VERBOSE) cout<<"TELLIE: trigger delay "<<trig_delay<<" ns, fibre delay "<<fibre_delay<<" ns, total offset "<<pca_offset<<" ns."<<endl;
    
    // Get fitted light position (from file)
    string fitfile = "../fibre_validation/TELLIE_FITRESULTS.txt";
    ifstream fit(fitfile.c_str());
    if (!fit) { cerr<<"Failed to open "<<fitfile<<endl; exit(1); }
    string fibrecheck;
    TVector3 fitpos(0,0,0);
    int runcheck;
    float dirx, diry, dirz, refx, refy, refz;
    for (int hdr=0; hdr<2; hdr++) {
      getline(fit,line);      // header
    }
    while (true) {
      fit >> runcheck >> fibrecheck >> dirx >> diry >> dirz >> refx >> refy >> refz;
      if (!fit.good()) break;
      if (runcheck != run) continue;
      if (fibrecheck != fibre) { cerr<<"Fibre mismatch: "<<fibrecheck<<" != "<<fibre<<" Check files."<<endl; exit(1); }
      fitpos.SetXYZ(dirx,diry,dirz);
      if (fitpos.Mag()!=0) break;
    }
    TVector3 fitdir = (fitpos-fibrepos).Unit();
    cout << "Loaded fit position: " << printVector(fitpos) << " mm, " << fitpos.Angle(lightpos)*180./pi << " deg deviation" << endl;
    
    // Initialise histograms
    const int NBINS = 24;
    const int MAXANG = 24;
    hcoarse = new TH1D("hcoarse",fibre.c_str(),300,-100,200);
    hpmtseg = new TH1D("hpmtseg",fibre.c_str(),NBINS,0,MAXANG);
    
    // First iteration: Get average hit time and number of PMTs in each segment
    int nevents=0, hitpmts=0;
    float allpmts[NPMTS], angpmts[NPMTS];
    memset( allpmts, 0, NPMTS*sizeof(float) );                // NPMTS only known at runtime
    memset( angpmts, 0, NPMTS*sizeof(float) );                // NPMTS only known at runtime
    TH2D *htime = new TH2D("htime","",NPMTS+1,0,NPMTS+1,600,-300,300);
    cout << "Looping over entries..." << endl;
    for(int iEntry=0; iEntry<dsreader.GetEntryCount(); iEntry++) {
      // Print progress
      if (iEntry % int(dsreader.GetEntryCount()/100.) == 0) {
        printProgress(iEntry, dsreader.GetEntryCount());
      }
      const RAT::DS::Entry& ds = dsreader.GetEntry(iEntry);
      for(int iEv=0; iEv<ds.GetEVCount(); iEv++) {            // mostly 1 event per entry
        const RAT::DS::EV& ev = ds.GetEV(iEv);
        int trig = ev.GetTrigType();
        if (!(trig & 0x8000)) continue;                       // EXT trigger only
        nevents++;
        const RAT::DS::CalPMTs& pmts = ev.GetCalPMTs();
        for(int iPMT=0; iPMT<pmts.GetNormalCount(); iPMT++) {
          RAT::DS::PMTCal pmt = pmts.GetNormalPMT(iPMT);
          int pmtID = pmt.GetID();
          if (!chs.IsTubeOnline(pmtID)) continue;             // test CHS
          if (pmt.GetCrossTalkFlag()) continue;               // remove crosstalk
          if (pmt.GetStatus().GetULong64_t(0) != 0) continue; // test PCA / ECA
          TVector3 pmtpos = pmtinfo.GetPosition(pmtID);
          TVector3 track = pmtpos-fibrepos;
          double theta = track.Angle(fitdir);                 // angle w.r.t. fibre [rad]
          double angle = theta*180./pi;                       // angle w.r.t. fibre [deg]
          double pmttime = pmt.GetTime();                     // hit time [ns]
          double corr = track.Mag()/c_water;                  // light travel time [ns]
          double offset = pmttime-corr+fibre_delay+trig_delay;
          if (angle > MAXANG) continue;
          hcoarse->Fill(offset-pca_offset);   // total offset minus PCA offset
          //htime->Fill(pmtID, pmttime-corr);   // corrected hit time vs PMT ID
          htime->Fill(pmtID, offset-pca_offset);  // corrected hit time vs PMT ID
          if (allpmts[pmtID]==0) {            // fill angles only once per PMT
            hpmtseg->Fill(angle);
            angpmts[pmtID] = angle;
          }
          allpmts[pmtID] += 1.;
          hitpmts++;
        } // pmt loop
      } // event loop
    } // entry loop
    cout << endl;

    // Turn array of PMTs into percentage (occupancy)
    cout << "Number of events was " << nevents << endl;
    for (int iPMT=0; iPMT<NPMTS; iPMT++) {
      allpmts[iPMT] /= nevents;
    }
    
    // Fit 1D Gaussian over PMT hit times (prompt peak)
    float fitgaus[4] = {0}; // fit results
    pmtfits = FitPromptPeaks(run, htime, NPMTS, allpmts, angpmts, fitgaus);
    if (!pmtfits) { cout << "Error fitting PMT prompt peaks!" << endl; return -999; }
    printf("RESULTS FOR FIT TIME: Constant ( %f +/- %f ) ns, Slope ( %f +/- %f ) ns/deg\n", 
            fitgaus[0], fitgaus[1], fitgaus[2], fitgaus[3]);
    //float tx=fitgaus[0], tex=fitgaus[1], ty=fitgaus[2], tey=fitgaus[3];
    //TGraphErrors *gfitgaus = new TGraphErrors(1,&tx,&ty,&tex,&tey);
    //return 0;
    
    // -----
    
    hmeantime = new TH1D("hmeantime","",NBINS,0,NBINS);
    hrmstime = new TH1D("hrmstime","",NBINS,0,NBINS);
    double *tx = pmtfits->GetX();
    double *ty = pmtfits->GetY();
    double *tex = pmtfits->GetEX();           // should always be zero
    double *tey = pmtfits->GetEY();
    int tb, tn[NBINS] = {0};
    double tavg[NBINS] = {0};
    double terr[NBINS] = {0};
    // Loop over all PMTs in the graph and bin the data
    for (int i=0; i<pmtfits->GetN(); i++) {
      tb = (int)tx[i];                        // bin index (1 degree bins)
      tn[tb]++;                               // number of PMTs in each bin
      //tavg[tb] += ty[i];                      // sum of hit times
      tavg[tb] += ty[i]/(tey[i]*tey[i]);      // sum of weighted hit times
      terr[tb] += 1./(tey[i]*tey[i]);         // sum of errors squared
    }
    // Loop over each bin to get the average and propagated error
    for (int j=0; j<NBINS; j++) {
      if (tn[j]==0) continue;                 // no PMTs in this bin
      tavg[j] = tavg[j]/terr[j];              // weighted mean hit time
      terr[j] = sqrt(1./terr[j]);             // error on weighted mean
      //else terr[j] = sqrt(terr[j]/(tn[j]-1)); // standard deviation
      //printf("Bin %2d, Hit time %8.3lf ns, error %8.3lf ns.\n",j,tavg[j],terr[j]);
      hmeantime->SetBinContent(j+1,tavg[j]);
      hmeantime->SetBinError(j+1,terr[j]);
      hrmstime->SetBinContent(j+1,terr[j]);
    }
    
    // TODO - Replace times below with these fit results!
    
    // -----
    
    // TODO - Update plots!
    
    
    // Central hit time for 2D plot
    //int commontime = hcoarse->GetXaxis()->GetBinLowEdge(hcoarse->GetMaximumBin());
    //int commontime = round(fitgaus[0] + fibre_delay + trig_delay - pca_offset);
    int commontime = round(fitgaus[0] + fitgaus[2]*12.0);
    if (VERBOSE) printf("Most common hit time (after fit): %d ns\n",commontime);
    commontime += 5/2;              // intermediate step
    commontime -= commontime % 5;   // rounded to the nearest multiple of 5
    
    // More histograms
    h1 = new TH1D("h1",fibre.c_str(),NBINS,0,MAXANG);
    h2 = new TH2D("h2",fibre.c_str(),NBINS,0,MAXANG,NBINS,commontime-15,commontime+15);
    herr = new TH1D("herr",fibre.c_str(),NBINS,0,MAXANG);
    
    // Second iteration: Fill histograms
    totalnhit=0;
    count=0;
    for(int iEntry=0; iEntry<dsreader.GetEntryCount(); iEntry++) {
      const RAT::DS::Entry& ds = dsreader.GetEntry(iEntry);
      RAT::DS::MC mc;
      if (isMC) mc = ds.GetMC();  // don't initialise this for real data (crashes)
      
      // Loop over triggered events in each entry
      for(int iEv=0; iEv<ds.GetEVCount(); iEv++) {              // mostly 1 event per entry
        const RAT::DS::EV& ev = ds.GetEV(iEv);
        
        // Trigger type
        int trig = ev.GetTrigType();
        if (!(trig & 0x8000)) continue;                         // EXT trigger only
        
        // Event observables
        int nhitscleaned = ev.GetNhitsCleaned();   // calibrated PMT hits with removed crosstalk
        totalnhit += nhitscleaned;
        count++;
        
        // PMT information
        const RAT::DS::CalPMTs& pmts = ev.GetCalPMTs();
        for(int iPMT=0; iPMT<pmts.GetNormalCount(); iPMT++) {
          RAT::DS::PMTCal pmt = pmts.GetNormalPMT(iPMT);
          int pmtID = pmt.GetID();
          if (!chs.IsTubeOnline(pmtID)) continue;               // test CHS
          if (pmt.GetCrossTalkFlag()) continue;                 // remove crosstalk
          if (pmt.GetStatus().GetULong64_t(0) != 0) continue;   // test PCA
          TVector3 pmtpos = pmtinfo.GetPosition(pmtID);
          TVector3 track = pmtpos-fibrepos;
          double theta = track.Angle(fitdir);                   // angle w.r.t. fibre [rad]
          double angle = theta*180./pi;                         // angle w.r.t. fibre [deg]
          double pmttime = pmt.GetTime();                       // hit time [ns]
          double corr = track.Mag()/c_water;                    // light travel time [ns]
          double offset = pmttime-corr+fibre_delay+trig_delay;  // total offset
          if (angle > MAXANG) continue;
          // TODO - only fill histogram if pmttime is within prompt peak?
          for (int iFit=0; iFit<pmtfits->GetN(); iFit++) {
            if (fabs(tx[iFit]-angle)>1e-6) continue;
            //if (tx[iFit] != angle) continue;
            //cout << "Compared " << angle << " to " << tx[iFit] << endl;
            if (fabs(offset-pca_offset - ty[iFit]) > 8.) continue; // prompt peak estimate (TODO - get proper Gaussian width somehow)
            h2->Fill(angle, offset-pca_offset);                   // angle [deg], time [ns]
            break;
          }
        } // pmt loop
      } // event loop
    } // entry loop
    
    // Fit angular slices with a gaussian distribution in time
    h2->FitSlicesY();
    h2_0 = (TH1D*)gDirectory->Get("h2_0");  // Constant
    h2_1 = (TH1D*)gDirectory->Get("h2_1");  // Mean
    h2_2 = (TH1D*)gDirectory->Get("h2_2");  // StdDev
    h2_3 = (TH1D*)gDirectory->Get("h2_chi2");  // chi^2/ndof  
    h2_py = NULL;
    
    // Put mean and RMS into 1D histogram
    int binscone=0;
    float avgmean=0., avgdev=0.;
    for (int b=0; b<NBINS+2; b++) {
      h1->SetBinContent(b, h2_1->GetBinContent(b));
      h1->SetBinError(b, h2_2->GetBinContent(b));
      h2->ProjectionY("_py",b,b);
      h2_py = (TH1D*)gDirectory->Get("h2_py");
      h2_0->SetBinContent(b, h2_py->GetSum());  // overwrite Gaussian constant with intensity
      h2_0->SetBinError(b, sqrt(h2_py->GetSum()));  // overwrite Gaussian constant error with sqrt(intensity)
      h2_1->SetBinContent(b, h2_py->GetMean()); // overwrite Gaussian mean with mean
      h2_1->SetBinError(b, h2_py->GetMeanError()); // overwrite Gaussian mean error with mean error
      h2_2->SetBinContent(b, h2_py->GetRMS());  // overwrite Gaussian sigma with RMS
      h2_2->SetBinError(b, h2_py->GetRMSError());  // overwrite Gaussian sigma error with RMS error
      herr->SetBinContent(b, h2_py->GetMeanError()); // not *quite* the same as error on Gaussian mean
      if (h1->GetBinCenter(b)>0 && h1->GetBinCenter(b)<MAXANG) {  // within 24 deg cone
        avgmean += h2_1->GetBinContent(b);
        avgdev += h2_2->GetBinContent(b);
        binscone++;
      }
    }
    avgmean/=binscone;
    avgdev/=binscone;
    if (VERBOSE) printf("Average hit time mean (< %d deg): %6.2f\n", MAXANG, avgmean);
    if (VERBOSE) printf("Average hit time RMS (< %d deg):  %6.2f\n", MAXANG, avgdev);
    
    // Normalise angular slices to number of PMTs in slice
    double tmp, pmts;
    for (int binx=0; binx<NBINS+2; binx++) {
      pmts = hpmtseg->GetBinContent(binx);
      tmp = h2_0->GetBinContent(binx);
      if (pmts<=0) { h2_0->SetBinContent(binx,0); h2_0->SetBinError(binx,0); }
      else { h2_0->SetBinContent(binx,tmp/pmts); h2_0->SetBinError(binx,sqrt(tmp)/pmts); }
      //if (VERBOSE) printf("Segment %2d contains %4d PMTs.\n",binx,(int)pmts);
      for (int biny=0; biny<NBINS+2; biny++) {
        tmp = h2->GetBinContent(binx,biny);
        if (pmts<=0) h2->SetBinContent(binx,biny,0);
        else h2->SetBinContent(binx,biny,tmp/pmts);
      }
    }
    
    // Write all objects to file
    pmtfits->Write("pmtfits");
    outfile->Write();
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
  
  // Mean
  pad0->cd()->SetGrid();
  pad0->SetLeftMargin(0.12);
  //pad1->DrawFrame(0,-5,24,5);
  pmtfits->SetTitle("Fitted PMT hit times (prompt peaks);Angle [deg];Time [ns]");
  pmtfits->SetMarkerColor(4);
  pmtfits->SetMarkerStyle(6);
  pmtfits->Draw("AP same");
  pmtfits->GetXaxis()->SetTitleOffset(1.2);
  pmtfits->GetYaxis()->SetTitleOffset(1.5);
  pmtfits->GetXaxis()->SetRangeUser(0,24);
  /*
  h2_1->SetTitle("Mean hit time;Angle [deg];Time [ns]");
  h2_1->SetLineWidth(2);
  h2_1->Draw();
  h2_1->GetXaxis()->SetTitleOffset(1.2);
  h2_1->GetYaxis()->SetTitleOffset(1.7);
  float minval = h2_1->GetMaximum();
  float maxval = h2_1->GetMinimum();
  float val;
  for (int b=1; b<=NBINS; b++) { 
    val = h2_1->GetBinContent(b);
    if (val==0.) continue;
    if (minval>val) minval=val;
    if (maxval<val) maxval=val;
  }
  h2_1->GetYaxis()->SetRangeUser(0.75*minval, 1.25*maxval);
  */
  
  // Normalised intensity profile as time vs. angle (2D)
  pad1->cd()->SetGrid();
  pad1->SetRightMargin(0.12);   // for TH2D color scale
  h2->SetTitle(Form("Intensity profile %s (norm.);Angle [deg];Time [ns]",fibre.c_str()));
  h2->Draw("colz");
  h2->GetXaxis()->SetTitleOffset(1.2);
  h2->GetYaxis()->SetTitleOffset(1.5);
  
  // Mean hit times (binned PMT fits)
  pad2->cd()->SetGrid();
  pad2->DrawFrame(0,0,24,5,"Binned PMT mean times (new);Angle [deg];Time [ns]");
  hmeantime->SetLineWidth(2);
  hmeantime->Draw("same");
  hmeantime->GetXaxis()->SetTitleOffset(1.2);
  hmeantime->GetYaxis()->SetTitleOffset(1.7);
  /*
  float minval = hmeantime->GetMaximum();
  float maxval = hmeantime->GetMinimum();
  float val;
  for (int b=1; b<=NBINS; b++) { 
    val = hmeantime->GetBinContent(b);
    if (val==0.) continue;
    if (minval>val) minval=val;
    if (maxval<val) maxval=val;
  }
  hmeantime->GetYaxis()->SetRangeUser(0.75*minval, 1.25*maxval);
  */
  /*
  h2_2->SetTitle("RMS hit time;Angle [deg];Time [ns]");
  h2_2->SetLineWidth(2);
  h2_2->Draw();
  h2_2->GetXaxis()->SetTitleOffset(1.2);
  h2_2->GetYaxis()->SetTitleOffset(1.4);
  h2_2->GetYaxis()->SetRangeUser(3.8, 6.2);
  */
  
  // Mean errors (fitted PMT widths added in quadrature)
  pad3->cd()->SetGrid();
  //pad3->DrawFrame(0,0,24,0.2,"Binned PMT errors;Angle [deg];Time [ns]");
  hrmstime->SetTitle("Binned PMT errors (new);Angle [deg];Time [ns]");
  hrmstime->SetLineWidth(2);
  hrmstime->Draw();
  hrmstime->GetXaxis()->SetTitleOffset(1.2);
  hrmstime->GetYaxis()->SetTitleOffset(1.7);
  hrmstime->GetYaxis()->SetRangeUser(0, 1.1*hrmstime->GetMaximum());
  /*
  float minerr = hrmstime->GetMaximum();
  float maxerr = hrmstime->GetMinimum();
  float err;
  for (int b=1; b<=NBINS; b++) { 
    err = hrmstime->GetBinContent(b);
    if (val==0.) continue;
    if (minerr>val) minerr=err;
    if (maxerr<val) maxerr=err;
  }
  hrmstime->GetYaxis()->SetRangeUser(0.75*minerr, 1.25*maxerr);
  */
  /*
  h2_1->SetTitle("Mean hit time;Angle [deg];Time [ns]");
  h2_1->SetLineWidth(2);
  h2_1->Draw();
  h2_1->GetXaxis()->SetTitleOffset(1.2);
  h2_1->GetYaxis()->SetTitleOffset(1.7);
  float minval1 = h2_1->GetMaximum();
  float maxval1 = h2_1->GetMinimum();
  for (int b=1; b<=NBINS; b++) { 
    val = h2_1->GetBinContent(b);
    if (val==0.) continue;
    if (minval1>val) minval1=val;
    if (maxval1<val) maxval1=val;
  }
  h2_1->GetYaxis()->SetRangeUser(0.75*minval1, 1.25*maxval1);
  */
  
  // Error on mean (old method)
  pad4->cd()->SetGrid();
  pad4->SetLeftMargin(0.2);   // for axis label
  herr->SetTitle("Error on mean hit time (old);Angle [deg];#Delta t_{mean} [ns]");
  herr->SetLineWidth(2);
  herr->Draw();
  herr->GetXaxis()->SetTitleOffset(1.2);
  herr->GetYaxis()->SetTitleOffset(1.7);
  herr->GetYaxis()->SetRangeUser(0, 0.1);

  // Summed intensities (projected profile)
  pad5->cd()->SetGrid();
  pad5->SetLeftMargin(0.2);   // for axis label
  h2_0->SetTitle("Intensity (norm.);Angle [deg];Avg. NHit/PMT [ ]");
  h2_0->SetLineWidth(2);
  h2_0->Draw();
  h2_0->GetXaxis()->SetTitleOffset(1.2);
  h2_0->GetYaxis()->SetTitleOffset(1.7);
  h2_0->GetYaxis()->SetRangeUser(0, 1e4);
  hpmtseg->Scale(10);
  hpmtseg->SetLineWidth(2);
  hpmtseg->SetLineColor(2);
  hpmtseg->SetTitle("NPMTs (#times10)");
  hpmtseg->Draw("same");
  pad5->BuildLegend();
  
  // Chisquare/NDOF
  /*
  pad5->cd()->SetGrid();
  pad5->SetLogy();
  h2_3->SetTitle("Gaussian fit #chi^{2}/ndof;Angle [deg];");
  h2_3->SetLineWidth(2);
  h2_3->Draw();
  h2_3->GetXaxis()->SetTitleOffset(1.2);
  h2_3->GetYaxis()->SetRangeUser(0.5, 500);
  */
  
  // Save canvas and close
  string imgfile = "angular";
  if (!TEST) imgfile = Form("images/angular_%s",fibre.c_str());
  c0->Print(Form("%s.png",imgfile.c_str()));
  c0->Print(Form("%s.pdf",imgfile.c_str()));
  c0->Close();
 
  // Close files and free memory
  outfile->Close();
  if (outfile) delete outfile;
  if (c0) delete c0;
  
  return (float)totalnhit/count;
}

