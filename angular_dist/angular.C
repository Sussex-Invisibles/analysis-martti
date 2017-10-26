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
const int RUN_CLUSTER = 1;  // whether running on cluster (0=local)
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

  string outstr = Form("output/angular_%s.root",fibre.c_str());
  TFile *outfile = NULL;
  int hitpmts, nevents;
  TH1D *hpeak=NULL, *hpeakerr=NULL;
  TH1D *hmean=NULL, *hmeanerr=NULL;
  TH1D *hwdth=NULL, *hwdtherr=NULL;
  TH1D *herr=NULL, *hpmtseg=NULL, *hpmtgood=NULL;
  TH1D *h2_0=NULL, *h2_1=NULL, *h2_2=NULL, *h2_3=NULL, *h2_py=NULL;
  TH2D *h2=NULL;
  TGraph2DErrors *pmtfits=NULL;
  TGraphErrors *gpmts=NULL;
  if (scanned_file) {
    outfile = new TFile(outstr.c_str(),"READ");
    h2 = (TH2D*)outfile->Get("h2"); // Time vs angle
    pmtfits = (TGraph2DErrors*)outfile->Get("pmtfits"); // Fitted prompt peaks (3D)
    gpmts = (TGraphErrors*)outfile->Get("gpmts"); // Fitted prompt peaks (time vs angle)
    hpeak = (TH1D*)outfile->Get("hpeak");       // Binned intensities
    hpeakerr = (TH1D*)outfile->Get("hpeakerr"); // Binned errors on intensities
    hmean = (TH1D*)outfile->Get("hmean");       // Binned mean times
    hmeanerr = (TH1D*)outfile->Get("hmeanerr"); // Binned errors on mean times
    hwdth = (TH1D*)outfile->Get("hwdth");       // Binned mean signal widths
    hwdtherr = (TH1D*)outfile->Get("hwdtherr"); // Binned errors on signal widths
    herr = (TH1D*)outfile->Get("herr");         // Error on mean time
    hpmtseg = (TH1D*)outfile->Get("hpmtseg");   // PMTs in angular bin
    hpmtgood = (TH1D*)outfile->Get("hpmtgood"); // Good PMTs in angular bin
    h2_0 = (TH1D*)outfile->Get("h2_0");         // Intensity (normalised)
    h2_1 = (TH1D*)outfile->Get("h2_1");         // Mean
    h2_2 = (TH1D*)outfile->Get("h2_2");         // RMS
    h2_3 = (TH1D*)outfile->Get("h2_chi2");      // Gaussian GOF
    hitpmts = -999; // h2_0->Integral(); // TODO - fix this!
    nevents = 1; // h2_0->GetEntries()*h2_0->GetNbinsX(); // histogram normalised
  } else {
    outfile = new TFile(outstr.c_str(),"RECREATE");

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
    hpmtseg = new TH1D("hpmtseg",fibre.c_str(),NBINS,0,MAXANG);
    hpmtgood = new TH1D("hpmtgood",fibre.c_str(),NBINS,0,MAXANG);
    
    // Loop over all entries
    nevents=0;
    hitpmts=0;
    float allpmts[NPMTS], angpmts[NPMTS];
    memset( allpmts, 0, NPMTS*sizeof(float) );                // NPMTS only known at runtime
    memset( angpmts, 0, NPMTS*sizeof(float) );                // NPMTS only known at runtime
    TH2D *htime = new TH2D("htime","",NPMTS,0,NPMTS,200,-100,100);
    cout << "Looping over entries..." << endl;
    for(int iEntry=0; iEntry<dsreader.GetEntryCount(); iEntry++) {
      const RAT::DS::Entry& ds = dsreader.GetEntry(iEntry);
      RAT::DS::MC mc;
      if (isMC) mc = ds.GetMC();  // don't initialise this for real data (crashes)
      
      // Print progress
      if (iEntry % int(dsreader.GetEntryCount()/100.) == 0) {
        printProgress(iEntry, dsreader.GetEntryCount());
      }
      
      // Loop over triggered events in each entry
      for(int iEv=0; iEv<ds.GetEVCount(); iEv++) {            // mostly 1 event per entry
        const RAT::DS::EV& ev = ds.GetEV(iEv);
        
        // Trigger type
        int trig = ev.GetTrigType();
        if (!(trig & 0x8000)) continue;                       // EXT trigger only
        nevents++;
        
        // PMT information
        const RAT::DS::CalPMTs& pmts = ev.GetCalPMTs();
        for(int iPMT=0; iPMT<pmts.GetNormalCount(); iPMT++) {
          RAT::DS::PMTCal pmt = pmts.GetNormalPMT(iPMT);
          int pmtID = pmt.GetID();
          hitpmts++;
          
          // Cuts on PMT status
          if (!chs.IsTubeOnline(pmtID)) continue;             // test CHS
          if (pmt.GetCrossTalkFlag()) continue;               // remove crosstalk
          if (pmt.GetStatus().GetULong64_t(0) != 0) continue; // test PCA / ECA
          
          // Get PMT info
          TVector3 pmtpos = pmtinfo.GetPosition(pmtID);
          TVector3 track = pmtpos-fibrepos;
          double theta = track.Angle(fitdir);                 // angle w.r.t. fibre [rad]
          double angle = theta*180./pi;                       // angle w.r.t. fibre [deg]
          double pmttime = pmt.GetTime();                     // hit time [ns]
          double corr = track.Mag()/c_water;                  // light travel time [ns]
          double offset = pmttime-corr+fibre_delay+trig_delay;
          
          // Fill histograms/arrays
          if (angle > MAXANG) continue;
          htime->Fill(pmtID, offset-pca_offset);    // corrected hit time vs PMT ID
          if (allpmts[pmtID]==0) {                  // fill angles only once per PMT
            hpmtseg->Fill(angle);
            angpmts[pmtID] = angle;
          }
          allpmts[pmtID] += 1.;
          
        } // pmt loop
      } // event loop
    } // entry loop
    cout << endl;

    // Turn array of PMTs into percentage (occupancy)
    cout << "Number of events was " << nevents << endl;
    for (int iPMT=0; iPMT<NPMTS; iPMT++) {
      allpmts[iPMT] /= nevents;
    }
    
    // Fit PMT hit times (Gaussian prompt peak)
    pmtfits = new TGraph2DErrors();
    FitPromptPeaks(htime, NPMTS, allpmts, angpmts, pmtfits);
    if (!pmtfits) { cout << "*** ERROR *** Could not fit PMT prompt peaks!" << endl; return -999; }
    if (!pmtfits->GetN()) { cout << "*** ERROR *** Graph contains zero points!" << endl; return -999; }
    
    double *tx = pmtfits->GetX();   // Constant (A)
    double *ty = pmtfits->GetY();   // Mean value (mu)
    double *tz = pmtfits->GetZ();   // Uncertainty (sigma)
    double *tex = pmtfits->GetEX();
    double *tey = pmtfits->GetEY();
    double *tez = pmtfits->GetEZ();
  
    // Fill 1D graph with time vs angle (use only PMTs of interest)
    gpmts = new TGraphErrors();
    int npmts=0;
    for (int iPMT=0; iPMT<NPMTS; iPMT++) {
      if (tx[iPMT]+ty[iPMT]+tz[iPMT]==0) continue;
      gpmts->SetPoint(npmts, angpmts[iPMT], ty[iPMT]);
      gpmts->SetPointError(npmts, 0, tey[iPMT]);
      npmts++;
    }

    // Investigate PMTs with unusual offsets w.r.t. mean hit time
    string outlogstr = Form("logs/unusual_timing_%d.log", run);
    FILE *outlog = fopen(outlogstr.c_str(),"w");
    double meanhittime = gpmts->GetMean(2);
    double rmshittime = gpmts->GetRMS(2);
    fprintf(outlog,"# Mean %.1lf ns, RMS %.1lf ns. PMTs with unusually high offsets (>3*RMS):\n",meanhittime,rmshittime);
    fprintf(outlog,"# PmtId Offset[ns]\n");
    for (int iPMT=0; iPMT<NPMTS; iPMT++) {
      if (tx[iPMT]+ty[iPMT]+tz[iPMT]==0) continue;
      // Check for unusual offsets (seen e.g. for FT019A, PMTs 4393-4416)
      if (fabs(ty[iPMT]-meanhittime)/fabs(rmshittime) > 3.) { // more than 3x RMS deviation
        fprintf(outlog,"%4d %5.1lf\n",iPMT,ty[iPMT]-meanhittime);
      }
    }
    fclose(outlog);
    
    // Fill fit results into histograms
    hpeak = new TH1D("hpeak","",NBINS,0,NBINS);         // mean intensity
    hmean = new TH1D("hmean","",NBINS,0,NBINS);         // mean hit time
    hwdth = new TH1D("hwdth","",NBINS,0,NBINS);         // mean signal width
    hpeakerr = new TH1D("hpeakerr","",NBINS,0,NBINS);   // error on mean intensity
    hmeanerr = new TH1D("hmeanerr","",NBINS,0,NBINS);   // error on mean hit time
    hwdtherr = new TH1D("hwdtherr","",NBINS,0,NBINS);   // error on mean signal width
    int tb, tn[NBINS] = {0};
    double tavgx[NBINS] = {0};
    double tavgy[NBINS] = {0};
    double tavgz[NBINS] = {0};
    double terrx[NBINS] = {0};
    double terry[NBINS] = {0};
    double terrz[NBINS] = {0};
    // Loop over all PMTs in the graph and bin the data
    for (int i=0; i<NPMTS; i++) {
      if (tx[i]+ty[i]+tz[i]==0) continue;     // not good PMT
      tb = (int)angpmts[i];                   // bin index (1 degree bins)
      tn[tb]++;                               // number of PMTs in each bin
      tavgx[tb] += tx[i]/(tex[i]*tex[i]);     // sum of weighted intensities
      terrx[tb] += 1./(tex[i]*tex[i]);        // sum of errors squared
      tavgy[tb] += ty[i]/(tey[i]*tey[i]);     // sum of weighted hit times
      terry[tb] += 1./(tey[i]*tey[i]);        // sum of errors squared
      tavgz[tb] += tz[i]/(tez[i]*tez[i]);     // sum of weighted signal widths
      terrz[tb] += 1./(tez[i]*tez[i]);        // sum of errors squared
    }
    // Loop over each bin to get the average and propagated error
    for (int j=0; j<NBINS; j++) {
      if (tn[j]==0) continue;                 // no PMTs in this bin
      // Mean intensity
      tavgx[j] = tavgx[j]/terrx[j];           // weighted mean intensity
      terrx[j] = sqrt(1./terrx[j]);           // error on weighted mean
      hpeak->SetBinContent(j+1,tavgx[j]);
      hpeak->SetBinError(j+1,terrx[j]);
      hpeakerr->SetBinContent(j+1,terrx[j]);
      // Mean hit time
      tavgy[j] = tavgy[j]/terry[j];           // weighted mean hit time
      terry[j] = sqrt(1./terry[j]);           // error on weighted mean
      hmean->SetBinContent(j+1,tavgy[j]);
      hmean->SetBinError(j+1,terry[j]);
      hmeanerr->SetBinContent(j+1,terry[j]);
      // Mean signal width
      tavgz[j] = tavgz[j]/terrz[j];           // weighted mean signal width
      terrz[j] = sqrt(1./terrz[j]);           // error on weighted mean
      hwdth->SetBinContent(j+1,tavgz[j]);
      hwdth->SetBinError(j+1,terrz[j]);
      hwdtherr->SetBinContent(j+1,terrz[j]);
    }
    
    // Central hit time for 2D plot
    //int centraltime = round(fitSyst->GetParameter(0) + fitSyst->GetParameter(1)*12.0);
    if (VERBOSE) printf("Mean hit time: %.2lf ns\n",meanhittime);
    int centraltime = (int)round(meanhittime + 5/2.);
    centraltime -= centraltime % 5;             // rounded to the nearest multiple of 5
    
    // More histograms
    h2 = new TH2D("h2",fibre.c_str(),NBINS,0,MAXANG,30,centraltime-15,centraltime+15);
    herr = new TH1D("herr",fibre.c_str(),NBINS,0,MAXANG);
    
    for (int i=0; i<NPMTS; i++) {
      if (tx[i]+ty[i]+tz[i]==0) continue;     // not good PMT
      hpmtgood->Fill(angpmts[i]);
      // Only fill histogram if pmttime is within prompt peak (90% area ~ 1.645 sigma)
      for (int j=0; j<=htime->GetNbinsY()+1; j++) {
        float time = htime->GetYaxis()->GetBinCenter(j+1);
        if (fabs(time - ty[i]) > 1.645*tz[i]) continue;     // not in prompt peak
        //printf("!!! PMT #%4d (bin %4d) at time %8.3f (bin %3d) has content %3d\n",i,i+1,time,j+1,(int)htime->GetBinContent(i+1,j+1));
        for (int k=0; k<(int)htime->GetBinContent(i+1,j+1); k++)
          h2->Fill(angpmts[i], time);         // Angle of PMT w.r.t. fitted fibre direction [deg], time [ns]
      }
    }
    
    // Fit angular slices with a gaussian distribution in time
    h2->FitSlicesY();
    h2_0 = (TH1D*)gDirectory->Get("h2_0");  // Constant
    h2_1 = (TH1D*)gDirectory->Get("h2_1");  // Mean
    h2_2 = (TH1D*)gDirectory->Get("h2_2");  // StdDev
    h2_3 = (TH1D*)gDirectory->Get("h2_chi2");  // chi^2/ndof  
    h2_py = NULL;
    
    // Put mean and RMS into 1D histogram
    for (int b=0; b<=NBINS+1; b++) {
      h2->ProjectionY("_py",b,b);
      h2_py = (TH1D*)gDirectory->Get("h2_py");
      h2_0->SetBinContent(b, h2_py->GetSum());  // overwrite Gaussian constant with intensity
      h2_0->SetBinError(b, sqrt(h2_py->GetSum()));  // overwrite Gaussian constant error with sqrt(intensity)
      h2_1->SetBinContent(b, h2_py->GetMean()); // overwrite Gaussian mean with mean
      h2_1->SetBinError(b, h2_py->GetMeanError()); // overwrite Gaussian mean error with mean error
      h2_2->SetBinContent(b, h2_py->GetRMS());  // overwrite Gaussian sigma with RMS
      h2_2->SetBinError(b, h2_py->GetRMSError());  // overwrite Gaussian sigma error with RMS error
      herr->SetBinContent(b, h2_py->GetMeanError()); // not *quite* the same as error on Gaussian mean
    }
    
    // Normalise angular slices to number of PMTs in slice
    double tmp, goodpmts;
    for (int binx=0; binx<=h2->GetNbinsX()+1; binx++) {     // include underflow/overflow
      goodpmts = hpmtgood->GetBinContent(binx);
      tmp = h2_0->GetBinContent(binx);
      //printf("Angular bin #%2d has %3d good PMTs with intensity %3d.\n",binx,(int)goodpmts,(int)tmp);
      if (goodpmts<=0) { h2_0->SetBinContent(binx,0); h2_0->SetBinError(binx,0); }
      else { h2_0->SetBinContent(binx,tmp/goodpmts); h2_0->SetBinError(binx,sqrt(tmp)/goodpmts); }
      //if (VERBOSE) printf("Bin #2d contains %3d good PMTs.\n",binx,(int)goodpmts);
      for (int biny=0; biny<=h2->GetNbinsY()+1; biny++) {   // include underflow/overflow
        tmp = h2->GetBinContent(binx,biny);
        if (goodpmts<=0) h2->SetBinContent(binx,biny,0);
        else h2->SetBinContent(binx,biny,tmp/goodpmts);
      }
    }
    
    // Write all objects to file (histograms included automatically, graphs not)
    pmtfits->Write("pmtfits");
    gpmts->Write("gpmts");
    outfile->Write();
  }

  // Parametrise angular systematic: y = a + b/cos(x)
  TCanvas *c = new TCanvas("c","",800,600);
  TF1 *fitSyst = new TF1("fitSyst", "[0] + [1]/cos(x/180.*pi)", 0, 24);
  fitSyst->SetParameter(0, gpmts->GetMean(2));  // overall mean hit time
  fitSyst->SetParameter(1, 0);                  // assume flat line a priori
  gpmts->Fit("fitSyst", "R,q");                 // force range, quiet mode
  printf("ANGULAR SYSTEMATIC: y = a + b/cos(x)\n --> a = ( %f +/- %f ) ns\n --> b = ( %f +/- %f ) ns\n", fitSyst->GetParameter(0), fitSyst->GetParError(0), fitSyst->GetParameter(1), fitSyst->GetParError(1));
  c->Close();
  if (c) delete c;
  
  
  // ******************
  //  PLOTTING SECTION
  // ******************
  // Plotting options
  gStyle->SetOptStat(0);
  gStyle->SetPadLeftMargin(0.1);
  gStyle->SetPadRightMargin(0.1);
  gStyle->SetPadTopMargin(0.1);
  gStyle->SetPadBottomMargin(0.1);
  gStyle->SetPadBorderSize(0);
  
  // *****
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
  
  // *****
  // Projected mean hit time (fit results) vs angle
  pad0->cd()->SetGrid();
  pad0->SetLeftMargin(0.12);    // for label
  string tstr = Form("TELLIE angular systematic (run %d); Angle of PMT w.r.t. fitted fibre direction [deg]; Hit time [ns]",run);
  gpmts->SetTitle(tstr.c_str());
  gpmts->SetMarkerColor(4);
  gpmts->SetMarkerStyle(6);
  gpmts->Draw("AP");
  gpmts->GetXaxis()->SetTitleOffset(1.2);
  gpmts->GetYaxis()->SetTitleOffset(1.5);
  gpmts->GetXaxis()->SetRangeUser(0,24);
  
  // *****
  // Normalised intensity profile, hit time vs angle
  pad1->cd()->SetGrid();
  pad1->SetRightMargin(0.15);   // for TH2D color scale
  h2->SetTitle(Form("Normalised intensity profile (%s);Angle of PMT w.r.t. fitted fibre direction [deg];Time [ns]",fibre.c_str()));
  h2->Draw("colz");
  h2->GetXaxis()->SetTitleOffset(1.2);
  h2->GetYaxis()->SetTitleOffset(1.5);
  
  // *****
  // Mean hit times (binned PMT fits)
  pad2->cd()->SetGrid();
  hmean->SetTitle("Mean PMT hit time;Angle of PMT w.r.t. fitted fibre direction [deg];Mean time [ns]");
  hmean->SetLineWidth(2);
  hmean->Draw();
  hmean->GetXaxis()->SetTitleOffset(1.2);
  hmean->GetYaxis()->SetTitleOffset(1.2);
  hmean->GetYaxis()->SetRangeUser(0,8);
  
  // *****
  // Mean errors (fitted PMT widths added in quadrature)
  pad3->cd()->SetGrid();
  hmeanerr->SetTitle("Error on mean PMT hit time;Angle of PMT w.r.t. fitted fibre direction [deg];Error on mean time [ns]");
  hmeanerr->SetLineWidth(2);
  hmeanerr->Draw();
  hmeanerr->GetXaxis()->SetTitleOffset(1.2);
  hmeanerr->GetYaxis()->SetTitleOffset(1.8);
  hmeanerr->GetYaxis()->SetRangeUser(0,0.16);
  
  // *****
  // Mean signal widths
  pad4->cd()->SetGrid();
  hwdth->SetTitle("Binned PMT signal width;Angle of PMT w.r.t. fitted fibre direction [deg];Mean signal width [ns]");
  hwdth->SetLineWidth(2);
  hwdth->Draw();
  hwdth->GetXaxis()->SetTitleOffset(1.2);
  hwdth->GetYaxis()->SetTitleOffset(1.5);
  hwdth->GetYaxis()->SetRangeUser(3.5,5.5);
  
  // *****
  // Mean intensities and number of PMTs in angular segment
  pad5->cd()->SetGrid();
  hpeak->SetTitle("Mean PMT intensity;Angle of PMT w.r.t. fitted fibre direction [deg];Mean peak intensity [a.u.]");
  hpeak->SetLineWidth(2);
  hpeak->Draw();
  hpeak->GetXaxis()->SetTitleOffset(1.2);
  hpeak->GetYaxis()->SetTitleOffset(1.7);
  hpeak->GetYaxis()->SetRangeUser(0, 900);
  // PMTs in that angular segment
  hpmtseg->SetLineWidth(2);
  hpmtseg->SetLineColor(2);
  hpmtseg->SetTitle("# All PMTs");
  hpmtseg->Draw("same");
  // PMTs with occupancy >1%
  hpmtgood->SetLineWidth(2);
  hpmtgood->SetLineColor(3);
  hpmtgood->SetTitle("# Good PMTs");
  hpmtgood->Draw("same");
  // Legend
  TLegend *leg = pad5->BuildLegend();
  
  // *****
  // Save canvas and close
  string imgfile = "angular";
  if (!TEST) imgfile = Form("images/angular_%s",fibre.c_str());
  c0->Print(Form("%s.png",imgfile.c_str()));
  c0->Print(Form("%s.pdf",imgfile.c_str()));
  c0->Close();
  
  // Close files and free memory
  outfile->Close();
  if (outfile) delete outfile;
  if (leg) delete leg;
  if (c0) delete c0;
  
  // Return averaged nhit/event
  return (float)hitpmts/nevents;
}

