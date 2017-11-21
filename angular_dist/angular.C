// ---------------------------------------------------------
// Goal:            Evaluate angular systematic of TELLIE
// Author:          Martti Nirkko, 21/11/2017
// Compile & run:   clear && g++ -o angular.exe angular.C `root-config --cflags --libs` -I$RATROOT/include/libpq -I$RATROOT/include -L$RATROOT/lib -lRATEvent_Linux && ./angular.exe
// ---------------------------------------------------------

// Helper functions (includes everything else)
#include "../HelperFunc.C"

// Run time parameters
const int RUN_CLUSTER = 1;      // whether running on cluster (0=local)
const int USE_RATDB = 1;        // whether to use RATDB to get fibre positions
const int VERBOSE = 1;          // verbosity flag
const int IS_MC = 0;            // Monte-Carlo flag 
const double LOCALITY = 10.0;   // accepted tolerance [mm] for LightPathCalculator

// Initialise functions
int angular(string, int, TF1*, bool, bool);

// Main program
int main(int argc, char** argv) {
  
  // Test flag (0 = process all runs, else run number)
  const int TEST = (RUN_CLUSTER) ? 0 : 102157;

  // Initialise input file
  string input = "../pca_runs/TELLIE_PCA.txt";
  ifstream in(input.c_str());
  if (!in) { cerr<<"Failed to open "<<input<<endl; exit(1); }
  string line;
  for (int hdr=0; hdr<2; hdr++) {
    getline(in,line);      // skip header
  }
  
  // Initialise fit result
  TF1 *fitResult = new TF1("fitResult", "[0] - [1] + [1]/cos(x/180.*pi)", 0, 24);
  double reset[2] = {-999, -999};
  
  // Initialise output file
  string output = "ANGULAR_FITRESULTS.txt";
  FILE *out = fopen(output.c_str(),"w");
  fprintf(out,"Fibre  Run    Val_a Err_a Val_b  Err_b Chi_sq NDF\n");
  fprintf(out,"-------------------------------------------------\n");
  
  // Loop over all fibres in list
  string fibre;
  int node, channel, run, ipw, photons, pin, rms;
  float nhit;
  int nfiles=0;
  while (true) {
    
    // Check for correct run
    in >> node >> fibre >> channel >> run >> ipw >> photons >> pin >> rms >> nhit;
    if (!in.good()) break;
    if (TEST && TEST!=run) continue; // only want specified run
    
    // Reset parameters
    fitResult->SetParameters(reset);
    fitResult->SetParErrors(reset);
    
    // Fit angular systematic
    int errors = angular(fibre, run, fitResult, IS_MC, TEST);
    if (errors) cerr<<"*** WARNING *** Run "<<run<<" was not processed correctly."<<endl;
    else {
      // Print fit results
      double a = fitResult->GetParameter(0);
      double b = fitResult->GetParameter(1);
      double sa = fitResult->GetParError(0);
      double sb = fitResult->GetParError(1);
      int chi2 = round(fitResult->GetChisquare());
      int ndof = fitResult->GetNDF();
      fprintf(out,"%s %d %.3lf %.3lf %.3lf %.3lf %d %d\n",fibre.c_str(),run,a,sa,b,sb,chi2,ndof);
      if (VERBOSE) printf("ANGULAR SYSTEMATIC: y = a - b + b/cos(x)\n --> a = ( %f +/- %f ) ns\n --> b = ( %f +/- %f ) ns\n --> chi^2 / NDF = %d / %d\n", a, sa, b, sb, chi2, ndof);
      nfiles++;
    }
  }
  
  // Close output file and exit
  fclose(out);
  printf("Ran over %d files.\n",nfiles);
  if (nfiles==0) {
    cerr<<"*** ERROR *** Did not process any files."<<endl;
    return 1; 
  }
  return 0;
}

// Define macro
int angular(string fibre, int run, TF1 *fitResult, bool isMC=false, bool TEST=false) {

  // ********************************************************************
  // Initialisation
  // ********************************************************************
  const int NBINS = 24;
  const int MAXANG = 24;
  
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
    return 1;
  } else if(g.good()) {  // file extracted, but not processed
    printf("not processed! Generating plots for fibre %s.\n",fibre.c_str());
    scanned_file = 1;
  } else if(!f.good()) {          // file not downloaded
    printf("not downloaded! Skipping fibre %s.\n",fibre.c_str());
    return 2;
  } else {                        // file downloaded, but not processed
    printf("OK. Extracting data for fibre %s.\n",fibre.c_str());
  }

  string outstr = Form("output/angular_%s.root",fibre.c_str());
  TFile *rootfile = NULL;
  TH1D *hpeak=NULL, *hpeakerr=NULL;
  TH1D *hmean=NULL, *hmeanerr=NULL;
  TH1D *hwdth=NULL, *hwdtherr=NULL;
  TH1D *hpmtseg=NULL, *hpmtgood=NULL;
  TH2D *hprofile=NULL;
  TGraph2DErrors *pmtfits=NULL;
  TGraphErrors *gpmts=NULL;
  if (scanned_file) {
    rootfile = new TFile(outstr.c_str(),"READ");
    hprofile = (TH2D*)rootfile->Get("hprofile"); // Time vs angle
    pmtfits = (TGraph2DErrors*)rootfile->Get("pmtfits"); // Fitted prompt peaks (3D)
    gpmts = (TGraphErrors*)rootfile->Get("gpmts"); // Fitted prompt peaks (time vs angle)
    hpeak = (TH1D*)rootfile->Get("hpeak");       // Binned intensities
    hpeakerr = (TH1D*)rootfile->Get("hpeakerr"); // Binned errors on intensities
    hmean = (TH1D*)rootfile->Get("hmean");       // Binned mean times
    hmeanerr = (TH1D*)rootfile->Get("hmeanerr"); // Binned errors on mean times
    hwdth = (TH1D*)rootfile->Get("hwdth");       // Binned mean signal widths
    hwdtherr = (TH1D*)rootfile->Get("hwdtherr"); // Binned errors on signal widths
    hpmtseg = (TH1D*)rootfile->Get("hpmtseg");   // PMTs in angular bin
    hpmtgood = (TH1D*)rootfile->Get("hpmtgood"); // Good PMTs in angular bin
    cout << "Read from file " << outstr << endl;
  } else {
    rootfile = new TFile(outstr.c_str(),"RECREATE");

    // Initialise RAT
    RAT::DU::DSReader dsreader(fname);
    //RAT::DU::Utility::Get()->BeginOfRun();
    RAT::DU::LightPathCalculator lpc = RAT::DU::Utility::Get()->GetLightPathCalculator();
    //lpc.BeginOfRun(); // TODO - find out if this is needed too!
    RAT::DU::GroupVelocity gv = RAT::DU::Utility::Get()->GetGroupVelocity();
    const RAT::DU::ChanHWStatus& chs = RAT::DU::Utility::Get()->GetChanHWStatus();
    //gv.BeginOfRun(); // TODO - find out if this is needed too!
    const RAT::DU::PMTInfo& pmtinfo = RAT::DU::Utility::Get()->GetPMTInfo();
    const int NPMTS = pmtinfo.GetCount();
    
    TVector3 fibrePos(0,0,0);
    TVector3 fibreDir(0,0,0);
    TVector3 lightPos(0,0,0);
    if (USE_RATDB) {
      // Get fibre info (from RATDB) 
      RAT::DB *db = RAT::DB::Get();
      //db->LoadDefaults();	  // Already done when calling DU::Utility::Get()
      RAT::DBLinkPtr entry = db->GetLink("FIBRE",fibre);
      fibrePos.SetXYZ(entry->GetD("x"), entry->GetD("y"), entry->GetD("z")); // position of fibre [mm]
      fibreDir.SetXYZ(entry->GetD("u"), entry->GetD("v"), entry->GetD("w")); // direction of fibre
      lightPos = fibrePos + 2*fibrePos.Mag()*fibreDir;                       // projected light spot centre
      if (VERBOSE)  cout << "RATDB: fibre " << fibre << ", pos " << printVector(fibrePos) << ", dir " << printVector(fibreDir) << endl;
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
      fibrePos.SetXYZ(10*fx,10*fy,10*fz);                // position of fibre [mm]
      fibreDir.SetXYZ(fu,fv,fw);                         // direction of fibre
      lightPos = fibrePos + 2*fibrePos.Mag()*fibreDir;   // projected light spot centre 
      if (VERBOSE) cout << "DocDB: fibre " << fibre << ", pos " << printVector(fibrePos) << ", dir " << printVector(fibreDir) << endl;
    }

    // Get fitted light position (from file)
    string fitfile = "../fibre_validation/TELLIE_FITRESULTS.txt";
    ifstream fit(fitfile.c_str());
    if (!fit) { cerr<<"ERROR - Failed to open "<<fitfile<<endl; exit(1); }
    string line, fibrecheck;
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
    TVector3 fitdir = (fitpos-fibrePos).Unit();
    if (fitpos.Mag()==0) {
      cerr << "*** ERROR - Could not load fit position!" << endl;
      exit(1);
    } else {
      cout << "Loaded fit position: " << printVector(fitpos) << " mm, " << fitpos.Angle(lightPos)*180./pi << " deg deviation" << endl;
    }
    
    // Initialise histograms
    hpmtseg = new TH1D("hpmtseg",fibre.c_str(),NBINS,0,MAXANG);
    hpmtgood = new TH1D("hpmtgood",fibre.c_str(),NBINS,0,MAXANG);
    TH2D *htime = new TH2D("htime","",NPMTS,0,NPMTS,500,0,500);
    
    // Loop over all entries
    int nevents=0;
    int hitpmts=0;
    float pmtOccup[NPMTS], pmtAngle[NPMTS];                   // occupancy and angle
    memset( pmtOccup, 0, NPMTS*sizeof(float) );               // NPMTS only known at runtime
    memset( pmtAngle, 0, NPMTS*sizeof(float) );               // NPMTS only known at runtime
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
          
          // Get info for this PMT
          TVector3 pmtPos = pmtinfo.GetPosition(pmtID);       // position [mm]
          TVector3 pmtDir = pmtinfo.GetDirection(pmtID);      // direction
          double pmtTime = pmt.GetTime();                     // hit time [ns]
          
          // Get light travel time (as done in BiPoLikelihoodDiff.cc)
          lpc.CalcByPosition(fibrePos, pmtPos, ENERGY, LOCALITY);
          double distInInnerAV = lpc.GetDistInInnerAV();
          double distInAV = lpc.GetDistInAV();
          double distInWater = lpc.GetDistInWater();
          double lightTravelTime = gv.CalcByDistance(distInInnerAV, distInAV, distInWater, ENERGY);
          
          // Get light bucket time (as done in DQLaserBallProc.cc)
          TVector3 endDir = lpc.GetIncidentVecOnPMT();        // end direction at PMT
          double thetaAtPMT = endDir.Angle(pmtDir)*180./pi;   // incident angle with bucket face
          double lightBucketTime = gv.PMTBucketTime(thetaAtPMT);  // DocDB 3138
          
          // Get residual time after correcting for all offsets
          double emissionTime = pmtTime - lightTravelTime - lightBucketTime;
          //double totalOffset = emissionTime + fibreDelay + triggerDelay; // deprecated
          //double residualTime = emissionTime - pcaOffset;     // residual w.r.t. Mark's analysis
          
          // Get light emission angle
          TVector3 startDir = lpc.GetInitialLightVec();       // start direction at fibre
          double theta = startDir.Angle(fitdir)*180./pi;      // emission angle at fibre
          
          /*
          // DEPRECATED - used this before light path calculator
          TVector3 track = pmtPos-fibrePos;                   // straight line assumption
          double theta_old = track.Angle(fitdir)*180./pi;     // angle w.r.t. fibre [deg]
          double lightTravelTime_old = track.Mag()/C_WATER;   // light travel time [ns]
          if (VERBOSE) printf("Photon at PMT #%04d has differences %6.3lf ns (travel time) and %6.3lf deg (initial angle)\n", pmtID, lightTravelTime-lightTravelTime_old, theta-theta_old);
          */
          
          // Fill histograms/arrays
          if (theta > MAXANG) continue;
          htime->Fill(pmtID, emissionTime);                   // time residual vs PMT ID
          if (pmtOccup[pmtID]==0) {                           // fill angles only once per PMT
            hpmtseg->Fill(theta);
            pmtAngle[pmtID] = theta;
          }
          pmtOccup[pmtID]++;
          
        } // pmt loop
      } // event loop
    } // entry loop

    // Turn array of PMTs into percentage (occupancy)
    cout << "Number of events was " << nevents << endl;
    for (int iPMT=0; iPMT<NPMTS; iPMT++) {
      pmtOccup[iPMT] /= (float)nevents;
    }
    
    // Fit PMT hit times (Gaussian prompt peak)
    pmtfits = new TGraph2DErrors();
    FitPromptPeaks(htime, NPMTS, pmtOccup, pmtAngle, pmtfits);
    if (!pmtfits) { cout << "*** ERROR *** Could not fit PMT prompt peaks!" << endl; return 4; }
    if (!pmtfits->GetN()) { cout << "*** ERROR *** Graph contains zero points!" << endl; return 4; }
    
    double *tx = pmtfits->GetX();   // Constant (A)
    double *ty = pmtfits->GetY();   // Mean value (mu)
    double *tz = pmtfits->GetZ();   // Uncertainty (sigma)
    double *tex = pmtfits->GetEX();
    double *tey = pmtfits->GetEY();
    double *tez = pmtfits->GetEZ();
    
    // Fill graph with time vs angle (use only PMTs of interest)
    gpmts = new TGraphErrors();
    int npmts=0;
    for (int iPMT=0; iPMT<NPMTS; iPMT++) {
      if (tx[iPMT]+ty[iPMT]+tz[iPMT]==0) continue;
      gpmts->SetPoint(npmts, pmtAngle[iPMT], ty[iPMT]);
      gpmts->SetPointError(npmts, 0, tey[iPMT]);
      npmts++;
    }
    
    // Mean and RMS of entire graph
    double meanhittime = gpmts->GetMean(2);
    double rmshittime = gpmts->GetRMS(2);
    
    // Investigate PMTs with unusual offsets w.r.t. mean hit time
    string outlogstr = Form("logs/unusual_timing_%d.log", run);
    FILE *outlog = fopen(outlogstr.c_str(),"w");
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
      tb = (int)pmtAngle[i];                  // bin index (1 degree bins)
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
    
    // Get central hit time and initialise 2D histogram
    if (VERBOSE) printf("Mean hit time: %.2lf ns\n",meanhittime);
    int centraltime = (int)round(meanhittime + 5/2.);
    centraltime -= centraltime % 5;           // rounded to the nearest multiple of 5
    hprofile = new TH2D("hprofile",fibre.c_str(),NBINS,0,MAXANG,30,centraltime-15,centraltime+15);
    
    for (int i=0; i<NPMTS; i++) {
      if (tx[i]+ty[i]+tz[i]==0) continue;     // not good PMT
      hpmtgood->Fill(pmtAngle[i]);
      // Only fill histogram if pmttime is within prompt peak (90% area ~ 1.645 sigma)
      for (int j=0; j<=htime->GetNbinsY()+1; j++) {
        float time = htime->GetYaxis()->GetBinCenter(j);
        if (fabs(time - ty[i]) > 1.645*tz[i]) continue; // not in prompt peak
        for (int k=0; k<(int)htime->GetBinContent(i+1,j+1); k++)
          hprofile->Fill(pmtAngle[i], time);  // Angle of PMT w.r.t. fitted fibre direction [deg], time [ns]
      }
    }
    
    // Normalise angular slices to number of PMTs in slice
    double tmp, goodpmts;
    for (int binx=0; binx<=hprofile->GetNbinsX()+1; binx++) {     // include underflow/overflow
      goodpmts = hpmtgood->GetBinContent(binx);
      for (int biny=0; biny<=hprofile->GetNbinsY()+1; biny++) {   // include underflow/overflow
        tmp = hprofile->GetBinContent(binx,biny);
        if (goodpmts<=0) hprofile->SetBinContent(binx,biny,0);
        else hprofile->SetBinContent(binx,biny,tmp/goodpmts);
      }
    }
    
    // Write all objects to file (histograms included automatically, graphs not)
    pmtfits->Write("pmtfits");
    gpmts->Write("gpmts");
    rootfile->Write();
    cout << "Wrote output to file " << outstr << endl;
  }

  // TELLIE fibre/trigger delays & Mark's PCA offsets
  string delayfile = "TELLIE_delays.txt";
  ifstream del(delayfile.c_str());
  if (!del) { cerr<<"ERROR - Failed to open "<<delayfile<<endl; exit(1); }
  int num, triggerDelay;
  float fibreDelay, pcaOffset;
  string line;
  for (int hdr=0; hdr<2; hdr++) {
    getline(del,line);      // header
  }
  while (true) {
    del >> num >> fibreDelay >> triggerDelay >> pcaOffset;
    if (!del.good()) break;
    if (num == run) break;
  }
  if (VERBOSE) cout<<"TELLIE: trigger delay "<<triggerDelay<<" ns, fibre delay "<<fibreDelay<<" ns, PCA offset "<<pcaOffset<<" ns."<<endl;
    
  // Apply global offset to histograms
  double x,t;
  for (int i=0; i<gpmts->GetN(); i++) {
    gpmts->GetPoint(i,x,t);
    gpmts->SetPoint(i,x,t-pcaOffset);
  }
  double meanhittime2 = gpmts->GetMean(2);
  int centraltime2 = (int)round(meanhittime2 + 5/2.);
  centraltime2 -= centraltime2 % 5;           // rounded to the nearest multiple of 5
  TH2D *hprofile2 = new TH2D("hprofile2",fibre.c_str(),NBINS,0,MAXANG,30,centraltime2-15,centraltime2+15);
  for (int i=0; i<=hprofile->GetNbinsX()+1; i++) {
    for (int j=0; j<=hprofile->GetNbinsY()+1; j++) {
      hprofile2->SetBinContent(i,j,hprofile->GetBinContent(i,j));
    }
  } 
  for (int k=0; k<=hmean->GetNbinsX()+1; k++) {
    t = hmean->GetBinContent(k);
    hmean->SetBinContent(k,t-pcaOffset);
  }
  
  // Fit angular systematic: y = a - b + b/cos(x)
  TCanvas *c1 = new TCanvas("c1","",800,600);
  fitResult->SetParameter(0, gpmts->GetMean(2));  // value at zero angle
  fitResult->SetParameter(1, 0);                  // assume flat line a priori
  gpmts->Fit("fitResult","R,q");                  // force range, quiet mode
  c1->Close(); delete c1;
  
  // Fill histograms with time residual
  TH1D *hresid = new TH1D("hresid","",30,-3,3);
  //TH1D *hpulls = new TH1D("hpulls","",40,-20,20);
  double y,ey,y0;
  for (int i=0; i<gpmts->GetN(); i++) {
    gpmts->GetPoint(i,x,y);
    //ey = gpmts->GetErrorY(i);
    y0 = fitResult->Eval(x);
    hresid->Fill(y-y0);
    //hpulls->Fill((y-y0)/ey);
  }
  c1 = new TCanvas("c1","",800,600);
  TF1 *fitresid = new TF1("fitresid","gaus",-3,3);
  fitresid->SetParameters(hresid->GetMaximum(),hresid->GetMean(),hresid->GetRMS());
  hresid->Fit("fitresid","R,q");
  double sigma_resid = fitresid->GetParameter(2);
  c1->Close(); delete c1;
  
  // Find out if direct light from fibre is affected by belly plates
  ifstream belly("../fibre_validation/BELLY_FIBRES.list");
  bool IS_BELLY_FIBRE = false;
  string thisfibre;
  while (belly.good()) {
    getline(belly,thisfibre);
    if(!belly.good()) break;
    if(thisfibre==fibre) IS_BELLY_FIBRE = true;
  }
  
  // Retrieve rotated detector view from fibre validation document
  const int NCOL=20;
  TFile *extfile=NULL;
  TH2D *hfineD=NULL;
  //TGraph2D *gDir2D=NULL;
  TGraph *gDir[NCOL+2]={NULL};
  TGraph *pFibDir=NULL, *pWgtDir=NULL, *pFitDir=NULL;
  TVector3 *fitDir=NULL, *maxima=NULL;
  float nearmax=-1;
  fpath = Form("%s/Software/SNOP/work/analysis/fibre_validation/output",getenv("HOME"));
  fname = Form("%s/PCA_%s.root",fpath.c_str(),fibre.c_str());
  ifstream ext(fname.c_str());
  if (ext.good()) {
    extfile = new TFile(fname.c_str(),"READ");
    hfineD = (TH2D*)extfile->Get("hfineD");
    //gDir2D = (TGraph2D*)extfile->Get("gDir2D");
    for(int s=0; s<NCOL+2; s++) {
      string ss = Form("[%d]",s);
      gDir[s] = (TGraph*)extfile->Get(("gDir"+ss).c_str());
    }
    pFibDir = (TGraph*)extfile->Get("pFibDir");
    pWgtDir = (TGraph*)extfile->Get("pWgtDir");
    pFitDir = (TGraph*)extfile->Get("pFitDir");
    fitDir = (TVector3*)extfile->Get("fitDir");
    maxima = (TVector3*)extfile->Get("maxima");
    nearmax = maxima->X();
  } else {
    cout << "*** WARNING: File containing 2D graphs not found!" << endl;
  }
  
  // ******************
  //  PLOTTING SECTION
  // ******************
  int minvalY = (int)round(fitResult->GetParameter(0)); // lower limit for times to display
  minvalY -= minvalY % 5;                               // round down to next multiple of 5
  
  // Plotting options
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(1);
  gStyle->SetStatX(0.6);
  gStyle->SetStatY(0.88);
  gStyle->SetStatW(0.25);
  gStyle->SetStatH(0.18);
  gStyle->SetPadLeftMargin(0.1);
  gStyle->SetPadRightMargin(0.1);
  gStyle->SetPadTopMargin(0.1);
  gStyle->SetPadBottomMargin(0.1);
  gStyle->SetPadBorderSize(0);
  
  // *****
  // Define canvas and pads
  TCanvas *c0 = new TCanvas("","",1200,1400);
  TPad *pad0 = new TPad("pad0","",0.01,0.58,0.48,0.99);
  TPad *pad1 = new TPad("pad1","",0.52,0.58,1.00,0.99);
  TPad *pad2 = new TPad("pad2","",0.01,0.295,0.33,0.57);
  TPad *pad3 = new TPad("pad3","",0.34,0.295,0.66,0.57);
  TPad *pad4 = new TPad("pad4","",0.67,0.295,0.99,0.57);
  TPad *pad5 = new TPad("pad5","",0.01,0.01,0.33,0.285);
  TPad *pad6 = new TPad("pad6","",0.34,0.01,0.66,0.285);
  TPad *pad7 = new TPad("pad7","",0.67,0.01,0.99,0.285);
  pad0->Draw();
  pad1->Draw();
  pad2->Draw();
  pad3->Draw();
  pad4->Draw();
  pad5->Draw();
  pad6->Draw();
  pad7->Draw();
  
  // *****
  // Text boxes to display fit results
  TBox *tbox = new TBox(0.6,minvalY+4.45,15.5,minvalY+6.85);
  tbox->SetLineColor(2);
  tbox->SetFillColor(kYellow-9);
  TLatex *tfit[4] = {NULL};
  tfit[0] = new TLatex(1, minvalY+6.7, "Fit function: y = a #plus b#left(#frac{1}{cos(x)} #minus 1#right)");
  tfit[1] = new TLatex(1, minvalY+5.9, Form(" #Rightarrow a = ( %.3lf #pm %.3lf ) ns",
                                            fitResult->GetParameter(0), fitResult->GetParError(0)));
  tfit[2] = new TLatex(1, minvalY+5.4, Form(" #Rightarrow b = ( %.3lf #pm %.3lf ) ns",
                                            fitResult->GetParameter(1), fitResult->GetParError(1)));
  tfit[3] = new TLatex(1, minvalY+4.95, Form(" #Rightarrow #chi^{2}/ndf = %d / %d",
                                             (int)round(fitResult->GetChisquare()), fitResult->GetNDF()));
  for (int l=0; l<4; l++) {
    tfit[l]->SetTextAlign(13);
    tfit[l]->SetTextFont(62);
    tfit[l]->SetTextSize(0.03);
  }
  
  // *****
  // PMT hit time offsets and parametrised angular systematic
  pad0->cd()->SetGrid();
  pad0->SetLeftMargin(0.12);    // for label
  string tstr = Form("TELLIE angular systematic (run %d);Angle of PMT w.r.t. fitted fibre direction [deg];Offset in PMT hit time [ns]",run);
  gpmts->SetTitle(tstr.c_str());
  gpmts->SetMarkerColor(4);
  gpmts->SetMarkerStyle(6);
  gpmts->Draw("AP");
  gpmts->GetXaxis()->SetTitleOffset(1.2);
  gpmts->GetYaxis()->SetTitleOffset(1.5);
  gpmts->GetXaxis()->SetLimits(0,24);
  gpmts->GetYaxis()->SetRangeUser(minvalY-3,minvalY+7); // suppresses outliers!
  //tbox->Draw("L same");
  //for (int l=0; l<4; l++) tfit[l]->Draw("same");

  //Create a histogram to hold the confidence intervals
  TH1D *hint = new TH1D("hint", "Fitted function with error band", 10*NBINS, 0, MAXANG);
  TH1D *hint2 = new TH1D("hint2", "Fitted function with error band", 10*NBINS, 0, MAXANG);
  TH1D *hint3 = new TH1D("hint3", "Fitted function with error band", 10*NBINS, 0, MAXANG);
  //(TVirtualFitter::GetFitter())->GetConfidenceIntervals(hint);
  //Now the "hint" histogram has the fitted function values as the
  //bin contents and the confidence intervals as bin errors
  for(int i=0;i<=hint->GetNbinsX()+1;i++) {
    float x = (float)MAXANG/hint->GetNbinsX()*(i-1);
    hint->SetBinContent(i,fitResult->Eval(x));
    hint->SetBinError(i,sigma_resid);
    hint2->SetBinContent(i,fitResult->Eval(x));
    hint2->SetBinError(i,2*sigma_resid);
    hint3->SetBinContent(i,fitResult->Eval(x));
    hint3->SetBinError(i,3*sigma_resid);
  }
  hint->SetStats(kFALSE);
  hint->SetLineColor(0);
  hint->SetFillColorAlpha(2, 0.2);
  hint->Draw("e3 same");
  hint2->SetStats(kFALSE);
  hint2->SetLineColor(0);
  hint2->SetFillColorAlpha(2, 0.2);
  hint2->Draw("e3 same");
  hint3->SetStats(kFALSE);
  hint3->SetLineColor(0);
  hint3->SetFillColorAlpha(2, 0.2);
  hint3->Draw("e3 same");
  gpmts->Draw("P same");
 
  // *****
  // Normalised intensity profile (time vs angle)
  pad1->cd()->SetGrid();
  pad1->SetRightMargin(0.15);   // for TH2D color scale
  hprofile2->SetTitle(Form("Normalised intensity profile (%s);Angle of PMT w.r.t. fitted fibre direction [deg];Offset in PMT hit time [ns]",fibre.c_str()));
  hprofile2->SetStats(0);
  hprofile2->Draw("colz");
  hprofile2->GetXaxis()->SetTitleOffset(1.2);
  hprofile2->GetYaxis()->SetTitleOffset(1.5);
  
  // *****
  // Mean hit times (binned PMTs)
  pad2->cd()->SetGrid();
  hmean->SetTitle("Signal mean (#mu);Angle of PMT w.r.t. fitted fibre direction [deg];Mean PMT signal offset [ns]");
  hmean->SetLineWidth(2);
  hmean->Draw();
  hmean->GetXaxis()->SetTitleOffset(1.2);
  hmean->GetYaxis()->SetTitleOffset(1.2);
  hmean->GetYaxis()->SetRangeUser(minvalY-3,minvalY+7);
  hmeanerr->GetYaxis()->SetRangeUser(0,0.16);
  
  // *****
  // Mean signal widths
  pad3->cd()->SetGrid();
  hwdth->SetTitle("Signal width (#sigma);Angle of PMT w.r.t. fitted fibre direction [deg];Mean PMT signal width [ns]");
  hwdth->SetLineWidth(2);
  hwdth->Draw();
  hwdth->GetXaxis()->SetTitleOffset(1.2);
  hwdth->GetYaxis()->SetTitleOffset(1.5);
  hwdth->GetYaxis()->SetRangeUser(3.8,5.2);
  
  // *****
  // Mean intensities and number of PMTs in angular segment
  pad4->cd()->SetGrid();
  pad4->SetLeftMargin(0.12);
  hpeak->SetTitle("Signal amplitude (A);Angle of PMT w.r.t. fitted fibre direction [deg];Mean PMT signal amplitude [a.u.]");
  hpeak->SetLineWidth(2);
  hpeak->Draw();
  hpeak->GetXaxis()->SetTitleOffset(1.2);
  hpeak->GetYaxis()->SetTitleOffset(1.7);
  hpeak->GetYaxis()->SetRangeUser(0, 900);
  // PMTs in that bin
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
  TLegend *leg = pad4->BuildLegend();
  
  // *****
  // Mean errors (calculated with weighted arithmetic mean method)
  pad5->cd()->SetGrid();
  pad5->SetLeftMargin(0.12);
  hmeanerr->SetTitle("Error on mean (s_{#mu});Angle of PMT w.r.t. fitted fibre direction [deg];Error on mean PMT signal offset [ns]");
  hmeanerr->SetLineWidth(2);
  hmeanerr->Draw();
  hmeanerr->GetXaxis()->SetTitleOffset(1.2);
  hmeanerr->GetYaxis()->SetTitleOffset(1.8);
  
  // *****
  // Time residuals
  pad6->cd()->SetGrid();
  hresid->SetTitle("PMT time residuals;Time residual [ns];#PMTs/bin");
  hresid->SetLineWidth(2);
  hresid->Draw();
  hresid->GetXaxis()->SetTitleOffset(1.2);
  hresid->GetYaxis()->SetTitleOffset(1.5);
  hresid->GetYaxis()->SetRangeUser(0,500);
  
  // *****
  /*
  // Pull variable
  pad7->cd()->SetGrid();
  hpulls->SetTitle("PMT pull variable;Pull [ ];#PMTs/bin");
  hpulls->SetLineWidth(2);
  hpulls->Draw();
  hpulls->GetXaxis()->SetTitleOffset(1.2);
  hpulls->GetYaxis()->SetTitleOffset(1.5);
  hpulls->GetYaxis()->SetRangeUser(0,200);
  */
  // View from direct light spot (fitted)
  pad7->cd()->SetGrid();
  if (ext.good()) {
    /*
    TH3F *hempty = new TH3F("hempty","",10,-10,10,10,-10,10,10,0,1000+1);
    hempty->Draw("");             // empty histogram for plot range
    pad7->SetTheta(90-0.001);      // view from above
    pad7->SetPhi(0+0.001);         // no x-y rotation
    gDir2D->Draw("pcol,same");
    */
    hfineD->SetTitle("Direct light spot;X' [m];Y' [m]");
    hfineD->GetXaxis()->SetTitleOffset(1.3);
    hfineD->GetYaxis()->SetTitleOffset(1.4);
    hfineD->Draw("scat");
    hfineD->SetStats(0);
    for(int s=0;s<NCOL+2;s++) {
      if(!gDir[s]) continue;
      if(gDir[s]->GetN()==0) continue;
      gDir[s]->SetMarkerStyle(6);
      gDir[s]->Draw("P same");
    }
    pFibDir->Draw("P same");
    pWgtDir->Draw("P same");
    TGraph *pFitDir2 = (TGraph*)pFitDir->Clone();
    double fitD_rotX = fitDir->X()/1e3;
    double fitD_rotY = fitDir->Y()/1e3;
    if (IS_BELLY_FIBRE) pFitDir2->SetMarkerStyle(28); // open cross to indicate Gaussian fit is not used
    pFitDir2->DrawGraph(1,&fitD_rotX,&fitD_rotY,"P same");
    // Indicate possible belly plate effect
    TLatex *txtB = new TLatex(-9.25,9.5,"BELLY PLATE");
    txtB->SetTextAlign(13);
    txtB->SetTextFont(102);
    txtB->SetTextSize(0.04);
    if (IS_BELLY_FIBRE) txtB->Draw();
    // Indicate colour scale
    TLatex *txtD = new TLatex(9.5,9.5,Form("Occup. #leq %.2f%%",100.*nearmax));
    txtD->SetTextAlign(33);
    txtD->SetTextFont(82);
    txtD->SetTextSize(0.04);
    txtD->Draw();
  }
  
  // *****
  // Save canvas and close
  string imgfile = "angular";
  if (!TEST) imgfile = Form("images/angular_%s",fibre.c_str());
  c0->Print(Form("%s.png",imgfile.c_str()));
  c0->Print(Form("%s.pdf",imgfile.c_str()));
  c0->Close();
  
  // Close files and free memory
  rootfile->Close();
  if (rootfile) delete rootfile;
  if (leg) delete leg;
  if (c0) delete c0;
  
  // Return no errors
  return 0;
}

