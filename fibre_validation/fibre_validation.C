// ---------------------------------------------------------
// Goal:          Validate fibre installation positions/directions listed in RATDB
// Author:        Martti Nirkko, 14/11/2017
// Compile & run: clear && g++ -g -o fibre_validation.exe fibre_validation.C `root-config --cflags --libs` -I$RATROOT/include/libpq -I$RATROOT/include -L$RATROOT/lib -lRATEvent_Linux && ./fibre_validation.exe
// ---------------------------------------------------------

// Helper functions (includes everything else)
#include "../HelperFunc.C"

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
void fibre_validation(string, int, int, int, int, float, TVector3*, TVector3*, bool, bool);

// Main program
int main(int argc, char** argv) {
  
  // Test flag (0 = process all runs, else run number)
  const int TEST = (RUN_CLUSTER) ? 0 : 102157;
  
  // Fit results for direct and reflected light
  TVector3 *dirfit = new TVector3(0,0,0);
  TVector3 *reffit = new TVector3(0,0,0);
  
  // Run over Monte Carlo file
  if (IS_MC) {
    string fibre = "FT040A";	// FT040A, FT079A
    fibre_validation(fibre, -999, -999, -999, -999, -999., dirfit, reffit, (bool)IS_MC, false);
    cout << "Results for fibre " << fibre << ":" << endl;
    cout << "- Direct light fit (x,y,z) [mm] = " << printVector(*dirfit) << endl;
    cout << "- Reflected light fit (x,y,z) [mm] = " << printVector(*reffit) << endl;
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
  FILE *fitresult = NULL;
  if (!TEST) {    // don't overwrite file unless running over all fibres
    fitresult = fopen("TELLIE_FITRESULTS.txt","w");
    fprintf(fitresult,"#Run Fibre Direct_light(xyz) Reflected_light(xyz)\n");
    fprintf(fitresult,"#------------------------------------------------\n");
  }
  while (true) {
    in >> node >> fibre >> channel >> run >> ipw >> photons >> pin >> rms >> nhit;
    if (!in.good()) break;
    if (TEST && TEST!=run) continue; // only want specified run
    if (VERBOSE) printf("%6s %2d %6d %5d %6d %5d %5d %.2f\n", fibre.c_str(), channel, run, ipw, photons, pin, rms, nhit);
    fibre_validation(fibre, channel, run, ipw, photons, nhit, dirfit, reffit, (bool)IS_MC, (bool)TEST);
    if (dirfit->Mag()==0) continue;
    if (!TEST) fprintf(fitresult, "%6d %6s %.3f %.3f %.3f %.3f %.3f %.3f\n", run, fibre.c_str(), dirfit->X(), dirfit->Y(), dirfit->Z(), reffit->X(), reffit->Y(), reffit->Z());
    cout << "- Direct light fit    (x,y,z) [mm] = " << printVector(*dirfit) << endl;
    cout << "- Reflected light fit (x,y,z) [mm] = " << printVector(*reffit) << endl;
    cout << endl << "* * * * * * * * * * * * * * * * * * * * " << endl << endl;
    nfiles++;
  }
  if (!TEST) fclose(fitresult);
  printf("Ran over %d files.\n",nfiles);
  if (nfiles==0) { 
    cerr<<"*** ERROR *** No input files found, or nothing to do!"<<endl;
    return 1;
  }
  return 0;
}

// Calculates the fitted light and fibre positions for a given fibre/run
void fibre_validation(string fibre, int channel, int run, int ipw, int photons, float nhit, TVector3 *dirfit, TVector3 *reffit, bool isMC=false, bool TEST=false) {

  // ********************************************************************
  // Initialisation
  // ********************************************************************
  
  // Check files for given run
  if(!TEST && VERBOSE) printf("*****\n");
  printf("Checking files for run %d... ", run);
  string fpath = Form("%s/Software/SNOP/work/data",getenv("HOME"));
  string fname, out, root, img;
  string mcopt = "";
  ifstream ff;
  if (isMC) {
    mcopt = "centre"; // centre, Z+10cm (shifted AV position)
    fpath = Form("%s/Software/SNOP/work/analysis/tellie_jobs/markstringer",getenv("HOME"));
    fname = Form("%s/MC_%s_AV_%s_processed.root",fpath.c_str(),fibre.c_str(),mcopt.c_str());
    ff.open(fname.c_str());
    run = 0; // TODO: figure out how to get this from file
    out  = Form("./mc/MC_%s_%s.out",fibre.c_str(),mcopt.c_str());
    root = Form("./mc/MC_%s_%s.root",fibre.c_str(),mcopt.c_str());
    img  = Form("./mc/MC_%s.pdf",fibre.c_str());
  } else {
    for (int pass=3;pass>=0;pass--) {
      fname = Form("%s/Analysis_r0000%d_s000_p00%d.root",fpath.c_str(),run,pass);
      ff.open(fname.c_str());
      if (ff.good()) break;
    }
    out  = Form("./output/PCA_%s.out",fibre.c_str());
    root = Form("./output/PCA_%s.root",fibre.c_str());
    img  = Form("./images/PCA_%s.pdf",fibre.c_str());
  }
  ifstream gg(out.c_str());
  ifstream hh(root.c_str());
  ifstream ii(img.c_str());
  int scanned_file = 0;
  if (!TEST && ii.good()) {  // file downloaded and processed
    printf("already processed! Skipping fibre %s.\n",fibre.c_str());
    return;
  } else if (hh.good()) {    // graphs created, but not plotted
    printf("already processed! Plotting graphs for fibre %s.\n",fibre.c_str());
    scanned_file = 2;
  } else if (gg.good()) {    // file extracted, but not processed
    printf("not processed! Generating graphs for fibre %s.\n",fibre.c_str());
    scanned_file = 1;
  } else if (!ff.good()) {   // file not downloaded
    printf("not downloaded! Skipping fibre %s.\n",fibre.c_str());
    return;
  } else {                   // file downloaded, but not processed
    printf("OK. Extracting data for fibre %s.\n",fibre.c_str());
  }

  // Find out if direct light from fibre is affected by belly plates
  ifstream belly("BELLY_FIBRES.list");
  bool IS_BELLY_FIBRE = false;
  string thisfibre;
  while (belly.good()) {
    getline(belly,thisfibre);
    if(!belly.good()) break;
    if(thisfibre==fibre) IS_BELLY_FIBRE = true;
  }
  
  // If file was scanned at least once, read number of PMT hits from text file
  int totalnhit=0, directnhit=0, reflectednhit=0, count=0, NPMTS=0;
  int *pmtlightcone=NULL; // is PMT in direct/reflected light cone (1/2), or not (0)
  float *occupancy=NULL; // occupancy
  if (scanned_file>0) {
    int checkrun;
    string dummy;
    gg >> dummy >> checkrun;
    if(!gg.good()) printf("*** ERROR *** Bad input - str=%s int=%d\n",dummy.c_str(),checkrun);
    if(checkrun != run && !isMC) { printf("*** ERROR *** Bad run number %d\n",checkrun); return; }
    gg >> dummy >> count;
    if(!gg.good()) { printf("*** ERROR *** Bad input - str=%s int=%d\n",dummy.c_str(),count); return; }
    gg >> dummy >> totalnhit;
    if(!gg.good()) { printf("*** ERROR *** Bad input - str=%s int=%d\n",dummy.c_str(),totalnhit); return; }
    gg >> dummy >> NPMTS;
    if(!gg.good()) { printf("*** ERROR *** Bad input - str=%s int=%d\n",dummy.c_str(),NPMTS); return; }
    printf("*** INFO *** Run %d has %d EXTA events with %d total NHits on %d PMTs.\n",checkrun,count,totalnhit,NPMTS);
    occupancy = new float[NPMTS];
    pmtlightcone = new int[NPMTS];
    memset( occupancy, 0., NPMTS*sizeof(float) ); // NPMTS only known at runtime
    memset( pmtlightcone, 0, NPMTS*sizeof(int) ); // NPMTS only known at runtime
    int pmtid, pmthits, onspot;
    while (gg.good()) {
      gg >> pmtid >> pmthits >> onspot;
      occupancy[pmtid]=pmthits;
      pmtlightcone[pmtid]=onspot;
      if(onspot==1) directnhit += pmthits;
      else if(onspot==2) reflectednhit += pmthits;
    }
  }
  
  // Initialise graphs and histograms
  TFile *rootfile = NULL;
  TH2D *hicos=NULL, *hcoarse=NULL;
  TH2D *hfineD=NULL, *hfineR=NULL;
  TH1F *hoccup=NULL, *hoccupoff=NULL, *hoccuplo=NULL, *hoccuphi=NULL;
  TGraph2D *gDir2D=NULL, *gRef2D=NULL;
  TGraph *icos[NCOL+2]={NULL}, *gDir[NCOL+2]={NULL}, *gRef[NCOL+2]={NULL};  // array of graphs (for 2D view)
  TGraph *pcontD=NULL, *pcontR=NULL;
  TGraph *pFibPos=NULL, *pFibDir=NULL, *pWgtDir=NULL, *pWgtRef=NULL;
  TGraph *pFitDir=NULL, *pFitRef=NULL;
  TGraph *pcircD1=NULL, *pcircD2=NULL, *pcircR1=NULL, *pcircR2=NULL;
  float nearmax=-1, farmax=-1, HOTLIMIT=-1, COLDLIMIT=-1;
  double DIRANG=-1, REFANG=-1;
  TVector2 *angles=NULL, *maxima=NULL, *limits=NULL;
  TVector3 *fitDir=NULL, *fitRef=NULL, *fitResD=NULL, *fitResR=NULL;
  
  // If file was scanned twice, read graphs from root file
  if (scanned_file==2) {
    rootfile = new TFile(root.c_str(),"READ");
    angles = (TVector2*)rootfile->Get("angles");
    DIRANG = angles->X();
    REFANG = angles->Y();
    maxima = (TVector2*)rootfile->Get("maxima");
    nearmax = maxima->X();
    farmax = maxima->Y();
    limits = (TVector2*)rootfile->Get("limits");
    HOTLIMIT = limits->X();
    COLDLIMIT = limits->Y();
    fitDir = (TVector3*)rootfile->Get("fitDir");  // rotated view
    fitRef = (TVector3*)rootfile->Get("fitRef");  // rotated view
    fitResD = (TVector3*)rootfile->Get("dirfit");  // detector coordinates
    fitResR = (TVector3*)rootfile->Get("reffit");  // detector coordinates
    hicos = (TH2D*)rootfile->Get("hicos");
    hfineD = (TH2D*)rootfile->Get("hfineD");
    hfineR = (TH2D*)rootfile->Get("hfineR");
    hoccup = (TH1F*)rootfile->Get("hoccup");
    hoccupoff = (TH1F*)rootfile->Get("hoccupoff");
    hoccuplo = (TH1F*)rootfile->Get("hoccuplo");
    hoccuphi = (TH1F*)rootfile->Get("hoccuphi");
    gDir2D = (TGraph2D*)rootfile->Get("gDir2D");
    gRef2D = (TGraph2D*)rootfile->Get("gRef2D");
    for(int s=0; s<NCOL+2; s++) {
      string ss = Form("[%d]",s);
      icos[s] = (TGraph*)rootfile->Get(("icos"+ss).c_str());
      gDir[s] = (TGraph*)rootfile->Get(("gDir"+ss).c_str());
      gRef[s] = (TGraph*)rootfile->Get(("gRef"+ss).c_str());
    }
    pcontD = (TGraph*)rootfile->Get("pcontD");
    pcontR = (TGraph*)rootfile->Get("pcontR");
    pFibPos = (TGraph*)rootfile->Get("pFibPos");
    pFibDir = (TGraph*)rootfile->Get("pFibDir");
    pWgtDir = (TGraph*)rootfile->Get("pWgtDir");
    pWgtRef = (TGraph*)rootfile->Get("pWgtRef");
    pFitDir = (TGraph*)rootfile->Get("pFitDir");
    pFitRef = (TGraph*)rootfile->Get("pFitRef");
    pcircD1 = (TGraph*)rootfile->Get("pcircD1");
    pcircD2 = (TGraph*)rootfile->Get("pcircD2");
    pcircR1 = (TGraph*)rootfile->Get("pcircR1");
    pcircR2 = (TGraph*)rootfile->Get("pcircR2");
    
    // Important to set output objects using these functions (ROOT is weird)
    dirfit->SetXYZ(fitResD->X(), fitResD->Y(), fitResD->Z());
    reffit->SetXYZ(fitResR->X(), fitResR->Y(), fitResR->Z());
    cout << "Read from file " << root << endl;
    
  } else {  // need to do the dirty work
  
    rootfile = new TFile(root.c_str(),"RECREATE");
  
    // Initialise RAT
    RAT::DU::DSReader dsreader(fname);
    const RAT::DU::ChanHWStatus& chs = RAT::DU::Utility::Get()->GetChanHWStatus();
    const RAT::DU::PMTInfo& pmtinfo = RAT::DU::Utility::Get()->GetPMTInfo();
    NPMTS = pmtinfo.GetCount();
    cout << "Initialised RAT. Number of PMTs is " << NPMTS << endl;
    
    TVector3 fibrepos, fibredir, lightpos;
    if (USE_RATDB) {
      // Get fibre info (from RATDB) 
      RAT::DB *db = RAT::DB::Get();
      //db->LoadDefaults();	  // Already done when calling DU::Utility::Get()
      RAT::DBLinkPtr entry = db->GetLink("FIBRE",fibre);
      fibrepos.SetXYZ(entry->GetD("x"), entry->GetD("y"), entry->GetD("z")); // position
      fibredir.SetXYZ(entry->GetD("u"), entry->GetD("v"), entry->GetD("w")); // direction
      lightpos = fibrepos + 2*fibrepos.Mag()*fibredir; // projected light spot centre
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
      fibrepos.SetXYZ(10*fx,10*fy,10*fz); // position of fibre [mm]
      fibredir.SetXYZ(fu,fv,fw); // direction of fibre
      lightpos = fibrepos + 2*fibrepos.Mag()*fibredir; // projected light spot centre 
      cout << "DocDB: fibre " << fibre << ", pos " << printVector(fibrepos) << ", dir " << printVector(fibredir) << endl;
    }
    
    if (scanned_file==0) {
    
      // ********************************************************************
      // Sum PMT hit counts for entire run
      // ********************************************************************
      occupancy = new float[NPMTS];
      pmtlightcone = new int[NPMTS];
      memset( occupancy, 0., NPMTS*sizeof(float) ); // NPMTS only known at runtime
      memset( pmtlightcone, 0, NPMTS*sizeof(int) ); // NPMTS only known at runtime
      // Loop over all entries in file
      cout << "Looping over entries..." << endl;
      for(int iEntry=0; iEntry<dsreader.GetEntryCount(); iEntry++) {
        const RAT::DS::Entry& ds = dsreader.GetEntry(iEntry);

        // Print progress
        if (iEntry % int(dsreader.GetEntryCount()/100.) == 0) {
          printProgress(iEntry, dsreader.GetEntryCount());
        }      
        
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
            occupancy[pmtID]++;
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
      fprintf(outFile, "NPMTS: %d\n", NPMTS);         // number of PMTs
      for(int id=0; id<NPMTS; id++) {
        pmtpos = pmtinfo.GetPosition(id);
        int in_spot = 0;
        if (pmtpos.Angle(lightpos) < DIR_CONE/180.*pi) in_spot = 1;
        else if (pmtpos.Angle(fibrepos) < REF_CONE/180.*pi) in_spot = 2;
        fprintf(outFile, "%d %d %d\n", id, (int)occupancy[id], in_spot);
      }
      fclose(outFile);
      
      scanned_file = 1;
    
    } // PMT hit counts extracted
    
    // Convert hit counts to occupancy
    for(int id=0; id<NPMTS; id++) {
      occupancy[id] /= count;
    }

    // Initialise histograms
    hicos = new TH2D("hicos","PMT positions",200,0,1,200,0,1); // icosahedral
    hcoarse = new TH2D("hcoarse","PMT positions",20,-1,1.001,10,-1.001,0); // units of pi
    hfineD = new TH2D("hfineD","PMT positions",1000,-10,10,1000,-10,10); // fine grained
    hfineR = new TH2D("hfineR","PMT positions",1000,-10,10,1000,-10,10); // fine grained
    hoccup = new TH1F("hoccup","PMT hit count",60,0,1); // occupancy
    BinLog(hoccup->GetXaxis(),1e-6); // minimum for log axis
    
    // Find threshold for "screamers"
    GetLimits(occupancy, NPMTS, HOTLIMIT, COLDLIMIT);
    for(int id=0; id<NPMTS; id++) {
      if(occupancy[id]==0) hoccup->Fill(2e-6); // visible on log scale
      else hoccup->Fill(occupancy[id]);
    }
   
    // ********************************************************************
    // Get icosahedral projection of detector
    // ********************************************************************
    // Use icosahedral projection functions from HelperFunc.C
    using namespace func;
    int icosN[NCOL+2]={0};  // off PMTs (0), N colours (1-NCOL), hot PMTs (NCOL+1)
    int pmtface[NPMTS];
    double icosX[NCOL+2][NPMTS], icosY[NCOL+2][NPMTS];
    memset( icosX, 0, (NCOL+2)*NPMTS*sizeof(double) ); // NPMTS only known at runtime
    memset( icosY, 0, (NCOL+2)*NPMTS*sizeof(double) ); // NPMTS only known at runtime
    memset( pmtface, 0, NPMTS*sizeof(int) ); // NPMTS only known at runtime
    
    TVector2 icospos;
    TVector3 pmtpos;
    int nfacepmts[20] = {0}, nfacegoodpmts[20] = {0};
    TVector3 facecenter[20], faceweight[20];
    cout << "Creating icosahedral projection..." << endl;
    for(int id=0; id<NPMTS; id++) {
      pmtpos = pmtinfo.GetPosition(id);
      if (pmtpos.Mag()==0) continue; // not a valid PMT
      if (pmtinfo.GetType(id) != 1) continue; // not a normal PMT
      
      // Obtain PSUP face
      int face; // side of icosahedron containing PMT
      icospos = func::IcosProject(pmtpos,face);
      pmtface[id] = face; // face of PMT (1-20)
      facecenter[face-1] += pmtpos;
      nfacepmts[face-1]++;
      //for(int j=0; j<occupancy[id]; j++) hicos->Fill(icospos.X(),icospos.Y());
      
      // Get correct bin for colour scale
      int step;
      if (occupancy[id] == 0) step=0;                  // off PMT
      else if (occupancy[id] < COLDLIMIT)  step=1;     // cold PMT
      else if (occupancy[id] > HOTLIMIT) step=NCOL+1;  // hot PMT
      // linear colour scale
      //else step = (int)TMath::Ceil(occupancy[id]/(1.*HOTLIMIT/NCOL))+1;
      // logarithmic colour scale
      else step = (int)TMath::Ceil((log10(occupancy[id])-log10(COLDLIMIT)) / ((log10(HOTLIMIT)-log10(COLDLIMIT))/(NCOL-1))) + 1;
      //if (step>1) printf("PMT #%d has occupancy %6.2f%% ==> step=%d\n",id,100.*occupancy[id],step);
      
      // Put PMT in correct graph
      icosX[step][icosN[step]] = icospos.X();
      icosY[step][icosN[step]] = icospos.Y();
      icosN[step]++;
      if (step<2) continue; // cold or off PMT
      if (step>NCOL) continue; // hot PMT
      
      // Calculate PSUP face heat
      faceweight[face-1] += pmtpos*occupancy[id];
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
    for (int s=0; s<NCOL+2; s++) {
      int col;
      if      (s==0) col = 16; // off PMTs (grey)
      else if (s==1) col = 51; // cold PMTs (purple)
      else if (s==21) col = 1; // hot PMTs (black)
      else col = (int)(50.+s*(50./NCOL)); // scale from 55-100 (s=2-20)
      if (icosN[s]==0) {icos[s]=NULL; continue;}
      icos[s] = new TGraph(icosN[s],icosX[s],icosY[s]);
      icos[s]->SetMarkerStyle(7);
      icos[s]->SetMarkerColor(col);
    }
    
    
    // ********************************************************************
    // First iteration: Fill coarse histogram with PMT hit positions
    // ********************************************************************
    float hitsum=0;
    double pmtrad = 0;
    TVector3 pmtsum(0,0,0);
    cout << "Filling coarse histogram..." << endl;
    for(int id=0; id<NPMTS; id++) {
      pmtpos = pmtinfo.GetPosition(id);
      if (pmtpos.Mag()==0) continue;            // not a valid PMT
      if (pmtinfo.GetType(id) != 1) continue;   // not a normal PMT
      if (occupancy[id] == 0) continue;         // off PMTs
      if (occupancy[id] > HOTLIMIT) continue;   // hot PMTs
      if (occupancy[id] < COLDLIMIT) continue;  // cold PMTs
      for (int j=0; j<occupancy[id]; j++) {
        hcoarse->Fill(pmtpos.Phi()/pi, -pmtpos.Theta()/pi);  // negative theta (neck on top)
      }
      pmtrad += pmtpos.Mag()*occupancy[id];
      pmtsum += pmtpos*occupancy[id];
      hitsum += occupancy[id];
    }
    pmtrad /= hitsum;
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
    float weight, empty1, empty2;
    GetMaxColVal(guess_dir, occupancy, NPMTS, nearmax, empty1, pmtinfo);
    GetMaxColVal(guess_ref, occupancy, NPMTS, farmax, empty2, pmtinfo);
    
    if(nearmax==0) cout << "*** WARNING *** No good PMTs in direct light cone!" << endl;
    if(farmax==0) cout << "*** WARNING *** No good PMTs in reflected light cone!" << endl;
    
    cout << "Calculating weighted average..." << endl;
    for(int id=0; id<NPMTS; id++) {
      //if (id % int(NPMTS/100.) == 0) printProgress(id, NPMTS);
      pmtpos = pmtinfo.GetPosition(id);
      if (pmtpos.Mag()==0) continue;                  // not a valid PMT position
      if (pmtinfo.GetType(id) != 1) continue;         // not a normal PMT
      weight = occupancy[id];
      if (weight == 0) continue;                      // off PMTs
      if (weight > HOTLIMIT) continue;                // hot PMTs
      if (weight < COLDLIMIT) continue;               // cold PMTs
      if (pmtpos.Angle(guess_dir) < pi*DIR_CONE/180.) // within cone of direct light spot
        direct += pmtpos*(1.*weight/hitsum);
      if (pmtpos.Angle(guess_ref) < pi*REF_CONE/180.) // within cone of reflected light spot
        reflected += pmtpos*(1.*weight/hitsum);
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
    gDir2D = new TGraph2D();    // 2D graph (for 3D view and fit)
    gRef2D = new TGraph2D();
    FillHemisphere(direct, occupancy, NPMTS, gDir, gDir2D, NCOL, nearmax, pmtinfo);
    FillHemisphere(reflected, occupancy, NPMTS, gRef, gRef2D, NCOL, farmax, pmtinfo);
    
    // TODO - use Gaussian fit amplitude as maximum for colour scale
    //nearmax = (float)totalnhit/count/8.;
    //farmax = (float)totalnhit/count/40.;
    
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
    cout << "Performing 2D Gaussian fits..." << endl;
    FitLightSpot(gDir2D,pmtrad/1e3,DIR_CONE,parDir); // units: dist [m], ang [deg]
    FitLightSpot(gRef2D,pmtrad/1e3,REF_CONE,parRef);
    if (VERBOSE) {
      cout << "FIT RESULTS:" << endl;
      cout << "- Direct light";
      for (int p=0; p<NPARS; p++) printf(" %.3lf", parDir[p]);
      cout << endl;
      cout << "- Reflected light";
      for (int p=0; p<NPARS; p++) printf(" %.3lf", parRef[p]);
      cout << endl;
    }
    // Translate to point on PSUP sphere
    fitDir = new TVector3(1e3*parDir[1],1e3*parDir[2],sqrt(pmtrad*pmtrad-1e6*parDir[1]*parDir[1]-1e6*parDir[2]*parDir[2]));
    fitRef = new TVector3(1e3*parRef[1],1e3*parRef[2],sqrt(pmtrad*pmtrad-1e6*parRef[1]*parRef[1]-1e6*parRef[2]*parRef[2]));
    // Standard deviation (Gaussian sigma) translated to angle [deg]
    double sigmaAngDir = asin(1e3*parDir[3]/pmtrad)*180./pi;
    double sigmaAngRef = asin(1e3*parRef[3]/pmtrad)*180./pi;
    
    // Output values (Gaussian fit)
    dirfit->SetXYZ(fitDir->X(), fitDir->Y(), fitDir->Z());
    dirfit->RotateY(rot_Z1);
    dirfit->RotateZ(rot_X1);
    reffit->SetXYZ(fitRef->X(), fitRef->Y(), fitRef->Z());
    reffit->RotateY(rot_Z2);
    reffit->RotateZ(rot_X2);
    DIRANG = lightpos.Angle(*dirfit);   // angle between expected and fitted light spot
    REFANG = fibrepos.Angle(*reffit);   // angle between fibre position and fitted light spot
    
    // Overwrite value for direct light if affected by belly plates
    if (IS_BELLY_FIBRE) {
      cout << "Direct light is affected by belly plates. Using weighted average instead of Gaussian fit." << endl;
      DIRANG = lightpos.Angle(direct);         // angle between expected and weighted light spot
    }
    
    
    // ********************************************************************
    // Create histograms and graphs for output
    // ********************************************************************
    
    // Get contours around fitted light spots in (phi,theta)
    TVector3 *dotsD[NDOTS], *dotsR[NDOTS];
    if (IS_BELLY_FIBRE) DrawCircle(direct, sigmaAngDir, dotsD, NDOTS); // weighted method (TODO - use fitted sigma?)
    else DrawCircle(*dirfit, sigmaAngDir, dotsD, NDOTS);               // Gaussian fit
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
    pcontD = new TGraph(NDOTS,pcontDx,pcontDy);
    pcontR = new TGraph(NDOTS,pcontRx,pcontRy);
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
    int weightMarker = 2;
    int fitMarker = 34;
    int trueMarker = 24;
    
    // Create markers for relevant points (icosahedral view)
    int nope;
    double ptX, ptY;
    if (IS_BELLY_FIBRE) icospos = func::IcosProject(direct,nope); // use weighted average
    else icospos = func::IcosProject(*dirfit,nope);               // use Gaussian fit
    ptX = icospos.X();
    ptY = icospos.Y();
    pFitDir = new TGraph(1,&ptX,&ptY);          // direct light (fitted)
    pFitDir->SetMarkerStyle(fitMarker);
    pFitDir->SetMarkerColor(1);
    pFitDir->SetMarkerSize(2);
    
    icospos = func::IcosProject(*reffit,nope);
    ptX = icospos.X();
    ptY = icospos.Y();
    pFitRef = new TGraph(1,&ptX,&ptY);          // reflected light (fitted)
    pFitRef->SetMarkerStyle(fitMarker);
    pFitRef->SetMarkerColor(1);
    pFitRef->SetMarkerSize(1.5);
    
    // Create markers for relevant points (rotated view)
    lightpos.RotateZ(-rot_X1);
    lightpos.RotateY(-rot_Z1);
    ptX = lightpos.X()/1e3;
    ptY = lightpos.Y()/1e3;
    pFibDir = new TGraph(1,&ptX,&ptY);      // expected light position
    pFibDir->SetMarkerStyle(trueMarker);
    pFibDir->SetMarkerColor(1);
    pFibDir->SetMarkerSize(2);
    
    fibrepos.RotateZ(-rot_X2);
    fibrepos.RotateY(-rot_Z2);
    ptX = fibrepos.X()/1e3;
    ptY = fibrepos.Y()/1e3;
    pFibPos = new TGraph(1,&ptX,&ptY);      // expected fibre position
    pFibPos->SetMarkerStyle(trueMarker);
    pFibPos->SetMarkerColor(1);
    pFibPos->SetMarkerSize(1.5);
    
    // View is centered over weighted light spot
    const int nil = 0;
    pWgtDir = new TGraph(1,&nil,&nil);         // direct light (weighted)
    if (IS_BELLY_FIBRE) pWgtDir->SetMarkerStyle(fitMarker);
    else pWgtDir->SetMarkerStyle(weightMarker);
    pWgtDir->SetMarkerColor(1);
    pWgtDir->SetMarkerSize(2);
    
    pWgtRef = new TGraph(1,&nil,&nil);         // reflected light (weighted)
    pWgtRef->SetMarkerStyle(weightMarker);
    pWgtRef->SetMarkerColor(1);
    pWgtRef->SetMarkerSize(1.5);

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
    if(nfrontD>0) { pcircD1 = new TGraph(nfrontD,pcircDx[0],pcircDy[0]);
                    pcircD1->SetMarkerStyle(6); }
    if(nbackD>0)  { pcircD2 = new TGraph(nbackD,pcircDx[1],pcircDy[1]);
                    pcircD2->SetMarkerStyle(1); }
    if(nfrontR>0) { pcircR1 = new TGraph(nfrontR,pcircRx[0],pcircRy[0]);
                    pcircR1->SetMarkerStyle(6); }
    if(nbackR>0)  { pcircR2 = new TGraph(nbackR,pcircRx[1],pcircRy[1]);
                    pcircR2->SetMarkerStyle(1); }
    
    // OVERWRITE DIRECT LIGHT FIT RESULT FOR FIBRES AFFECTED BY BELLY PLATES
    if (IS_BELLY_FIBRE) dirfit->SetXYZ(direct.X(), direct.Y(), direct.Z());
    // --> dirfit is one of the output variables, written to file in main()
    
    // Write all objects to file (histograms included automatically, graphs not)
    angles = new TVector2(DIRANG,REFANG);
    angles->Write("angles");
    maxima = new TVector2(nearmax,farmax);
    maxima->Write("maxima");
    limits = new TVector2(HOTLIMIT,COLDLIMIT);
    limits->Write("limits");
    fitDir->Write("fitDir");
    fitRef->Write("fitRef");
    dirfit->Write("dirfit");
    reffit->Write("reffit");
    gDir2D->Write("gDir2D");
    for (int s=0; s<NCOL+2; s++) {
      string ss = Form("[%d]",s);
      if(icos[s]) icos[s]->Write(("icos"+ss).c_str());
      if(gDir[s]) gDir[s]->Write(("gDir"+ss).c_str());
      if(gRef[s]) gRef[s]->Write(("gRef"+ss).c_str());
    }
    pcontD->Write("pcontD");
    pcontR->Write("pcontR");
    pFibPos->Write("pFibPos");
    pFibDir->Write("pFibDir");
    pWgtDir->Write("pWgtDir");
    pWgtRef->Write("pWgtRef");
    pFitDir->Write("pFitDir");
    pFitRef->Write("pFitRef");
    if(nfrontD>0) pcircD1->Write("pcircD1");
    if(nbackD>0)  pcircD2->Write("pcircD2");
    if(nfrontR>0) pcircR1->Write("pcircR1");
    if(nbackR>0)  pcircR2->Write("pcircR2");
    rootfile->Write();
    cout << "Wrote output to file " << root << endl;
    
    scanned_file = 2;
  }
  
  
  // ********************************************************************
  // Plotting section
  // ********************************************************************
    
  // No rotation required, fit was performed in rotated view
  double fitD_rotX = fitDir->X()/1e3;
  double fitD_rotY = fitDir->Y()/1e3;
  double fitR_rotX = fitRef->X()/1e3;
  double fitR_rotY = fitRef->Y()/1e3;
  
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

  // Indicate logarithmic scale
  TLatex *txtL = new TLatex(0.01,0.9,"LOG SCALE");
  txtL->SetTextAlign(13);
  txtL->SetTextFont(102);
  txtL->SetTextSize(0.045);

  // Indicate possible belly plate effect
  TLatex *txtB = new TLatex(-9.25,9.5,"BELLY PLATE");
  txtB->SetTextAlign(13);
  txtB->SetTextFont(102);
  txtB->SetTextSize(0.04);
  
  // Indicate colour scale
  TLatex *txtD = new TLatex(9.5,9.5,Form("Occup. #leq %.2f%%",100.*nearmax));
  txtD->SetTextAlign(33);
  txtD->SetTextFont(82);
  txtD->SetTextSize(0.04);
  TLatex *txtR = new TLatex(9.5,9.5,Form("Occup. #leq %.2f%%",100.*farmax));
  txtR->SetTextAlign(33);
  txtR->SetTextFont(82);
  txtR->SetTextSize(0.04);
  
  // ----------
  // PMT hit count histogram
  pad1->cd()->SetGrid();
  pad1->SetLogx();
  pad1->SetLogy();
  hoccup->SetMinimum(0.5);
  hoccup->SetMaximum(5e3);
  hoccup->SetTitle("PMT occupancy;Occupancy;NPMTs");
  hoccup->SetLineWidth(2);
  hoccup->SetLineColor(1);
  hoccup->GetXaxis()->SetTitleOffset(1.3);
  //hoccup->GetYaxis()->SetTitleOffset(1.3);
  hoccup->Draw();
  
  // Draw off, cold and hot PMT bins in their respective colour
  hoccupoff = (TH1F*)hoccup->Clone("hoccupoff");
  hoccupoff->SetFillColor(16);
  hoccupoff->GetXaxis()->SetRangeUser(1e-6,3e-6);
  hoccupoff->Draw("same");
  hoccuplo = (TH1F*)hoccup->Clone("hoccuplo");
  hoccuplo->SetFillColor(51);
  double coldlimedge;
  for (int b=hoccuplo->GetNbinsX(); b>0; b--) {
    if(hoccuplo->GetBinCenter(b)>COLDLIMIT) continue;
    coldlimedge = hoccuplo->GetBinCenter(b);
    break;
  }
  hoccuplo->GetXaxis()->SetRangeUser(3e-6,coldlimedge);
  hoccuplo->Draw("same");
  hoccuphi = (TH1F*)hoccup->Clone("hoccuphi");
  hoccuphi->SetFillColor(1);
  double hotlimedge;
  for (int b=0; b<hoccuphi->GetNbinsX(); b++) {
    if(hoccuphi->GetBinCenter(b)<HOTLIMIT) continue;
    hotlimedge = hoccuphi->GetBinCenter(b);
    break;
  }
  hoccuphi->GetXaxis()->SetRangeUser(hotlimedge,1.);
  hoccuphi->Draw("same");
  
  // Draw line indicating "hot PMT" limit
  TLine *lhot = new TLine(HOTLIMIT,0,HOTLIMIT,5e3);
  lhot->SetLineWidth(2);
  lhot->SetLineColor(2);
  lhot->Draw("same");
  
  // Draw line indicating "cold PMT" limit
  TLine *lcold = new TLine(COLDLIMIT,0,COLDLIMIT,5e3);
  lcold->SetLineWidth(2);
  lcold->SetLineColor(4);
  lcold->Draw("same");
  
  // ----------
  // Icosahedral projection of detector display
  pad2->cd()->SetGrid();
  //hicos->SetTitle("Detector display (PMT hit sum)");
  //hicos->Draw();
  for(int s=0;s<NCOL+2;s++) { 
    if(!icos[s]) continue;
    if(icos[s]->GetN()==0) continue;
    icos[s]->Draw("P same");
  }
  //phot->Draw("P same");
  pFitDir->Draw("P same");
  pFitRef->Draw("P same");
  pcontD->Draw("P same");
  pcontR->Draw("P same");
  txtL->Draw();
  
  // ----------
  // View from direct light spot (fitted)
  pad3->cd()->SetGrid();
  hfineD->SetTitle("Direct light (PMT hit sum);X' [m];Y' [m]");
  hfineD->GetXaxis()->SetTitleOffset(1.3);
  hfineD->GetYaxis()->SetTitleOffset(1.4);
  hfineD->Draw("scat");
  hfineD->SetStats(0);
  for(int s=0;s<NCOL+2;s++) {
    if(!gDir[s]) continue;
    if(gDir[s]->GetN()==0) continue;
    gDir[s]->Draw("P same");
  }
  pFibDir->Draw("P same");
  pWgtDir->Draw("P same");
  TGraph *pFitDir2 = (TGraph*)pFitDir->Clone();
  if (IS_BELLY_FIBRE) pFitDir2->SetMarkerStyle(28); // open cross = Gaussian fit not used
  pFitDir2->DrawGraph(1,&fitD_rotX,&fitD_rotY,"P same");
  if(pcircD1) pcircD1->Draw("P same");
  if(pcircD2) pcircD2->Draw("P same");
  txtD->Draw();  
  if (IS_BELLY_FIBRE) txtB->Draw();  

  // ----------
  // View from reflected light spot (fitted)
  pad4->cd()->SetGrid();
  hfineR->SetTitle("Reflected light (PMT hit sum);X'' [m];Y'' [m]");
  hfineR->GetXaxis()->SetTitleOffset(1.3);
  hfineR->GetYaxis()->SetTitleOffset(1.4);
  hfineR->Draw("scat");
  hfineR->SetStats(0);
  for(int s=0;s<NCOL+2;s++) {
    if(!gRef[s]) continue;
    if(gRef[s]->GetN()==0) continue;
    gRef[s]->Draw("P same");
  }
  pFibPos->Draw("P same");
  pWgtRef->Draw("P same");
  pFitRef->DrawGraph(1,&fitR_rotX,&fitR_rotY,"P same");
  if(pcircR1) pcircR1->Draw("P same");
  if(pcircR2) pcircR2->Draw("P same");
  txtR->Draw();
  
  // ----------
  // Save canvas and close
  string outfile = "fibre_validation";
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
  if(hoccup) delete hoccup;
  if(hoccupoff) delete hoccupoff;
  if(hoccuplo) delete hoccuplo;
  if(hoccuphi) delete hoccuphi;
  //if(gDir2D) delete gDir2D;
  //if(gRef2D) delete gRef2D;
  if (icos) { for(int s=0; s<NCOL+2; s++) delete icos[s]; }
  if (gDir) { for(int s=0; s<NCOL+2; s++) delete gDir[s]; }
  if (gRef) { for(int s=0; s<NCOL+2; s++) delete gRef[s]; }
  if(pcontD) delete pcontD;
  if(pcontR) delete pcontR;
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
