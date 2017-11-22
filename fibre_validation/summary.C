// ---------------------------------------------------------
// Goal:            Summarise fibre validation fit results
// Author:          Martti Nirkko, 22/11/2017
// Compile & run:   clear && g++ -o summary.exe summary.C `root-config --cflags --libs` -I$RATROOT/include/libpq -I$RATROOT/include -L$RATROOT/lib -lRATEvent_Linux && ./summary.exe
// ---------------------------------------------------------

// Helper functions (includes everything else)
#include "../HelperFunc.C"

// Run time parameters
const int RUN_CLUSTER = 1;      // whether running on cluster (0=local)
const int USE_RATDB = 0;        // whether to use RATDB to get fibre positions
const int VERBOSE = 1;          // verbosity flag
const int IS_MC = 0;            // Monte-Carlo flag 
const int NBINS = 24;           // number of bins

// Initialise functions
int summary(string, int, float&, float&, bool, bool);

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
  
  // Initialise variables and histograms
  string fibre;
  int node, channel, run, ipw, photons, pin, rms;
  float nhit;
  TH1F *hDir = new TH1F("hDir","",NBINS,0,12);
  TH1F *hRef = new TH1F("hRef","",NBINS,0,12);
  
  // Loop over all files in list
  int nfiles=0;
  while (true) {
    
    // Check for correct run
    in >> node >> fibre >> channel >> run >> ipw >> photons >> pin >> rms >> nhit;
    if (!in.good()) break;
    if (TEST && TEST!=run) continue; // only want specified run
    
    // Get fit results for fibre
    float dirAng, refAng;
    int errors = summary(fibre, run, dirAng, refAng, IS_MC, TEST);
    if (errors) {
      cerr<<"*** WARNING *** Run "<<run<<" was not processed correctly."<<endl;
    } else {
      hDir->Fill(dirAng);
      hRef->Fill(refAng);
      nfiles++;
    }
  }
  
  // Plot histograms
  gStyle->SetOptStat(0);
  TCanvas *c0 = new TCanvas("c0","",800,600);
  c0->SetGrid();
  hDir->SetTitle("Fibre validation summary;Fit result deviation [deg];Fibres");
  hDir->SetLineColor(2);
  hRef->SetLineColor(4);
  hDir->SetLineWidth(2);
  hRef->SetLineWidth(2);
  hDir->Draw();
  hRef->Draw("same");
  hDir->GetYaxis()->SetRangeUser(0,1.1*std::max(hDir->GetMaximum(),hRef->GetMaximum()));
  TLegend *leg = new TLegend(0.6,0.7,0.88,0.88);
  leg->AddEntry(hDir,"Direct light","L");
  leg->AddEntry(hRef,"Reflected light","L");
  leg->Draw("same");
  c0->Print("summary.png");
  c0->Print("summary.pdf");
  c0->Close();
  
  // Summarize, delete pointers and exit
  printf("Ran over %d files.\n",nfiles);
  if (c0) delete c0;
  if (leg) delete leg;
  if (hDir) delete hDir;
  if (hRef) delete hRef;
  if (nfiles==0) {
    cerr<<"*** ERROR *** Did not process any files."<<endl;
    return 1; 
  }
  return 0;
}

// Define macro
int summary(string fibre, int run, float& dirAng, float& refAng, bool isMC=false, bool TEST=false) {
  
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
  TVector3 fitDirPos(0,0,0);
  TVector3 fitRefPos(0,0,0);
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
    fitDirPos.SetXYZ(dirx,diry,dirz);
    fitRefPos.SetXYZ(refx,refy,refz);
    if (fitDirPos.Mag()!=0) break;
    if (fitRefPos.Mag()!=0) break;
  }
  TVector3 fitdir = (fitDirPos-fibrePos).Unit();
  if (fitDirPos.Mag()==0) {
    cerr << "*** ERROR - Could not load fit position!" << endl;
    exit(1);
  } else {
    cout << "Loaded fit position: " << printVector(fitDirPos) << " mm, " << fitDirPos.Angle(lightPos)*180./pi << " deg deviation" << endl;
  }
  cout << "Deviation seen from fibre: " << fitdir.Angle(fibreDir)*180./pi << " degrees (CHECK: 2.0 == " << fitDirPos.Angle(lightPos)/fitdir.Angle(fibreDir) << ")." << endl << endl;
  
  // Output
  dirAng = fitDirPos.Angle(lightPos)*180./pi;
  refAng = fitRefPos.Angle(fibrePos)*180./pi;
  return 0;
}

