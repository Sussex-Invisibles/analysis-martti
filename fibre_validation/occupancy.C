// ---------------------------------------------------------
// Goal:            Plot PMT coverage for full PCA dataset
// Author:          Martti Nirkko, 21/11/2017
// Compile & run:   clear && g++ -g -o occupancy.exe occupancy.C `root-config --cflags --libs` -I$RATROOT/include/libpq -I$RATROOT/include -L$RATROOT/lib -lRATEvent_Linux && ./occupancy.exe
// ---------------------------------------------------------

// Helper functions, header files etc.
#include "../include/HelperFunc.C"

// Global constants
const int RUN_CLUSTER = 1;  // whether running on cluster (0=local)
const int IS_MC = 0;        // Monte-Carlo flag
const int MAXCOVERAGE = 10; // maximal coverage for colour scale

// Initialise functions
void occupancy(string, int, int*, int, bool, bool);
void GetProjection(TGraph2D*, TGraph2D*, int*, const int, const RAT::DU::PMTInfo&);
using namespace std;

// Main program
int main(int argc, char** argv) {
  
  // Test flag (0 = process all runs, else run number)
  const int TEST = (RUN_CLUSTER) ? 0 : 101834;

  // Loop over all fibres in list
  string input = "../pca_runs/TELLIE_PCA_Sep2018.txt";
  ifstream in(input.c_str());
  if (!in) { cerr<<"Failed to open "<<input<<endl; exit(1); }
  string line, fibre;
  int node, channel, run, ipw, photons, pin, rms;
  float nhit;
  for (int hdr=0; hdr<2; hdr++) {
    getline(in,line);      // header
  }
  
  // Initialise RAT
  string example = Form("%s/Software/SNOP/work/data/Analysis_r0000110264_s000_p000.root",getenv("HOME"));
  RAT::DU::DSReader dsreader(example);
  const RAT::DU::PMTInfo& pmtinfo = RAT::DU::Utility::Get()->GetPMTInfo();
  const int NPMTS = pmtinfo.GetCount();
  printf("Initialised DSReader & PMTInfo (%d PMTs).\n",NPMTS);
  
  int cvrg[NPMTS];
  memset( cvrg, 0, NPMTS*sizeof(int) ); // NPMTS only known at runtime
  for (int id=0; id<NPMTS; id++) cvrg[id]=-999;
  
  int nfiles=0; 
  while (true) {
    in >> node >> fibre >> channel >> run >> ipw >> photons >> pin >> rms >> nhit;
    if (!in.good()) break;
    if (TEST && TEST!=run) continue; // only want specified run
    //printf("%6s %2d %6d %5d %6d\n", fibre.c_str(), channel, run, ipw, photons);
    occupancy(fibre, run, cvrg, NPMTS, IS_MC, TEST);
    nfiles++;
  }
  printf("Ran over %d files.\n",nfiles);
  if (nfiles==0) { cerr<<"*** ERROR *** No input files found."<<endl; return 1; }

  TGraph2D *coverage = new TGraph2D();
  TGraph2D *offpmts = new TGraph2D();
  GetProjection(coverage, offpmts, cvrg, NPMTS, pmtinfo);
  printf("Graph has %d good PMTs and %d bad PMTs. Coverage was capped to %d.\n", coverage->GetN(), offpmts->GetN(), MAXCOVERAGE);
  
  // Plotting section
  TCanvas *c0 = new TCanvas("","",1600,800);
  gStyle->SetOptStat(0);
  gStyle->SetLabelOffset(999,"XY");
  gStyle->SetCanvasBorderSize(0);
  gStyle->SetFrameBorderSize(0);
  
  // Configure pad
  c0->SetGrid();
  c0->SetLeftMargin(0.005);
  c0->SetRightMargin(0.075); // for TH2 color scale
  c0->SetTopMargin(0.02);
  c0->SetBottomMargin(0.02);
  c0->SetBorderSize(0);
  
  // Draw empty 3D histogram as boundary for 2D graphs
  TH3F *hicos = new TH3F("hicos","",10,0,1,10,0,1,10,0,MAXCOVERAGE); // flat map
  hicos->Draw("a,fb,bb");     // suppress axis, front box, back box
  c0->SetTheta(90-0.001);     // view from above
  c0->SetPhi(0+0.001);        // no x-y rotation
  
  // Draw active PMT coverage
  coverage->SetMarkerStyle(7);
  coverage->SetMinimum(0);                            // TODO - takes >10min (!?)
  coverage->SetMaximum(MAXCOVERAGE);                  // TODO - takes >10min (!?)
  //coverage->GetZaxis()->SetLimits(0,MAXCOVERAGE);     // TODO - doesn't work
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

/// Get the PMT occupancies for a given fibre/run, and fill the coverage array
void occupancy(string fibre, int run, int* cvrg, int NPMTS, bool isMC=false, bool TEST=false) {

  // Check files for given run
  if(!TEST) printf("-----\n");
  printf("Checking files for run %d... ", run);
  string out = Form("./output/PCA_%s.out",fibre.c_str());
  ifstream g(out.c_str());
  if(g.good()) {  // file extracted
    printf("OK! Adding direct light for fibre %s.\n",fibre.c_str());
  } else {        // file not available
    printf("not available! Skipping fibre %s.\n",fibre.c_str());
    return;
  }
  
  // 
  float occupancy[NPMTS];
  int pmtlightcone[NPMTS];
  memset( occupancy, 0., NPMTS*sizeof(float) ); // NPMTS only known at runtime
  memset( pmtlightcone, 0, NPMTS*sizeof(int) ); // NPMTS only known at runtime
  int checkrun, count, totalnhit, npmts;
  int pmtid, pmthits, onspot;
  string dummy;
  g >> dummy >> checkrun;
  if(!g.good()) printf("*** ERROR *** Bad input - str=%s int=%d\n",dummy.c_str(),checkrun);
  if(checkrun != run) { printf("*** ERROR *** Bad run number %d\n",checkrun); return; }
  g >> dummy >> count;
  if(!g.good()) { printf("*** ERROR *** Bad input - str=%s int=%d\n",dummy.c_str(),count); return; }
  g >> dummy >> totalnhit;
  if(!g.good()) { printf("*** ERROR *** Bad input - str=%s int=%d\n",dummy.c_str(),totalnhit); return; }
  g >> dummy >> npmts;
  if(!g.good()) { printf("*** ERROR *** Bad input - str=%s int=%d\n",dummy.c_str(),npmts); return; }
  // Print run info here
  float avgnhit = (float)totalnhit/count;  
  printf("*** INFO *** Run %d has %d EXTA events with an average %.2f nhits.\n",checkrun,count,avgnhit);
  
  // Get occupancy of each PMT
  while (g.good()) {
    g >> pmtid >> pmthits >> onspot;
    if(!g.good()) break;
    occupancy[pmtid]=(float)pmthits/count;
    pmtlightcone[pmtid]=onspot;
  }
  
  // Find thresholds for hot/cold PMTs
  float HOTLIMIT, COLDLIMIT;
  GetLimits(occupancy, NPMTS, HOTLIMIT, COLDLIMIT);
  printf("Found occupancy limits: [ %f, %f ]\n", COLDLIMIT, HOTLIMIT);
  
  // Increase counter for PMTs with more than 1% occupancy
  for(int id=0; id<NPMTS; id++) {
    if (cvrg[id]<0) cvrg[id]=0;                   // PMT ID is valid, reset counter
    if (pmtlightcone[id] != 1) continue;          // not in direct light cone
    if (occupancy[id] == 0) continue;             // broken/disabled tube
    if (occupancy[id] < COLDLIMIT) continue;      // low occupancy tube
    if (occupancy[id] > HOTLIMIT) continue;       // high occupancy tube
    if (occupancy[id] >= 0.01) cvrg[id]++;        // PMT is "covered" by fibre
  }
  
}

void GetProjection(TGraph2D *coverage, TGraph2D* offpmts, int* cvrg, const int NPMTS, const RAT::DU::PMTInfo& pmtinfo) {
  
  // Use icosahedral projection function from HelperFunc.C
  using namespace func;
  TVector2 icospos;
  TVector3 pmtpos;
  printf("Generating icosahedral projection... ");
  int goodpmts=0, badpmts=0;
  for(int id=0; id<NPMTS; id++) {
    pmtpos = pmtinfo.GetPosition(id);
    if (pmtpos.Mag()==0) continue;                // not a valid PMT
    if (pmtinfo.GetType(id) != 1) continue;       // not a normal PMT
    int face;                                     // side of PSUP icosahedron
    icospos = func::IcosProject(pmtpos,face);     // PMT position on flatmap
    if(cvrg[id] <= 0) {                           // inactive PMT
      offpmts->SetPoint(badpmts,icospos.X(),icospos.Y(),0.001); // must be non-zero!
      badpmts++;
    } else {
      if(cvrg[id] > MAXCOVERAGE) {                // high coverage
        cvrg[id] = MAXCOVERAGE;
      } else if (cvrg[id] < 3) {                  // low coverage
        printf("*** WARNING *** PMT #%d is covered by only %d good fibres.\n",id,cvrg[id]);
      }
      coverage->SetPoint(goodpmts,icospos.X(),icospos.Y(),cvrg[id]);
      goodpmts++;
    }
  } // pmt loop
  printf("done.\n");
  
}
