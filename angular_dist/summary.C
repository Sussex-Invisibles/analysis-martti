// ---------------------------------------------------------
// Goal:            Plot angular response of TELLIE fibres
// Author:          Martti Nirkko, 26/04/2017
// Compile & run:   clear && g++ -o summary.exe summary.C `root-config --cflags --libs` -I$RATROOT/include/libpq -I$RATROOT/include -L$RATROOT/lib -lRATEvent_Linux && ./summary.exe
// ---------------------------------------------------------

// Helper functions
#include "../include/HelperFunc.C"

// Global constants
const int RUN_CLUSTER = 1;  // whether running on cluster (0=local)
const int VERBOSE = 1;      // verbosity flag
const int IS_MC = 0;        // Monte-Carlo flag 

// Initialise functions
TH1D* get_histo(string, int, bool, bool);

// Main program
int main(int argc, char** argv) {
  
  // Test flag (0 = process all runs, else run number)
  const int TEST = (RUN_CLUSTER) ? 0 : 102315;// 101410;

  // Initialise input file
  string input = "ANGULAR_FITRESULTS.txt";
  ifstream in(input.c_str());
  if (!in) { cerr<<"Failed to open "<<input<<endl; exit(1); }
  string line;
  for (int hdr=0; hdr<2; hdr++) {
    getline(in,line);      // header
  }

  // Initialise graph
  TGraph* graph[95] = {NULL};
  TF1* func[95] = {NULL}; 
 
  // Loop over all fibres in list
  string fibre;
  int run, chisq, ndof;
  float a, sa, b, sb;
  float xmin=1e6, ymin=1e6, amin=1e6, bmin=1e6;
  float xmax=-1e6, ymax=-1e6, amax=-1e6, bmax=-1e6;
  float xavg=0, yavg=0, aavg=0, bavg=0;
  float f0, f1, f, sf;
  int nfiles=0;
  while (true) {
    in >> fibre >> run >> a >> sa >> b >> sb >> chisq >> ndof;
    if (!in.good()) break;
    if (TEST && TEST!=run) continue; // only want specified run
    graph[nfiles] = new TGraph();
    graph[nfiles]->SetTitle(fibre.c_str());
    graph[nfiles]->SetPoint(0,a,b);
    func[nfiles] = new TF1(fibre.c_str(), "[0] - [1] + [1]/cos(x/180.*pi)", 0, 24);
    func[nfiles]->SetParameters(a,b);
    f0 = func[nfiles]->GetMinimum();
    f1 = func[nfiles]->GetMaximum();
    f  = (f1+f0)/2.;
    //f = func[nfiles]->Mean(0,24);
    //sf = func[nfiles]->GetRMS();
    nfiles++;

    // Get max/min/avg
    if (a<amin) amin=a;
    if (b<bmin) bmin=b;
    if (f0<ymin) ymin=f0;
    if (a>amax) amax=a;
    if (b>bmax) bmax=b;
    if (f1>ymax) ymax=f1;
    //if (fibre=="FT047A" || fibre=="FT067A") continue; // skip bad fibres for average
    aavg+=a;
    bavg+=b;
    yavg+=f;
  }
  printf("Ran over %d files.\n",nfiles);
  if (nfiles==0) { 
    cerr<<"*** ERROR *** No input files found."<<endl;
    return 1; 
  }
  aavg/=nfiles;
  bavg/=nfiles;
  yavg/=nfiles;
  
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
  c0->cd(1)->DrawFrame(0,ymin-5,24,ymax+5,"TELLIE angular systematic fits;Angle of PMT w.r.t. fitted fibre direction [deg];Mean hit time offset [ns]");
  for (int fib=0; fib<95; fib++) {
    if (!func[fib]) continue;
    func[fib]->SetLineWidth(1);
    func[fib]->SetLineColor(100-fib);
    func[fib]->Draw("L same");
   
    // Highlight unusual fibres
    float yval = func[fib]->GetParameter(0);
    if(yval > 5) {
      tname->DrawLatex(10,yval-1.5,func[fib]->GetName());
    }
  }
  // Box for text
  TBox *tbox = new TBox(0.6,ymax+2,14.8,ymax+4.5);
  tbox->SetLineColor(2);
  tbox->SetFillColor(kYellow-9);
  tbox->Draw("L same");
  // Fit function as text
  TLatex *text = new TLatex(1,ymax+4,"Fit function: y = a #plus b#left(#frac{1}{cos(x)} #minus 1#right)");
  text->SetTextAlign(13);
  text->SetTextFont(62);
  text->SetTextSize(0.04);
  text->Draw("same");

  // Fit lines (zoomed in)
  c0->cd(2)->SetGrid();
  c0->cd(2)->DrawFrame(0,yavg-3,16,yavg+1,"TELLIE angular systematic fits (zoom);Angle of PMT w.r.t. fitted fibre direction [deg];Mean hit time offset [ns]");
  for (int fib=0; fib<95; fib++) {
    if (!func[fib]) continue;
    func[fib]->SetLineWidth(1);
    func[fib]->SetLineColor(100-fib);
    func[fib]->Draw("L same");
    
    // Highlight unusual fibres
    float yval = func[fib]->GetParameter(0);
    if(yval > 0) {
      tname->DrawLatex(2,yval+0.1,func[fib]->GetName());
    }
  }
  
  // Fit parameters (full range)
  c0->cd(3)->SetGrid();
  c0->cd(3)->DrawFrame(amin-2,bmin-10,amax+2,bmax+10,"Fit parameters;Fit parameter a [ns];Fit parameter b [ns]");
  for (int fib=0; fib<95; fib++) {
    if (!graph[fib]) continue;
    graph[fib]->SetMarkerStyle(7);
    graph[fib]->SetMarkerColor(100-fib);
    graph[fib]->Draw("P same");
    
    // Highlight unusual fibres
    double xval,yval;
    graph[fib]->GetPoint(0,xval,yval);
    if(yval < 10) {
      tname->SetTextAlign(23);
      tname->DrawLatex(xval-0.1,yval-2,graph[fib]->GetTitle());
    }
  }
  
  // Fit parameters (zoomed in)
  c0->cd(4)->SetGrid();
  c0->cd(4)->DrawFrame(aavg-1,bavg-30,aavg+1,bavg+30,"Fit parameters (zoom);Fit parameter a [ns];Fit parameter b [ns]");
  for (int fib=0; fib<95; fib++) {
    if (!graph[fib]) continue;
    graph[fib]->SetMarkerStyle(7);
    graph[fib]->SetMarkerColor(100-fib);
    graph[fib]->Draw("P same");
    /*
    // Highlight unusual fibres
    double xval,yval;
    graph[fib]->GetPoint(0,xval,yval);
    if(xval > aavg+0.7 || yval < bavg-20) {
      tname->SetTextAlign(12);
      tname->DrawLatex(xval+0.03,yval,graph[fib]->GetTitle());
    }
    */
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
