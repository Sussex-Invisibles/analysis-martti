// ---------------------------------------------------------
// Goal:          Evaluate gamma peaks in simulated AmBe events
// Author:        Martti Nirkko, 10/01/2018
// Compile & run: clear && g++ -g -o testMC.exe testMC.C `root-config --cflags --libs` -I$RATROOT/include/libpq -I$RATROOT/include -L$RATROOT/lib -lRATEvent_Linux -lGeom && ./testMC.exe
// ---------------------------------------------------------

// Helper functions (includes everything else)
#include "../HelperFunc.C"

// ROOT functions to get elements
#include <TGeoElement.h>
#include <TGeoManager.h>

using std::setw;
using std::setprecision;
using std::fixed;

// Global parameters
const int TEST = 0;                    // test flag
const int VERBOSE = 0;                 // verbosity flag
const int isMC = 1;                    // Monte-Carlo flag
const int NBINS = 120;                 // number of bins for most histograms
const double EMAX = 12;                // maximal energy for histograms
const int OLD_GEOM = 0;                // old geometry (filled stem instead of rod)
const double Z_OFFSET = -17.5;         // mm offset of the AmBe source

// Turn points into a graph
void FillGraphs(TGraph *g, int n, double *x, double *y, double *z, int col) {
  new (g) TGraph();
  g->SetPoint(0,y[0],z[0]+Z_OFFSET);
  for (int i=0; i<n; i++) {
    g->SetPoint(i+1,x[i],z[i]+Z_OFFSET);
    g->SetPoint(n+i+1,y[n-i-1],z[n-i-1]+Z_OFFSET);
  }
  g->SetLineWidth(1);
  g->SetLineColor(col);
};

// Get table of element names in ROOT
std::vector<std::string> InitialiseElements() {
  new TGeoManager("world", "the simplest geometry");
  TGeoElementTable *table = gGeoManager->GetElementTable();
  std::vector<std::string> elements;
  for (int i=0; i<112; i++) {
    elements.push_back(table->GetElement(i+1)->GetName());
    if (elements[i].length()>1) elements[i][1] = std::tolower( elements[i][1] );
    if (elements[i].length()>2) elements[i][2] = std::tolower( elements[i][2] );
    //std::cout << i+1 << elements[i] << std::endl;
  }
  std::cout << "Initialised " << (int)elements.size() << " elements." << std::endl;
  return elements;
};

// Get isotope name from PDG code
std::string GetIsotope(int pdg, std::vector<std::string>& elements) {
  std::string isotope = "NULL";
  // Check that elements exist
  if (elements.size()==0) {
    std::cout << "ERROR - Table of elements not initialised!" << std::endl;
    return isotope;
  }
  // PDG code must be of nucleus
  if (fabs(pdg) < 1e9) {
    std::cout << "ERROR - Can only convert nuclear isotopes!" << std::endl;
    return isotope;
  }
  // Nuclear codes are given as 10-digit numbers ±10LZZZAAAI.
  int Z = (pdg % 10000000) / 10000; // number of protons
  int A = (pdg %    10000) /    10; // number of nucleons
  isotope = Form("%s%d",elements[Z-1].c_str(),A); // Z-1 because of vector index
  return isotope;
};

// Main program
int main(int argc, char** argv) {

  // Initialise elements
  std::vector<std::string> elements = InitialiseElements();
  
  // Materials used as neutron absorber
  string MATERIALS[] = {"Chromium","Gadolinium","Iron","Lead","Nickel","Titanium","Zinc"};
  //string MATERIAL = "Nickel_OLD";
  if (TEST) MATERIALS[0] = "Test";
  int NMAT = sizeof(MATERIALS)/sizeof(string);
  
  for (int mat=0; mat<NMAT; mat++) {
  
    if (TEST && mat>0) continue;
    string MATERIAL = MATERIALS[mat];
    
    // Input/Output files
    string fpath = ((TEST) ? "test7" : MATERIAL);
    string inroot = fpath;
    if (!TEST) inroot = inroot + "/output";
    inroot = inroot + "/AmBe_*.root";

    string outlog = fpath;
    if (!TEST) outlog = outlog + "/results";
    outlog = outlog + "/output.log";
   
    // Plotting options 
    gStyle->SetOptStat(0);
    gStyle->SetTitleBorderSize(0);
    gStyle->SetHistLineWidth(2);
    
    // Initialise stuff
    TFile *outfile=NULL;
    TH1D *hGam=NULL, *hEle=NULL, *hNeu=NULL, *hPro=NULL, *hOth=NULL;
    TH1D *hGamAV=NULL, *hGamAVn=NULL, *hGamFrac=NULL, *hGamFracn=NULL;
    TH2D *hGamTime=NULL, *hnCapPos=NULL, *hnCapPosZ=NULL;
    TH2I *hnCapIso=NULL;
    
    string outroot = fpath;
    if (!TEST) outroot = outroot + "/results";
    outroot = outroot + "/histograms.root";
    std::ifstream chk(outroot.c_str());
    bool processed = chk.good();
    if (processed) {  // if already processed, read output file
      outfile = new TFile(outroot.c_str(),"READ");
      hGam      = (TH1D*)outfile->Get("hGam");
      hEle      = (TH1D*)outfile->Get("hEle");
      hNeu      = (TH1D*)outfile->Get("hNeu");
      hPro      = (TH1D*)outfile->Get("hPro");
      hOth      = (TH1D*)outfile->Get("hOth");
      hGamAV    = (TH1D*)outfile->Get("hGamAV");
      hGamAVn   = (TH1D*)outfile->Get("hGamAVn");
      hGamFrac  = (TH1D*)outfile->Get("hGamFrac");
      hGamFracn = (TH1D*)outfile->Get("hGamFracn");
      hGamTime  = (TH2D*)outfile->Get("hGamTime");
      hnCapPos  = (TH2D*)outfile->Get("hnCapPos");
      hnCapPosZ = (TH2D*)outfile->Get("hnCapPosZ");
      hnCapIso  = (TH2I*)outfile->Get("hnCapIso");
    } else {          // if not processed, read input file & generate output
      outfile = new TFile(outroot.c_str(),"NEW");
      std::ofstream out(outlog.c_str());
      
      // Initialise RAT
      RAT::DU::DSReader dsreader(inroot);
      const RAT::DU::PMTInfo& pmtinfo = RAT::DU::Utility::Get()->GetPMTInfo();
      const int NPMTS = pmtinfo.GetCount();

      // Initialise histograms
      hGam      = new TH1D("hGam","Gamma",NBINS,0,EMAX);
      hEle      = new TH1D("hEle","Electron",NBINS,0,EMAX);
      hNeu      = new TH1D("hNeu","Neutron",NBINS,0,EMAX);
      hPro      = new TH1D("hPro","Proton",NBINS,0,EMAX);
      hOth      = new TH1D("hOth","Other",NBINS,0,EMAX);
      hGamAV    = new TH1D("hGamAV","Gamma (AV)",NBINS,0,EMAX);
      hGamAVn   = new TH1D("hGamAVn","Gamma (AV) from nCapture",NBINS,0,EMAX);
      hGamFrac  = new TH1D("hGamFrac","Gamma (AV frac.)",NBINS,0,EMAX);
      hGamFracn = new TH1D("hGamFracn","Gamma (AV frac.) from nCapture",NBINS,0,EMAX);
      hGamTime  = new TH2D("hGamTime","Gamma (AV) timing",NBINS,0,EMAX,NBINS,0,4e6);
      hnCapPos  = new TH2D("hnCapPos","Neutron capture position",200,0,100,400,-500,500);
      hnCapPosZ = new TH2D("hnCapPosZ","Neutron capture position",1000,0,100,2000,-500,500);
      hnCapIso  = new TH2I("hnCapIso","Neutron capture isotope",100,0,100,100,0,100);
      BinLog(hGamTime->GetYaxis(), 0.04); // log axis for timing (0.04 ns - 4 ms)
      
      // Loop over all entries in file
      int nevents = 0, ntracks = 0;
      for(int iEntry=0; iEntry<dsreader.GetEntryCount(); iEntry++) {
        const RAT::DS::Entry& ds = dsreader.GetEntry(iEntry);
        nevents++;
        
        // Print progress
        printProgress(iEntry, dsreader.GetEntryCount());
        
        // Initialise MC for this event
        RAT::DS::MC mc;
        if (isMC) mc = ds.GetMC();  // don't initialise this for real data (crashes)
        
        /*
        // Loop over all MC *primary* particles
        for(int iPart=0; iPart<mc.GetMCParticleCount(); iPart++) {
          const RAT::DS::MCParticle& part = mc.GetMCParticle(iPart);
          
          int pdg = part.GetPDGCode();
          double eTrk = part.GetKineticEnergy();
          if (pdg == 2112) hNeu->Fill(eTrk);     // neutron
          else if (pdg == 22) hGam->Fill(eTrk);  // gamma
          else cout << "Particle " << pdg << " with energy " << eTrk << endl;
        }
        */
        
        // Loop over all MC tracks
        out << std::right;
        for(int iTrk=0; iTrk<mc.GetMCTrackCount(); iTrk++) {
          const RAT::DS::MCTrack& trk = mc.GetMCTrackFromIndex(iTrk);
          ntracks++;      
          
          int pdg = trk.GetPDGCode();
          double eTrk = trk.GetFirstMCTrackStep().GetKineticEnergy();
          if (pdg > 1e9) {           // nucleus
            // Test for neutron capture isotopes
            if (trk.GetFirstMCTrackStep().GetProcess() == "nCapture") {
	      if (fabs(pdg) > 1e9) { // Nuclear codes are given as 10-digit numbers ±10LZZZAAAI.
                int Z = (pdg % 10000000) / 10000; // atomic number
                int A = (pdg % 10000) / 10;       // number of nucleons
                hnCapIso->Fill(A-1-Z,Z);          // parent nucleus that captured the neutron
                if (VERBOSE) {
                  out << "Event " << setw(2) << iEntry << ", track " << setw(2) << iTrk;
                  out << " is " << setw(8) << trk.GetParticleName();
                  out << " (" << setw(5) << setprecision(3) << fixed << eTrk << " MeV)";
                  out << ", parent " << setw(7) << mc.GetMCTrack(trk.GetParentID()).GetParticleName();
	          string iso = GetIsotope(pdg-10, elements);  // same nucleus with 1 less neutron
	          out << " captured on " << iso << endl;
	          // Names of particles lighter than nuclei can be obtained as follows:
	          //if (nucpdg==1000010010) nucpdg = 2212; // identify H1 as proton
	          //const TParticlePDG *part = pdgdb->GetParticle(nucpdg);
                }
              }
            }
          } else if (pdg == 2212) {  // proton
            hPro->Fill(eTrk);
          } else if (pdg == 2112) {  // neutron
            hNeu->Fill(eTrk);
          } else if (pdg == 11) {    // electrons
            hEle->Fill(eTrk);
          } else if (pdg ==   22) {  // gamma
            hGam->Fill(eTrk);
          } else {                   // other
            hOth->Fill(eTrk);
            if (pdg == -11) {        // positrons (don't care)
              continue;
            }
            if (fabs(pdg) == 12) {   // neutrinos (really don't care)
              continue;
            }
            out << "*** INFO - Event " << setw(2) << iEntry << ",";
            out << " track " << setw(2) << iTrk << ":";
            out << " Unusual particle (" << setw(6) << setprecision(3) << fixed << eTrk << " MeV)";
            out << " with PDG " << setw(4) << pdg;
            out << " (" << setw(4) << trk.GetParticleName() << ")" << endl;
          }
          
          // Print information for high-energy gammas to file - from Phys. Rev. 89, 375:
          // "The nickel spectrum contains an intense gamma-ray with an energy of (8.997 ±  0.005) MeV, 
          // which is produced in 50 percent of the captures in Ni-58. Another prominent gamma-ray at 
          // (8.532 ± 0.008) MeV may represent the transition to the ground state in Ni-61; if so, it 
          // accounts for some 80 percent of captures in Ni-60. From considerations of intensity, five 
          // of the remaining nickel y-rays can be ascribed to transitions to excited states in Ni-59."
          if (pdg!=22) continue;
          for (int iStep=0; iStep<trk.GetMCTrackStepCount(); iStep++) {
            const RAT::DS::MCTrackStep& stp = trk.GetMCTrackStep(iStep);
            double eStep = stp.GetKineticEnergy();
            if (VERBOSE) { 
              if(iStep==0) {
                out << "-----" << endl;
                out << "Event " << iEntry << ", track " << iTrk << " is from gamma (" << setprecision(3) << eTrk << " MeV):" << endl;
              }
              out << "Step #" << setw(2) << iStep;
              out << " (" << setw(6) << setprecision(3) << fixed << eStep << " MeV)";
              out << ": " << setw(20) << stp.GetStartVolume();
              out << " --> " << setw(20) << std::left << stp.GetEndVolume();
              out << " (" << setw(6) << std::right << setprecision(1) << fixed << stp.GetLength() << " mm)";
              out << " process = " << stp.GetProcess() << endl;
            }
            if (stp.GetEndVolume() != "inner_av") continue; // get first step that ends in inner AV
            hGamAV->Fill(eStep); // gamma energy when entering AV (visible energy)
            hGamTime->Fill(eStep,stp.GetGlobalTime());
            // Gammas originating from neutron capture
            if (trk.GetFirstMCTrackStep().GetProcess() == "nCapture") {
              hGamAVn->Fill(eStep);
              TVector3 pos = trk.GetFirstMCTrackStep().GetPosition();
              hnCapPos->Fill(pos.Perp(),pos.Z());
              hnCapPosZ->Fill(pos.Perp(),pos.Z());
            }
            break;
          } // steps
        } // tracks
      } // entries
      out << "-----" << endl;
      out << "Found " << nevents << " MC events with " << ntracks << " tracks, consisting of:" << endl;
      out << "- Protons:   " << setw(6) << (int)hPro->GetEntries() << " (" << hPro->GetEntries()/nevents << " / event)" << endl;
      out << "- Neutrons:  " << setw(6) << (int)hNeu->GetEntries() << " (" << hNeu->GetEntries()/nevents << " / event)" << endl;
      out << "- Electrons: " << setw(6) << (int)hEle->GetEntries() << " (" << hEle->GetEntries()/nevents << " / event)" << endl;
      out << "- Gammas:    " << setw(6) << (int)hGam->GetEntries() << " (" << hGam->GetEntries()/nevents << " / event)" << endl;
      out << "- Other:     " << setw(6) << (int)hOth->GetEntries() << " (" << hOth->GetEntries()/nevents << " / event)" << endl;
      out.close();
      
      // Save output file
      outfile->Write();
    }

    // Plotting section
    string imgname = fpath;
    if (!TEST) imgname = imgname + "/results";
    imgname = imgname + "/AmBe_" + MATERIAL;
    
    // Energy of various particles 
    double maxval = std::max(hPro->GetMaximum(),hGam->GetMaximum());
    TCanvas *c = new TCanvas("c","",1000,600);
    c->SetLogy();
    c->SetGrid();
    TH1F *pad = c->DrawFrame(0,0.5,EMAX,1.5*maxval);
    pad->SetTitle("AmBe simulation with neutron absorber;Kinetic energy [MeV];MC tracks");
    pad->Draw("");
    hPro->SetLineColor(3); hPro->Draw("same");
    hNeu->SetLineColor(4); hNeu->Draw("same");
    hGam->SetLineColor(2); hGam->Draw("same");
    hOth->SetLineColor(1); hOth->Draw("same");
    TLegend *leg = new TLegend(0.7,0.7,0.89,0.89);
    leg->AddEntry(hPro);
    leg->AddEntry(hNeu);
    leg->AddEntry(hGam);
    leg->AddEntry(hOth);
    leg->Draw();
    c->Print((imgname+"_1.png").c_str());
    c->Close();
    
    // Gamma rays that reach inner AV (absolute)
    c = new TCanvas("c","",1000,600);
    c->SetLogy();
    c->SetGrid();
    pad = c->DrawFrame(0,0.5,EMAX,1.5*maxval);
    pad->SetTitle("AmBe #gamma-rays entering inner AV;Photon energy [MeV];MC tracks");
    pad->Draw("");
    hGamAVn->SetLineColorAlpha(2,0.0);
    hGamAVn->SetFillColorAlpha(2,0.5);
    hGamAVn->Draw("same");
    hGamAV->SetLineColor(2);
    hGamAV->Draw("same");
    leg = new TLegend(0.59,0.74,0.89,0.89);
    leg->AddEntry(hGamAV);
    leg->AddEntry(hGamAVn);
    leg->Draw();
    c->Print((imgname+"_2.png").c_str());
    c->Close();
    
    // Gamma rays that reach inner AV (relative)
    c = new TCanvas("c","",1000,600);
    c->SetLogy();
    c->SetGrid();
    pad = c->DrawFrame(0,0.01,EMAX,10.);
    pad->SetTitle("Fraction of #gamma-rays entering inner AV;Photon energy [MeV];Fraction");
    pad->Draw("");
    for(int i=0;i<=NBINS+1;i++) {
      if (hGam->GetBinContent(i)==0) {
        hGamFrac->SetBinContent(i,0);
        hGamFracn->SetBinContent(i,0);
      } else {
        hGamFrac->SetBinContent(i,(float)hGamAV->GetBinContent(i)/(float)hGam->GetBinContent(i));
        hGamFracn->SetBinContent(i,(float)hGamAVn->GetBinContent(i)/(float)hGam->GetBinContent(i));
      }
    }
    hGamFracn->SetLineColorAlpha(2,0.0);
    hGamFracn->SetFillColorAlpha(2,0.5);
    hGamFracn->Draw("same");
    hGamFrac->SetLineColor(2);
    hGamFrac->Draw("same");
    leg = new TLegend(0.59,0.74,0.89,0.89);
    leg->AddEntry(hGamAV);
    leg->AddEntry(hGamAVn);
    leg->Draw();
    c->Print((imgname+"_3.png").c_str());
    c->Close();
    
    // Gamma rays that reach inner AV (time vs energy)
    c = new TCanvas("c","",1000,600);
    c->SetLogy();
    c->SetLogz();
    c->SetGrid();
    hGamTime->SetTitle("#gamma-rays entering inner AV;Photon energy [MeV];Global time [ns]");
    hGamTime->Draw("colz");
    c->Print((imgname+"_4.png").c_str());
    c->Close();  

    // Envelope volume (Air)
    double z0[] = { 0.4, 60.7, 60.7, 73.8, 73.8, 173.5, 369.3, 369.3, 378.5};
    double y0[] = { 0.0,  0.0,  0.0,  0.0,  0.0,   0.0,   0.0,   0.0,   0.0};
    double x0[] = {22.2, 22.2, 14.8, 14.8,  7.2,  14.4,  14.4,  31.9,  31.9};
    int n0 = sizeof(z0)/sizeof(double);
    TGraph *g0 = new TGraph();
    FillGraphs(g0,n0,x0,y0,z0,16);
    
    // Source Can (Steel)
    double z1[] = { 0.5,  2.5,  2.5, 57.6, 57.6, 60.6, 60.6, 66.95};
    double y1[] = { 0.0,  0.0, 21.0, 21.0,  0.0,  0.0,  0.0,  0.0};
    double x1[] = {22.13, 22.13, 22.13, 22.13, 22.13, 22.13, 14.72, 14.72};
    int n1 = sizeof(x1)/sizeof(double);
    TGraph *g1 = new TGraph();
    FillGraphs(g1,n1,x1,y1,z1,4);
    
    // Gamma Shield (Lead)
    double z2[] = { 2.5,  4.5,  4.5, 55.6, 55.6, 57.6};
    double y2[] = { 0.0,  0.0, 19.0, 19.0,  0.0,  0.0};
    double x2[] = {21.0, 21.0, 21.0, 21.0, 21.0, 21.0};
    int n2 = sizeof(x2)/sizeof(double);
    TGraph *g2 = new TGraph();
    FillGraphs(g2,n2,x2,y2,z2,3);
    
    // Source Body (Steel?)
    double z3[] = { 5.0, 39.8, 39.8, 49.8, 49.8, 55.1};
    double y3[] = { 0.0,  0.0,  0.0,  0.0,  0.0,  0.0};
    double x3[] = {14.6, 14.6, 16.2, 16.2, 14.6, 14.6};
    int n3 = sizeof(x3)/sizeof(double);
    TGraph *g3 = new TGraph();
    FillGraphs(g3,n3,x3,y3,z3,2);
    
    // Teflon Stem (PTFE)
    double z4[] = {67.35, 73.75, 73.75, 73.75, 173.45, 369.35, 369.35, 378.45};
    double y4[] = {  0.0,   0.0,   0.0,  5.05,   5.05,   5.05,    0.0,    0.0};
    double x4[] = {14.72, 14.72,  7.15,  7.15,  14.30,  14.30,  31.81,  31.81};
    if (OLD_GEOM) { y4[4]=12.3; y4[5]=12.3; } // wider rod
    int n4 = sizeof(x4)/sizeof(double);
    TGraph *g4 = new TGraph();
    FillGraphs(g4,n4,x4,y4,z4,1);
    
    // Absorber Rod (neutron absorber material)
    double z5[] = {74.15, 173.45, 368.95};
    double y5[] = { 0.0,    0.0,    0.0};
    double x5[] = { 5.0,    5.0,    5.0};
    if (OLD_GEOM) { x5[1]=12.25; x5[2]=12.25; } // wider rod
    int n5 = sizeof(x5)/sizeof(double);
    TGraph *g5 = new TGraph();
    FillGraphs(g5,n5,x5,y5,z5,6);
   
    // Neutron capture positions (z vs r)
    c = new TCanvas("c","",1000,600);
    c->Divide(2,1);
    // Scatterplot
    c->cd(1)->SetGrid();
    c->cd(1)->SetLeftMargin(0.12);
    c->cd(1)->SetRightMargin(0.14); // same pad width as right pad
    c->cd(1)->DrawFrame(0,-100,100,400,"Position of neutron captures;r [mm];z [mm]");
    hnCapPosZ->SetMarkerStyle(1);
    hnCapPosZ->Draw("scat same");
    // Source geometry
    g1->Draw("L same");
    g2->Draw("L same");
    g3->Draw("L same");
    g4->Draw("L same");
    g5->Draw("L same");
    // Normalised by volume
    c->cd(2)->SetGrid();
    c->cd(2)->SetLogz();
    c->cd(2)->SetLeftMargin(0.12);
    c->cd(2)->SetRightMargin(0.14); // colz
    c->cd(2)->DrawFrame(0,-100,100,400,"Neutron captures by volume [mm^{#minus3}];r [mm];z [mm]");
    // Normalise neutron capture histograms (radial dependence)
    for (int i=0; i<=hnCapPos->GetNbinsX(); i++) {
      double dx = hnCapPos->GetXaxis()->GetBinWidth(i);
      double rad = hnCapPos->GetXaxis()->GetBinCenter(i);
      for (int j=0; j<=hnCapPos->GetNbinsY(); j++) {
        double dy = hnCapPos->GetYaxis()->GetBinWidth(j);
        double val = hnCapPos->GetBinContent(i,j);
        hnCapPos->SetBinContent(i,j,val/(2.*pi*rad*dx*dy));
      }
    }
    TH2D *hnCapPos2 = (TH2D*)hnCapPos->Rebin2D(3,3,"hnCapPos2");
    hnCapPos2->Draw("colz same");
    // Source geometry
    g1->Draw("L same");
    g2->Draw("L same");
    g3->Draw("L same");
    g4->Draw("L same");
    g5->Draw("L same");
    c->Print((imgname+"_5.png").c_str());
    c->Close();
    
    // Neutron capture positions (z vs r) - zoomed in on source
    c = new TCanvas("c","",1000,600);
    c->Divide(2,1);
    // Scatterplot
    c->cd(1)->SetGrid();
    c->cd(1)->SetLeftMargin(0.12);
    c->cd(1)->SetRightMargin(0.14); // same pad width as right pad
    c->cd(1)->DrawFrame(0,-30,40,170,"Position of neutron captures;r [mm];z [mm]");
    hnCapPosZ->SetMarkerStyle(6);
    hnCapPosZ->Draw("scat same");
    // Source geometry
    g1->Draw("L same");
    g2->Draw("L same");
    g3->Draw("L same");
    g4->Draw("L same");
    g5->Draw("L same");
    // Normalised by volume
    c->cd(2)->SetGrid();
    c->cd(2)->SetLogz();
    c->cd(2)->SetLeftMargin(0.12);
    c->cd(2)->SetRightMargin(0.14); // colz
    c->cd(2)->DrawFrame(0,-30,40,170,"Neutron captures by volume [mm^{#minus3}];r [mm];z [mm]");
    // Neutron capture histograms already normalised
    hnCapPos->Draw("colz same");
    // Source geometry
    g1->Draw("L same");
    g2->Draw("L same");
    g3->Draw("L same");
    g4->Draw("L same");
    g5->Draw("L same");
    c->Print((imgname+"_5_zoom.png").c_str());
    c->Close();
 
    // Neutron capture isotopes (A vs Z)
    c = new TCanvas("c","",1000,800);
    c->SetLogz();
    c->SetGrid();
    c->SetRightMargin(0.14);
    hnCapIso->SetTitle("Neutron capture isotopes;N;Z");
    hnCapIso->Draw("colz");
    c->Print((imgname+"_6.png").c_str());
    c->Close();
    
    // Free memory
    if (hPro) delete hPro;
    if (hNeu) delete hNeu;
    if (hEle) delete hEle;
    if (hGam) delete hGam;
    if (hOth) delete hOth;
    if (hGamAV) delete hGamAV;
    if (hGamAVn) delete hGamAVn;
    if (leg) delete leg;
    if (c) delete c;
    
    // Close files and exit
    outfile->Close();

  }

  return 0;
  
}
