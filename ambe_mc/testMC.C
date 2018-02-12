// ---------------------------------------------------------
// Goal:          Evaluate gamma peaks in simulated AmBe events
// Author:        Martti Nirkko, 10/01/2018
// Compile & run: clear && g++ -g -o testMC.exe testMC.C `root-config --cflags --libs` -I$RATROOT/include/libpq -I$RATROOT/include -L$RATROOT/lib -lRATEvent_Linux && ./testMC.exe
// ---------------------------------------------------------

// Helper functions (includes everything else)
#include "../HelperFunc.C"
using std::setw;
using std::setprecision;
using std::fixed;

// Global parameters
const int TEST = 0;
const int VERBOSE = 0;
const int isMC = 1;
const int NBINS = 120;
const double EMAX = 12;
const int OLD_GEOM = 1;
const double Z_OFFSET = -17.5; // mm offset of the AmBe source

// Turn points into a graph
void FillGraphs(TGraph *g, int n, double *x, double *y, double *z, int col) {
  new (g) TGraph();
  g->SetPoint(0,y[0],z[0]+Z_OFFSET);
  for (int i=0; i<n; i++) {
    g->SetPoint(i+1,x[i],z[i]+Z_OFFSET);
    g->SetPoint(n+i+1,y[n-i-1],z[n-i-1]+Z_OFFSET);
  }
  g->SetLineWidth(2);
  g->SetLineColor(col);
};

// Main program
int main(int argc, char** argv) {
  
  // Input/Output files
  string inroot;
  if (TEST) inroot = "test7/AmBe_*.root";
  else      inroot = "nickel/output/AmBe_*.root";
  std::ofstream out("testMC.out");
 
  // Plotting options 
  gStyle->SetOptStat(0);
  gStyle->SetTitleBorderSize(0);
  gStyle->SetHistLineWidth(2);
  
  // Initialise stuff
  TFile *outroot=NULL;
  TH1D *hGam=NULL, *hNeu=NULL, *hPro=NULL, *hOth=NULL;
  TH1D *hGamAV=NULL, *hGamAVn=NULL, *hGamFrac=NULL, *hGamFracn=NULL;
  TH2D *hGamTime=NULL, *hnCapPos=NULL;
  
  string outfile = "testMC.root";
  if (!TEST) outfile = "nickel/results/histograms.root";
  std::ifstream chk(outfile.c_str());
  bool processed = chk.good();
  if (processed) {  // if already processed, read output file
    outroot = new TFile(outfile.c_str(),"READ");
    hGam = (TH1D*)outroot->Get("hGam");
    hNeu = (TH1D*)outroot->Get("hNeu");
    hPro = (TH1D*)outroot->Get("hPro");
    hOth = (TH1D*)outroot->Get("hOth");
    hGamAV = (TH1D*)outroot->Get("hGamAV");
    hGamAVn = (TH1D*)outroot->Get("hGamAVn");
    hGamFrac = (TH1D*)outroot->Get("hGamFrac");
    hGamFracn = (TH1D*)outroot->Get("hGamFracn");
    hGamTime = (TH2D*)outroot->Get("hGamTime");
    hnCapPos = (TH2D*)outroot->Get("hnCapPos");
  } else {          // if not processed, read input file & generate output
    outroot = new TFile(outfile.c_str(),"NEW");
    
    // Initialise RAT
    RAT::DU::DSReader dsreader(inroot);
    const RAT::DU::PMTInfo& pmtinfo = RAT::DU::Utility::Get()->GetPMTInfo();
    const int NPMTS = pmtinfo.GetCount();

    // Initialise histograms
    hGam = new TH1D("hGam","Gamma",NBINS,0,EMAX);
    hNeu = new TH1D("hNeu","Neutron",NBINS,0,EMAX);
    hPro = new TH1D("hPro","Proton",NBINS,0,EMAX);
    hOth = new TH1D("hOth","Other",NBINS,0,EMAX);
    hGamAV = new TH1D("hGamAV","Gamma (AV)",NBINS,0,EMAX);
    hGamAVn = new TH1D("hGamAVn","Gamma (AV) from nCapture",NBINS,0,EMAX);
    hGamFrac = new TH1D("hGamFrac","Gamma (AV frac.)",NBINS,0,EMAX);
    hGamFracn = new TH1D("hGamFracn","Gamma (AV frac.) from nCapture",NBINS,0,EMAX);
    hGamTime = new TH2D("hGamTime","Gamma (AV) timing",NBINS,0,EMAX,NBINS,0,4e6);
    BinLog(hGamTime->GetYaxis(), 0.04); // log axis for timing (0.04 ns - 4 ms)
    hnCapPos = new TH2D("hnCapPos","Neutron capture position",1e3,0,1e3,2e3,-1e3,1e3);
    
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
          continue;
        } else if (pdg == 2212) {  // proton
          hPro->Fill(eTrk);
        } else if (pdg == 2112) {  // neutron
          hNeu->Fill(eTrk);
        } else if (pdg ==   22) {  // gamma
          hGam->Fill(eTrk);
        } else {                   // other
          hOth->Fill(eTrk);
          if (pdg == -11) {        // positron (don't care, already omitted electrons)
            continue;
          }
          if (fabs(pdg) == 12) {   // (anti-)nu_e (really don't care)
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
        if (pdg != 22) continue;    // gamma
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
          if (stp.GetEndVolume() != "inner_av") continue; // step ends in inner AV
          hGamAV->Fill(eStep); // gamma energy when entering AV (visible energy)
          hGamTime->Fill(eStep,stp.GetGlobalTime());
          //if (trk.GetFirstMCTrackStep().GetStartVolume() != "AmBeAbsorber") continue; // particle produced on Nickel
          if (trk.GetFirstMCTrackStep().GetProcess() != "nCapture") break; // gamma from neutron capture
          hGamAVn->Fill(eStep);
          TVector3 pos = trk.GetFirstMCTrackStep().GetPosition();
          hnCapPos->Fill(pos.Perp(),pos.Z());
          break;
        } // steps
      } // tracks
    } // entries
    out << "-----" << endl;
    out << "Found " << nevents << " MC events with " << ntracks << " tracks, consisting of:" << endl;
    out << "- Protons:  " << setw(6) << (int)hPro->GetEntries() << " (" << hPro->GetEntries()/nevents << " / event)" << endl;
    out << "- Neutrons: " << setw(6) << (int)hNeu->GetEntries() << " (" << hNeu->GetEntries()/nevents << " / event)" << endl;
    out << "- Gammas:   " << setw(6) << (int)hGam->GetEntries() << " (" << hGam->GetEntries()/nevents << " / event)" << endl;
    out << "- Other:    " << setw(6) << (int)hOth->GetEntries() << " (" << hOth->GetEntries()/nevents << " / event)" << endl;
    out.close();
    
    // Save output file
    outroot->Write();
  }
  
  // Plotting section
  string imgname = "testMC";
  if (!TEST) imgname = "nickel/results/AmBe_plot";
  double maxval = std::max(hPro->GetMaximum(),hGam->GetMaximum());
  TCanvas *c = new TCanvas("c","",800,600);
  c->SetLogy();
  c->SetGrid();
  TH1F *pad = c->DrawFrame(0,0.5,EMAX,1.5*maxval);
  pad->SetTitle("AmBe simulation with Nickel absorber;Kinetic energy [MeV];MC tracks");
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
  
  c = new TCanvas("c","",800,600);
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
  
  c = new TCanvas("c","",800,600);
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
  
  c = new TCanvas("c","",800,600);
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
  
  // Draw geometry
  c = new TCanvas("c","",800,600);
  c->Divide(2,1);
  c->cd(1)->SetGrid();
  c->cd(1)->DrawFrame(0,-400,600,400,"Position of neutron capture;r [mm];z [mm]");
  // Neutron captures
  hnCapPos->SetMarkerStyle(6);
  hnCapPos->Draw("scat same");
  // Source geometry
  //g0->Draw("L same");
  g1->Draw("L same");
  g2->Draw("L same");
  g3->Draw("L same");
  g4->Draw("L same");
  g5->Draw("L same");
  c->cd(2)->SetGrid();
  c->cd(2)->DrawFrame(0,-30,25,170,"Position of neutron capture;r [mm];z [mm]");
  // Neutron captures
  hnCapPos->SetMarkerStyle(6);
  hnCapPos->Draw("scat same");
  // Source geometry
  //g0->Draw("L same");
  g1->Draw("L same");
  g2->Draw("L same");
  g3->Draw("L same");
  g4->Draw("L same");
  g5->Draw("L same");
  c->Print((imgname+"_5.png").c_str());
  c->Close();
  
  if (hPro) delete hPro;
  if (hNeu) delete hNeu;
  if (hGam) delete hGam;
  if (hOth) delete hOth;
  if (hGamAV) delete hGamAV;
  if (hGamAVn) delete hGamAVn;
  if (leg) delete leg;
  if (c) delete c;
  return 0;
  
}
