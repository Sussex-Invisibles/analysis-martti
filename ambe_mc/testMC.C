// ---------------------------------------------------------
// Goal:          Evaluate gamma peaks in simulated AmBe events
// Author:        Martti Nirkko, 10/01/2018
// Compile & run: clear && g++ -g -o testMC.exe testMC.C `root-config --cflags --libs` -I$RATROOT/include/libpq -I$RATROOT/include -L$RATROOT/lib -lRATEvent_Linux && ./testMC.exe
// ---------------------------------------------------------

// Helper functions (includes everything else)
#include "../analysis/HelperFunc.C"
using std::setw;
using std::setprecision;
using std::fixed;

// Global parameters
const int TEST = 0;
const int VERBOSE = 0;
const int isMC = 1;
const int NBINS = 120;
const double EMAX = 12;

// Main program
int main(int argc, char** argv) {
  
  // Input/Output files
  string file;
  if (TEST) file = "test2/AmBeTest.root";
  else      file = "cluster/output/AmBe_Nickel_*.root";
  ofstream out("testMC.out");
  
  // Initialise RAT
  RAT::DU::DSReader dsreader(file);
  const RAT::DU::PMTInfo& pmtinfo = RAT::DU::Utility::Get()->GetPMTInfo();
  const int NPMTS = pmtinfo.GetCount();
 
  // Plotting options 
  gStyle->SetOptStat(0);
  gStyle->SetTitleBorderSize(0);
  gStyle->SetHistLineWidth(2);  

  // Initialise histograms
  TH1D *hGam = new TH1D("hGam","Gamma",NBINS,0,EMAX);
  TH1D *hNeu = new TH1D("hNeu","Neutron",NBINS,0,EMAX);
  TH1D *hPro = new TH1D("hPro","Proton",NBINS,0,EMAX);
  TH1D *hOth = new TH1D("hOth","Other",NBINS,0,EMAX);
  TH1D *hGamAV = new TH1D("hGamAV","Gamma (AV)",NBINS,0,EMAX);
  TH1D *hGamAVn = new TH1D("hGamAVn","Gamma (AV) from nCapture",NBINS,0,EMAX);
  TH1D *hGamFrac = new TH1D("hGamFrac","Gamma (AV frac.)",NBINS,0,EMAX);
  TH1D *hGamFracn = new TH1D("hGamFracn","Gamma (AV frac.) from nCapture",NBINS,0,EMAX);
  TH2D *hGamTime = new TH2D("hGamTime","Gamma (AV) timing",NBINS,0,EMAX,NBINS,0,4e6);
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
  
  // Plotting section
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
  c->Print("testMC.png");
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
  c->Print("testMC_2.png");
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
  c->Print("testMC_3.png");
  c->Close();
  
  c = new TCanvas("c","",800,600);
  c->SetLogy();
  c->SetLogz();
  c->SetGrid();
  hGamTime->SetTitle("#gamma-rays entering inner AV;Photon energy [MeV];Global time [ns]");
  hGamTime->Draw("colz");
  c->Print("testMC_4.png");
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
