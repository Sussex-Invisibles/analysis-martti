// ---------------------------------------------------------
// Goal:          Evaluate gamma peaks in simulated dinucleon decay events
// Author:        Martti Nirkko, 31/08/2018
// Compile & run: clear && g++ -g -o Dinucleon.exe Dinucleon.C `root-config --cflags --libs` -I$RATROOT/include/libpq -I$RATROOT/include -L$RATROOT/lib -lRATEvent_Linux -lGeom && ./Dinucleon.exe
// ---------------------------------------------------------

// Helper functions (includes everything else)
#include "../HelperFunc.C"

// Global parameters
const int TEST = 0;                    // test flag
const int VERBOSE = 0;                 // verbosity flag
const int isMC = 1;                    // MC flag
const int NBINS = 100;                 // number of bins for most histograms
const double EMAX = 10;                // maximal energy for histograms

// Main program
int main(int argc, char** argv) {

  const int NDEC=3;
  const string DECAYS[NDEC] = {"PP","PN","NN"};

  for (int dec=0; dec<NDEC; dec++) {
  
    if (TEST && dec>0) continue;
    string DECAY = DECAYS[dec];
    
    // Input files
    string fpath = DECAY;
    string inroot = fpath + "/output";
    if (!TEST) inroot = inroot + "/" + DECAY + "_decay_*.root";
    else       inroot = inroot + "/" + DECAY + "_decay_1.root";
    
    // Plotting options 
    gStyle->SetOptStat(0);
    gStyle->SetTitleBorderSize(0);
    gStyle->SetHistLineWidth(2);
    
    // Initialise stuff
    TFile *outfile=NULL;
    TH1D *hEtrue=NULL, *hEreco=NULL;
    TH2D *hEnergy=NULL, *hDeltaE=NULL;
    
    // Output file
    string outroot = "histograms.root";
    if (!TEST) outroot = fpath + "/results/" + outroot;
    std::ifstream chk(outroot.c_str());
    bool processed = chk.good();
    if (processed) {  // if already processed, read output file
      outfile = new TFile(outroot.c_str(),"READ");
      hEtrue    = (TH1D*)outfile->Get("hEtrue");
      hEreco    = (TH1D*)outfile->Get("hEreco");
      hEnergy   = (TH2D*)outfile->Get("hEnergy");
      hDeltaE   = (TH2D*)outfile->Get("hDeltaE");
    } else {          // if not processed, read input file & generate output
      outfile = new TFile(outroot.c_str(),"NEW");
      
      // Aggro offline mode
      RAT::DB *db = RAT::DB::Get();
      db->SetAirplaneModeStatus(true);
      db->SetDefaultPlaneLockStatus(false);
      
      // Initialise RAT
      RAT::DU::DSReader dsreader(inroot);
      const RAT::DU::PMTInfo& pmtinfo = RAT::DU::Utility::Get()->GetPMTInfo();
      const int NPMTS = pmtinfo.GetCount();
      
      // Initialise histograms
      hEtrue    = new TH1D("hEtrue","MC true energy",NBINS,0,2*EMAX);
      hEreco    = new TH1D("hEreco","MC reconstructed energy",NBINS,0,2*EMAX);
      hEnergy   = new TH2D("hEnergy","Reconstructed vs true energy",NBINS,0,2*EMAX,NBINS,0,2*EMAX);
      hDeltaE   = new TH2D("hDeltaE","Difference vs true energy",NBINS,0,2*EMAX,NBINS,-2,2);
      
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
        
        double mcEnergy = 0;
        for (int ipart=0; ipart<mc.GetMCParticleCount(); ipart++) {
          const RAT::DS::MCParticle& part = mc.GetMCParticle(ipart);
          if (part.GetPDGCode() != 22) continue; // gamma
          if (ipart>0) { 
            cout << "WARNING: more than one gamma in event! Ignoring..." << endl;
            continue;
          }
          mcEnergy += part.GetKineticEnergy();
          //cout << "MC gamma " << ipart << " has true energy " << mcEnergy << " MeV." << endl;
        }

        // Loop over events
        double fitEnergy = 0;
        for (int iEvent=0; iEvent<ds.GetEVCount(); iEvent++) {
          RAT::DS::EV evt = ds.GetEV(iEvent);
          
          int nhits_cleaned = evt.GetNhitsCleaned();
          if (nhits_cleaned == 0) continue;
          
          string lFit = evt.GetDefaultFitName();
          try {
            RAT::DS::FitVertex fitVertex = evt.GetDefaultFitVertex();
            if (fitVertex.ContainsEnergy() && fitVertex.ValidEnergy())
              fitEnergy += fitVertex.GetEnergy();
          } catch (RAT::DS::FitResult::NoVertexError& e) {
            cout << lFit << " has not reconstructed a vertex. Continuing..." << endl;
            continue;
          } catch (std::runtime_error& e) {
            cout << lFit << " failed for event " << iEvent << ". Continuing..." << endl;
            continue;
          }
        }
        
        hEtrue->Fill(mcEnergy);
        //cout << "True E = " << mcEnergy << "MeV, fit E = " << fitEnergy << " MeV." << endl;
        if(fitEnergy==0) continue;
        hEreco->Fill(fitEnergy);
        hEnergy->Fill(mcEnergy,fitEnergy);
        hDeltaE->Fill(mcEnergy,(fitEnergy-mcEnergy)/mcEnergy);
        
      } // entries
      
      // Save output file
      outfile->Write();
    }

    // Plotting section
    string imgname = "Dinucleon";
    if (!TEST) imgname = fpath + "/results/" + imgname;
    
    // Plotting options
    int col = dec+2;
    hEtrue->SetFillColorAlpha(col,0.25);
    hEtrue->SetLineColorAlpha(col,0);
    hEtrue->SetMarkerColorAlpha(col,0);
    hEreco->SetLineColor(col);
    hEreco->SetLineWidth(2);
    TCanvas *c = NULL;
    TH1F *pad = NULL;
    string tstr;
    
    // Energy of various particles 
    double maxval = std::max(hEtrue->GetMaximum(),hEreco->GetMaximum());
    c = new TCanvas("c","",800,600);
    c->SetLogy();
    c->SetGrid();
    pad = c->DrawFrame(0,0.5,2*EMAX,1.5*maxval);
    tstr = Form("MC truth vs reconstructed energy (%s decay);E_{#gamma} [MeV];Events",DECAY.c_str());
    pad->SetTitle(tstr.c_str());
    pad->Draw("");
    hEtrue->Draw("same");
    hEreco->Draw("same");
    TLegend *leg = new TLegend(0.6,0.8,0.89,0.89);
    leg->AddEntry(hEtrue);
    leg->AddEntry(hEreco);
    leg->Draw();
    c->Print((imgname+"_1.png").c_str());
    c->Print((imgname+"_1.pdf").c_str());
    c->Close();
    
    c = new TCanvas("c","",800,600);
    c->SetGrid();
    pad = c->DrawFrame(0,0,2*EMAX,2*EMAX);
    tstr = Form("MC truth vs reconstructed energy (%s decay);E_true [MeV];E_reco [MeV]",DECAY.c_str());
    pad->SetTitle(tstr.c_str());
    pad->Draw("");
    hEnergy->Draw("colz same");
    TLine *line = new TLine(0,0,2*EMAX,2*EMAX);
    line->SetLineColor(1);
    line->SetLineWidth(2);
    line->Draw("same");
    c->Print((imgname+"_2.png").c_str());
    c->Print((imgname+"_2.pdf").c_str());
    c->Close();
    
    c = new TCanvas("c","",800,600);
    c->SetGrid();
    pad = c->DrawFrame(0,-2,2*EMAX,2);
    tstr = Form("Energy resolution (%s decay);E_true [MeV];(E_reco #minus E_true) / E_true [MeV]",DECAY.c_str());
    pad->SetTitle(tstr.c_str());
    pad->Draw("");
    hDeltaE->Draw("colz same");
    line->DrawLine(0,0,2*EMAX,0);
    c->Print((imgname+"_3.png").c_str());
    c->Print((imgname+"_3.pdf").c_str());
    c->Close();
    
    // Free memory
    if (leg) delete leg;
    if (c) delete c;
    
    // Close files and exit
    outfile->Close();

  }

  return 0;
  
}
