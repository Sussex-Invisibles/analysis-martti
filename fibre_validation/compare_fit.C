// ---------------------------------------------------------
// Goal:          Check variation of fit results
// Author:        Martti Nirkko, 21/05/2019
// Compile & run: clear && g++ -g -o compare_fit.exe compare_fit.C `root-config --cflags --libs` -I$RATROOT/include/libpq -I$RATROOT/include -L$RATROOT/lib -lRATEvent_Linux && ./compare_fit.exe
// ---------------------------------------------------------

// Helper functions (includes everything else)
#include "../include/HelperFunc.C"

// Main program
int main(int argc, char** argv) {
  
  // Global parameters
  const int NPCAS = 6;
  const int NFIBRES = 95;
  
  // Files to load
  string fpath = Form("%s/pCloud/PostDoc/Calibration/TELLIE/Fibre_validation/",getenv("HOME"));
  string fname = "TELLIE_FITRESULTS_";
  string set[NPCAS] = {"2017-11","2018-03","2018-06","2018-09","2018-12","2019-03"};
  
  // Graph showing all PMTs
  RAT::DU::DSReader dsreader("../../data/Analysis_r0000102315_s000_p000.root");
  const RAT::DU::PMTCalStatus& pmtStatus = RAT::DU::Utility::Get()->GetPMTCalStatus();
  const RAT::DU::PMTInfo& pmtinfo = RAT::DU::Utility::Get()->GetPMTInfo();
  int NPMTS = pmtinfo.GetCount();
  cout << NPMTS << " PMTs found." << endl;
  TGraph *psup = new TGraph();
  TVector3 pmtpos(0,0,0);
  TVector2 icospos;
  int npmt=0;
  for(int id=0; id<NPMTS; id++) {
    pmtpos = pmtinfo.GetPosition(id);
    if (pmtpos.Mag()==0) continue; // not a valid PMT
    if (pmtinfo.GetType(id) != 1) continue; // not a normal PMT
    int face; // side of icosahedron containing PMT
    icospos = func::IcosProject(pmtpos,face);
    psup->SetPoint(npmt,icospos.X(),icospos.Y());
    npmt++;
  }
  
  // Fibre positions and directions
  TVector3 fibrepos[NFIBRES];
  TVector3 fibredir[NFIBRES];
  TVector3 lightpos[NFIBRES];
  string fibre_table = "Fibre_Positions_DocDB1730.csv";
  ifstream tab(fibre_table.c_str());
  if (!tab) { cerr<<"Failed to open "<<fibre_table<<endl; exit(1); }
  string line, ff, fn;
  double fx,fy,fz,fu,fv,fw;
  getline(tab,line);      // header
  while (tab.good()) {
    tab >> ff >> fx >> fy >> fz >> fu >> fv >> fw >> fn;
    if (!tab.good()) break;
    if (ff.substr(5,5)=="B") continue;
    int fib = std::atoi(ff.substr(2,4).c_str());
    if (fib==101) fib=0; // set to zero
    if (fib>NFIBRES) continue;
    fibrepos[fib].SetXYZ(10*fx,10*fy,10*fz); // position of fibre [mm]
    fibredir[fib].SetXYZ(fu,fv,fw); // direction of fibre
    lightpos[fib] = fibrepos[fib] + 2*fibrepos[fib].Mag()*fibredir[fib]; // projected light spot centre (bad approximation)
  }
  tab.close();
  
  // Vectors and graphs for direct and reflected light
  TVector3 dirfit(0,0,0);
  TVector3 reffit(0,0,0);
  TGraph *gfib[NFIBRES] = {NULL};
  TGraph *gexp[NFIBRES] = {NULL};
  TGraph *gdir[NFIBRES] = {NULL};
  TGraph *gref[NFIBRES] = {NULL};
  
  // Loop over all fibres in list
  for (int s=0; s<NPCAS; s++) {
  
    // Load file
    string input = fpath+fname+set[s]+".txt";
    ifstream in(input.c_str());
    if (!in) { cerr<<"Failed to open "<<input<<endl; exit(1); }
    string line, fibre;
    int run, face, fib;
    float x,y,z,u,v,w;
    for (int hdr=0; hdr<2; hdr++) {
      getline(in,line);      // header
    }
    
    // Loop over lines in file
    int nfiles = 0;
    while (true) {
    
      // Read fit positions
      in >> run >> fibre >> x >> y >> z >> u >> v >> w;
      if (!in.good()) break;
      dirfit.SetXYZ(x,y,z);
      reffit.SetXYZ(u,v,w);
      
      // Get fibre number
      fib = std::atoi(fibre.substr(2,4).c_str());
      if (fib==101) fib=0; // set to zero
      if (s==0) {
        gdir[fib] = new TGraph();
        gref[fib] = new TGraph();
        gfib[fib] = new TGraph();
        icospos = func::IcosProject(fibrepos[fib],face);
        gfib[fib]->SetPoint(s,icospos.X(),icospos.Y());
        gexp[fib] = new TGraph();
        icospos = func::IcosProject(lightpos[fib],face);
        gexp[fib]->SetPoint(s,icospos.X(),icospos.Y());
      }
      
      // Project fit positions
      icospos = func::IcosProject(dirfit,face);
      gdir[fib]->SetPoint(s,icospos.X(),icospos.Y());
      icospos = func::IcosProject(reffit,face);
      gref[fib]->SetPoint(s,icospos.X(),icospos.Y());
      //cout << "Run " << run << ", fibre " << fibre << " has fit results:" << endl;
      //cout << "- Direct light fit    (x,y,z) [mm] = " << printVector(*dirfit) << endl;
      //cout << "- Reflected light fit (x,y,z) [mm] = " << printVector(*reffit) << endl;
      
      // Obtain fitted direction as follows
      /*
      const double ENERGY   = RAT::util::WavelengthToEnergy(5e-4); // photon energy [MeV]
      const double LOCALITY = 10.0;   // accepted tolerance [mm] for LightPathCalculator
      RAT::DU::LightPathCalculator lpc = RAT::DU::Utility::Get()->GetLightPathCalculator();
      //lpc.BeginOfRun(); // TODO - find out if this is required
      lpc.CalcByPosition(fibrePos, fitPos, ENERGY, LOCALITY);
      TVector3 fitDir = lpc.GetInitialLightVec(); // fitted direction of fibre
      double angDiff = fibreDir.Angle(fitDir)*180./pi; // difference to RATDB (should be small)
      */
      
      nfiles++;
    }
    cout << "PCA dataset " << set[s] << " had " << nfiles << " runs." << endl;
  }
  
  // Plotting section
  TCanvas c("c","",1200,1000);
  c.Divide(1,2);
  c.cd(1);
  psup->SetMarkerStyle(6);
  psup->SetMarkerColor(16);
  psup->Draw("P");
  for (int r=0;r<NFIBRES;r++) {
    // Direct light
    gexp[r]->SetMarkerStyle(4);
    gexp[r]->SetMarkerColor(1);
    gexp[r]->Draw("p same");
    gdir[r]->SetMarkerStyle(5);
    gdir[r]->SetMarkerColor(r+1);
    gdir[r]->Draw("p same");
  }
  c.cd(2);
  psup->Draw("P");
  for (int r=0;r<NFIBRES;r++) {
    // Reflected light
    gfib[r]->SetMarkerStyle(4);
    gfib[r]->SetMarkerColor(1);
    gfib[r]->Draw("p same");
    gref[r]->SetMarkerStyle(5);
    gref[r]->SetMarkerColor(r+1);
    gref[r]->Draw("p same");
  }
  c.Print("compare_fit.png");
  c.Print("compare_fit.pdf");
  c.Close();
  
  return 0;
}
