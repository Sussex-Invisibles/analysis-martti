// ---------------------------------------------------------
// Goal:            Test LightPathCalculator for different settings
// Author:          Martti Nirkko, 04/07/2019
// Compile & run:   clear && g++ -o lightpath.exe lightpath.C `root-config --cflags --libs` -I$RATROOT/include/libpq -I$RATROOT/include -L$RATROOT/lib -lRATEvent_Linux && ./lightpath.exe
// ---------------------------------------------------------

// Helper functions (includes everything else)
#include "../include/HelperFunc.C"
#include "../include/CommonStyle.H"

// Global parameters
const int NLOC = 5;
const int NBINS = 1000;
const int locality[NLOC] = {0,5,10,15,20};
const int col[NLOC] = {12,50,75,63,92};
string lightPathType[8] = { "Scint/InnerAV->AV->Water",
                            "AV->Water",
                            "AV->Scint/InnerAV->AV->Water",
                            "Water->AV->Scint/InnerAV->AV->Water",
                            "Water->AV->Water",
                            "Water",
                            "Water->Reflection->Water",
                            "Light Path Uninitialised"
                          };
  
int main() {

    // Select a TELLIE file/fibre
    string fname = "../../data/Calibration_r0000204420_s000_p000.root";
    string fibre = "FT008A";
    
    // Initialise RAT
    RAT::DU::DSReader dsreader(fname);
    RAT::DU::Utility::Get()->BeginOfRun();
    RAT::DU::GroupVelocity gv = RAT::DU::Utility::Get()->GetGroupVelocity();
    RAT::DU::LightPathCalculator lpc = RAT::DU::Utility::Get()->GetLightPathCalculator();
    lpc.SetELLIEEvent(true); // event originates outside AV (see PR #2621)
    const RAT::DU::PMTInfo& pmtinfo = RAT::DU::Utility::Get()->GetPMTInfo();
    const int NPMTS = pmtinfo.GetCount();
    cout << NPMTS << " PMTs found" << endl;
    
    // Get fibre info (from RATDB) 
    RAT::DB *db = RAT::DB::Get();
    //db->LoadDefaults();	  // Already done when calling DU::Utility::Get()
    RAT::DBLinkPtr entry = db->GetLink("FIBRE",fibre);
    TVector3 fibrePos(entry->GetD("x"), entry->GetD("y"), entry->GetD("z")); // position of fibre [mm]
    TVector3 fibreDir(entry->GetD("u"), entry->GetD("v"), entry->GetD("w")); // direction of fibre
    //TVector3 lightPos = fibrePos + 2*fibrePos.Mag()*fibreDir;                // projected light spot centre
    cout << "RATDB: fibre " << fibre << ", pos " << printVector(fibrePos) << ", dir " << printVector(fibreDir) << endl;
    
    // Declare variables
    int pmtcount[NLOC] = {0};
    int goodLP[NLOC] = {0};
    
    TH2D *hang[NLOC] = {NULL};
    TH2D *htyp[NLOC] = {NULL};
    TH2D *hdis[NLOC] = {NULL};
    TH2D *htof[NLOC] = {NULL};
    TH2D *hbin[NLOC] = {NULL};
    TH2D *hbkt[NLOC] = {NULL};
    
    TH2D *hinner[NLOC] = {NULL};
    TH2D *hinAV[NLOC] = {NULL};
    TH2D *hwater[NLOC] = {NULL};
    
    TH1D *hdisSL[NLOC] = {NULL};
    TH1D *hdisLP[NLOC] = {NULL};
    TH1D *hangAV[NLOC] = {NULL};
    TH1D *hangLP[NLOC] = {NULL};
    TH1D *hcosAV[NLOC] = {NULL};
    TH1D *hcosLP[NLOC] = {NULL};
    
    string hname = "";
    
    // Loop over different values for locality
    for (int l=0; l<NLOC; l++) {
    
      // Initialise histograms
      hname = Form("hang_%d",l);
      hang[l] = new TH2D(hname.c_str(),"",NBINS,0,100,NBINS,0,100);
      hname = Form("htyp_%d",l);
      htyp[l] = new TH2D(hname.c_str(),"",NBINS,0,100,8,0,8);
      hname = Form("hdis_%d",l);
      hdis[l] = new TH2D(hname.c_str(),"",NBINS,0,100,NBINS,0,20);
      hname = Form("htof_%d",l);
      htof[l] = new TH2D(hname.c_str(),"",NBINS,0,100,NBINS,0,100);
      hname = Form("hbin_%d",l);
      hbin[l] = new TH2D(hname.c_str(),"",NBINS,0,100,NBINS,0,100);
      hname = Form("hbkt_%d",l);
      hbkt[l] = new TH2D(hname.c_str(),"",NBINS,0,100,NBINS,0,1);
      
      hname = Form("hinner_%d",l);
      hinner[l] = new TH2D(hname.c_str(),"",NBINS,0,100,NBINS,0,20);
      hname = Form("hinAV_%d",l);
      hinAV[l] = new TH2D(hname.c_str(),"",NBINS,0,100,NBINS,0,20);
      hname = Form("hwater_%d",l);
      hwater[l] = new TH2D(hname.c_str(),"",NBINS,0,100,NBINS,0,20);
      
      hname = Form("hdisSL_%d",l);
      hdisSL[l] = new TH1D(hname.c_str(),"",100,0,20);
      hname = Form("hdisLP_%d",l);
      hdisLP[l] = new TH1D(hname.c_str(),"",100,0,20);
      hname = Form("hangAV_%d",l);
      hangAV[l] = new TH1D(hname.c_str(),"",100,0,180);
      hname = Form("hangLP_%d",l);
      hangLP[l] = new TH1D(hname.c_str(),"",100,0,180);
      hname = Form("hcosAV_%d",l);
      hcosAV[l] = new TH1D(hname.c_str(),"",100,-1,1);
      hname = Form("hcosLP_%d",l);
      hcosLP[l] = new TH1D(hname.c_str(),"",100,-1,1);
      
      // Loop over PMTs
      for (int it=0; it<NPMTS; it++) {
      
        // Use normal or HQE PMTs only
        if (pmtinfo.GetType(it) != 1 && pmtinfo.GetType(it) != 7) continue;
        pmtcount[l]++;
        
        // Get PMT information
        TVector3 pmtPos = pmtinfo.GetPosition(it);       // position [mm]
        TVector3 pmtDir = pmtinfo.GetDirection(it);      // direction
        double pmtAng = fibrePos.Angle(pmtPos)*180./pi;  // angle w.r.t. fibre position
        
        // Get straight line for reference
        TVector3 lineDir = (pmtPos-fibrePos).Unit();      // straight line assumption
        double thetaLine = fibreDir.Angle(lineDir)*180./pi; // angle w.r.t. fibre direction
        
        // Calculate light path
        lpc.CalcByPosition(fibrePos,pmtPos,ENERGY,locality[l]);
        if (!lpc.GetPathValid()) continue;
        goodLP[l]++;
        
        // Get light path type
        string type = lpc.GetLightPathType();
        int k;
        for (k=0; k<8; k++) {
          if (lightPathType[k]==type) break;
        }
        
        // Get light travel time (as done in BiPoLikelihoodDiff.cc)
        double distInInnerAV = lpc.GetDistInInnerAV();
        double distInAV = lpc.GetDistInAV();
        double distInWater = lpc.GetDistInWater();
        double transitTime = gv.CalcByDistance(distInInnerAV, distInAV, distInWater, ENERGY);
        double totalDist = distInInnerAV + distInAV + distInWater;
        
        // Get light bucket time (as done in DQLaserBallProc.cc)
        TVector3 endDir = lpc.GetIncidentVecOnPMT();        // end direction at PMT
        double thetaAtPMT = endDir.Angle(pmtDir)*180./pi;   // incident angle with bucket face
        double bucketTime = gv.PMTBucketTime(thetaAtPMT);   // DocDB 3138
        
        // Get start direction at fibre
        TVector3 startDir = lpc.GetInitialLightVec();
        double thetaLPC = fibreDir.Angle(startDir)*180./pi; // angle w.r.t. fibre
        
        /*
        if (fabs(thetaLPC-thetaLine)>10) {
          cout << "PMT #" << it << " has angles " << thetaLine << " - " << thetaLPC << " = " << thetaLPC-thetaLine << " degrees!" << endl;
        }
        */
        
        // Fill histograms
        hang[l]->Fill(thetaLine,thetaLPC);
        hdis[l]->Fill(thetaLine,totalDist/1e3);
        htof[l]->Fill(thetaLine,transitTime);
        hbin[l]->Fill(thetaLine,thetaAtPMT);
        hbkt[l]->Fill(thetaLine,bucketTime);
        htyp[l]->Fill(thetaLine,k);
        hinner[l]->Fill(thetaLine,distInInnerAV/1e3);
        hinAV[l]->Fill(thetaLine,distInAV/1e3);
        hwater[l]->Fill(thetaLine,distInWater/1e3);
        hdisSL[l]->Fill((pmtPos-fibrePos).Mag()/1e3);
        hdisLP[l]->Fill(totalDist/1e3);
        hangAV[l]->Fill(pmtAng);
        hangLP[l]->Fill(thetaLine);
        hcosAV[l]->Fill(cos(fibrePos.Angle(pmtPos)));
        hcosLP[l]->Fill(cos(fibreDir.Angle(lineDir)));
      }
      
      // Output LPC success rate
      cout << "Counted " << pmtcount[l] << " PMTs";
      cout << " with " << goodLP[l] << " good LightPaths";
      cout << " (" << 100.*goodLP[l]/pmtcount[l] << "%)";
      cout << " for locality " << locality[l] << endl;
    }
    
    // Plotting section
    CommonStyle();
    gStyle->SetTitleOffset(1.2, "y");
    TCanvas *c = NULL;
    
    TText lpType(0,0,"");
    lpType.SetTextSize(0.048);
    lpType.SetTextColor(kRed+2);
    
    c = new TCanvas("c","",1600,1200);
    c->Divide(3,2);
    // Angle comparison
    c->cd(1)->SetGrid();
    c->cd(1)->DrawFrame(0,0,100,100,"Angle comparison;(StraightLine) PMT angle from fibre [deg];LightPath initial direction from fibre [deg]");
    for (int l=NLOC-1; l>=0; l--) {
      hang[l]->SetMarkerStyle(6);
      hang[l]->SetMarkerColor(col[l]);
      hang[l]->SetLineColor(col[l]);
      hang[l]->Draw("scat same");
    }
    // Legend on first plot only
    TLegend leg(0.16,0.62,0.56,0.95);
    for (int l=0; l<NLOC; l++) leg.AddEntry(hang[l],Form("Locality %d mm",locality[l]),"lp");
    leg.Draw();
    // Angle comparison
    c->cd(2)->SetGrid();
    c->cd(2)->DrawFrame(0,0,100,100,"Angle comparison;(StraightLine) PMT angle from fibre [deg];Incoming angle at PMT [deg]");
    for (int l=NLOC-1; l>=0; l--) {
      hbin[l]->SetMarkerStyle(6);
      hbin[l]->SetMarkerColor(col[l]);
      hbin[l]->SetLineColor(col[l]);
      hbin[l]->Draw("scat same");
    }
    // Distance comparison
    c->cd(3)->SetGrid();
    c->cd(3)->DrawFrame(0,0,100,20,"Distance comparison;(StraightLine) PMT angle from fibre [deg];Total light path length [m]");
    for (int l=NLOC-1; l>=0; l--) {
      hdis[l]->SetMarkerStyle(6);
      hdis[l]->SetMarkerColor(col[l]);
      hdis[l]->SetLineColor(col[l]);
      hdis[l]->Draw("scat same");
    }
    // Light path type
    c->cd(4)->SetGrid();
    c->cd(4)->DrawFrame(0,0,100,8,"Type comparison;(StraightLine) PMT angle from fibre [deg];Light path type [ ]");
    for (int l=NLOC-1; l>=0; l--) {
      htyp[l]->SetMarkerStyle(6);
      htyp[l]->SetMarkerColor(col[l]);
      htyp[l]->SetLineColor(col[l]);
      htyp[l]->Draw("scat same");
    }
    for (int m=0; m<8; m++) lpType.DrawText(5,m+0.35,lightPathType[m].c_str());
    // Bucket time comparison
    c->cd(5)->SetGrid();
    c->cd(5)->DrawFrame(0,0.4,100,0.75,"Bucket time comparison;(StraightLine) PMT angle from fibre [deg];PMT bucket time [ns]");
    for (int l=NLOC-1; l>=0; l--) {
      hbkt[l]->SetMarkerStyle(6);
      hbkt[l]->SetMarkerColor(col[l]);
      hbkt[l]->SetLineColor(col[l]);
      hbkt[l]->Draw("scat same");
    }
    // Transit time comparison
    c->cd(6)->SetGrid();
    c->cd(6)->DrawFrame(0,0,100,100,"Transit time comparison;(StraightLine) PMT angle from fibre [deg];Transit time [ns]");
    for (int l=NLOC-1; l>=0; l--) {
      htof[l]->SetMarkerStyle(6);
      htof[l]->SetMarkerColor(col[l]);
      htof[l]->SetLineColor(col[l]);
      htof[l]->Draw("scat same");
    }
    // Save and close
    c->Print("lightpath.png");
    //c->Print("lightpath.pdf");
    c->Close();
    
    // Zoomed in version of same thing
    c = new TCanvas("c","",1600,1200);
    c->Divide(3,2);
    // Angle comparison
    c->cd(1)->SetGrid();
    c->cd(1)->DrawFrame(30,30,60,50,"Angle comparison;(StraightLine) PMT angle from fibre [deg];LightPath initial direction from fibre [deg]");
    for (int l=NLOC-1; l>=0; l--) {
      hang[l]->SetMarkerStyle(6);
      hang[l]->SetMarkerColor(col[l]);
      hang[l]->SetLineColor(col[l]);
      hang[l]->Draw("scat same");
    }
    //leg.Draw();
    // PMT angle
    c->cd(2)->SetGrid();
    c->cd(2)->DrawFrame(30,20,60,60,"Angle comparison;(StraightLine) PMT angle from fibre [deg];Incoming angle at PMT [deg]");
    for (int l=NLOC-1; l>=0; l--) {
      hbin[l]->SetMarkerStyle(6);
      hbin[l]->SetMarkerColor(col[l]);
      hbin[l]->SetLineColor(col[l]);
      hbin[l]->Draw("scat same");
    }
    // Total light path
    c->cd(3)->SetGrid();
    c->cd(3)->DrawFrame(30,10,60,16,"Distance comparison;(StraightLine) PMT angle from fibre [deg];Total light path length [m]");
    for (int l=NLOC-1; l>=0; l--) {
      hdis[l]->SetMarkerStyle(6);
      hdis[l]->SetMarkerColor(col[l]);
      hdis[l]->SetLineColor(col[l]);
      hdis[l]->Draw("scat same");
    }
    // Light path type
    c->cd(4)->SetGrid();
    c->cd(4)->DrawFrame(30,0,60,8,"Type comparison;(StraightLine) PMT angle from fibre [deg];Light path type [ ]");
    for (int l=NLOC-1; l>=0; l--) {
      htyp[l]->SetMarkerStyle(6);
      htyp[l]->SetMarkerColor(col[l]);
      htyp[l]->SetLineColor(col[l]);
      htyp[l]->Draw("scat same");
    }
    // Bucket time
    c->cd(5)->SetGrid();
    c->cd(5)->DrawFrame(30,0.45,60,0.75,"Bucket time comparison;(StraightLine) PMT angle from fibre [deg];PMT bucket time [ns]");
    for (int l=NLOC-1; l>=0; l--) {
      hbkt[l]->SetMarkerStyle(6);
      hbkt[l]->SetMarkerColor(col[l]);
      hbkt[l]->SetLineColor(col[l]);
      hbkt[l]->Draw("scat same");
    }
    // Transit time
    c->cd(6)->SetGrid();
    c->cd(6)->DrawFrame(30,40,60,70,"Transit time comparison;(StraightLine) PMT angle from fibre [deg];Transit time [ns]");
    for (int l=NLOC-1; l>=0; l--) {
      htof[l]->SetMarkerStyle(6);
      htof[l]->SetMarkerColor(col[l]);
      htof[l]->SetLineColor(col[l]);
      htof[l]->Draw("scat same");
    }
    // Save and close
    c->Print("lightpath_zoom.png");
    //c->Print("lightpath_zoom.pdf");
    c->Close();
    
    // Comparison of different path lengths
    c = new TCanvas("c","",1200,800);
    c->Divide(2,2);
    // Angle comparison
    c->cd(1)->SetGrid();
    c->cd(1)->DrawFrame(0,0,100,13,"Angle comparison;(StraightLine) PMT angle from fibre [deg];Path in inner AV [m]");
    for (int l=0; l<NLOC; l++) {
      hinner[l]->SetMarkerStyle(6);
      hinner[l]->SetMarkerColor(col[l]);
      hinner[l]->SetLineColor(col[l]);
      hinner[l]->Draw("scat same");
    }
    //leg.DrawLegend(0.58,0.88,0.58,0.88);
    // Angle comparison
    c->cd(2)->SetGrid();
    c->cd(2)->DrawFrame(0,0,100,2,"Angle comparison;(StraightLine) PMT angle from fibre [deg];Path in AV [m]");
    for (int l=0; l<NLOC; l++) {
      hinAV[l]->SetMarkerStyle(6);
      hinAV[l]->SetMarkerColor(col[l]);
      hinAV[l]->SetLineColor(col[l]);
      hinAV[l]->Draw("scat same");
    }
    // Angle comparison
    c->cd(3)->SetGrid();
    c->cd(3)->DrawFrame(0,0,100,13,"Angle comparison;(StraightLine) PMT angle from fibre [deg];Path in ext. water [m]");
    for (int l=0; l<NLOC; l++) {
      hwater[l]->SetMarkerStyle(6);
      hwater[l]->SetMarkerColor(col[l]);
      hwater[l]->SetLineColor(col[l]);
      hwater[l]->Draw("scat same");
    }
    // Angle comparison
    c->cd(4)->SetGrid();
    c->cd(4)->DrawFrame(0,0,100,20,"Angle comparison;(StraightLine) PMT angle from fibre [deg];Total path length [m]");
    for (int l=0; l<NLOC; l++) {
      hdis[l]->SetMarkerStyle(6);
      hdis[l]->SetMarkerColor(col[l]);
      hdis[l]->SetLineColor(col[l]);
      hdis[l]->Draw("scat same");
    }
    // Save and close
    c->Print("lightpath_dist.png");
    //c->Print("lightpath_dist.pdf");
    c->Close();
    
    // Distance and angular distributions (1D)
    c = new TCanvas("c","",1600,1200);
    c->Divide(3,2);
    // Distance from straight line
    c->cd(1)->SetGrid();
    c->cd(1)->DrawFrame(0,0,20,1.5*hdisLP[0]->GetMaximum(),";Absolute distance (Fibre #leftrightarrow PMT) d_{PMT} [m];");
    for (int l=NLOC-1; l>=0; l--) {
      hdisSL[l]->SetLineColorAlpha(col[l],0.5);
      hdisSL[l]->SetLineWidth(2);
      hdisSL[l]->Draw("same");
    }
    leg.Draw();
    // Angle seen from AV [deg]
    c->cd(2)->SetGrid();
    c->cd(2)->DrawFrame(0,0,180,1.5*hangAV[0]->GetMaximum(),";Angle seen from AV centre (Fibre #leftrightarrow PMT) #theta [deg];");
    for (int l=NLOC-1; l>=0; l--) {
      hangAV[l]->SetLineColorAlpha(col[l],0.5);
      hangAV[l]->SetLineWidth(2);
      hangAV[l]->Draw("same");
    }
    // Angle seen from fibre [deg]
    c->cd(3)->SetGrid();
    c->cd(3)->DrawFrame(0,0,180,1.5*hangLP[0]->GetMaximum(),";Angle seen from fibre (Fibre #leftrightarrow PMT) #theta [deg];");
    for (int l=NLOC-1; l>=0; l--) {
      hangLP[l]->SetLineColorAlpha(col[l],0.5);
      hangLP[l]->SetLineWidth(2);
      hangLP[l]->Draw("same");
    }
    // Distance from LightPathCalculator
    c->cd(4)->SetGrid();
    c->cd(4)->DrawFrame(0,0,20,1.5*hdisLP[0]->GetMaximum(),";LightPath distance (Fibre #leftrightarrow PMT) l_{LPC} [m];");
    for (int l=NLOC-1; l>=0; l--) {
      hdisLP[l]->SetLineColorAlpha(col[l],0.5);
      hdisLP[l]->SetLineWidth(2);
      hdisLP[l]->Draw("same");
    }
    // Angle seen from AV cos(theta)
    c->cd(5)->SetGrid();
    c->cd(5)->DrawFrame(-1,0,1,1.5*hcosAV[0]->GetMaximum(),";Angle seen from AV centre (Fibre #leftrightarrow PMT) cos(#theta) [];");
    for (int l=NLOC-1; l>=0; l--) {
      hcosAV[l]->SetLineColorAlpha(col[l],0.5);
      hcosAV[l]->SetLineWidth(2);
      hcosAV[l]->Draw("same");
    }
    // Angle seen from fibre cos(theta)
    c->cd(6)->SetGrid();
    c->cd(6)->DrawFrame(-1,0,1,1.5*hcosLP[0]->GetMaximum(),";Angle seen from fibre (Fibre #leftrightarrow PMT) cos(#theta) [];");
    for (int l=NLOC-1; l>=0; l--) {
      hcosLP[l]->SetLineColorAlpha(col[l],0.5);
      hcosLP[l]->SetLineWidth(2);
      hcosLP[l]->Draw("same");
    }
    // Save and close
    c->Print("lightpath_param.png");
    //c->Print("lightpath_param.pdf");
    c->Close();
    
    // Free memory
    if (c) delete c;
    for (int l=0; l<NLOC; l++) {
      if (hang[l]) delete hang[l];
      if (hbin[l]) delete hbin[l];
      if (hdis[l]) delete hdis[l];
      if (htyp[l]) delete htyp[l];
      if (hbkt[l]) delete hbkt[l];
      if (htof[l]) delete htof[l];
      if (hinner[l]) delete hinner[l];
      if (hinAV[l]) delete hinAV[l];
      if (hwater[l]) delete hwater[l];
      if (hdisLP[l]) delete hdisLP[l];
      if (hangAV[l]) delete hangAV[l];
      if (hangLP[l]) delete hangLP[l];
      if (hcosAV[l]) delete hcosAV[l];
      if (hcosLP[l]) delete hcosLP[l];
    }
    return 0;
}
