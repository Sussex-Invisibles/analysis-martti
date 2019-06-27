// ---------------------------------------------------
//  Compile:  g++ -g -O2 -o 0vbb_phasespace.exe 0vbb_phasespace.C `root-config --cflags --libs` -I$RATROOT/include/libpq -I$RATROOT/include -L$RATROOT/lib -lRATEvent_Linux -lGeom
//  Run:      ./0vbb_phasespace.exe
//  Info:     Will take 2-3 min to generate distributions. Run a second time to plot.
//  Author:   Martti Nirkko, April 2019 (based on work done at INSS 2015)
// ---------------------------------------------------
// Include helper functions / plotting style
#include "../include/HelperFunc.C"
#include "../include/CommonStyle.H"
#include <TArrow.h>

// Global constants
const int NSTATS = 1e5;           // statistics for contour limits
const int NTHROWS = 1e9;          // full phase space
const int NTHROWSLOWMASS = 1e6;   // low m_eff phase space
const int NBINS = 4e2;            // number of bins

// Global plot options
const int PALETTE = 51;           // ROOT color palette
const int DEFCOL  = kGray+3;      // text color for Majorana phase definitions
const int HIERCOL = kBlue+3;      // text color for hierarchy labels
const int BANDCOL = kOrange+1;    // fill color for sensitivity bands
const int LINECOL = kOrange+2;    // line color for sensitivity bands
const int TEXTCOL = kOrange+4;    // text color for sensitivity bands
const int CONTCOL = kViolet;      // line color for 0-3 sigma contours
const int CONTWDT = 1.5;          // line width for 0-3 sigma contours
//const int CONTSTY[4] = {1,10,6,4};    // line styles for contours

// Neutrino mixing parameters (PDG 2018)
const double S12 = 0.307;
const double ES12 = 0.013;
const double DM21 = 7.53e-5;
const double EDM21 = 0.18e-5;
const double DM32NH = 2.51e-3; // NH
const double EDM32NH = 0.05e-3;
const double DM32IH = -2.56e-3; // IH
const double EDM32IH = 0.04e-3;
const double S13 = 2.12e-2;
const double ES13 = 0.08e-2;

// Global functions and generators
double meff(double m0, double& msum, double& mtrit, double alpha, double beta, float sigma=0., int inverted=0);
double mefflimits(double m0, double& msumlim, float sigma=1., int index=0);

// Random number generator
TRandom3* rng = new TRandom3();

// *****************************************************************************
int main() {

  // I/O file (either plot from ROOT file, or generate it)
  std::string fname = "0vbb_phasespace";
  std::string output = fname+".root";
  ifstream io(output.c_str());
  
  // Initialise objects
  TFile *file = NULL;
  TH2D *hpsNH=NULL, *hpsIH=NULL, *hpsNHlo=NULL, *hangNHlo=NULL;
  TH2D *hpsNHc=NULL, *hpsIHc=NULL, *hpsNHs=NULL, *hpsIHs=NULL, *hpsNHt=NULL, *hpsIHt=NULL;
  TGraph *contours[4][4] = {NULL}, *contours2[4][4] = {NULL}, *contours3[4][4] = {NULL};
  TGraph *lines[4] = {NULL};
  TGraph *f0=NULL, *f1=NULL;
  TCanvas *c=NULL;
  
  if (!io.good()) {  // File does not exist
    
    // Create ROOT output file
    file = new TFile(output.c_str(),"RECREATE");
    
    // Histograms for normal and inverted hierarchy
    hpsNH = new TH2D("hpsNH","Normal hierarchy",NBINS,1e-4,1e0,NBINS,1e-4,1e0);
    BinLog(hpsNH->GetXaxis()); BinLog(hpsNH->GetYaxis());
    hpsIH = new TH2D("hpsIH","Inverted hierarchy",NBINS,1e-4,1e0,NBINS,1e-4,1e0);
    BinLog(hpsIH->GetXaxis()); BinLog(hpsIH->GetYaxis());
    
    // Histograms for very small effective masses (NH)
    hpsNHlo = new TH2D("hpsNHlo","Small effective mass",NBINS,1e-4,1e0,NBINS,1e-8,1e-4);
    BinLog(hpsNHlo->GetXaxis()); BinLog(hpsNHlo->GetYaxis());
    hangNHlo = new TH2D("hangNHlo","Majorana hangNHlo",NBINS,0,1,NBINS,0,1);
    
    // Histograms for plotting sum of neutrino masses as x-axis
    hpsNHc = new TH2D("hpsNHc","Normal hierarchy",NBINS,3e-2,3e0,NBINS,1e-4,1e0);
    BinLog(hpsNHc->GetXaxis()); BinLog(hpsNHc->GetYaxis());
    hpsIHc = new TH2D("hpsIHc","Inverted hierarchy",NBINS,3e-2,3e0,NBINS,1e-4,1e0);
    BinLog(hpsIHc->GetXaxis()); BinLog(hpsIHc->GetYaxis());
    
    // Histograms for plotting sum of neutrino masses as y-axis
    hpsNHs = new TH2D("hpsNHs","Normal hierarchy",NBINS,1e-4,1e0,10*NBINS,3e-2,3e0);
    BinLog(hpsNHs->GetXaxis()); BinLog(hpsNHs->GetYaxis());
    hpsIHs = new TH2D("hpsIHs","Inverted hierarchy",NBINS,1e-4,1e0,10*NBINS,3e-2,3e0);
    BinLog(hpsIHs->GetXaxis()); BinLog(hpsIHs->GetYaxis());
    
    // Histograms for effective neutrino mass for tritium decay as y-axis
    hpsNHt = new TH2D("hpsNHt","Normal hierarchy",NBINS,1e-4,1e0,10*NBINS,1e-3,1e0);
    BinLog(hpsNHt->GetXaxis()); BinLog(hpsNHt->GetYaxis());
    hpsIHt = new TH2D("hpsIHt","Inverted hierarchy",NBINS,1e-4,1e0,10*NBINS,1e-3,1e0);
    BinLog(hpsIHt->GetXaxis()); BinLog(hpsIHt->GetYaxis());
    
    // lightest and effective neutrino masses
    double m0, meffNH, meffIH;
    
    // Get maximum lines for allowed phase space region
    std::cout << "Generating contour limits..." << std::endl;
    std::string cname;
    double mefflim, msumlim;
    for (int i=0; i<4; i++) { // min/max for NH/IH
      for (int s=0; s<4; s++) { // sigmas
        contours[i][s] = new TGraph();
        contours2[i][s] = new TGraph();
        contours3[i][s] = new TGraph();
        for (int b=0; b<=NBINS; b++) {
          m0 = hpsNH->GetXaxis()->GetBinCenter(b);
          mefflim = mefflimits(m0,msumlim,s,i);
          contours[i][s]->SetPoint(b,m0,mefflim);
          contours2[i][s]->SetPoint(b,msumlim,mefflim);
          contours3[i][s]->SetPoint(b,m0,msumlim);
          printProgress(i*4*NBINS+s*NBINS+b,16*NBINS);
        }
        cname = Form("contours_%d_%d",i,s);
        contours[i][s]->SetName(cname.c_str());
        cname = Form("contours2_%d_%d",i,s);
        contours2[i][s]->SetName(cname.c_str());
        cname = Form("contours3_%d_%d",i,s);
        contours3[i][s]->SetName(cname.c_str());
      }
    }
    
    // Get interesting lines for specific Majorana angles (NH only)
    double dummy;
    for (int i=0; i<4; i++) lines[i] = new TGraph();
    for (int b=0; b<=1e4; b++) {
      m0 = hpsNH->GetXaxis()->GetBinCenter(b);
      lines[0]->SetPoint(b,m0,meff(m0,dummy,dummy,0,0,0,0));
      lines[1]->SetPoint(b,m0,meff(m0,dummy,dummy,0,pi/2,0,0));
      lines[2]->SetPoint(b,m0,meff(m0,dummy,dummy,pi/2,0,0,0));
      lines[3]->SetPoint(b,m0,meff(m0,dummy,dummy,pi/2,pi/2,0,0));
    }
    for (int i=0; i<4; i++) {
      cname = Form("lines_%d",i);
      lines[i]->SetName(cname.c_str());
    }
    
    // Majorana angles and global minima (actually using phi2, phi3 convention)
    double alpha, beta;
    double alpha0=0., beta0=0., m0min=1.e6, meffmin=1.e6;
    double msumNH, msumIH, mtritNH, mtritIH;
  
    // Random number seed
    rng->SetSeed(0);
    
    // Generate contour plots
    long throws=0, count=0;
    std::cout << "Filling allowed phase space regions..." << std::endl;
    while (throws<NTHROWS) {
      m0 = pow(10.,4.*(rng->Rndm()-1.)); // log prior (100% of phase space)
      alpha = pi*rng->Rndm(); // 0-pi (100% of phase space)
      beta = pi*rng->Rndm(); // 0-pi (100% of phase space)
      meffNH = meff(m0,msumNH,mtritNH,alpha,beta,-1,0);
      meffIH = meff(m0,msumIH,mtritIH,alpha,beta,-1,1);
      throws++;
      hpsNH->Fill(m0, meffNH);
      hpsIH->Fill(m0, meffIH);
      hpsNHc->Fill(msumNH, meffNH);
      hpsIHc->Fill(msumIH, meffIH);
      hpsNHs->Fill(m0, msumNH);
      hpsIHs->Fill(m0, msumIH);
      hpsNHt->Fill(m0, mtritNH);
      hpsIHt->Fill(m0, mtritIH);
      printProgress(throws,NTHROWS);
      if (meffNH>1e-4) continue;
      count++;
      /*
      hpsNHlo->Fill(m0, meffNH); 
      hangNHlo->Fill(alpha/pi,beta/pi);
      */
    }
    hpsNH->Scale(1./throws);
    hpsIH->Scale(1./throws);
    hpsNHc->Scale(1./throws);
    hpsIHc->Scale(1./throws);
    hpsNHs->Scale(1./throws);
    hpsIHs->Scale(1./throws);
    hpsNHt->Scale(1./throws);
    hpsIHt->Scale(1./throws);
    //hpsNHlo->Scale(1./throws);
    //hangNHlo->Scale(1./throws);
    std::cout<<"Generated " << (float)throws << " throws, of which " << throws-count;
    std::cout<<" (" << 100.*(throws-count)/throws << "%) were used for histogram." << std::endl;
    
    long throws2=0, count2=0;
    std::cout << "Filling phase space region for low effective mass..." << std::endl;
    // Region in which values for alpha are generated
    f0 = new TGraph();
    f1 = new TGraph();
    for (int i=0; i<NBINS; i++) {
      double beta = (double)i/NBINS;
      double alpha = 0.5 + 0.065*sin(2.*beta*pi);  // 5% of phase space
      f1->SetPoint(i,alpha+0.025,beta);
      f0->SetPoint(i,alpha-0.025,beta);
    }
    f0->SetName("alpha_lo");
    f1->SetName("alpha_hi");
    while (count2<NTHROWSLOWMASS) {
      m0 = pow(10.,(rng->Rndm()-3.)); // log prior (25% of phase space)
      beta = pi*rng->Rndm();
      alpha = pi*(0.5+0.065*sin(2*beta)+0.05*(rng->Rndm()-0.25)); // 5% of phase space
      meffNH = meff(m0,msumNH,mtritNH,alpha,beta,-1,0);
      throws2++;
      if (meffNH>1e-4) continue;
      hpsNHlo->Fill(m0, meffNH); 
      hangNHlo->Fill(alpha/pi,beta/pi);
      count2++;
      printProgress(count2,NTHROWSLOWMASS);
      if (meffNH<meffmin) { alpha0=alpha; beta0=beta; m0min=m0; meffmin=meffNH; }
    }
    hpsNHlo->Scale(1./throws2*0.25*0.05);
    hangNHlo->Scale(1./throws2*0.25*0.05);
    std::cout<<"Generated " << (float)throws2 << " throws, of which " << hpsNHlo->GetEntries();
    std::cout<<" (" << 100.*hpsNHlo->GetEntries()/throws2 << "%) were used for histogram." << std::endl;
    std::cout<<"MINIMUM: meff="<<1.e3*meffmin<<" meV for m0="<<1.e3*m0min<<" meV, alpha="<<alpha0/pi<<" pi, beta="<<beta0/pi<<" pi"<<std::endl;
    
    file->Write();
    for (int i=0; i<4; i++) { // min/max for NH/IH
      for (int s=0; s<4; s++) { // sigmas
        contours[i][s]->Write();
        contours2[i][s]->Write();
        contours3[i][s]->Write();
      }
    }
    for (int i=0; i<4; i++) lines[i]->Write();
    f0->Write();
    f1->Write();
        
  } else {  // File already exists, plot results
    
    // Use T2K style
    CommonStyle();
    gStyle->SetPadBottomMargin(0.14);
    gStyle->SetPadLeftMargin(0.14);//to include both large/small font options
    gStyle->SetPadRightMargin(0.14);//to include both large/small font options
    gStyle->SetLabelSize(0.05,"xyz");
    gStyle->SetTitleSize(0.05,"xyz");
    gStyle->SetLabelOffset(0.01,"xyz");
    gStyle->SetTitleOffset(1.4,"x");
    gStyle->SetTitleOffset(1.3,"y");
    gStyle->SetMarkerStyle(7);
    gStyle->SetPalette(PALETTE);
    gROOT->ForceStyle();  // force histograms read in below to use this style
  
    // Read from ROOT output file
    file = new TFile(output.c_str(),"READ");
    hpsNH = (TH2D*)file->Get("hpsNH");
    hpsIH = (TH2D*)file->Get("hpsIH");
    hpsNHc = (TH2D*)file->Get("hpsNHc");
    hpsIHc = (TH2D*)file->Get("hpsIHc");
    hpsNHs = (TH2D*)file->Get("hpsNHs");
    hpsIHs = (TH2D*)file->Get("hpsIHs");
    hpsNHt = (TH2D*)file->Get("hpsNHt");
    hpsIHt = (TH2D*)file->Get("hpsIHt");
    hpsNHlo = (TH2D*)file->Get("hpsNHlo");
    hangNHlo = (TH2D*)file->Get("hangNHlo");
    std::string cname;
    TGraph *temp = NULL;
    for (int i=0; i<4; i++) { // min/max for NH/IH
    
      cname = Form("lines_%d",i);
      lines[i] = (TGraph*)file->Get(cname.c_str());
      //lines[i] = new TGraph();
      //Smooth(lines[i],temp,5);
      
      for (int s=0; s<4; s++) { // sigmas
      
        cname = Form("contours_%d_%d",i,s);
        temp = (TGraph*)file->Get(cname.c_str());
        contours[i][s] = new TGraph();
        Smooth(contours[i][s],temp,5);
        
        cname = Form("contours2_%d_%d",i,s);
        temp = (TGraph*)file->Get(cname.c_str());
        contours2[i][s] = new TGraph();
        Smooth(contours2[i][s],temp,5);
        
        cname = Form("contours3_%d_%d",i,s);
        temp = (TGraph*)file->Get(cname.c_str());
        contours3[i][s] = new TGraph();
        Smooth(contours3[i][s],temp,5);
        
        if(temp) delete temp;
      }
    }
    f0 = (TGraph*)file->Get("alpha_lo");
    f1 = (TGraph*)file->Get("alpha_hi");
    
    // ------------------
    //  PLOTTING SECTION
    // ------------------
    
    // Text in figures
    TLatex *Tl = new TLatex();
    Tl->SetTextSize(0.041);
    
    // SNO+ sensitivity
    double fac = sqrt(10./2.1); // phase II/I lifetime sensitivity (1e26y)
    double sx[4] = {1e-5,1e1,1e1,1e-5};
    double sy1[4] = {89e-3,89e-3,37e-3,37e-3}; // eV (phase I)
    double sy2[4] = {0};
    for (int k=0; k<4; k++) sy2[k] = sy1[k]/fac; // eV (phase II)
    TGraph *gs1 = new TGraph(4,sx,sy1);
    TGraph *gs2 = new TGraph(4,sx,sy2);
    
    // DBD results (eV)
    double rgerda[4] = {0.2,0.2,0.4,0.4};
    double rklzen[4] = {0.061,0.061,0.165,0.165};
    double rcuore[4] = {0.27,0.27,0.76,0.76};
    double rnemo[4] = {0.33,0.33,0.62,0.62};
    double rtotal[4] = {0.061,0.061,0.76,0.76};
    
    // DBD sensitivities
    double sklzen[4] = {0.04,0.04,0.108,0.108};
    double ssnemo[4] = {0.2,0.2,0.4,0.4};
    double ssnop[4] = {0.037,0.037,0.089,0.089};
    double snext[4] = {0.08,0.08,0.16,0.16};
    double stotal[4] = {0.037,0.037,0.4,0.4};
    
    // Future sensitivities
    double fnexo[4] = {0.0057,0.0057,0.017,0.017};
    double fsnop[4] = {0.017,0.017,0.041,0.041};
    double ftotal[4] = {0.0057,0.0057,0.017,0.017};
    
    // Graphs
    TGraph *gres = new TGraph(4,sx,rklzen);
    TGraph *gsen = new TGraph(4,sx,ssnop);
    TGraph *gfut = new TGraph(4,sx,ftotal);
    gres->SetLineStyle(5);
    gres->SetLineWidth(2);
    gres->SetLineColorAlpha(kGreen+2,0.75);
    gres->SetFillColorAlpha(kGreen+2,0.5);
    gsen->SetLineStyle(5);
    gsen->SetLineWidth(2);
    gsen->SetLineColorAlpha(kOrange,0.75);
    gsen->SetFillColorAlpha(kOrange,0.5);
    gfut->SetLineStyle(5);
    gfut->SetLineWidth(2);
    gfut->SetLineColorAlpha(kRed,0.75);
    gfut->SetFillColorAlpha(kRed,0.5);
    
    // Cosmology limits
    double cx[4] = {0.12,0.12};
    double cy[4] = {1e-5,1e1};
    TGraph *gc = new TGraph(2,cx,cy);
    TGraph *gs = new TGraph(2,cy,cx);
    
    // KATRIN sensitivity
    double tx[4] = {1e-5,1e1};
    double ty[4] = {0.2,0.2};
    TGraph *gt = new TGraph(2,tx,ty);
    
    // Plot full phase space for double beta decay (NH+IH)
    // ---------------------------------------------------
    c = new TCanvas("c","Phase space plot",1200,900);
    c->SetLogx(); c->SetLogy(); c->SetLogz(); c->SetGrid();
    //hpsIH->SetTitle("#bf{0#nu#beta#beta decay - Type-I seesaw}");
    c->DrawFrame(1e-4,1e-4,1e0,1e0,";m_{0} (eV);#LTm_{#beta#beta}#GT (eV)");
    // Draw allowed phase space for IH
    hpsIH->Draw("colz same");
    hpsIH->GetZaxis()->SetRangeUser(2e-10,2e-4);
    for (int i=2; i<4; i++) { // contours for IH
      for (int s=0; s<4; s++) { // contours for IH
        //contours[i][s]->SetLineStyle(CONTSTY[s]);
        contours[i][s]->SetLineColor(CONTCOL+s);
        contours[i][s]->SetLineWidth(CONTWDT);
        contours[i][s]->Draw("l");
      }
    }
    // Draw allowed phase space for NH
    hpsNH->Draw("colz same");
    hpsNH->GetZaxis()->SetRangeUser(2e-10,2e-4);
    for (int i=0; i<2; i++) { // contours for NH
      for (int s=0; s<4; s++) { // contours for NH
        //contours[i][s]->SetLineStyle(CONTSTY[s]);
        contours[i][s]->SetLineColor(CONTCOL+s);
        contours[i][s]->SetLineWidth(CONTWDT);
        contours[i][s]->Draw("l");
      }
    }
    // Draw legend for contours
    TLegend *leg = new TLegend(0.17,0.78,0.34,0.94,"#bf{PDG 2018}");
    for (int s=0; s<4; s++) {
      std::string legstr = Form("%d #sigma",s);
      leg->AddEntry(contours[0][s],legstr.c_str(),"l");
    }
    leg->SetFillStyle(0);
    leg->Draw();
    // Draw text describing plots
    Tl->SetTextSize(0.035);
    Tl->SetTextColor(HIERCOL);
    Tl->DrawLatex(1.5e-4,2.5e-2,"IH");
    Tl->DrawLatex(1.5e-4,2.0e-3,"NH");
    // Evince seems to have trouble displaying uppercase delta (font not embedded)
    //Tl->DrawLatex(1.5e-4,2.5e-2,"#Deltam_{32}^{2} < 0");
    //Tl->DrawLatex(1.5e-4,2.0e-3,"#Deltam_{32}^{2} > 0");
    // Save and close
    c->Update();
    c->RedrawAxis();
    c->Print((fname+".png").c_str());
    c->Print((fname+".pdf").c_str());
    c->Close();
    
    // Same plot with SNO+ sensitivity bands
    // ---------------------------------------------------
    c = new TCanvas("c","Phase space plot",1200,900);
    c->SetLogx(); c->SetLogy(); c->SetLogz(); c->SetGrid();
    //hpsIH->SetTitle("#bf{0#nu#beta#beta decay - Type-I seesaw}");
    c->DrawFrame(1e-4,1e-4,1e0,1e0,";m_{0} (eV);#LTm_{#beta#beta}#GT (eV)");
    // Draw SNO+ sensitivity bands (areas only)
    gs1->SetFillColorAlpha(BANDCOL,0.5);
    gs2->SetFillColorAlpha(BANDCOL,0.5);
    gs1->Draw("f");
    gs2->Draw("f");
    // Draw allowed phase space for IH
    hpsIH->Draw("colz same");
    hpsIH->GetZaxis()->SetRangeUser(2e-10,2e-4);
    for (int i=2; i<4; i++) { // contours for IH
      for (int s=0; s<4; s++) { // contours for IH
        //contours[i][s]->SetLineStyle(CONTSTY[s]);
        contours[i][s]->SetLineColor(CONTCOL+s);
        contours[i][s]->SetLineWidth(CONTWDT);
        contours[i][s]->Draw("l");
      }
    }
    // Draw allowed phase space for NH
    hpsNH->Draw("colz same");
    hpsNH->GetZaxis()->SetRangeUser(2e-10,2e-4);
    for (int i=0; i<2; i++) { // contours for NH
      for (int s=0; s<4; s++) { // contours for NH
        //contours[i][s]->SetLineStyle(CONTSTY[s]);
        contours[i][s]->SetLineColor(CONTCOL+s);
        contours[i][s]->SetLineWidth(CONTWDT);
        contours[i][s]->Draw("l");
      }
    }
    // Draw legend for contours
    leg->Draw();
    // Draw text describing plots
    Tl->SetTextSize(0.035);
    Tl->SetTextColor(HIERCOL);
    Tl->DrawLatex(1.5e-4,2.5e-2,"IH");
    Tl->DrawLatex(1.5e-4,2.0e-3,"NH");
    // Evince seems to have trouble displaying uppercase delta (font not embedded)
    //Tl->DrawLatex(1.5e-4,2.5e-2,"#Deltam_{32}^{2} < 0");
    //Tl->DrawLatex(1.5e-4,2.0e-3,"#Deltam_{32}^{2} > 0");
    Tl->SetTextSize(0.025);
    Tl->SetTextColor(TEXTCOL);
    Tl->DrawLatex(2e-4,6.0e-2,"SNO+ Phase I (T_{1/2} > 2.1 #times 10^{26} y)");
    Tl->DrawLatex(1.5e-1,2.4e-2,"SNO+ Phase II");
    // Redraw sensitivity bands (lines only)
    gs1->SetLineColorAlpha(LINECOL,0.75);
    gs2->SetLineColorAlpha(LINECOL,0.75);
    gs1->SetLineStyle(5);
    gs2->SetLineStyle(5);
    gs1->SetLineWidth(2);
    gs2->SetLineWidth(2);
    gs1->Draw("l");
    gs2->Draw("l");
    // Save and close
    c->Update();
    c->RedrawAxis();
    c->Print((fname+"_snoplus.png").c_str());
    c->Print((fname+"_snoplus.pdf").c_str());
    c->Close();
    
    // Plot full phase space for double beta decay (NH+IH)
    // ---------------------------------------------------
    c = new TCanvas("c","Phase space plot",1200,900);
    c->SetLogx(); c->SetLogy(); c->SetLogz(); c->SetGrid();
    //hpsIH->SetTitle("#bf{0#nu#beta#beta decay - Type-I seesaw}");
    c->DrawFrame(1e-4,1e-4,1e0,1e0,";m_{0} (eV);#LTm_{#beta#beta}#GT (eV)");
    // Draw allowed phase space for IH
    hpsIH->Draw("colz same");
    hpsIH->GetZaxis()->SetRangeUser(2e-10,2e-4);
    for (int i=2; i<4; i++) { // contours for IH
      for (int s=0; s<4; s++) { // contours for IH
        //contours[i][s]->SetLineStyle(CONTSTY[s]);
        contours[i][s]->SetLineColor(CONTCOL+s);
        contours[i][s]->SetLineWidth(CONTWDT);
        contours[i][s]->Draw("l");
      }
    }
    // Draw allowed phase space for NH
    hpsNH->Draw("colz same");
    hpsNH->GetZaxis()->SetRangeUser(2e-10,2e-4);
    for (int i=0; i<2; i++) { // contours for NH
      for (int s=0; s<4; s++) { // contours for NH
        //contours[i][s]->SetLineStyle(CONTSTY[s]);
        contours[i][s]->SetLineColor(CONTCOL+s);
        contours[i][s]->SetLineWidth(CONTWDT);
        contours[i][s]->Draw("l");
      }
    }
    // Draw text describing plots
    Tl->SetTextSize(0.035);
    Tl->SetTextColor(HIERCOL);
    Tl->DrawLatex(1.5e-4,2.5e-2,"IH");
    Tl->DrawLatex(1.5e-4,2.0e-3,"NH");
    // Sensitivity arrows
    TArrow *arr1 = new TArrow(0,0,1,1,0.01,"|->");
    TArrow *arr2 = new TArrow(0,0,1,1,0.01,"|->");
    TArrow *arr3 = new TArrow(0,0,1,1,0.01,"|->");
    arr1->SetLineWidth(2);
    arr2->SetLineWidth(2);
    arr3->SetLineWidth(2);
    arr1->SetLineColor(kGreen+2);   // completed
    arr2->SetLineColor(kOrange+2);  // construction
    arr3->SetLineColor(kRed+1);     // proposed
    arr1->DrawArrow(3.125e-1,0.76,3.125e-1,0.11); // EXO-200, GERDA I, CUORE-0, NEMO3
    arr1->DrawArrow(6.25e-2,0.165,6.25e-2,0.061); // Kamland-Zen I-II
    arr2->DrawArrow(1.25e-2,0.1,1.25e-2,0.037); // SNO+ I, Kamland-Zen 800
    arr3->DrawArrow(2.5e-3,0.06,2.5e-3,0.017); // SNO+ II, PandaX-III, CUPID, AMoRE
    arr3->DrawArrow(5e-4,0.0177,5e-4,0.0057); // nEXO, LEGEND-1000
    // Describe bands
    Tl->SetTextSize(0.022);
    Tl->SetTextColor(kGreen+2);
    Tl->DrawLatex(1.7e-4,0.6,"Completed");
    Tl->DrawLatex(1.15e-2,0.6,"CUORE-0, EXO-200, GERDA, NEMO3");
    Tl->DrawLatex(2e-2,0.19,"KamLAND-Zen I-II");
    Tl->SetTextColor(kOrange+2);
    Tl->DrawLatex(1.7e-4,0.4,"Upcoming");
    Tl->DrawLatex(1.2e-3,0.11,"KamLAND-Zen 800, SNO+ I, SuperNEMO");
    Tl->SetTextColor(kRed+1);
    Tl->DrawLatex(1.7e-4,0.27,"Proposed");
    Tl->DrawLatex(4e-4,0.073,"AMoRE, CUPID, KamLAND2-Zen,");
    Tl->DrawLatex(4e-4,0.057,"PandaX-III, SNO+ II");
    Tl->DrawLatex(1.2e-4,7.5e-3,"LEGEND-1000, nEXO");
    // Save and close
    c->Update();
    c->RedrawAxis();
    c->Print((fname+"_sensitivities.png").c_str());
    c->Print((fname+"_sensitivities.pdf").c_str());
    c->Close();
    
    // Plot phase space for double beta decay (minimalistic)
    // -----------------------------------------------------
    c = new TCanvas("c","Phase space plot",1200,900);
    c->SetLogx(); c->SetLogy(); c->SetLogz(); c->SetGrid();
    //hpsIH->SetTitle("#bf{0#nu#beta#beta decay - Type-I seesaw}");
    c->DrawFrame(1e-4,1e-4,1e0,1e0,";m_{0} (eV);#LTm_{#beta#beta}#GT (eV)");
    // Draw allowed phase space for IH
    hpsIH->Draw("colz same");
    hpsIH->GetZaxis()->SetRangeUser(2e-10,2e-4);
    // Draw allowed phase space for NH
    hpsNH->Draw("colz same");
    hpsNH->GetZaxis()->SetRangeUser(2e-10,2e-4);
    Tl->SetTextSize(0.035);
    Tl->SetTextColor(HIERCOL);
    Tl->DrawLatex(1.5e-4,2.5e-2,"IH");
    Tl->DrawLatex(1.5e-4,2.0e-3,"NH");
    // Save a first time
    c->Update();
    c->RedrawAxis();
    c->Print((fname+"_minimal.png").c_str());
    c->Print((fname+"_minimal.pdf").c_str());
    // Draw lines of interest for NH
    for (int i=0; i<4; i++) {
      lines[i]->SetLineWidth(2);
      lines[i]->SetLineColor(kOrange+i);
      lines[i]->Draw();
    }
    TLegend *leglines = new TLegend(0.17,0.79,0.31,0.94,"#bf{Majorana phases}");
    leglines->AddEntry(lines[0],"#varphi_{2}=0#lower[-0.2]{#circ}, #varphi_{3}=0#lower[-0.2]{#circ}","l");
    leglines->AddEntry(lines[1],"#varphi_{2}=0#lower[-0.2]{#circ}, #varphi_{3}=90#lower[-0.2]{#circ}","l");
    leglines->AddEntry(lines[2],"#varphi_{2}=90#lower[-0.2]{#circ}, #varphi_{3}=0#lower[-0.2]{#circ}","l");
    leglines->AddEntry(lines[3],"#varphi_{2}=90#lower[-0.2]{#circ}, #varphi_{3}=90#lower[-0.2]{#circ}","l");
    leglines->Draw();
    // Save again and close
    c->Update();
    c->RedrawAxis();
    c->Print((fname+"_phenom.png").c_str());
    c->Print((fname+"_phenom.pdf").c_str());
    c->Close();
    
    // Plot effective mass vs sum of masses (cosmology)
    // ------------------------------------------------
    c = new TCanvas("c","Phase space plot",1200,900);
    c->SetLogx(); c->SetLogy(); c->SetLogz(); c->SetGrid();
    //hpsNHlo->SetTitle("#bf{0#nu#beta#beta decay - allowed phase space}");
    c->DrawFrame(3e-2,1e-4,3e0,1e0,";#Sigma m_{#nu} (eV);#LTm_{#beta#beta}#GT (eV)");
    // Draw sensitivity bands (areas only)
    gres->Draw("f");
    // Draw allowed phase space
    hpsNHc->Draw("colz same");
    hpsNHc->GetZaxis()->SetRangeUser(2e-10,2e-2);
    for (int i=0; i<2; i++) { // contours for NH
      for (int s=0; s<4; s++) { // contours for NH
        //contours2[i][s]->SetLineStyle(CONTSTY[s]);
        contours2[i][s]->SetLineColor(CONTCOL+s);
        contours2[i][s]->SetLineWidth(CONTWDT);
        contours2[i][s]->Draw("l");
      }
    }
    hpsIHc->Draw("colz same");
    hpsIHc->GetZaxis()->SetRangeUser(2e-10,2e-2);
    for (int i=2; i<4; i++) { // contours for IH
      for (int s=0; s<4; s++) { // contours for IH
        //contours2[i][s]->SetLineStyle(CONTSTY[s]);
        contours2[i][s]->SetLineColor(CONTCOL+s);
        contours2[i][s]->SetLineWidth(CONTWDT);
        contours2[i][s]->Draw("l");
      }
    }
    // Close contours at x-cutoff by joining min/max graphs
    TLine join(0,0,1,1); // dummy constructor
    for (int i=0; i<2; i++) {
      for (int s=0; s<4; s++) {
        join.SetLineColor(CONTCOL+s);
        join.SetLineWidth(CONTWDT);
        join.DrawLine(contours2[2*i][s]->GetX()[1],contours2[2*i][s]->GetY()[1],
                      contours2[2*i+1][s]->GetX()[1],contours2[2*i][s]->GetY()[1]);
        join.DrawLine(contours2[2*i+1][s]->GetX()[1],contours2[2*i][s]->GetY()[1],
                      contours2[2*i+1][s]->GetX()[1],contours2[2*i+1][s]->GetY()[1]);
      }
    }
    // Draw legend for contours
    leg->Draw();
    // Write SNO+ sensitivities and cosmology limits
    Tl->SetTextSize(0.025);
    Tl->SetTextColor(kGreen+4);
    Tl->DrawLatex(4.2e-2,9e-2,"KamLAND-Zen");
    Tl->DrawLatex(1.3e-1,1.3e-4,"#color[13]{PLANCK 2018}");
    // Write mass ordering
    Tl->SetTextSize(0.035);
    Tl->SetTextColor(0);
    Tl->DrawLatex(7.0e-2,5.0e-3,"NH");
    Tl->DrawLatex(1.2e-1,2.7e-2,"IH");
    // Redraw sensitivity bands (lines only)
    gres->Draw("l");
    gc->SetLineColorAlpha(13,0.8);
    gc->SetLineStyle(5);
    gc->SetLineWidth(2);
    gc->Draw("l");
    // Save and close
    c->Print((fname+"_cosmology.png").c_str());
    c->Print((fname+"_cosmology.pdf").c_str());
    c->Close();
    
    // Plot lightest mass vs sum of masses (cosmology)
    // ------------------------------------------------
    c = new TCanvas("c","Phase space plot",1200,900);
    c->SetLogx(); c->SetLogy(); c->SetLogz(); c->SetGrid();
    //hpsNHlo->SetTitle("#bf{0#nu#beta#beta decay - allowed phase space}");
    c->DrawFrame(1e-4,3e-2,1e0,3e0,";m_{0} (eV);#Sigma m_{#nu} (eV)");
    // Draw degenerate region (line)
    TLine deg(1e-2,3e-2,1e0,3e0);
    deg.SetLineWidth(3);
    deg.SetLineStyle(7);
    deg.SetLineColorAlpha(2,0.5);
    deg.Draw("l");
    Tl->SetTextSize(0.025);
    Tl->SetTextColorAlpha(2,0.8);
    Tl->DrawLatex(1.5e-2,3.5e-2,"#Sigma m_{#nu} = 3m_{0}");
    // Draw allowed phase space
    hpsNHs->Draw("colz same");
    hpsNHs->GetZaxis()->SetRangeUser(2e-7,2e-4);
    hpsIHs->Draw("colz same");
    hpsIHs->GetZaxis()->SetRangeUser(2e-7,2e-4);
    // Write cosmology limits
    Tl->SetTextSize(0.025);
    Tl->DrawLatex(1.4e-4,1.3e-1,"#color[13]{PLANCK 2018}");
    // Write mass ordering
    Tl->SetTextSize(0.035);
    Tl->SetTextColor(HIERCOL);
    Tl->DrawLatex(1.4e-4,4.8e-2,"NH");
    Tl->DrawLatex(1.4e-4,8.0e-2,"IH");
    // Draw sensitivity line
    gs->SetLineColorAlpha(13,0.8);
    gs->SetLineStyle(5);
    gs->SetLineWidth(2);
    gs->Draw("l");
    // Save and close
    c->Print((fname+"_sum.png").c_str());
    c->Print((fname+"_sum.pdf").c_str());
    c->Close();
    
    // Plot lightest mass vs effective mass (tritium)
    // ----------------------------------------------
    c = new TCanvas("c","Phase space plot",1200,900);
    c->SetLogx(); c->SetLogy(); c->SetLogz(); c->SetGrid();
    //hpsNHlo->SetTitle("#bf{0#nu#beta#beta decay - allowed phase space}");
    c->DrawFrame(1e-4,1e-3,1e0,1e0,";m_{0} (eV);#LTm_{#nu_{e}}#GT (eV)");
    // Draw degenerate region (line)
    deg.DrawLine(1e-3,1e-3,1e0,1e0);
    Tl->SetTextSize(0.025);
    Tl->SetTextColorAlpha(2,0.8);
    Tl->DrawLatex(2e-3,1.4e-3,"#LTm_{#nu_{e}}#GT = m_{0}");
    // Draw allowed phase space
    hpsNHt->Draw("colz same");
    hpsNHt->GetZaxis()->SetRangeUser(3.2e-7,3.2e-4);
    hpsIHt->Draw("colz same");
    hpsIHt->GetZaxis()->SetRangeUser(3.2e-7,3.2e-4);
    // Write KATRIN sensitivity
    Tl->SetTextSize(0.025);
    Tl->SetTextColor(TEXTCOL);
    Tl->DrawLatex(1.4e-4,2.3e-1,"KATRIN (3 years)");
    // Write mass ordering
    Tl->SetTextSize(0.035);
    Tl->SetTextColor(HIERCOL);
    Tl->DrawLatex(1.4e-4,1.1e-2,"NH");
    Tl->DrawLatex(1.4e-4,6.0e-2,"IH");
    // Draw sensitivity line
    gt->SetLineColorAlpha(LINECOL,0.75);
    gt->SetLineStyle(5);
    gt->SetLineWidth(2);
    gt->Draw("l");
    // Save and close
    c->Print((fname+"_tritium.png").c_str());
    c->Print((fname+"_tritium.pdf").c_str());
    c->Close();
    
    // Plot low m_eff region for double beta decay (NH)
    // ------------------------------------------------
    c = new TCanvas("c","Phase space plot",1200,900);
    c->SetLogx(); c->SetLogy(); c->SetLogz(); c->SetGrid();
    //hpsNHlo->SetTitle("#bf{0#nu#beta#beta decay - allowed phase space}");
    c->DrawFrame(1e-4,1e-7,1e0,1e-4,";m_{0} (eV);#LTm_{#beta#beta}#GT (eV)");
    hpsNHlo->Draw("colz same");
    hpsNHlo->GetZaxis()->SetRangeUser(5e-12,2e-7);
    // Contours are rather jaggy (stats limited) here
    for (int s=0; s<4; s++) { // contours for NH
      contours[0][s]->Draw("l");
      contours[1][s]->Draw("l");
    }
    // Draw legend (defined above)
    leg->Draw();
    // Save and close
    c->Print((fname+"_lowmass.png").c_str());
    c->Print((fname+"_lowmass.pdf").c_str());
    c->Close();
    
    // Plot Majorana angles for low m_eff region (NH)
    // ----------------------------------------------
    c = new TCanvas("c","Majorana angle plot",1200,900);
    c->SetLogz(); c->SetGrid();
    //hangNHlo->SetTitle("#bf{Majorana hangNHlo that lead to m_{#beta#beta} < 10^{-4} eV}");
    c->DrawFrame(0,0,1,1,";#varphi_{2} [#pi];#varphi_{3} [#pi]");
    hangNHlo->Draw("colz same");
    //hangNHlo->GetZaxis()->SetRangeUser(8e-12,4e-8);
    // Draw region in which points were generated
    f0->SetLineColor(LINECOL);
    f1->SetLineColor(LINECOL);
    f0->Draw("l");
    f1->Draw("l");
    // Save and close
    c->Print((fname+"_lowmass_angles.png").c_str());
    c->Print((fname+"_lowmass_angles.pdf").c_str());
    c->Close();
    
  }
  
  // Free memory
  if (hpsNH) delete hpsNH;
  if (hpsIH) delete hpsIH;
  if (hpsNHc) delete hpsNHc;
  if (hpsIHc) delete hpsIHc;
  if (hpsNHs) delete hpsNHs;
  if (hpsIHs) delete hpsIHs;
  if (hpsNHt) delete hpsNHt;
  if (hpsIHt) delete hpsIHt;
  if (hpsNHlo) delete hpsNHlo;
  if (hangNHlo) delete hangNHlo;
  for (int i=0; i<4; i++) { // min/max for NH/IH
    for (int s=0; s<4; s++) { // sigmas
      if (contours[i][s]) delete contours[i][s];
      if (contours2[i][s]) delete contours2[i][s];
    }
  }
  if (f0) delete f0;
  if (f1) delete f1;
  if (c) delete c;
  
  return 0;
  
}

// *****************************************************************************
double meff(double m0, double& msum, double& mtrit, double alpha, double beta, 
            float sigma, int inverted) {
  
  // Throw neutrino mixing parameters
  double s12, s13, dm21, dm32;
  if (sigma >= 0) { // fixed sigma, flat priors
    s12 = S12 + sigma*(2*rng->Rndm()-1.)*ES12;
    s13 = S13 + sigma*(2*rng->Rndm()-1.)*ES13;
    dm21 = DM21 + sigma*(2*rng->Rndm()-1.)*EDM21;
    dm32;
    if (inverted) dm32 = DM32IH + sigma*(2*rng->Rndm()-1.)*EDM32IH;
    else          dm32 = DM32NH + sigma*(2*rng->Rndm()-1.)*EDM32NH;
  } else { // simulate Gaussian throws
    s12 = rng->Gaus(S12,ES12);
    s13 = rng->Gaus(S13,ES13);
    dm21 = rng->Gaus(DM21,EDM21);
    if (inverted) dm32 = rng->Gaus(DM32IH,EDM32IH);
    else          dm32 = rng->Gaus(DM32NH,EDM32NH);
  }
  double c12 = 1.-s12;
  double c13 = 1.-s13;
  
  // Define squared mass eigenstates
  double m1sq, m2sq, m3sq;
  if (inverted) {
    m1sq = m0*m0-dm32-dm21;   // dm21 > 0
    m2sq = m0*m0-dm32;        // dm32 < 0
    m3sq = m0*m0;             // m3 = m0
  } else {
    m1sq = m0*m0;             // m1 = m0
    m2sq = m0*m0+dm21;        // dm21 > 0
    m3sq = m0*m0+dm21+dm32;   // dm32 > 0
  }
  
  // Calculate sum of nu masses and effective nu mass from tritium decay
  msum = sqrt(m1sq) + sqrt(m2sq) + sqrt(m3sq);
  mtrit = sqrt(c12*c13*m1sq + s12*c13*m2sq + s13*m3sq);
  
  // Define complex numbers for 0vbb decay (in polar form)
  TComplex fac1(c12*c13*sqrt(m1sq),0.,true);
  TComplex fac2(s12*c13*sqrt(m2sq),2.*alpha,true);
  TComplex fac3(s13*sqrt(m3sq),2.*beta,true);
  
  // Add the three components and return absolute value
  TComplex res(fac1+fac2+fac3);
  return res.Abs(res);
}

// *****************************************************************************
double mefflimits(double m0, double& msumlim, float s, int index) {

  // Initialise very large limits
  double limit = 1e6;
  msumlim = 1e6;
  if (index%2) limit*=-1;
  else msumlim*=-1;
  if (index==0 && m0<2.3e-3) msumlim = 1e6;
  
  // Find min/max values for NH/IH contours
  double mass, msum, mtrit;
  for (int i=0; i<NSTATS; i++) {
    switch (index) {
      case 1 :  // NH max
        mass = meff(m0,msum,mtrit,0,0,s,0);
        if (mass>limit) limit=mass;
        if (msum<msumlim) msumlim=msum;
        break;
      case 2 :  // IH min
        mass = meff(m0,msum,mtrit,pi/2,pi/2,s,1);
        if (mass<limit) limit=mass;
        if (msum>msumlim) msumlim=msum;
        break;
      case 3 :  // IH max
        mass = meff(m0,msum,mtrit,0,0,s,1);
        if (mass>limit) limit=mass;
        if (msum<msumlim) msumlim=msum;
        break;
      default : // NH min
        if (m0<2.3e-3) {
          mass = meff(m0,msum,mtrit,pi/2,0,s,0);
          if (msum<msumlim) msumlim=msum;
        } else if (m0<6.6e-3) {
          mass = 0;
          msumlim = 6.5e-2;
        } else {
          mass = meff(m0,msum,mtrit,pi/2,pi/2,s,0);
          if (msum>msumlim) msumlim=msum;
        }
        if (mass<limit) limit=mass;
        break;
    }
  }
  
  return limit;
}

