// ---------------------------------------------------
//  Compile and run with: root -l -b -q 0vbb_phasespace.C+g
//  Author: Martti Nirkko, 22.08.2015 (updated 11.04.2019)
// ---------------------------------------------------
// Include helper functions / plotting style
#include "../include/HelperFunc.C"
#include "../include/CommonStyle.H"

// Global constants
const int NSTATS = 1e4;           // statistics for contour limits
const int NTHROWS = 1e8;          // full phase space
const int NTHROWSLOWMASS = 1e5;   // low m_eff phase space
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
double mefflimits(double m0, float sigma=1., int index=0);
double meff(double m0, double alpha, double beta, float sigma=0., int inverted=0);
  
// Random number generator
TRandom3* rng = new TRandom3();

int main() {

  // I/O file (either plot from ROOT file, or generate it)
  std::string fname = "0vbb_phasespace";
  std::string output = fname+".root";
  ifstream io(output.c_str());
  
  // Initialise objects
  TFile *file = NULL;
  TH2D *hpsNH=NULL, *hpsIH=NULL, *hpsNHlo=NULL, *hangNHlo=NULL;
  TGraph *contours[4][4] = {NULL};
  TGraph *f0=NULL, *f1=NULL;
  TCanvas *c1=NULL, *c2=NULL, *c3=NULL;
  
  if (!io.good()) {  // File does not exist
    
    // Create ROOT output file
    file = new TFile(output.c_str(),"RECREATE");
    
    // Histograms for normal and inverted hierarchy
    TH2D* hpsNH = new TH2D("hpsNH","Normal hierarchy",NBINS,1e-4,1e0,NBINS,1e-4,1e0);
    BinLog(hpsNH->GetXaxis()); BinLog(hpsNH->GetYaxis());
    TH2D* hpsIH = new TH2D("hpsIH","Inverted hierarchy",NBINS,1e-4,1e0,NBINS,1e-4,1e0);
    BinLog(hpsIH->GetXaxis()); BinLog(hpsIH->GetYaxis());
    
    // Histograms for very small effective masses (NH)
    TH2D* hpsNHlo = new TH2D("hpsNHlo","Small effective mass",NBINS,1e-4,1e0,NBINS,1e-8,1e-4);
    BinLog(hpsNHlo->GetXaxis()); BinLog(hpsNHlo->GetYaxis());
    TH2D* hangNHlo = new TH2D("hangNHlo","Majorana hangNHlo",NBINS,0.8,1.2,NBINS,0,2);
    
    // lightest and effective neutrino masses
    double m0, meffNH, meffIH;
    
    // Get maximum lines for allowed phase space region
    std::cout << "Generating countour limits..." << std::endl;
    std::string cname;
    for (int i=0; i<4; i++) { // min/max for NH/IH
      for (int s=0; s<4; s++) { // sigmas
        contours[i][s] = new TGraph();
        for (int b=0; b<=NBINS; b++) {
          m0 = hpsNH->GetXaxis()->GetBinCenter(b);
          contours[i][s]->SetPoint(b,m0,mefflimits(m0,s,i));
          printProgress(i*4*NBINS+s*NBINS+b,16*NBINS);
        }
        cname = Form("contours_%d_%d",i,s);
        contours[i][s]->SetName(cname.c_str());
      }
    }
    
    // Majorana angles and global minima
    double alpha, beta;
    double alpha0=0., beta0=0., m0min=1.e6, meffmin=1.e6;
  
    // Random number seed
    rng->SetSeed(0);
    
    // Generate contour plots
    long throws=0, count=0;
    std::cout << "Filling allowed phase space regions..." << std::endl;
    while (throws<NTHROWS) {
      m0 = pow(10.,4.*(rng->Rndm()-1.)); // log prior (100% of phase space)
      alpha = 2*pi*rng->Rndm(); // 0-2 pi (100% of phase space)
      beta = 2*pi*rng->Rndm(); // 0-2 pi (100% of phase space)
      meffNH = meff(m0,alpha,beta,-1,0);
      meffIH = meff(m0,alpha,beta,-1,1);
      throws++;
      hpsNH->Fill(m0, meffNH);
      hpsIH->Fill(m0, meffIH);
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
      double beta = 2.*i/NBINS;
      double alpha = 1. + 0.13*sin(beta*pi);  // 5% of phase space
      f1->SetPoint(i,alpha+0.05,beta);
      f0->SetPoint(i,alpha-0.05,beta);
    }
    f0->SetName("alpha_lo");
    f1->SetName("alpha_hi");
    while (count2<NTHROWSLOWMASS) {
      m0 = pow(10.,(rng->Rndm()-3.)); // log prior (25% of phase space)
      beta = 2*pi*rng->Rndm(); // 0-2 pi (100% of phase space)
      alpha = pi*(1.+0.13*sin(beta)+0.1*(rng->Rndm()-0.5)); // 5% of phase space
      //alpha = pi*(1.+0.4*(rng->Rndm()-0.5)); // 0.8-1.2 pi (20% of phase space)
      meffNH = meff(m0,alpha,beta,-1,0);
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
      }
    }
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
    hpsIH = (TH2D*)file->Get("hpsIH");
    hpsNH = (TH2D*)file->Get("hpsNH");
    hpsNHlo = (TH2D*)file->Get("hpsNHlo");
    hangNHlo = (TH2D*)file->Get("hangNHlo");
    std::string cname;
    for (int i=0; i<4; i++) { // min/max for NH/IH
      for (int s=0; s<4; s++) { // sigmas
        cname = Form("contours_%d_%d",i,s);
        TGraph *temp = (TGraph*)file->Get(cname.c_str());
        contours[i][s] = new TGraph();
        Smooth(contours[i][s],temp,5);
        if(temp) delete temp;
      }
    }
    f0 = (TGraph*)file->Get("alpha_lo");
    f1 = (TGraph*)file->Get("alpha_hi");
    
    // Text in figures
    TLatex *Tl = new TLatex();
    Tl->SetTextSize(0.041);
    
    // SNO+ sensitivity
    double fac=10./1.9; // phase II /I lifetime sensitivity (1e26y)
    double sx1[4] = {1e-4,1e0,1e0,1e-4};
    double sy1[4] = {99e-3,99e-3,41e-3,41e-3}; // meV (phase I)
    double sy2[4] = {99e-3/fac,99e-3/fac,41e-3/fac,41e-3/fac}; // meV (phase II)
    TGraph *gs1 = new TGraph(4,sx1,sy1);
    TGraph *gs2 = new TGraph(4,sx1,sy2);
    
    // Plot full phase space for double beta decay (NH+IH)
    // ---------------------------------------------------
    c1 = new TCanvas("c1","Phase space plot",1200,900);
    c1->SetLogx(); c1->SetLogy(); c1->SetLogz(); c1->SetGrid();
    //hpsIH->SetTitle("#bf{0#nu#beta#beta decay - Type-I seesaw}");
    c1->DrawFrame(1e-4,1e-4,1e0,1e0,";m_{0} (eV);#LTm_{#beta#beta}#GT (eV)");
    // Draw SNO+ sensitivity bands (areas only)
    gs1->SetFillColorAlpha(BANDCOL,0.5);
    gs2->SetFillColorAlpha(BANDCOL,0.5);
    gs1->Draw("f");
    gs2->Draw("f");
    // Draw allowed phase space for IH
    hpsIH->Draw("colz same");
    hpsIH->GetZaxis()->SetRangeUser(2e-9,2e-4);
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
    hpsNH->GetZaxis()->SetRangeUser(2e-9,2e-4);
    for (int i=0; i<2; i++) { // contours for NH
      for (int s=0; s<4; s++) { // contours for NH
        //contours[i][s]->SetLineStyle(CONTSTY[s]);
        contours[i][s]->SetLineColor(CONTCOL+s);
        contours[i][s]->SetLineWidth(CONTWDT);
        contours[i][s]->Draw("l");
      }
    }
    // Draw legend for contours
    TLegend *leg = new TLegend(0.16,0.16,0.33,0.32,"#bf{PDG 2018}");
    for (int s=0; s<4; s++) {
      std::string legstr = Form("%d #sigma",s);
      leg->AddEntry(contours[0][s],legstr.c_str(),"l");
    }
    leg->SetFillStyle(0);
    leg->Draw();
    // Draw text describing plots
    Tl->SetTextColor(DEFCOL);
    Tl->SetTextSize(0.03);
    Tl->DrawLatex(1.5e-4,6.4e-1,"#alpha #equiv #alpha_{1}");
    Tl->DrawLatex(1.5e-4,4.0e-1,"#beta #equiv #alpha_{2}#font[122]{-} 2#delta_{CP}");
    Tl->DrawLatex(1.5e-4,2.5e-1,"#alpha, #beta #in [0, 2#pi]");
    Tl->SetTextSize(0.035);
    Tl->SetTextColor(HIERCOL);
    Tl->DrawLatex(1.5e-4,2.5e-2,"IH");
    Tl->DrawLatex(1.5e-4,2.0e-3,"NH");
    // Evince seems to have trouble displaying uppercase delta (font not embedded)
    //Tl->DrawLatex(1.5e-4,2.5e-2,"#Deltam_{32}^{2} < 0");
    //Tl->DrawLatex(1.5e-4,2.0e-3,"#Deltam_{32}^{2} > 0");
    Tl->SetTextSize(0.025);
    Tl->SetTextColor(TEXTCOL);
    Tl->DrawLatex(2e-4,6.0e-2,"SNO+ Phase I (T_{1/2} > 1.9 #times 10^{26} y)");
    Tl->DrawLatex(2e-4,1.0e-2,"SNO+ Phase II (T_{1/2} > 10^{27} y)");
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
    c1->Update();
    c1->RedrawAxis();
    c1->Print("0vbb_phasespace.png");
    c1->Print("0vbb_phasespace.pdf");
    c1->Close();
    
    // Plot low m_eff region for double beta decay (NH)
    // ------------------------------------------------
    c2 = new TCanvas("c2","Phase space plot",1200,900);
    c2->SetLogx(); c2->SetLogy(); c2->SetLogz(); c2->SetGrid();
    //hpsNHlo->SetTitle("#bf{0#nu#beta#beta decay - allowed phase space}");
    c2->DrawFrame(1e-4,1e-7,1e0,1e-4,";m_{0} (eV);#LTm_{#beta#beta}#GT (eV)");
    hpsNHlo->Draw("colz same");
    hpsNHlo->GetZaxis()->SetRangeUser(5e-10,2e-7);
    // Contours are rather jaggy (stats limited) here
    for (int s=0; s<4; s++) { // contours for NH
      contours[0][s]->Draw("l");
      contours[1][s]->Draw("l");
    }
    // Draw text describing plots
    Tl->SetTextSize(0.03);
    Tl->SetTextColor(1);
    Tl->DrawLatex(4e-2,7.00e-5,"#alpha #equiv #alpha_{1}");
    Tl->DrawLatex(4e-2,5.00e-5,"#beta #equiv #alpha_{2}#font[122]{-} 2#delta_{CP}");
    Tl->DrawLatex(4e-2,3.57e-5,"#beta #in [0, 2]");
    Tl->DrawLatex(4e-2,2.55e-5,"#alpha #in [1 #plus 0.13 sin#beta #minus 0.05,");
    Tl->DrawLatex(6.8e-2,1.82e-5,"1 #plus 0.13 sin#beta #plus 0.05]");
    // Draw legend (defined above)
    leg->Draw();
    // Save and close
    c2->Print("0vbb_phasespace_lowmass.png");
    c2->Print("0vbb_phasespace_lowmass.pdf");
    c2->Close();
    
    // Plot Majorana angles for low m_eff region (NH)
    // ----------------------------------------------
    c3 = new TCanvas("c3","Majorana angle plot",1200,900);
    c3->SetLogz(); c3->SetGrid();
    //hangNHlo->SetTitle("#bf{Majorana hangNHlo that lead to m_{#beta#beta} < 10^{-4} eV}");
    c3->DrawFrame(0.8,0,1.2,2,";#alpha #equiv #alpha_{1} [#pi];#beta #equiv #alpha_{2}#font[122]{-} 2#delta_{CP} [#pi]");
    hangNHlo->Draw("colz same");
    hangNHlo->GetZaxis()->SetRangeUser(8e-10,4e-8);
    // Draw region in which points were generated
    f0->SetLineColor(LINECOL);
    f1->SetLineColor(LINECOL);
    f0->Draw("l");
    f1->Draw("l");
    // Draw text describing plots
    Tl->SetTextSize(0.03);
    Tl->SetTextColor(1);
    Tl->DrawLatex(1.05,1.72,"#beta #in [0, 2]");
    Tl->DrawLatex(1.05,1.62,"#alpha #in [1 #plus 0.13 sin#beta #minus 0.05,");
    Tl->DrawLatex(1.073,1.52,"1 #plus 0.13 sin#beta #plus 0.05]");
    // Save and close
    c3->Print("0vbb_phasespace_lowmass_2.png");
    c3->Print("0vbb_phasespace_lowmass_2.pdf");
    c3->Close();
    
  }
  
  // Free memory
  if (hpsNH) delete hpsNH;
  if (hpsIH) delete hpsIH;
  if (hpsNHlo) delete hpsNHlo;
  if (hangNHlo) delete hangNHlo;
  for (int i=0; i<4; i++) { // min/max for NH/IH
    for (int s=0; s<4; s++) { // sigmas
      if (contours[i][s]) delete contours[i][s];
    }
  }
  if (f0) delete f0;
  if (f1) delete f1;
  if (c1) delete c1;
  if (c2) delete c2;
  if (c3) delete c3;
  
  return 0;
  
}

double meff(double m0, double alpha, double beta, float sigma, int inverted) {

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
  
  // Define complex numbers in polar form
  TComplex *fac1 = new TComplex(1,0,true);
  TComplex *fac3 = new TComplex(1,alpha,true);
  TComplex *fac2 = new TComplex(1,beta,true);
  
  // Multiply by PMNS matrix element and neutrino mass eigenstate
  if (inverted) {
    *fac1 *= c12*c13*sqrt(m0*m0-dm32-dm21);   // dm21>0
    *fac3 *= s12*c13*sqrt(m0*m0-dm32);        // dm32<0
    *fac2 *= s13*m0;                          // m3 = m0
  } else {
    *fac1 *= c12*c13*m0;                      // m1 = m0
    *fac3 *= s12*c13*sqrt(m0*m0+dm21);        // dm21>0
    *fac2 *= s13*sqrt(m0*m0+dm21+dm32);       // dm32>0
  }
  
  // Add the three components and return absolute value
  TComplex *res = new TComplex(*fac1 + *fac3 + *fac2);
  return res->Abs(*res);
}

double mefflimits(double m0, float s, int index) {

  // Initialise very large limits
  double limit = 1e6;
  if (index%2) limit*=-1;
  
  // Find min/max values for NH/IH contours
  double mass;
  for (int i=0; i<NSTATS; i++) {
    switch (index) {
      case 1 :  // NH max
        mass = meff(m0,0,0,s,0);
        if (mass>limit) limit=mass;
        break;
      case 2 :  // IH min
        mass = meff(m0,pi,pi,s,1);
        if (mass<limit) limit=mass;
        break;
      case 3 :  // IH max
        mass = meff(m0,0,0,s,1);
        if (mass>limit) limit=mass;
        break;
      default : // NH min
        if      (m0<2.3e-3) mass = meff(m0,pi,0,s,0);
        else if (m0<6.6e-3) mass = 0;
        else                mass = meff(m0,pi,pi,s,0);
        if (mass<limit) limit=mass;
        break;
    }
  }
  
  return limit;
}

