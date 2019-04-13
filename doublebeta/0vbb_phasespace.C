// ---------------------------------------------------
//  Compile and run with: root -l -b -q 0vbb_phasespace.C+g
//  Author: Martti Nirkko, 22.08.2015 (updated 11.04.2019)
// ---------------------------------------------------
#include <TColor.h>
#include <TROOT.h>

// Include helper functions
#include "../HelperFunc.C"

// Use official T2K plotting style
#include "../CommonStyle.H"

// global parameters and constants
const int NSTATS = 2e4; // statistics for limits
const int NTHROWS = 1e8;
const int NTHROWSLOWMASS = 1e5;
const int NBINS = 4e2;

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

// Get cos^2(theta) values
const double C12 = 1.-S12;
const double C13 = 1.-S13;

// Global functions and generators
double mefflimits(double m0, float sigma=1., int index=0);
double meff(double m0, double alpha, double beta, float sigma=0., int inverted=0);
  
// Random number generator
TRandom3* duran = new TRandom3();

int main() {
  
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
  gStyle->SetPalette(55);
  
  // Random number seed
  duran->SetSeed(0);
  
  // lightest and effective neutrino masses
  double m0, meffNH, meffIH;
  
  // Histograms for normal and inverted hierarchy
  TH2D* hpsNH = new TH2D("hpsNH","Normal hierarchy",NBINS,1e-4,1e0,NBINS,1e-4,1e0);
  BinLog(hpsNH->GetXaxis()); BinLog(hpsNH->GetYaxis());
  TH2D* hpsIH = new TH2D("hpsIH","Inverted hierarchy",NBINS,1e-4,1e0,NBINS,1e-4,1e0);
  BinLog(hpsIH->GetXaxis()); BinLog(hpsIH->GetYaxis());
  
  // Histograms for very small effective masses (NH)
  TH2D* hpsNHlo = new TH2D("hpsNHlo","Small effective mass",NBINS,1e-4,1e0,NBINS,1e-8,1e-4);
  BinLog(hpsNHlo->GetXaxis()); BinLog(hpsNHlo->GetYaxis());
  TH2D* hangNHlo = new TH2D("hangNHlo","Majorana hangNHlo",NBINS,0.8,1.2,NBINS,0,2);
  
  // Text in figures
  TLatex *Tl = new TLatex();
  Tl->SetTextSize(0.041);
  
  // Get maximum lines for allowed phase space region
  int COL = 1;
  int WDT = 1.5;
  //int sty[4] = {1,10,6,4};
  TGraph *contours[4][4] = {NULL};
  std::cout << "Generating countour limits..." << std::endl;
  for (int i=0; i<4; i++) { // min/max for NH/IH
    for (int s=0; s<4; s++) { // sigmas
      contours[i][s] = new TGraph();
      for (int b=0; b<=NBINS; b++) {
        m0 = hpsNH->GetXaxis()->GetBinCenter(b);
        contours[i][s]->SetPoint(b,m0,mefflimits(m0,s,i));
        printProgress(i*4*NBINS+s*NBINS+b,16*NBINS);
      }
      contours[i][s]->SetLineColor(COL);
      contours[i][s]->SetLineWidth(WDT);
      //contours[i][s]->SetLineStyle(sty[s]);
    }
  }
  
  // Define parameters
  double alpha0=0., beta0=0., m0min=1.e6, meffmin=1.e6;
  double alpha, beta;
  
  // Generate contour plots
  long throws=0, count=0;
  std::cout << "Filling allowed phase space regions..." << std::endl;
  while (throws<NTHROWS) {
    m0 = pow(10.,4.*(duran->Rndm()-1.)); // log prior (100% of phase space)
    alpha = 2*pi*duran->Rndm(); // 0-2 pi (100% of phase space)
    beta = 2*pi*duran->Rndm(); // 0-2 pi (100% of phase space)
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
  
  TCanvas *c1 = new TCanvas("c1","Phase space plot",1200,900);
  c1->SetLogx(); c1->SetLogy(); c1->SetLogz(); c1->SetGrid();
  //hpsIH->SetTitle("#bf{0#nu#beta#beta decay - Type-I seesaw}");
  hpsIH->GetXaxis()->SetTitle("m_{0} (eV)");
  hpsIH->GetYaxis()->SetTitle("#LTm_{#beta#beta}#GT (eV)");
  hpsIH->Draw("colz");
  hpsIH->GetZaxis()->SetRangeUser(2e-9,2e-4);
  for (int s=0; s<4; s++) { // contours for IH
    contours[2][s]->Draw("c");
    contours[3][s]->Draw("c");
  }
  hpsNH->Draw("colz same");
  hpsNH->GetZaxis()->SetRangeUser(2e-9,2e-4);
  for (int s=0; s<4; s++) { // contours for NH
    contours[0][s]->Draw("c");
    contours[1][s]->Draw("c");
  }
  Tl->DrawLatex(1.4e-4,5.6e-1,"#alpha #equiv #alpha_{1}");
  Tl->DrawLatex(1.4e-4,2.8e-1,"#beta #equiv #alpha_{2}#font[122]{-} 2#delta_{CP}");
  Tl->DrawLatex(1.4e-4,1.4e-1,"#alpha, #beta #in [0, 2#pi]");
  Tl->DrawLatex(1.4e-4,2.6e-2,"IH");
  Tl->DrawLatex(1.4e-4,2.0e-3,"NH");
  c1->Print("0vbb_phasespace.png");
  c1->Print("0vbb_phasespace.pdf");
  c1->Close();
  
  if (hpsNH) delete hpsNH;
  if (hpsIH) delete hpsIH;
  
  long throws2=0, count2=0;
  std::cout << "Filling phase space region for low effective mass..." << std::endl;
  while (count2<NTHROWSLOWMASS) {
    m0 = pow(10.,(duran->Rndm()-3.)); // log prior (25% of phase space)
    beta = 2*pi*duran->Rndm(); // 0-2 pi (100% of phase space)
    alpha = pi*(1.+0.13*sin(beta)+0.1*(duran->Rndm()-0.5)); // 5% of phase space
    //alpha = pi*(1.+0.4*(duran->Rndm()-0.5)); // 0.8-1.2 pi (20% of phase space)
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
  
  TCanvas *c2 = new TCanvas("c2","Phase space plot",1200,900);
  c2->SetLogx(); c2->SetLogy(); c2->SetLogz(); c2->SetGrid();
  //hpsNHlo->SetTitle("#bf{0#nu#beta#beta decay - allowed phase space}");
  hpsNHlo->GetXaxis()->SetTitle("m_{0} (eV)");
  hpsNHlo->GetYaxis()->SetTitle("#LTm_{#beta#beta}#GT (eV)");
  hpsNHlo->Draw("colz");
  hpsNHlo->GetZaxis()->SetRangeUser(5e-10,2e-7);
  Tl->DrawLatex(2e-4,0.4,"#alpha #equiv #alpha_{1} #in [0, 2#pi]");
  Tl->DrawLatex(2e-4,0.2,"#beta #equiv #alpha_{2}#font[122]{-} 2#delta_{CP} #in [0, 2#pi]");
  /*
  // Contours are too jaggy (stats limited) here
  for (int i=0; i<2; i++) { // min/max for NH only
    for (int s=0; s<4; s++) { // sigmas
      contours[i][s]->Draw("c");
    }
  }
  */
  c2->Print("0vbb_phasespace_lowmass.png");
  c2->Print("0vbb_phasespace_lowmass.pdf");
  c2->Close();
  
  gStyle->SetOptStat(1111); gStyle->SetStatX(0.885); gStyle->SetStatY(0.88);
  TCanvas *c3 = new TCanvas("c3","Majorana angle plot",1200,900);
  c3->SetLogz(); c3->SetGrid();
  //hangNHlo->SetTitle("#bf{Majorana hangNHlo that lead to m_{#beta#beta} < 10^{-4} eV}");
  hangNHlo->GetXaxis()->SetTitle("#alpha #equiv #alpha_{1} [#pi]");
  hangNHlo->GetYaxis()->SetTitle("#beta #equiv #alpha_{2}#font[122]{-} 2#delta_{CP} [#pi]");
  hangNHlo->Draw("colz");
  hangNHlo->GetXaxis()->SetRangeUser(0.8,1.2);
  hangNHlo->GetZaxis()->SetRangeUser(8e-10,4e-8);
  
  TGraph *f1 = new TGraph();
  TGraph *f0 = new TGraph();
  for (int i=0; i<NBINS; i++) {
    beta = 2.*i/NBINS;
    alpha = 1. + 0.13*sin(beta*pi);
    f1->SetPoint(i,alpha+0.05,beta);
    f0->SetPoint(i,alpha-0.05,beta);
  }
  f1->SetLineColor(1);
  f1->Draw("l");
  f0->SetLineColor(1);
  f0->Draw("l");
  
  c3->Print("0vbb_phasespace_lowmass_2.png");
  c3->Print("0vbb_phasespace_lowmass_2.pdf");
  c3->Close();
  
  if (hpsNHlo) delete hpsNHlo;
  if (hangNHlo) delete hangNHlo;
  
  if (c1) delete c1;
  if (c2) delete c2;
  if (c3) delete c3;
  
  return 0;
}

double meff(double m0, double alpha, double beta, float sigma, int inverted) {

  // Throw neutrino mixing parameters
  double s12, s13, dm21, dm32;
  if (sigma >= 0) { // fixed sigma, flat priors
    s12 = S12 + sigma*(2*duran->Rndm()-1.)*ES12;
    s13 = S13 + sigma*(2*duran->Rndm()-1.)*ES13;
    dm21 = DM21 + sigma*(2*duran->Rndm()-1.)*EDM21;
    dm32;
    if (inverted) dm32 = DM32IH + sigma*(2*duran->Rndm()-1.)*EDM32IH;
    else          dm32 = DM32NH + sigma*(2*duran->Rndm()-1.)*EDM32NH;
  } else { // simulate Gaussian throws
    s12 = duran->Gaus(S12,ES12);
    s13 = duran->Gaus(S13,ES13);
    dm21 = duran->Gaus(DM21,EDM21);
    if (inverted) dm32 = duran->Gaus(DM32IH,EDM32IH);
    else          dm32 = duran->Gaus(DM32NH,EDM32NH);
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
        if      (m0<2.5e-3) mass = meff(m0,pi,0,s,0);
        else if (m0<6.5e-3) mass = 0;
        else                mass = meff(m0,pi,pi,s,0);
        if (mass<limit) limit=mass;
        break;
    }
  }
  
  return limit;
}

