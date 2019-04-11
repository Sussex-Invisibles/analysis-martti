// ---------------------------------------------------
//  Compile and run with: root -l -b -q 0vbb_phasespace.C+g
//  Author: Martti Nirkko, 22.08.2015 (updated 11.04.2019)
// ---------------------------------------------------
#include <TColor.h>
#include <TROOT.h>

// Include code by Dr. Xianguo Lu (Oxford) for logarithmically binned axes
#include "../Xianguo.C"

// Use official T2K plotting style
#include "/home/nirkko/pCloud/PhD/ND280/CommonStyle.H"

// global parameters and constants
const int NUMTHROWS = 1e7; // 1e6 throws take about 2 seconds
const int NBINS = 4e2;
const double pi = TMath::Pi();

// Neutrino mixing parameters (PDG 2018)
const double S12 = 0.307;
const double ES12 = 0.013;
const double DM21 = 7.53e-5;
const double EDM21 = 0.18e-5;
const double DM32 = 2.51e-3; // NH
const double EDM32 = 0.05e-3;
//const double DM32 = -2.56e-3; // IH
//const double EDM32 = 0.04e-3;
const double S13 = 2.12e-2;
const double ES13 = 0.08e-2;

// Get cos^2(theta) values
const double C12 = 1.-S12;
const double C13 = 1.-S13;

// global functions and generators
//double mbeta(double m0, int inverted=0);
double meff(double m0, double alpha, double beta, int inverted=0);
  
// Random number generator
TRandom3* duran = new TRandom3();

int main() {
  
  // Use T2K style
  CommonStyle();
  gStyle->SetPadBottomMargin(0.14);
  gStyle->SetPadLeftMargin(0.14);//to include both large/small font options
  gStyle->SetPadRightMargin(0.13);//to include both large/small font options
  gStyle->SetLabelSize(0.05,"xyz");
  gStyle->SetTitleSize(0.05,"xyz");
  gStyle->SetLabelOffset(0.015,"xyz");
  gStyle->SetTitleOffset(1.4,"x");
  gStyle->SetTitleOffset(1.3,"y");
  gStyle->SetMarkerStyle(7);
  gStyle->SetPalette(55);
  
  // lightest and effective neutrino masses
  double m0, meff1, meff2, meff3;
  
  // Histograms for neutrinoless double beta decay
  TH2D* hbb1 = new TH2D("hbb1","Normal hierarchy",NBINS,1e-4,1e0,NBINS,1e-4,1e0);
  TH2D* hbb2 = new TH2D("hbb2","Inverted hierarchy",NBINS,1e-4,1e0,NBINS,1e-4,1e0);
  BinLog(hbb1->GetXaxis()); BinLog(hbb1->GetYaxis());
  BinLog(hbb2->GetXaxis()); BinLog(hbb2->GetYaxis());
  
  // Histogram for Majorana angles resulting in small effective masses
  TH2D* hist3 = new TH2D("hist3","Small effective mass",NBINS,1e-4,1e0,NBINS,1e-8,1e-4);
  BinLog(hist3->GetXaxis()); BinLog(hist3->GetYaxis());
  
  TH2D* angles = new TH2D("angles","Majorana angles",NBINS,0,2,NBINS,0,2);
  double alpha0=0., beta0=0., m00=1.e6, meff0=1.e6;
  
  // Define Majorana angles
  double alpha = 0.;
  double beta = 0.;
  
  // Produce N throws for random values of m0 (with random Majorana phases)
  for (int i=0; i<NUMTHROWS; i++) {
    m0 = pow(10.,4.*(duran->Rndm()-1.)); // log prior
    alpha = 2.*pi*duran->Rndm();    // alpha_1
    beta  = 2.*pi*duran->Rndm();    // alpha_2 - 2*delta_CP
    meff1 = meff(m0, alpha, beta, 0);
    meff2 = meff(m0, alpha, beta, 1);
    hbb1->Fill(m0, meff1);
    hbb2->Fill(m0, meff2);
    if (meff1<meff0) { alpha0=alpha; beta0=beta; m00=m0; meff0=meff1; }
  }
  
  // Output values for minimal effective mass
  m00=1.e6; meff0=1.e6;
  int throws=0, count=0;
  while (count<NUMTHROWS/100) {
    m0 = pow(10.,(duran->Rndm()-3.)); // log prior (25% of phase space)
    alpha = pi*(1.+0.4*(duran->Rndm()-0.5)); // 0.8-1.2 pi (20% of phase space)
    beta = 2*pi*duran->Rndm(); // 0-2 pi (100% of phase space)
    meff3 = meff(m0,alpha,beta,0);
    throws++;
    if (meff3>1e-4) continue;
    angles->Fill(alpha/pi,beta/pi);
    hist3->Fill(m0, meff3); 
    count++;
    if (meff3<meff0) { alpha0=alpha; beta0=beta; m00=m0; meff0=meff3; }
  }
  angles->Scale(1./throws*0.25*0.2);
  hist3->Scale(1./throws*0.25*0.2);
  std::cout<<"MINIMUM: meff="<<1.e3*meff0<<" meV for m0="<<1.e3*m00<<" meV, alpha="<<alpha0/pi<<" pi, beta="<<beta0/pi<<" pi"<<std::endl;
  
  // Text in figures
  TLatex *Tl = new TLatex();
  Tl->SetTextSize(0.041);
  TBox *white = new TBox(1.3e-4, 1.2e-1, 9.0e-4, 8.0e-1);
  white->SetFillColor(0); white->SetLineColor(0);
  
  // -------------------
  //  PLOTTING SECTION
  // ------------------
  
  TCanvas *c1 = new TCanvas("c1","Phase space plot",1200,900);
  c1->SetLogx(); c1->SetLogy(); c1->SetLogz(); c1->SetGrid();
  //hbb1->SetTitle("#bf{0#nu#beta#beta decay - Type-I seesaw}");
  hbb1->GetXaxis()->SetTitle("m_{0} (eV)");
  hbb1->GetYaxis()->SetTitle("#LTm_{#beta#beta}#GT (eV)");
  //hbb1->GetXaxis()->SetTitleOffset(1.4);
  //hbb1->SetMarkerColorAlpha(kRed, 0.5);
  //hbb2->SetMarkerColorAlpha(kGreen, 0.5);
  //hbb1->Draw("scat");
  //hbb2->Draw("scat same");
  hbb1->Scale(1./NUMTHROWS);
  hbb1->Draw("colz");
  Tl->DrawLatex(1.5e-4,5.6e-1,"#alpha #equiv #alpha_{1}");
  Tl->DrawLatex(1.5e-4,2.8e-1,"#beta #equiv #alpha_{2}#font[122]{-} 2#delta_{CP}");
  Tl->DrawLatex(1.5e-4,1.4e-1,"#alpha, #beta #in [0, 2#pi]");
  c1->Print("0vbb_phasespace.png");
  c1->Print("0vbb_phasespace.pdf");
  c1->Close();
  
  TCanvas *c3 = new TCanvas("c3","Phase space plot",1200,900);
  c3->SetLogx(); c3->SetLogy(); c3->SetLogz(); c3->SetGrid();
  //hist3->SetTitle("#bf{0#nu#beta#beta decay - allowed phase space}");
  hist3->GetXaxis()->SetTitle("lightest #nu mass m_{0} (eV)");
  hist3->GetYaxis()->SetTitle("effective #nu_{e} mass m_{#beta#beta} (eV)");
  hist3->Draw("colz");
  Tl->DrawLatex(2e-4,0.4,"#alpha #equiv #alpha_{1} #in [0, 2#pi]");
  Tl->DrawLatex(2e-4,0.2,"#beta #equiv #alpha_{2}#font[122]{-} 2#delta_{CP} #in [0, 2#pi]");
  c3->Print("0vbb_phasespace_lowmass.png");
  c3->Print("0vbb_phasespace_lowmass.pdf");
  c3->Close();
  
  gStyle->SetOptStat(1111); gStyle->SetStatX(0.885); gStyle->SetStatY(0.88);
  TCanvas *c2 = new TCanvas("c2","Majorana angle plot",1200,900);
  c2->SetLogz(); c2->SetGrid();
  //angles->SetTitle("#bf{Majorana angles that lead to m_{#beta#beta} < 10^{-4} eV}");
  angles->GetXaxis()->SetTitle("#alpha #equiv #alpha_{1} [#pi]");
  angles->GetYaxis()->SetTitle("#beta #equiv #alpha_{2}#font[122]{-} 2#delta_{CP} [#pi]");
  angles->Draw("colz");
  c2->Print("0vbb_phasespace_lowmass_2.png");
  c2->Print("0vbb_phasespace_lowmass_2.pdf");
  c2->Close();
  
  return 0;
}

double meff(double m0, double alpha, double beta, int inverted) {

  // Throw neutrino mixing parameters
  double s12 = S12 + (2*duran->Rndm()-1.)*ES12;
  double c12 = 1.-s12;
  double s13 = S13 + (2*duran->Rndm()-1.)*ES13;
  double c13 = 1.-s13;
  double dm21 = DM21 + (2*duran->Rndm()-1.)*EDM21;
  double dm32 = DM32 + (2*duran->Rndm()-1.)*EDM32;
  
  // Define complex numbers in polar form
  TComplex *fac1 = new TComplex(1,0,true);
  TComplex *fac2 = new TComplex(1,alpha,true);
  TComplex *fac3 = new TComplex(1,beta,true);
  
  // Multiply by PMNS matrix element and neutrino mass eigenstate
  if (inverted) {
    *fac1 *= c12*c13*sqrt(m0*m0+dm32-dm21);
    *fac2 *= s12*c13*sqrt(m0*m0+dm32);
    *fac3 *= s13*m0;      // m3 = m0
  } else {
    *fac1 *= c12*c13*m0;  // m1 = m0
    *fac2 *= s12*c13*sqrt(m0*m0+dm21);
    *fac3 *= s13*sqrt(m0*m0+dm21+dm32);
  }
  
  // Add the three components and return absolute value
  TComplex *res = new TComplex(*fac1 + *fac2 + *fac3);
  return res->Abs(*res);

}

