// ---------------------------------------------------
//  Compile and run with: root -l -b -q 0vbb_phasespace.C+g
//  Author: Martti Nirkko, 22.08.2015 (updated 11.04.2019)
// ---------------------------------------------------
#include <TColor.h>
#include <TROOT.h>

// Include code by Dr. Xianguo Lu (Oxford) for logarithmically binned axes
#include "../HelperFunc.C"

// Use official T2K plotting style
#include "/home/nirkko/pCloud/PhD/ND280/CommonStyle.H"

// global parameters and constants
const int NTHROWS = 1e8;
const int NTHROWSLOWMASS = 1e5;
const int NSTATS = 1e3; // statistics for each bin
const int NBINS = 4e2;
const int NBINSCOARSE = NBINS/100;

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
  gStyle->SetPadRightMargin(0.14);//to include both large/small font options
  gStyle->SetLabelSize(0.05,"xyz");
  gStyle->SetTitleSize(0.05,"xyz");
  gStyle->SetLabelOffset(0.01,"xyz");
  gStyle->SetTitleOffset(1.4,"x");
  gStyle->SetTitleOffset(1.3,"y");
  gStyle->SetMarkerStyle(7);
  gStyle->SetPalette(55);
  
  // lightest and effective neutrino masses
  double m0, meffIH, meff2, meffNH;
  
  /*
  // Histograms for neutrinoless double beta decay
  TH2D* hbb1 = new TH2D("hbb1","Normal hierarchy",NBINS,1e-4,1e0,NBINS,1e-4,1e0);
  TH2D* hbb2 = new TH2D("hbb2","Inverted hierarchy",NBINS,1e-4,1e0,NBINS,1e-4,1e0);
  BinLog(hbb1->GetXaxis()); BinLog(hbb1->GetYaxis());
  BinLog(hbb2->GetXaxis()); BinLog(hbb2->GetYaxis());
  */
  
  TH2D* hpsNH = new TH2D("hpsNH","Normal hierarchy",NBINS,1e-4,1e0,NBINS,1e-4,1e0);
  BinLog(hpsNH->GetXaxis()); BinLog(hpsNH->GetYaxis());
  TH2D* hpsIH = new TH2D("hpsIH","Inverted hierarchy",NBINS,1e-4,1e0,NBINS,1e-4,1e0);
  BinLog(hpsIH->GetXaxis()); BinLog(hpsIH->GetYaxis());
  
  TH2D* htemp = new TH2D("htemp","Temporary",NBINS,1e-4,1e0,NBINS,1e-4,1e0);
  BinLog(htemp->GetXaxis()); BinLog(htemp->GetYaxis());
  
  TH2I* hcnNH = new TH2I("hcnNH","Counter",NBINSCOARSE,1e-4,1e0,NBINSCOARSE,1e-4,1e0);
  BinLog(hcnNH->GetXaxis()); BinLog(hcnNH->GetYaxis());
  double minalpha[NBINSCOARSE][NBINSCOARSE] = {1e6};
  double maxalpha[NBINSCOARSE][NBINSCOARSE] = {-1e6};
  double minbeta[NBINSCOARSE][NBINSCOARSE] = {1e6};
  double maxbeta[NBINSCOARSE][NBINSCOARSE] = {-1e6};
  
  // Histograms for very small effective masses (NH)
  TH2D* hpsNHlo = new TH2D("hpsNHlo","Small effective mass",NBINS,1e-4,1e0,NBINS,1e-8,1e-4);
  BinLog(hpsNHlo->GetXaxis()); BinLog(hpsNHlo->GetYaxis());
  TH2D* hangNHlo = new TH2D("hangNHlo","Majorana hangNHlo",NBINS,0.8,1.2,NBINS,0,2);
  
  // Define parameters
  double alpha0=0., beta0=0., m0min=1.e6, meffmin=1.e6;
  double alpha = 0.;
  double beta = 0.;
  long throws, count;
  
  // Generate inverted hierarchy plot (easy)
  throws=0; count=0;
  while (count<NTHROWS) {
    m0 = pow(10.,4.*(duran->Rndm()-1.)); // log prior (100% of phase space)
    alpha = 2*pi*duran->Rndm(); // 0-2 pi (100% of phase space)
    beta = 2*pi*duran->Rndm(); // 0-2 pi (100% of phase space)
    meffNH = meff(m0,alpha,beta,0);
    meffIH = meff(m0,alpha,beta,1);
    throws++;
    hpsNH->Fill(m0, meffNH);
    hpsIH->Fill(m0, meffIH);
    count++;
    printProgress(count,NTHROWS);
    
    // Fill coarse NH plot too (useful later?)
    hcnNH->Fill(m0, meffNH);
    int i = hcnNH->GetXaxis()->FindBin(m0);
    int j = hcnNH->GetYaxis()->FindBin(meffNH);
    if (alpha < minalpha[i][j]) minalpha[i][j] = alpha;
    if (alpha > maxalpha[i][j]) maxalpha[i][j] = alpha;
    if (beta < minbeta[i][j]) minbeta[i][j] = beta;
    if (beta > maxbeta[i][j]) maxbeta[i][j] = beta;
  }
  hpsNH->Scale(1./throws);
  hpsIH->Scale(1./throws);
  std::cout<<"Generated " << (float)throws << " throws, of which " << hpsIH->GetEntries();
  std::cout<<" (" << 100.*hpsIH->GetEntries()/throws << "%) were used for histogram." << std::endl;
  
  /*
  // How many bins in course histogram have entries?
  int nbinsx=hcnNH->GetXaxis()->GetNbins();
  int nbinsy=hcnNH->GetYaxis()->GetNbins();
  int nbinsPS=0;
  for (int i=0; i<=nbinsx; i++) {
    for (int j=0; j<=nbinsy; j++) {
      if (hcnNH->GetBinContent(i,j)==0) continue;
      if (minalpha[i][j]<0) minalpha[i][j]=0;
      if (maxalpha[i][j]>2*pi) maxalpha[i][j]=2*pi;
      if (minbeta[i][j]<0) minbeta[i][j]=0;
      if (maxbeta[i][j]>2*pi) maxbeta[i][j]=2*pi;
      nbinsPS++;
    }
  }
  
  // Generate normal hierarchy plot (tricky since we want it to be smooth)
  throws=0; count=0;
  long sumthrows=0, sumentries=0, sumcount=0;
  double xmin, ymin, xmax, ymax, norm;
  for (int i=0; i<=nbinsx; i++) {
    for (int j=0; j<=nbinsy; j++) {
      if (hcnNH->GetBinContent(i,j)==0) continue;
      xmin = hcnNH->GetXaxis()->GetBinLowEdge(i);
      xmax = hcnNH->GetXaxis()->GetBinUpEdge(i);
      ymin = hcnNH->GetYaxis()->GetBinLowEdge(j);
      ymax = hcnNH->GetYaxis()->GetBinUpEdge(j);
      throws=0;
      count=0;
      htemp->Reset();
      while (count<NSTATS) {
        m0 = xmin*pow(xmax/xmin,duran->Rndm()); // log prior (TODO - check)
        alpha = minalpha[i][j]+(maxalpha[i][j]-minalpha[i][j])*duran->Rndm(); // restrict phase space
        beta = minbeta[i][j]+(maxbeta[i][j]-minbeta[i][j])*duran->Rndm(); // restrict phase space
        meffNH = meff(m0,alpha,beta,0);
        throws++;
        htemp->Fill(m0, meffNH);
        if (m0<xmin) continue;
        if (m0>=xmax) continue;
        if (meffNH<ymin) continue;
        if (meffNH>=ymax) continue;
        count++;
        printProgress(sumcount+count,nbinsPS*NSTATS);
      }
      norm = 1./NBINSCOARSE; // normalisation factor for m0
      norm *= (maxalpha[i][j]-minalpha[i][j]) / (2*pi); // normalisation factor for alpha
      norm *= (maxbeta[i][j]-minbeta[i][j]) / (2*pi); // normalisation factor for beta
      htemp->Scale(1./throws * norm);
      sumthrows += throws;
      sumcount += count;
      sumentries += htemp->GetEntries();
      hpsNH->Add(htemp);
    }
  }
  std::cout<<"Generated " << (float)sumthrows << " throws, of which " << (float)sumentries;
  std::cout<<" (" << 100.*sumentries/sumthrows << "%) were used for histogram." << std::endl;
  */
  
  long throws2=0, count2=0;
  while (count2<NTHROWSLOWMASS) {
    m0 = pow(10.,(duran->Rndm()-3.)); // log prior (25% of phase space)
    alpha = pi*(1.+0.4*(duran->Rndm()-0.5)); // 0.8-1.2 pi (20% of phase space)
    beta = 2*pi*duran->Rndm(); // 0-2 pi (100% of phase space)
    meffNH = meff(m0,alpha,beta,0);
    throws2++;
    if (meffNH>1e-4) continue;
    hangNHlo->Fill(alpha/pi,beta/pi);
    hpsNHlo->Fill(m0, meffNH); 
    count2++;
    printProgress(count2,NTHROWSLOWMASS);
    if (meffNH<meffmin) { alpha0=alpha; beta0=beta; m0min=m0; meffmin=meffNH; }
  }
  hpsNHlo->Scale(1./throws2*0.25*0.2);
  hangNHlo->Scale(1./throws2*0.25*0.2);
  std::cout<<"Generated " << (float)throws2 << " throws, of which " << hpsNHlo->GetEntries();
  std::cout<<" (" << 100.*hpsNHlo->GetEntries()/throws2 << "%) were used for histogram." << std::endl;
  std::cout<<"MINIMUM: meff="<<1.e3*meffmin<<" meV for m0="<<1.e3*m0min<<" meV, alpha="<<alpha0/pi<<" pi, beta="<<beta0/pi<<" pi"<<std::endl;
  
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
  //hpsIH->SetTitle("#bf{0#nu#beta#beta decay - Type-I seesaw}");
  hpsIH->GetXaxis()->SetTitle("m_{0} (eV)");
  hpsIH->GetYaxis()->SetTitle("#LTm_{#beta#beta}#GT (eV)");
  hpsIH->Draw("colz");
  hpsIH->GetZaxis()->SetRangeUser(2e-9,2e-4);
  hpsNH->Draw("colz same");
  hpsNH->GetZaxis()->SetRangeUser(2e-9,2e-4);
  Tl->DrawLatex(1.4e-4,5.6e-1,"#alpha #equiv #alpha_{1}");
  Tl->DrawLatex(1.4e-4,2.8e-1,"#beta #equiv #alpha_{2}#font[122]{-} 2#delta_{CP}");
  Tl->DrawLatex(1.4e-4,1.4e-1,"#alpha, #beta #in [0, 2#pi]");
  Tl->DrawLatex(1.4e-4,2.6e-2,"IH");
  Tl->DrawLatex(1.4e-4,2.0e-3,"NH");
  c1->Print("0vbb_phasespace.png");
  c1->Print("0vbb_phasespace.pdf");
  c1->Close();
  
  TCanvas *c3 = new TCanvas("c3","Phase space plot",1200,900);
  c3->SetLogx(); c3->SetLogy(); c3->SetLogz(); c3->SetGrid();
  //hpsNHlo->SetTitle("#bf{0#nu#beta#beta decay - allowed phase space}");
  hpsNHlo->GetXaxis()->SetTitle("m_{0} (eV)");
  hpsNHlo->GetYaxis()->SetTitle("#LTm_{#beta#beta}#GT (eV)");
  hpsNHlo->Draw("colz");
  hpsNHlo->GetZaxis()->SetRangeUser(3e-10,2e-7);
  Tl->DrawLatex(2e-4,0.4,"#alpha #equiv #alpha_{1} #in [0, 2#pi]");
  Tl->DrawLatex(2e-4,0.2,"#beta #equiv #alpha_{2}#font[122]{-} 2#delta_{CP} #in [0, 2#pi]");
  c3->Print("0vbb_phasespace_lowmass.png");
  c3->Print("0vbb_phasespace_lowmass.pdf");
  c3->Close();
  
  gStyle->SetOptStat(1111); gStyle->SetStatX(0.885); gStyle->SetStatY(0.88);
  TCanvas *c2 = new TCanvas("c2","Majorana angle plot",1200,900);
  c2->SetLogz(); c2->SetGrid();
  //hangNHlo->SetTitle("#bf{Majorana hangNHlo that lead to m_{#beta#beta} < 10^{-4} eV}");
  hangNHlo->GetXaxis()->SetTitle("#alpha #equiv #alpha_{1} [#pi]");
  hangNHlo->GetYaxis()->SetTitle("#beta #equiv #alpha_{2}#font[122]{-} 2#delta_{CP} [#pi]");
  hangNHlo->Draw("colz");
  hangNHlo->GetXaxis()->SetRangeUser(0.8,1.2);
  hangNHlo->GetZaxis()->SetRangeUser(3e-10,3e-8);
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

