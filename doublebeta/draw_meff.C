// ---------------------------------------------------
//  Compile:  g++ -g -O2 -o 0vbb_phasespace.exe 0vbb_phasespace.C `root-config --cflags --libs` -I$RATROOT/include/libpq -I$RATROOT/include -L$RATROOT/lib -lRATEvent_Linux -lGeom
//  Run:      ./0vbb_phasespace.exe
//  Info:     Will take 2-3 min to generate distributions. Run a second time to plot.
//  Author:   Martti Nirkko, April 2019 (based on work done at INSS 2015)
// ---------------------------------------------------
// Include helper functions / plotting style
#include "../include/HelperFunc.C"
#include "../include/CommonStyle.H"

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
double meff(double m0, TComplex* mcomp, double alpha, double beta, int inverted=0);

// Random number generator
TRandom3* rng = new TRandom3();

int main() {

  double alpha = 0.42*pi;
  double beta = 0.66*pi;
  double m0 = 5e-3;
  TComplex mcompNH[3], mcompIH[3];
  double meffNH = meff(m0,mcompNH,alpha,beta,0);
  double meffIH = meff(m0,mcompIH,alpha,beta,1);
  
  // Convert to meV
  meffNH *= 1e3;
  meffIH *= 1e3;
  for (int i=0; i<3; i++) {
    mcompNH[i] *= 1e3;
    mcompIH[i] *= 1e3;
  }
  
  // Print results
  cout << "***" << endl << "Normal hierarchy:" << endl;
  cout << "Effective mass = " << meffNH << " meV (NH)" << endl;
  cout << "z = " << mcompNH[0].Re() << " + " << mcompNH[0].Im() << "*i, |z| = " << mcompNH[0].Abs(mcompNH[0]) << endl;
  cout << "z = " << mcompNH[1].Re() << " + " << mcompNH[1].Im() << "*i, |z| = " << mcompNH[1].Abs(mcompNH[1]) << endl;
  cout << "z = " << mcompNH[2].Re() << " + " << mcompNH[2].Im() << "*i, |z| = " << mcompNH[2].Abs(mcompNH[2]) << endl;
  
  cout << "***" << endl << "Inverted hierarchy:" << endl;
  cout << "Effective mass = " << meffIH << " meV (IH)" << endl;
  cout << "z = " << mcompIH[0].Re() << " + " << mcompIH[0].Im() << "*i, |z| = " << mcompIH[0].Abs(mcompIH[0]) << endl;
  cout << "z = " << mcompIH[1].Re() << " + " << mcompIH[1].Im() << "*i, |z| = " << mcompIH[1].Abs(mcompIH[1]) << endl;
  cout << "z = " << mcompIH[2].Re() << " + " << mcompIH[2].Im() << "*i, |z| = " << mcompIH[2].Abs(mcompIH[2]) << endl;
  
  // Add complex numbers (vector addition)
  double x[4][2] = {0};
  double y[4][2] = {0};
  for (int i=0; i<3; i++) {
    for (int j=0; j<=i; j++) {
      x[i+1][0] += mcompNH[j].Re();
      y[i+1][0] += mcompNH[j].Im();
      x[i+1][1] += mcompIH[j].Re();
      y[i+1][1] += mcompIH[j].Im();
    }
    cout << "Point for NH is (" << x[i+1][0] << ", " << y[i+1][0] << ")" << endl;
  }
    
  // Use T2K style
  CommonStyle();
  gStyle->SetPadBottomMargin(0.14);
  gStyle->SetPadLeftMargin(0.14);//to include both large/small font options
  gStyle->SetPadRightMargin(0.03);//to include both large/small font options
  gStyle->SetLabelSize(0.05,"xyz");
  gStyle->SetTitleSize(0.05,"xyz");
  gStyle->SetLabelOffset(0.01,"xyz");
  gStyle->SetTitleOffset(1.4,"x");
  gStyle->SetTitleOffset(1.3,"y");
  gStyle->SetMarkerStyle(7);
  gROOT->ForceStyle();  // force histograms read in below to use this style
  
  // Canvas
  TCanvas c("c","",1300,600);
  c.Divide(2,1);
    
  // Text in figures
  TLatex Tl(0,0,"");
  Tl.SetTextSize(0.045);
  
  // Line in figures
  TLine mcomp(0,0,1,1);
  mcomp.SetLineWidth(2);
  
  // Normal hierarchy
  c.cd(1)->SetGrid();
  c.cd(1)->DrawFrame(0,0,4,4,"NH;Re(m) [meV];Im(m) [meV]");
  mcomp.SetLineColor(2);
  for (int i=0; i<3; i++) {
    mcomp.DrawLine(x[i][0],y[i][0],x[i+1][0],y[i+1][0]);
  }
  mcomp.SetLineColor(4);
  mcomp.DrawLine(x[0][0],y[0][0],x[3][0],y[3][0]);
  Tl.DrawLatex(1.3,0.15,"#color[2]{c_{12}c_{13}m_{1}}");
  Tl.DrawLatex(2.0,0.85,"#color[2]{s_{12}c_{13}m_{2}}");
  Tl.DrawLatex(0.5,0.85,"#color[2]{s_{13}m_{3}}");
  Tl.DrawLatex(0.15,0.2,"#color[4]{m_{#beta#beta}}");
  Tl.DrawLatex(0.25,3.65,"#bf{Normal hierarchy}");
  Tl.DrawLatex(0.25,3.35,"m_{0} = 5 meV");
  Tl.DrawLatex(0.25,3.05,"#varphi_{2} = 0.42 #pi");
  Tl.DrawLatex(0.25,2.75,"#varphi_{3} = 0.66 #pi");
  
  // Inverted hierarchy
  c.cd(2)->SetGrid();
  c.cd(2)->DrawFrame(0,0,40,40,"IH;Re(m) [meV];Im(m) [meV]");
  mcomp.SetLineColor(2);
  for (int i=0; i<3; i++) mcomp.DrawLine(x[i][1],y[i][1],x[i+1][1],y[i+1][1]);
  mcomp.SetLineColor(4);
  mcomp.DrawLine(x[0][1],y[0][1],x[3][1],y[3][1]);
  Tl.DrawLatex(16,1.5,"#color[2]{c_{12}c_{13}m_{1}}");
  Tl.DrawLatex(27,4.5,"#color[2]{s_{12}c_{13}m_{2}}");
  Tl.DrawLatex(18,8.3,"#color[2]{s_{13}m_{3}}");
  Tl.DrawLatex(8.5,5.2,"#color[4]{m_{#beta#beta}}");
  Tl.DrawLatex(2.5,36.5,"#bf{Inverted hierarchy}");
  Tl.DrawLatex(2.5,33.5,"m_{0} = 5 meV");
  Tl.DrawLatex(2.5,30.5,"#varphi_{2} = 0.42 #pi");
  Tl.DrawLatex(2.5,27.5,"#varphi_{3} = 0.66 #pi");
  c.Print("draw_meff.png");
  c.Print("draw_meff.pdf");
  c.Close();
  
  return 0;
}

double meff(double m0, TComplex* mcomp, double alpha, double beta, int inverted) {

  // Neutrino mixing parameters
  double s12, s13, dm21, dm32;
  s12 = S12;
  s13 = S13;
  dm21 = DM21;
  if (inverted) dm32 = DM32IH;
  else          dm32 = DM32NH;
  double c12 = 1.-s12;
  double c13 = 1.-s13;
  
  // Define complex numbers in polar form
  TComplex *fac1 = new TComplex(1.,0.,true);
  TComplex *fac2 = new TComplex(1.,2.*alpha,true);
  TComplex *fac3 = new TComplex(1.,2.*beta,true);
  
  // Multiply by PMNS matrix element and neutrino mass eigenstate
  if (inverted) {
    *fac1 *= c12*c13*sqrt(m0*m0-dm32-dm21);   // dm21>0
    *fac2 *= s12*c13*sqrt(m0*m0-dm32);        // dm32<0
    *fac3 *= s13*m0;                          // m3 = m0
  } else {
    *fac1 *= c12*c13*m0;                      // m1 = m0
    *fac2 *= s12*c13*sqrt(m0*m0+dm21);        // dm21>0
    *fac3 *= s13*sqrt(m0*m0+dm21+dm32);       // dm32>0
  }
  
  mcomp[0] = *fac1;
  mcomp[1] = *fac2;
  mcomp[2] = *fac3;
  
  // Add the three components and return absolute value
  TComplex *res = new TComplex(*fac1 + *fac2 + *fac3);
  return res->Abs(*res);
}
