#include "SetT2KStyle.H"

int CommonStyle() {

  // T2K style
  TString localStyleName = "T2K";
  // -- WhichStyle --
  // 1 = presentation large fonts
  // 2 = presentation small fonts
  // 3 = publication/paper
  Int_t localWhichStyle = 3;

  TStyle* t2kstyle = SetT2KStyle(localWhichStyle, localStyleName);
  gROOT->SetStyle(t2kstyle->GetName());

  // -- margin --
  //  gStyle->SetPadTopMargin(0.05);
  //  gStyle->SetPadBottomMargin(0.16);
  //  gStyle->SetPadRightMargin(0.075);//more space if no colz option 
  //  gStyle->SetPadLeftMargin(0.165); 

  // -- title/lable offset --
  //   gStyle->SetTitleOffset(1.5, "x");
  //   gStyle->SetTitleOffset(1.7, "y");
  //   gStyle->SetLabelOffset(0.02, "x");
  //   gStyle->SetLabelOffset(0.02, "y");

  // -- title/label size --
  //   gStyle->SetTitleSize(0.04, "x"); 
  //   gStyle->SetTitleSize(0.04, "y"); 
  //   gStyle->SetTitleSize(0.04, "z"); 
  //   gStyle->SetLabelSize(0.034,"x"); 
  //   gStyle->SetLabelSize(0.034,"y"); 
  //   gStyle->SetLabelSize(0.034,"z"); 

  // -- statistic and title info --
  //  gStyle->SetOptTitle(1);
  //  gStyle->SetOptStat(1111);

  // -- lines --
  //  gStyle->SetLineWidth(4);

  // -- fills/color --
  //  gStyle->SetFrameFillColor(0); // white color for backgroud
  //  gStyle->SetFillColor(1);

  // -- color scheme --
  // - "warm"/red-ish -
  //  const Int_t NRGBs = 3;
  //  const Int_t NCont = 500;
  //  Double_t red[]   = {1.00, 1.00, 0.25 };
  //  Double_t green[] = {1.00, 0.00, 0.00 };
  //  Double_t blue[]  = {0.00, 0.00, 0.00 };
  //  Double_t stops[] = { 0.25, 0.75, 1.00 };
  //  TColor::CreateGradientColorTable(NRGBs, stops, red, green, blue, NCont);
  //  gStyle->SetNumberContours(NCont);   
  //  - rainbow -
  //  gStyle->SetPalette(1);  // use the rainbow color set

  // -- horizontal error bars back --
  //  gStyle->SetErrorX(0.5);

  // -- transparent stuff --
  //  gStyle->SetFillStyle(0);

  // -- axis --
  //  t2kStyle->SetStripDecimals(kFALSE); // don't do 1.0 -> 1
  //  TGaxis::SetMaxDigits(3); // doesn't seem to work :<
  // supressed zeroes!
  //  gStyle->SetHistMinimumZero(kFALSE);


  return(0);
}
