// ---------------------------------------------------
//  Compile:  g++ -g -O2 -o isotopes.exe isotopes.C `root-config --cflags --libs` -I$RATROOT/include/libpq -I$RATROOT/include -L$RATROOT/lib -lRATEvent_Linux -lGeom
//  Run:      ./isotopes.exe
//  Info:     R
//  Author:   Martti Nirkko, May 2019
// ---------------------------------------------------
// Include helper functions / plotting style
#include "../include/HelperFunc.C"
#include "../include/CommonStyle.H"

int main() {

  // Number of DBD isotopes
  const int N = 16;
  
  // Read from file
  int a[N];
  string s[N], t[N];
  float x[N], y[N];
  int A, Q, it=0;
  float IA, G, H;
  string S, line;
  ifstream in("isotopes.dat");
  getline(in,line);
  while (in >> A >> S >> Q >> IA >> G >> H) {
    a[it] = A;
    s[it] = S;
    x[it] = IA; // %
    y[it] = Q/1e3; // MeV
    it++;
  }
  
  // Generate graph
  TGraph g(N,x,y);
  g.SetMarkerStyle(8);
  g.SetMarkerColor(2);
  
  // Use T2K style
  CommonStyle();
  gStyle->SetPadBottomMargin(0.14);
  gStyle->SetPadLeftMargin(0.1);//to include both large/small font options
  gStyle->SetPadRightMargin(0.03);//to include both large/small font options
  gStyle->SetLabelSize(0.05,"xyz");
  gStyle->SetTitleSize(0.05,"xyz");
  gStyle->SetLabelOffset(0.01,"xyz");
  gStyle->SetTitleOffset(1.3,"x");
  gStyle->SetTitleOffset(0.9,"y");
  gStyle->SetMarkerStyle(7);
  
  // Plotting section
  TLatex l(0,0,"");
  l.SetTextSize(0.04);
  TCanvas c("c","",800,600);
  c.SetGrid();
  //c.SetLogx();
  TH1F *p = c.DrawFrame(0,0,40,5,";Natural abundance (%); Q-value (MeV)");
  //p->GetXaxis()->SetNdivisions(-1004);
  p->GetYaxis()->SetNdivisions(-505);
  g.Draw("p");
  for (int i=0; i<N; i++) {
    //cout << a[i] << s[i] << endl;
    string tmp = Form("^{%d}%s",a[i],s[i].c_str());
    switch(a[i]) {
      case 48 :
        l.DrawLatex(x[i]+0.5,y[i]-0.1,tmp.c_str());
        break;
      case 76:
      case 96 :
        l.DrawLatex(x[i]-0.2,y[i],tmp.c_str());
        break;
      case 100 :
        l.DrawLatex(x[i]-0.5,y[i],tmp.c_str());
        break;
      case 82 :
        l.DrawLatex(x[i]-2.5,y[i],tmp.c_str());
        break;
      case 116 :
      case 124 :
        l.DrawLatex(x[i]-3.1,y[i]-0.05,tmp.c_str());
        break;
      case 148 :
        l.DrawLatex(x[i]-1.8,y[i]-0.25,tmp.c_str());
        break;
      default :
        l.DrawLatex(x[i]-0.5,y[i],tmp.c_str());
        break;
    }
  }
  c.Print("isotopes.png");
  c.Print("isotopes.pdf");
  c.Close();
  
  return 0;
}

