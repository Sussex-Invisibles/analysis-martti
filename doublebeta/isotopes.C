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

  // Max. number of isotopes
  const int N = 100;
  
  // Read from file
  int a[N];
  string s[N], t[N];
  float x[N], y[N];
  int A, Q, count=0;
  float IA, G, H;
  string S, line;
  ifstream in("isotopes2.dat");
  getline(in,line);
  while (in >> A >> S >> Q >> IA) {
    if (IA==0) continue; // isotope doesn't exist naturally
    a[count] = A;
    s[count] = S;
    x[count] = IA; // %
    y[count] = Q/1e3; // MeV
    count++;
  }
  cout << "Found " << count << " isotopes." << endl;
  
  // Generate graph (all isotopes)
  TGraph g(count,x,y);
  g.SetMarkerStyle(4);
  g.SetMarkerColor(2);
  
  // Graph with only relevant isotopes
  TGraph h(9);
  int it=0;
  for (int i=0; i<count; i++) {
    if (x[i]==0) continue; // isotope doesn't exist naturally
    if (a[i]!=48 && x[i]<2.0) continue; // low abundance (except 48Ca)
    if (a[i]!=76 && y[i]<2.4) continue; // low Q-value (except 76Ge)
    cout << a[i] << s[i] << endl;
    h.SetPoint(it,x[i],y[i]);
    it++;
  }
  h.SetMarkerStyle(8);
  h.SetMarkerColor(2);
  
  // Use T2K style
  CommonStyle();
  gStyle->SetPadBottomMargin(0.14);
  gStyle->SetPadLeftMargin(0.1);
  gStyle->SetPadRightMargin(0.03);
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
  
  // Plotting range
  int XMIN = 0;
  int XMAX = 50;
  int YMIN = 0;
  int YMAX = 5;
  
  // Draw frame
  TH1F *p = c.DrawFrame(XMIN,YMIN,XMAX,YMAX,";Natural abundance (%); Q-value (MeV)");
  p->GetXaxis()->SetNdivisions(-505);
  p->GetYaxis()->SetNdivisions(-505);
  
  // Draw markers
  g.Draw("p");
  h.Draw("p");
  
  // Draw background lines
  TLine b(0,0,0,0);
  b.SetLineStyle(7);
  b.SetLineColor(4);
  b.DrawLine(XMIN,2.615,XMAX,2.615);
  b.DrawLine(XMIN,3.27,XMAX,3.27);
  
  // Draw background text
  l.DrawLatex(35,2.715,"#color[4]{E_{#gamma}(^{208}Tl) = 2.615 MeV}");
  l.DrawLatex(35,3.37,"#color[4]{E_{#beta}(^{214}Bi) = 3.27 MeV}");
  
  // Draw isotope names near markers
  for (int i=0; i<count; i++) {
    if (x[i]==0) continue; // isotope doesn't exist naturally
    if (a[i]!=48 && x[i]<2.0) continue; // low abundance (except 48Ca)
    if (a[i]!=76 && y[i]<2.4) continue; // low Q-value (except 76Ge)
    //cout << a[i] << s[i] << endl;
    string tmp = Form("^{%d}%s",a[i],s[i].c_str());
    switch(a[i]) {
      case 48 :
        l.DrawLatex(x[i]+0.6,y[i]-0.1,tmp.c_str());
        break;
      case 76:
        l.DrawLatex(x[i]-0.2,y[i],tmp.c_str());
        break;
      case 82 :
        l.DrawLatex(x[i]-3.0,y[i],tmp.c_str());
        break;
      case 96 :
        l.DrawLatex(x[i]-1.0,y[i]+0.1,tmp.c_str());
        break;
      case 100 :
        l.DrawLatex(x[i]-0.5,y[i],tmp.c_str());
        break;
      case 116 :
      case 136 :
        l.DrawLatex(x[i]-3.8,y[i]-0.1,tmp.c_str());
        break;
      case 130 :
        l.DrawLatex(x[i]-2.0,y[i]-0.25,tmp.c_str());
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

