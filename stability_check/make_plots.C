// Includes and helper functions
#include "../HelperFunc.C"

// Main program
int main(int argc, char** argv) {

  string input = "TELLIE_Stability_runs.list";
  ifstream in(input.c_str());
  if (!in) { cerr<<"Failed to open "<<input<<endl; exit(1); }
  
  string line, fibre, out;
  int set, fib;
  int run, subrun, ipw, pin;
  float rms, nhit;
  
  const int NPLOTS = 3;
  const int NSETS = 4;
  const int NFIBRES = 5;
  const string FIBRES[NFIBRES] = {"FT010A","FT033A","FT048A","FT079A","FT090A"};
  
  TGraph *g[NPLOTS][NSETS][NFIBRES] = {{{NULL}}};
  int nsubruns[NSETS][NFIBRES] = {{0}};
  
  while (true) {
  
    in >> run;
    if (!in.good()) break;
    
    string fdata = Form("output/%d.out",run);
    ifstream data(fdata.c_str());
    if (!data) { cerr<<"Failed to open "<<fdata<<endl; break; }
    
    while (true) {
      
      data >> fibre >> run >> subrun >> ipw >> pin >> rms >> nhit;
      if (!data.good()) break;
          
      // Determine dataset index
      set = 0;
      if (run > 111000) set = 1;
      if (run > 111800) set = 2;
      if (run > 115000) set = 3;
        
      // Determine fibre index
      fib = 0;
      for (int ifib=0; ifib<NFIBRES; ifib++) {
        if (fibre != FIBRES[ifib]) continue;
        fib = ifib;
        break;
      }
      
      for (int i=0; i<NPLOTS; i++)
        if(!g[i][set][fib]) g[i][set][fib] = new TGraph();
      
      g[0][set][fib]->SetPoint(nsubruns[set][fib],ipw,pin);
      g[1][set][fib]->SetPoint(nsubruns[set][fib],ipw,nhit);
      g[2][set][fib]->SetPoint(nsubruns[set][fib],pin,nhit);
      nsubruns[set][fib]++;
      
    }
    
  }
  
  // Plot options
  gErrorIgnoreLevel = kWarning;
  gStyle->SetTitleOffset(1.3,"x");
  gStyle->SetTitleOffset(1.3,"y");
  
  // Plot graphs
  TCanvas *c = NULL;
  TVirtualPad *p;
  string title[NPLOTS] = {""};
  string img[NPLOTS] = {"stability_pin-vs-ipw","stability_nhit-vs-ipw","stability_nhit-vs-pin"};
  int col[NSETS] = {2,3,4,92};
  int MINIPW[NFIBRES] = {6950,5900,8950,6050,6250};
  int MINPIN[NFIBRES] = {600,1100,500,800,700};
  
  // Plot legend
  TLegend l(0.1,0.1,0.9,0.9);
  if(g[0][0][0]) l.AddEntry(g[0][0][0],"23 March 2018","p");
  if(g[0][1][0]) l.AddEntry(g[0][1][0],"25 March 2018","p");
  if(g[0][2][0]) l.AddEntry(g[0][2][0],"17 April 2018","p");
  if(g[0][3][0]) l.AddEntry(g[0][3][0],"10 July 2018","p");
  
  // Make plots
  for (int iplt=0; iplt<NPLOTS; iplt++) {
    c = new TCanvas("c","",1500,1000);
    c->Divide(3,2);
    for (int ifib=0; ifib<NFIBRES; ifib++) {
      p = c->cd(ifib+1);
      p->SetGrid();
      title[0] = Form("TELLIE fibre %s;IPW setpoint;PIN reading",FIBRES[ifib].c_str());
      title[1] = Form("TELLIE fibre %s;IPW setpoint;Average NHit",FIBRES[ifib].c_str());
      title[2] = Form("TELLIE fibre %s;PIN reading;Average NHit",FIBRES[ifib].c_str());
      if (iplt==0) p->DrawFrame(MINIPW[ifib],MINPIN[ifib],MINIPW[ifib]+300,MINPIN[ifib]+700)->SetTitle(title[iplt].c_str());
      if (iplt==1) p->DrawFrame(MINIPW[ifib],0,MINIPW[ifib]+300,200)->SetTitle(title[iplt].c_str());
      if (iplt==2) p->DrawFrame(MINPIN[ifib],0,MINPIN[ifib]+700,200)->SetTitle(title[iplt].c_str());
      //p->SetLogx();
      //p->SetLogy();
      //for (int iplt=0; iplt<NPLOTS; iplt++) {
        for (int iset=0; iset<NSETS; iset++) {
          if (!g[iplt][iset][ifib]) continue;
          if (g[iplt][iset][ifib]->GetN()==0) continue;
          g[iplt][iset][ifib]->SetMarkerStyle(7);
          g[iplt][iset][ifib]->SetMarkerColor(col[iset]);
          g[iplt][iset][ifib]->Draw("P same");
        }
      //}
    }
    p = c->cd(6);
    l.Draw();
    c->Print((img[iplt]+".pdf").c_str());
    c->Print((img[iplt]+".png").c_str());
    c->Close();
  }

  // End of program
  //outfile->Close();
  if (c) delete c;
  
}
