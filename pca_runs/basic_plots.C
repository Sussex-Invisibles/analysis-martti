#include <fstream>
#include <iostream>
#include <stdio.h>
#include <TH1F.h>
#include <TCanvas.h>

void basic_plots();

using namespace std;
int main() {
  basic_plots();
  return 0;
}

void basic_plots() {
  ifstream in("spreadsheet.txt");
  string fibre;
  int channel, run, ipw, photons;
  float nhit;
  //TFile *f = new TFile("pca.root","RECREATE");
  TH1F *h[4];
  h[0] = new TH1F("h1","Run number (10xxxx)",100,0,3000);
  h[1] = new TH1F("h2","Pulse width (IPW setpoint)",30,0,15000);
  h[2] = new TH1F("h3","Photons/pulse (from calibDB)",40,0,2e5);
  h[3] = new TH1F("h4","NHit/pulse (observed)",30,0,60);
  //TNtuple *ntuple = new TNtuple("ntuple","data from ascii file","channel:run:ipw:photons:nhit");
  int nlines=0;
  while (1) {
     in >> fibre >> channel >> run >> ipw >> photons >> nhit;
     if (!in.good()) break;
     if(run>0)     h[0]->Fill(run-1e5);
     if(ipw>0)     h[1]->Fill(ipw);
     if(photons>0) h[2]->Fill(photons);
     if(nhit>0)    h[3]->Fill(nhit);
     //ntuple->Fill(channel,run,ipw,photons,nhit);
     nlines++;
  }
  in.close();
  printf(" found %d points\n",nlines);

  TCanvas *c = new TCanvas("c","",1200,400); 
  c->Divide(3,1);
  for (int i=1; i<4; i++) {
    c->cd(i)->SetGrid();
    h[i]->SetLineWidth(2);
    h[i]->SetLineColor(4);
    h[i]->Draw();
  }
  c->Print("basic_plots.png");
  c->Close();
  
  //f->Write();
  //f->Close();
  if(h) { for (int j=0; j<4; j++) delete h[j]; }
  if(c) delete c;
}
