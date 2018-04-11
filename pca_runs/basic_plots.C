#include <fstream>
#include <iostream>
#include <stdio.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TCanvas.h>

void basic_plots();

using namespace std;
int main() {
  basic_plots();
  return 0;
}

void basic_plots() {
  ifstream in("TELLIE_PCA_Slave.txt");
  string fibre;
  int node, channel, run, ipw, photons, pin, rms;
  float nhit;
  //TFile *f = new TFile("pca.root","RECREATE");
  TH1F *h[6];
  h[0] = new TH1F("h0","Run number (11xxxx)",30,0,1500);
  h[1] = new TH1F("h1","Pulse width (IPW setpoint)",30,0,15000);
  h[2] = new TH1F("h2","Photons/pulse (from calibDB)",30,0,3e4);
  h[3] = new TH1F("h3","PIN",30,0,3000);
  h[4] = new TH1F("h4","RMS",50,0,250);
  h[5] = new TH1F("h5","EXTA-NHit (mean)",30,30,60);
  //TH2F *h2d = new TH2F("h2d","NHit vs PIN",100,0,3000,100,0,60);
  //TNtuple *ntuple = new TNtuple("ntuple","data from ascii file","channel:run:ipw:photons:nhit");
  string line;
  for (int l=0; l<2; l++) getline(in,line); // ignore header
  int nlines=0;
  while (1) {
     in >> node >> fibre >> channel >> run >> ipw >> photons >> pin >> rms >> nhit;
     if (!in.good()) break;
     if(run>0)     h[0]->Fill(run-110000);
     if(ipw>0)     h[1]->Fill(ipw);
     if(photons>0) h[2]->Fill(photons);
     if(pin>0)     h[3]->Fill(pin);
     if(rms>0)     h[4]->Fill(rms);
     if(nhit>0)    h[5]->Fill(nhit);
     //h2d->Fill(pin,nhit);
     //ntuple->Fill(channel,run,ipw,photons,nhit);
     nlines++;
  }
  in.close();
  printf(" found %d points\n",nlines);

  TCanvas *c = new TCanvas("c","",900,600); 
  c->Divide(3,2);
  for (int i=0; i<6; i++) {
    c->cd(i+1)->SetGrid();
    h[i]->SetLineWidth(2);
    h[i]->SetLineColor(4);
    h[i]->Draw();
  }
  c->Print("basic_plots.png");
  c->Print("basic_plots.pdf");
  c->Close();
  
  //TCanvas *d = new TCanvas("d","",600,600); 
  //h2d->SetMarkerStyle(7);
  //h2d->Draw("scat");
  //d->Print("plot2d.pdf");
  //d->Close();

  //f->Write();
  //f->Close();
  if(h) { for (int j=0; j<4; j++) delete h[j]; }
  if(c) delete c;
  //if(d) delete d;
}
