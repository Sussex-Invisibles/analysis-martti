{
  const int NDEC=3;
  string decay[NDEC]={"PP","PN","NN"};
  float runtime;
  string file, line, tstr;
  TH1F *h[NDEC] = {NULL};
  TCanvas *c = new TCanvas("c","",800,600);
  c->SetGrid();
  c->DrawFrame(0,0,10,60,"Runtimes");
  for (int i=0; i<NDEC; i++) {
    file = Form("%s_runtime.log",decay[i].c_str());
    ifstream in(file.c_str());
    if(!in.good()) continue;
    int count=0;
    tstr = Form("h_%d",i);
    h[i] = new TH1F(tstr.c_str(),decay[i].c_str(),100,0,10);
    while (!in.eof()) {
      in >> runtime;
      h[i]->Fill(runtime/3600.);
    }
    h[i]->SetLineColor(i+2);
    h[i]->SetLineWidth(2);
    h[i]->Draw("same");
  }
  c->BuildLegend();
  c->Print("runtime.png");
  c->Print("runtime.pdf");
  c->Close();
}
