#include <TAxis.h>
#include <TCanvas.h>
#include <TComplex.h>
#include <TF1.h>
#include <TH1.h>
#include <TH2.h>
#include <TLatex.h>
#include <TLegend.h>
#include <TList.h>
#include <TMath.h>
#include <TRandom3.h>
#include <TString.h>
#include <TStyle.h>

const Double_t EPSILON = 1e-10;

Double_t * GetAxisArray(TAxis * aa);
void FitSlicesY(const TH2D *hh, TH1D *&hnor, TH1D *&hmpv, TH1D *&hwid, TH1D *&hres, TH1D *&hchi, const TString formula, const Double_t thres, TList *ll=0x0);
void BinLog(TAxis *axis, const Double_t non0start=-999.);

// *****************************************************************************
void BinLog(TAxis *axis, const Double_t non0start)
// *****************************************************************************
{
  
  const Int_t bins = axis->GetNbins();

  const Double_t xmin = axis->GetXmin();
  const Double_t xmax = axis->GetXmax();

  Bool_t k0start = kFALSE;
  if (xmin<EPSILON){
    k0start = kTRUE;
    if(non0start<EPSILON){
      printf("NeutrinoTools::BinLog bad non0start %f\n", non0start); exit(1);
    }
  }
  
  Double_t *new_bins = new Double_t[bins + 1];

  const Double_t factor = k0start? (TMath::Power(xmax/non0start, 1./(bins-1))) : (TMath::Power(xmax/xmin, 1./bins)) ;

  new_bins[0] = xmin;
  new_bins[1] = k0start ? non0start : (new_bins[0]*factor);

  for (int i = 2; i <= bins; i++) {
    new_bins[i] = factor * new_bins[i-1];
  }
  axis->Set(bins, new_bins);
  delete [] new_bins;
}

// *****************************************************************************
void FitSlicesY(const TH2D *hh, TH1D *&hnor, TH1D *&hmpv, TH1D *&hwid, TH1D *&hres, TH1D *&hchi, const TString formula, const Double_t thres, TList *ll)
// *****************************************************************************
{
  const Int_t x0 = hh->GetXaxis()->GetFirst();
  const Int_t x1 = hh->GetXaxis()->GetLast();
  const Int_t y0 = hh->GetYaxis()->GetFirst();
  const Int_t y1 = hh->GetYaxis()->GetLast();

  const Int_t nx = hh->GetNbinsX();
  const Int_t ny = hh->GetNbinsY();
  const Double_t xmin = hh->GetXaxis()->GetXmin();
  const Double_t xmax = hh->GetXaxis()->GetXmax();
  const Double_t ymin = hh->GetYaxis()->GetXmin();
  const Double_t ymax = hh->GetYaxis()->GetXmax();

  hnor = new TH1D(Form("%s_%samp",hh->GetName(), formula.Data()), "", nx, xmin, xmax); if(ll){ll->Add(hnor);}
  hmpv = new TH1D(Form("%s_%smpv",hh->GetName(), formula.Data()), "", nx, xmin, xmax); if(ll){ll->Add(hmpv);}
  hwid = new TH1D(Form("%s_%swid",hh->GetName(), formula.Data()), "", nx, xmin, xmax); if(ll){ll->Add(hwid);}
  hres = new TH1D(Form("%s_%sres",hh->GetName(), formula.Data()), "", nx, xmin, xmax); if(ll){ll->Add(hres);}
  hchi = new TH1D(Form("%s_%schi",hh->GetName(), formula.Data()), "", nx, xmin, xmax); if(ll){ll->Add(hchi);}

  const Double_t *hxbins = GetAxisArray(hh->GetXaxis());
  
  hnor->GetXaxis()->Set(nx, hxbins);
  hmpv->GetXaxis()->Set(nx, hxbins);
  hwid->GetXaxis()->Set(nx, hxbins);
  hres->GetXaxis()->Set(nx, hxbins);
  hchi->GetXaxis()->Set(nx, hxbins);
  delete hxbins;

  for(Int_t ix=x0; ix<=x1; ix++){
    TH1D *htmp = new TH1D(Form("%s_%s%d", hh->GetName(), formula.Data(), ix),"",ny, ymin, ymax);
    //checked, ok
    const Double_t *hhybins = GetAxisArray(hh->GetYaxis());
    htmp->GetXaxis()->Set(ny, hhybins);
    delete hhybins;
    
    Double_t ntot = 0;
    for(Int_t iy=y0; iy<=y1; iy++){
      const Double_t be = hh->GetBinError(ix,iy);
      const Double_t bc = hh->GetBinContent(ix, iy);

      if(be<EPSILON){
        if(bc>EPSILON){
          printf("NeutrinoTools::FitSlicesY error %d %d %e %e\n", ix, iy, be, bc); exit(1);
        }
        continue;
      }

      htmp->SetBinContent(iy, bc);
      htmp->SetBinError(iy, be);

      ntot += (bc/be)*(bc/be);

      //if(be) printf("test %d %d : %f %f %f\n", ix, iy, bc, be, pow(bc/be,2));
    }


    hnor->SetBinContent(ix, ntot);
    hnor->SetBinError(  ix, 0);

    if(ntot<thres || htmp->GetRMS()<EPSILON){
      delete htmp;
      continue;
    }

    //test htmp->Draw();
    Double_t pars[10]={htmp->Integral(0,htmp->GetNbinsX()+1)*htmp->GetBinWidth(1), htmp->GetMean()*0.9, htmp->GetRMS()*0.5};
    Double_t errs[10]={0,0,0,0,0,0,0,0,0,0}, chi[10]={0,0,0,0,0,0,0,0,0,0};

    if(formula=="RMS"){
      pars[1]=htmp->GetMean();
      errs[1]=htmp->GetMeanError();

      //remember that GetRMS is affected by SetRangeUser!!
      pars[2]=htmp->GetRMS();
      errs[2]=htmp->GetRMSError();
    }
    else{
      TF1 * tmpf1 = 0x0;
      if(formula.Contains("Gaus")){
        tmpf1 = new TF1("tmpf1","TMath::Abs([0])*TMath::Gaus(x,[1],TMath::Abs([2]),1)", xmin, xmax);
      }
      else if(formula.Contains("Cauchy")){
        tmpf1 = new TF1("tmpf1","TMath::Abs([0])*TMath::CauchyDist(x,[1],TMath::Abs([2]))", xmin, xmax);
      }
      else{
        printf("NeutrinoTools::FitSlicesY not known formula %s\n", formula.Data());exit(1);
      }

      tmpf1->SetParameters(pars);
      htmp->Fit(tmpf1, Form("%sQ",formula.Contains("LS")?"":"L"));

      tmpf1->GetParameters(pars);
      for(Int_t ipar=0; ipar<tmpf1->GetNpar(); ipar++){
        errs[ipar]=tmpf1->GetParError(ipar);
      }

      chi[0]=tmpf1->GetChisquare();
      chi[1]=tmpf1->GetNDF();

      delete tmpf1;
    }

    pars[2]=TMath::Abs(pars[2]);
    //hnor->SetBinContent(ix, htmp->GetBinContent(htmp->GetMaximumBin()));//htmp->Integral(0,htmp->GetNbinsX()+1));
    hmpv->SetBinContent(ix, pars[1]);
    hmpv->SetBinError(  ix, errs[1]);

    hwid->SetBinContent(ix, pars[2]);
    hwid->SetBinError(  ix, errs[2]);

    hres->SetBinContent(ix, fabs(pars[1])>EPSILON? pars[2]/fabs(pars[1]):0);
    hres->SetBinError(  ix, fabs(pars[1])>EPSILON? errs[2]/fabs(pars[1]):0);

    hchi->SetBinContent(ix, chi[1]>=1 ? chi[0]/chi[1]: 0);
    hchi->SetBinError(ix, 0);

    if(ll){
      ll->Add(htmp);
    }
    else{
      delete htmp;
    }
  }

  TH1 *hhs[]={hnor, hmpv, hwid, hres, hchi};
  const TString yt[]={"N", "MPV", "#sigma", "#sigma/MPV", "#chi^{2}/NDOF"};
  const Int_t nh = sizeof(hhs)/sizeof(TH1*);
  for(Int_t ii=0; ii<nh; ii++){
    hhs[ii]->SetYTitle(Form("%s of %s", yt[ii].Data(), hh->GetYaxis()->GetTitle()));
    hhs[ii]->SetXTitle(hh->GetXaxis()->GetTitle());
    hhs[ii]->GetYaxis()->SetTitleOffset(hh->GetYaxis()->GetTitleOffset());
    hhs[ii]->SetTitle(hh->GetTitle());
  }
}

// *****************************************************************************
Double_t * GetAxisArray(TAxis * aa)
// *****************************************************************************
{
  const Int_t nbin=aa->GetNbins();
  Double_t *bins = new Double_t[nbin+1];
  
  for(Int_t ii=0; ii<=nbin; ii++){
    bins[ii] = aa->GetBinUpEdge(ii);
  }

  return bins;
}

