// ---------------------------------------------------------
// Content: Helper functions uses for analysing TELLIE data
// Author:  Martti Nirkko, University of Sussex (2017)
// ---------------------------------------------------------

// C++ stuff
#include <fstream>
#include <sstream>
#include <iostream>
#include <stdio.h>
#include <stdlib.h>

// ROOT stuff
#include <TCanvas.h>
#include <TF2.h>
#include <TEllipse.h>
#include <TError.h>
#include <TFile.h>
#include <TGaxis.h>
#include <TGraph2D.h>
#include <TGraph2DErrors.h>
#include <TGraphErrors.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TH3F.h>
#include <TMath.h>
#include <TPad.h>
#include <TROOT.h>
#include <TVector2.h>
#include <TVector3.h>

// RAT stuff
#include "RAT/DU/DSReader.hh"
#include "RAT/DU/PMTInfo.hh"
#include "RAT/DU/Utility.hh"
#include "RAT/DS/Run.hh"
#include "RAT/DS/Entry.hh"
#include "RAT/DS/MC.hh"
#include <RAT/PhysicsUtil.hh>

// Helpful functions
#include "Xianguo.C"

// Include frequently used functions from 'std' (entire namespace is bad practice)
using std::cerr;
using std::cout;
using std::endl;
using std::flush;

// -----------------------------------------------------------------------------
// Run time parameters
const int MORE_OUTPUT = 0;                  // additional plots for testing

// Global constants
const double pi = TMath::Pi();
const TVector3 e1(1,0,0);
const TVector3 e2(0,1,0);
const TVector3 e3(0,0,1);

// Physical constants
const double LAMBDA   = 5.0e-4;             // mm -> TELLIE wavelength
const double N_WATER  = 1.33772;            // at 500 nm -> http://www.philiplaven.com/p20.html
const double C_VACUUM = 299.792458;         // mm/ns (detector units)
const double C_WATER  = C_VACUUM/N_WATER;   // mm/ns
const double ENERGY   = RAT::util::WavelengthToEnergy(LAMBDA); // photon energy [MeV]

// -----------------------------------------------------------------------------
// Initialise functions
string printVector(const TVector3&);
void printProgress(int, int);
void GetRotationAngles(const TVector3&, double&, double&);
void GetLimits(float*, int, float&, float&);
void GetMaxColVal(const TVector3&, float*, int, float&, float&, const RAT::DU::PMTInfo&);
void FillHemisphere(const TVector3&, float*, int, TGraph**, TGraph2D*, float, float, const RAT::DU::PMTInfo&);
void DrawCircle(const TVector3&, double, TVector3**, int);
void FitPromptPeaks(TH2D*, int, float*, float*, TGraph2DErrors*);
void FitLightSpot(TGraph2D*, double, double, double*);
string TriggerToString(int);
double getFWHM(TH1*);

// -----------------------------------------------------------------------------
/// Display vector as a string
string printVector(const TVector3& v) {
  string out;
  if (v.Mag() < 10) out = Form("(%.3f, %.3f, %.3f)", v.X(),  v.Y(), v.Z());
  else              out = Form("(%.1f, %.1f, %.1f)", v.X(),  v.Y(), v.Z());
  return out.c_str();
}

// -----------------------------------------------------------------------------
/// Display progress bar within a loop (can't have any other output in loop!)
void printProgress(int it, int n) {
  int div = (n - (n % 100))/100;              // 1% division
  if (it % div != 0 && it != n-1) return;
  float prog = (float)it/n;
  int barWidth = 70;
  cout << "[";
  int pos = barWidth * prog;
  for (int i=0; i<barWidth; ++i) {
    if (i < pos)                cout << "=";  // processed
    else if (pos+1 == barWidth) cout << "=";  // reached 100%
    else if (i == pos)          cout << ">";  // processing
    else                        cout << " ";  // not yet processed
  }
  int perc = (int)round(100.*prog);
  cout << "] " << perc << "%\r" << flush;
  if (it == n-1) cout << endl;
}

// -----------------------------------------------------------------------------
/// Get rotation angles to rotate from given direction to xyz-frame
void GetRotationAngles(const TVector3& vec, double &rot_Z, double &rot_X) {

  // Get rotation angle from z-axis to central direction
  rot_Z = acos(e3*(vec.Unit()));         // rotate counter-clockwise, towards x-axis by [0,pi]
  
  // Get rotation angle from x-axis to central direction (projected to transverse plane)
  int sign = 0;
  if (vec[1] > 0) sign = +1;
  else            sign = -1;
  TVector3 vect(vec[0],vec[1],0);
  rot_X = sign*acos(e1*(vect.Unit()));   // rotate around z-axis by [-pi,pi]
}

// -----------------------------------------------------------------------------
/// Get intensity limit (hits/PMT) above which PMT is considered a screamer
void GetLimits(float* occupancy, int NPMTS, float &HOTLIMIT, float &COLDLIMIT) {
  const int NBINS = 60;
  float MAX_NHIT = *std::max_element(occupancy,occupancy+NPMTS);   // hottest PMT
  TH1F *hOccup = new TH1F("hOccup","",NBINS,0.,1.); // 0-100% occupancy
  BinLog(hOccup->GetXaxis(),1e-6); // minimum for log scale
  for(int id=0; id<NPMTS; id++) {
    if (occupancy[id]==0) continue;
    hOccup->Fill(occupancy[id]);
  }
  int i = hOccup->GetMaximumBin(); // start at max. bin
  int j = hOccup->GetMaximumBin();
  //cout << "Starting at bin=" << i << " with occupancy " << hOccup->GetBinCenter(i) << " and nevents " << hOccup->GetBinContent(i) << endl;
  while (hOccup->GetBinContent(i) >= 10) i++;
  while (hOccup->GetBinContent(j) >= 10) j--;
  HOTLIMIT = hOccup->GetBinLowEdge(i);
  COLDLIMIT = hOccup->GetBinLowEdge(j+1);
  //cout << "Stopping at bin=" << i << " with occupancy " << hOccup->GetBinCenter(i) << " and nevents " << hOccup->GetBinContent(i) << endl;
  if (hOccup) delete hOccup;
}

// -----------------------------------------------------------------------------
/// Get maximum nhit values within both hemispheres, excluding hot PMTs
void GetMaxColVal(const TVector3& center, float* occupancy, int NPMTS, float &nearval, float &farval, const RAT::DU::PMTInfo& pmtinfo) {
  const int NBINS = 60;
  float HOTLIMIT, COLDLIMIT;
  GetLimits(occupancy, NPMTS, HOTLIMIT, COLDLIMIT);
  TH1F *hNear = new TH1F("hNear","",NBINS,0,1);
  TH1F *hFar  = new TH1F("hFar","",NBINS,0,1);
  BinLog(hNear->GetXaxis(),1e-6);
  BinLog(hFar->GetXaxis(),1e-6);
  TVector3 pmtpos, newpos;
  for(int id=0; id<NPMTS; id++) {
    pmtpos = pmtinfo.GetPosition(id);
    if (pmtpos.Mag()==0) continue;                // not a valid PMT position
    if (pmtinfo.GetType(id) != 1) continue;       // not a normal PMT
    if (occupancy[id] > HOTLIMIT) continue;       // hot PMT
    if (occupancy[id] < COLDLIMIT) continue;      // cold PMT
    if (center.Angle(pmtpos)*180./pi <= 24)       // narrow cone around central point
      hNear->Fill(occupancy[id]);
    else if (center.Angle(pmtpos)*180./pi >= 156) // same around opposite side
      hFar->Fill(occupancy[id]);
  }
  int i=hNear->GetMaximumBin();
  int j=hFar->GetMaximumBin();
  //while (hNear->Integral(1,j)/hNear->Integral(1,NBINS) < 0.99) j++;
  //while (hFar->Integral(1,k)/hFar->Integral(1,NBINS) < 0.99) k++;
  while (hNear->GetBinContent(i) >= 1) i++;
  while (hFar->GetBinContent(j) >= 1) j++;
  nearval = hNear->GetXaxis()->GetBinLowEdge(i);
  farval = hFar->GetXaxis()->GetBinLowEdge(j);
  cout << "Nearval is " << nearval << endl; // TODO - replace this function with hard limits? e.g. nearval=NHit/8, farval=NHit/40
  if(hNear) delete hNear;
  if(hFar) delete hFar;
}

// -----------------------------------------------------------------------------
/// Create detector view from above central point
void FillHemisphere(const TVector3& center, float* occupancy, int NPMTS, TGraph** dots, TGraph2D* graph, int NCOL, float MAXVAL, const RAT::DU::PMTInfo& pmtinfo) {
  
  // Get maximum value for color scale
  float HOTLIMIT, COLDLIMIT;
  GetLimits(occupancy, NPMTS, HOTLIMIT, COLDLIMIT);
   
  int ndot[NCOL+2];
  double dotx[NCOL+2][NPMTS], doty[NCOL+2][NPMTS];
  // Constants only known at runtime
  memset( ndot, 0, (NCOL+2)*sizeof(int) );
  memset( dotx, 0, (NCOL+2)*NPMTS*sizeof(double) );
  memset( doty, 0, (NCOL+2)*NPMTS*sizeof(double) );

  // Find rotation angles for this frame
  double rot_X, rot_Z;
  GetRotationAngles(center,rot_Z,rot_X);
  
  // Loop over all PMTs
  int counter=0;
  TVector3 pmtpos, newpos;
  for(int id=0; id<NPMTS; id++) {
    
    pmtpos = pmtinfo.GetPosition(id);
    if (pmtpos.Mag()==0) continue;                 // not a valid PMT position
    if (center.Angle(pmtpos) > pi/2.) continue;    // not in same hemisphere as central point
    
    // Rotate particles into correct frame
    newpos = pmtpos;
    newpos.RotateZ(-rot_X);
    newpos.RotateY(-rot_Z);

    // Get correct bin for colour scale
    int step;
    if (occupancy[id] == 0) step = 0;                 // off PMT
    else if (occupancy[id] < COLDLIMIT) step = 1;     // cold PMT
    else if (occupancy[id] > HOTLIMIT) step = NCOL+1; // hot PMT
    else if (occupancy[id] > MAXVAL) step = NCOL;     // cap color range
    // linear colour scale
    else step = (int)TMath::Ceil(occupancy[id]/(1.*MAXVAL/(NCOL-1)))+1;
    // logarithmic colour scale
    //else step = (int)TMath::Ceil((log10(occupancy[id])-log10(COLDLIMIT)) / ((log10(MAXVAL)-log10(COLDLIMIT))/(NCOL-1))) + 1;
    
    // Fill array of 1D graphs (bad practice) - TODO: replace this (currently plotted)
    dotx[step][ndot[step]]=newpos.X()/1e3;
    doty[step][ndot[step]]=newpos.Y()/1e3;
    ndot[step]++;
    
    // Fill 2D graph - TODO: more effective? (currently used for fit)
    if (newpos.Z() <= 0) continue;                 // not in hemisphere (safety check)
    if (pmtinfo.GetType(id) != 1) continue;        // not a normal PMT (remove OWLEs)
    else if (occupancy[id] == 0) continue;         // off PMT
    else if (occupancy[id] < COLDLIMIT) continue;  // cold PMT
    else if (occupancy[id] > HOTLIMIT) continue;   // hot PMT
    else if (occupancy[id] > MAXVAL) occupancy[id]=MAXVAL; // cap range
    double xpt = newpos.X()/1e3;//*cos(pi/4)/cos(newpos.Theta());
    double ypt = newpos.Y()/1e3;//*cos(pi/4)/cos(newpos.Theta());
    graph->SetPoint(counter, xpt, ypt, occupancy[id]);
    counter++;
  }

  for (int s=0; s<NCOL+2; s++) {
    int col;
    if      (s==0)  col=16; // off PMTs (grey)
    else if (s==1)  col=51; // cold PMTs (violet)
    else if (s==21) col=1;  // hot PMTs (black)
    else            col=(int)(50.+s*(50./NCOL)); // scale from 55-100 (s=2-20)
    if (ndot[s]==0) {dots[s]=NULL; continue;}
    dots[s] = new TGraph(ndot[s],dotx[s],doty[s]);
    dots[s]->SetMarkerStyle(7);
    dots[s]->SetMarkerColor(col);
  }

}

// -----------------------------------------------------------------------------
/// Fit the prompt peaks (direct light) of the individual PMTs hit time
/// Input:  run number, time-vs-pmtID histogram, number of PMTs, PMT occupancy,
///         PMT angle w.r.t. fibre, array for overall fit parameters
/// Output: 3D graph containing fit results for each PMT
void FitPromptPeaks(TH2D *htime, int NPMTS, float *occupancy, float *pmtangs, TGraph2DErrors *result) {
  
  double x[NPMTS], y[NPMTS], z[NPMTS], ex[NPMTS], ey[NPMTS], ez[NPMTS];
  memset( x,  0, NPMTS*sizeof(double) );
  memset( y,  0, NPMTS*sizeof(double) );
  memset( z,  0, NPMTS*sizeof(double) );
  memset( ex, 0, NPMTS*sizeof(double) );
  memset( ey, 0, NPMTS*sizeof(double) );
  memset( ez, 0, NPMTS*sizeof(double) );

  int saved_pmt[24] = {0};
  string histname[24] = {""};
  string histtitle[24] = {""};
  string fitname[24] = {""};
  TF1 *plotted_fit[24] = {NULL};
  TH1D *plotted_hist[24] = {NULL};
  TCanvas *c1 = NULL;
  
  cout << "Fitting PMT prompt peaks..." << endl;
  for (int iPMT=0; iPMT<NPMTS; iPMT++) {
    if (iPMT % (int)round(NPMTS/100.) == 0) printProgress(iPMT, NPMTS);
    
    // Reject PMTs outside ROI
    if (occupancy[iPMT]<0.01) continue; // only consider PMTs with >=1% occupancy
    if (pmtangs[iPMT]>24.0) continue;   // only consider PMTs within 2x nominal aperture (12 deg)
    TH1D *temp = htime->ProjectionY("temp",iPMT+1,iPMT+1,""); // histogram bins in [1,N]

    // Define prompt peak range (>20% of max. intensity)
    int lobin = temp->GetMaximumBin();
    int hibin = temp->GetMaximumBin();
    while(temp->GetBinContent(lobin) > 0.2*temp->GetMaximum()) lobin--;
    while(temp->GetBinContent(hibin) > 0.2*temp->GetMaximum()) hibin++;
    
    // Fit 1D Gaussian to each slice to get prompt hit time
    c1 = new TCanvas("c1","NULL",800,600);  // suppress default Canvas creation
    TF1 *fitPMT = new TF1("fitPMT", "gaus", temp->GetBinLowEdge(lobin), temp->GetBinLowEdge(hibin+1));
    fitPMT->SetParameters(temp->GetMaximum(), temp->GetBinCenter(htime->GetMaximumBin()), temp->GetRMS());
    temp->Fit("fitPMT","R,q");  // force range, quiet mode
    c1->Close();
    c1 = NULL;
    
    // Fill fitted hit times vs. angle into graph
    x[iPMT]  = fitPMT->GetParameter(0);
    y[iPMT]  = fitPMT->GetParameter(1);
    z[iPMT]  = fitPMT->GetParameter(2);
    ex[iPMT] = fitPMT->GetParError(0);
    ey[iPMT] = fitPMT->GetParError(1);
    ez[iPMT] = fitPMT->GetParError(2);
    
    // Save fit results for individual PMT (for proof of concept)
    int thisbin = (int)floor(pmtangs[iPMT]);
    if (MORE_OUTPUT && saved_pmt[thisbin]==0) {
      // Histogram for this bin
      histname[thisbin] = Form("hangle_%d",thisbin);
      histtitle[thisbin] = Form("PMT #%d (%d deg)",iPMT,thisbin);
      plotted_hist[thisbin] = (TH1D*)temp->Clone();
      plotted_hist[thisbin]->SetDirectory(0); // required to avoid segfault (?)
      plotted_hist[thisbin]->SetName(histname[thisbin].c_str());
      // Fit result for this bin
      fitname[thisbin] = Form("hfit_%d",thisbin);
      plotted_fit[thisbin] = (TF1*)fitPMT->Clone();
      plotted_fit[thisbin]->SetName(fitname[thisbin].c_str());
      saved_pmt[thisbin] = 1;
    }
    
    if (fitPMT) delete fitPMT;
    if (temp) delete temp;
    fitPMT = NULL;
    temp = NULL;
    
  } // PMT loop
  
  if (c1) delete c1;
  
  // Plot fit results for individual PMT (for proof of concept)
  TCanvas *c0 = NULL;
  if (MORE_OUTPUT) {
    TCanvas *c0 = new TCanvas("c0","Single PMT fits",3000,2000);
    c0->Divide(6,4);
    for (int thisbin=0; thisbin<24; thisbin++) {
      if (saved_pmt[thisbin]==0 || plotted_hist[thisbin]==NULL) continue;
      //cout << "PMT at angle " << thisbin << " deg has " << plotted_hist[thisbin]->GetEntries() << " entries around " << plotted_hist[thisbin]->GetMean() << endl;
      c0->cd(thisbin+1)->SetGrid();
      //c0->cd(thisbin+1)->SetLogy();
      plotted_hist[thisbin]->SetStats(1);
      plotted_hist[thisbin]->Draw();
      plotted_hist[thisbin]->Fit(fitname[thisbin].c_str(),"R,q");
      plotted_hist[thisbin]->SetTitle(histtitle[thisbin].c_str());
      plotted_hist[thisbin]->GetXaxis()->SetTitle("Hit time offset [ns]");
      plotted_hist[thisbin]->GetYaxis()->SetTitle("Number of events");
      plotted_hist[thisbin]->GetXaxis()->SetRangeUser(170,220);
      plotted_hist[thisbin]->GetYaxis()->SetRangeUser(0,900); // linear scale
      //plotted_hist[thisbin]->GetYaxis()->SetRangeUser(0.5,900); // log scale
      //plotted_hist[thisbin]->GetXaxis()->SetRangeUser(y[iPMT]-5*z[iPMT],y[iPMT]+5*z[iPMT]);
      //plotted_hist[thisbin]->GetYaxis()->SetRangeUser(1.5e-3*x[iPMT],1.5*x[iPMT]);
      plotted_hist[thisbin]->GetXaxis()->SetTitleOffset(1.3);
      plotted_hist[thisbin]->GetYaxis()->SetTitleOffset(1.6);
    }
    string sglstr = Form("angular_singlePMTs");
    c0->Print((sglstr+".png").c_str());
    c0->Print((sglstr+".pdf").c_str());
    c0->Close();
  }
  if (c0) delete c0;
  
  // Return 2D graph with fit results (by reconstructing object at given address)
  new (result) TGraph2DErrors(NPMTS,x,y,z,ex,ey,ez);
  
}

// -----------------------------------------------------------------------------
/// Fit 2D Gaussian surface over intensities for PMTs in a plane
void FitLightSpot(TGraph2D* graph, double radius, double cone, double* params) {
  // Get input graph entries
  int npts = graph->GetN();
  double *xpts = graph->GetX();   // X projection
  double *ypts = graph->GetY();   // Y projection
  double *zpts = graph->GetZ();   // Intensity
  
  // First loop: Find maximum value in desired range
  double maxval=-1;
  for (int n=0; n<npts; n++) {
    double rpt = sqrt(xpts[n]*xpts[n] + ypts[n]*ypts[n]);
    if (rpt>radius) continue;
    double ang = asin(rpt/radius);
    if (ang/pi*180. > cone) continue;
    if (maxval<zpts[n]) maxval=zpts[n];
  }
  
  // Second loop: Fill values above 20% of maximum values into output graph
  TGraph2D *graf = new TGraph2D();
  int pts=0;
  for (int n=0; n<npts; n++) {
    if(zpts[n] < maxval*0.2) continue;
    double rpt = sqrt(xpts[n]*xpts[n] + ypts[n]*ypts[n]);
    if (rpt>radius) continue;
    double ang = asin(rpt/radius);
    if (ang/pi*180. > cone) continue;
    double scale = 1./cos(ang);
    graf->SetPoint(pts,scale*xpts[n],scale*ypts[n],zpts[n]);
    pts++;
  }
  
  // Fit 2D Gaussian surface to selected points
  TF2 *fit = new TF2("gaus2d","[0]*TMath::Gaus(x,[1],[3])*TMath::Gaus(y,[2],[3])",-10,10,-10,10);
  double aperture = radius*tan(cone/2. / 180.*pi);  // half input angle <=> 12 deg aperture
  fit->SetParameters(0.8*maxval,0.,0.,aperture);    // a priori: 80% max. intensity, centered, nominal aperture
  fit->SetParLimits(0,0.2*maxval,maxval);           // amplitude range [20%, 100%] of max. intensity
  fit->SetParLimits(1,-aperture,aperture);          // x-pos range [-12, +12] degrees from weighted centre
  fit->SetParLimits(2,-aperture,aperture);          // y-pos range [-12, +12] degrees from weighted centre
  fit->SetParLimits(3,0.5*aperture,1.5*aperture);   // sigma range [-50%, +50%] of nominal aperture
  graf->Fit("gaus2d","R,q");                        // fit gaussian (force range, quiet mode)
  
  // Raw fit parameters
  double amp = fit->GetParameter(0);
  double mux = fit->GetParameter(1);
  double muy = fit->GetParameter(2);
  double sig = fit->GetParameter(3);
  
  // Scale mean value back to sphere
  double rad = sqrt(mux*mux+muy*muy);
  double sx = mux*cos(atan(rad/radius));
  double sy = muy*cos(atan(rad/radius));
  
  // Scale (mean +- sigma) values back to sphere
  double xp = (mux+sig);
  double rxp = sqrt(xp*xp+muy*muy);
  double sxp = xp*cos(atan(rxp/radius));
  double xm = (mux-sig);
  double rxm = sqrt(xm*xm+muy*muy);
  double sxm = xm*cos(atan(rxm/radius));
  double yp = (muy+sig);
  double ryp = sqrt(yp*yp+mux*mux);
  double syp = yp*cos(atan(ryp/radius));
  double ym = (muy-sig);
  double rym = sqrt(ym*ym+mux*mux);
  double sym = ym*cos(atan(rym/radius));
  
  // Average over all scaled sigmas
  double sigma = (fabs(sxp-sx)+fabs(sxm-sx)+fabs(syp-sy)+fabs(sym-sy))/4.;
  //printf("Fit result scaled to sphere:\tp1 = %.3lf\tp2 = %.3lf\tp3 = %.3lf\n",sx,sy,sigma);
  
  // Set fit results
  params[0] = amp;    // amplitude
  params[1] = sx;     // mu_x
  params[2] = sy;     // mu_y
  params[3] = sigma;  // sigma
  
  // Plot fit results (for proof of concept)
  if (MORE_OUTPUT) {
    gStyle->SetOptStat(0);
    TCanvas *c = new TCanvas("","",800,800);
    c->SetGrid();
    TH3F *hempty = new TH3F("hempty","",10,-10,10,10,-10,10,10,0,maxval+1e-6);
    hempty->Draw("");             // empty histogram for plot range
    //c->SetTheta(90-0.001);      // view from above
    //c->SetPhi(0+0.001);         // no x-y rotation
    graf->SetMarkerStyle(8);
    graf->Draw("pcol,same");
    fit->Draw("surf,same");
    string name = Form("fit_%d",(int)cone);
    c->Print((name+".png").c_str());
    c->Print((name+".pdf").c_str());
    c->Close();
    if (hempty) delete hempty;
    if (c) delete c;
  }
  
  if (graf) delete graf;
  if (fit) delete fit;
}

// -----------------------------------------------------------------------------
/// Draw a circle around a point in a plane orthogonal to an angle with N dots
void DrawCircle(const TVector3& center, double angle, TVector3** dots, int NDOTS) {
  TVector3 dot(0,0,0);
  for (int i=0; i<NDOTS; i++) {
    dot.SetMagThetaPhi(center.Mag(), pi*angle/180., i*2.*pi/NDOTS);
    dot.RotateUz(center.Unit());   // rotate into frame of input vector
    dots[i]=new TVector3();
    dots[i]->SetXYZ(dot.X(),dot.Y(),dot.Z());
  }
}


// ***********************
// *** HERE BE DRAGONS ***
// ***********************
// This code lazily copied from $RATROOT/src/calib/quality/DataQualityProc.hh
// after failing to access it directly - do not modify!
// -----------------------------------------------------------------------------
string TriggerToString( int trigger ) {
  std::stringstream triggerStream;
  if( trigger == 0x0 )
    triggerStream << "ORPH";
  if( trigger & 0x01 )
    triggerStream << "100L,";
  if( trigger & 0x02 )
    triggerStream << "100M,";
  if( trigger & 0x04 )
    triggerStream << "10HL,";
  if( trigger & 0x08 )
    triggerStream << "20,";
  if( trigger & 0x10 )
    triggerStream << "20L,";
  if( trigger & 0x20 )
    triggerStream << "ESUML,";
  if( trigger & 0x40 )
    triggerStream << "ESUMH,";
  if( trigger & 0x80 )
    triggerStream << "OWL,";
  if( trigger & 0x100 )
    triggerStream << "OWLL,";
  if( trigger & 0x200 )
    triggerStream << "OWLH,";
  if( trigger & 0x400 )
    triggerStream << "PUL,";
  if( trigger & 0x800 )
    triggerStream << "PRE,";
  if( trigger & 0x1000 )
    triggerStream << "PED,";
  if( trigger & 0x2000 )
    triggerStream << "PONG,";
  if( trigger & 0x4000 )
    triggerStream << "SYNC,";
  if( trigger & 0x8000 )
    triggerStream << "EXT,";
  if( trigger & 0x10000 )
    triggerStream << "EXT2,";
  if( trigger & 0x20000 )
    triggerStream << "EXT3,";
  if( trigger & 0x40000 )
    triggerStream << "EXT4,";
  if( trigger & 0x80000 )
    triggerStream << "EXT5,";
  if( trigger & 0x100000 )
    triggerStream << "EXT6,";
  if( trigger & 0x200000 )
    triggerStream << "EXT7,";
  if( trigger & 0x400000 )
    triggerStream << "EXT8,";
  if( trigger & 0x800000 )
    triggerStream << "SPRaw,";
  if( trigger & 0x1000000 )
    triggerStream << "NCD,";
  if( trigger & 0x2000000 )
    triggerStream << "SoftGT,";
  if( trigger & 0x4000000 )
    triggerStream << "MISS";
  return triggerStream.str();
}

// Get FWHM of histogram
double getFWHM(TH1* hist) {
  int bin1 = hist->FindFirstBinAbove(hist->GetMaximum()/2.);
  int bin2 = hist->FindLastBinAbove(hist->GetMaximum()/2.);
  double fwhm = hist->GetXaxis()->GetBinUpEdge(bin2) - hist->GetXaxis()->GetBinLowEdge(bin1);
  return fwhm;
}

// -----------------------------------------------------------------------------
/// Namespace taken from DataQualityProc class to get flatmap detector view
namespace func{

  // Vectors For PSUP Projection
  double fa = 1.0 / 5.5;
  double fb = fa * sqrt( 3.0 ) / 2.0;
  TVector2 *fA12a = new TVector2( fa / 2.0, 0.0 );
  TVector2 *fA12b = new TVector2( 3.0 * fa / 2.0, 0.0 );
  TVector2 *fA12c = new TVector2( 5.0 * fa / 2.0, 0.0 );
  TVector2 *fA12d = new TVector2( 7.0 *fa / 2.0, 0.0 );
  TVector2 *fA12e = new TVector2( 9.0 * fa / 2.0, 0.0 );
  TVector2 *fA2a = new TVector2( 0.0, fb );
  TVector2 *fA2b = new TVector2( 5.0 * fa, fb );
  TVector2 *fA17a = new TVector2( fa / 2.0 , 2.0 * fb );
  TVector2 *fA17b = new TVector2( 11.0 * fa / 2.0 , 2.0 * fb );
  TVector2 *fA51a = new TVector2( fa, 3.0 * fb );
  TVector2 *fA51b = new TVector2( 2.0 * fa, 3.0 * fb );
  TVector2 *fA51c = new TVector2( 3.0 * fa, 3.0 * fb );
  TVector2 *fA51d = new TVector2( 4.0 * fa, 3.0 * fb );
  TVector2 *fA51e = new TVector2( 5.0 * fa, 3.0 * fb );
  TVector2 *fA27 = new TVector2( 4.0 * fa, fb );
  TVector2 *fA46 = new TVector2( 3.0 * fa, fb );
  TVector2 *fA31 = new TVector2( 2.0 * fa, fb );
  TVector2 *fA6 = new TVector2( fa, fb );
  TVector2 *fA37 = new TVector2( 9.0 * fa / 2.0 , 2.0 * fb );
  TVector2 *fA33 = new TVector2( 3.0 * fa / 2.0 , 2.0 * fb );
  TVector2 *fA58 = new TVector2( 5.0 * fa / 2.0 , 2.0 * fb );
  TVector2 *fA54 = new TVector2( 7.0 * fa / 2.0 , 2.0 * fb );

  // ---------------------------------------------------------------------------
  TVector2 TransformCoord( const TVector3& V1, const TVector3& V2, const TVector3& V3, const TVector2& A1, const TVector2& A2, const TVector2& A3,const TVector3& P ) {
    TVector3 xV = V2 - V1;
    TVector3 yV = ( ( V3 - V1 ) + ( V3 - V2 ) ) * 0.5;
    TVector3 zV = xV.Cross( yV ).Unit();

    double planeD = V1.Dot( zV );

    double t = planeD / P.Dot( zV );

    TVector3 localP = t*P - V1;

    TVector2 xA = A2 - A1;
    TVector2 yA = ( ( A3 - A1 ) +( A3 - A2 ) ) * 0.5;

    double convUnits = xA.Mod() / xV.Mag();

    TVector2 result;
    result = localP.Dot( xV.Unit() ) * xA.Unit() * convUnits;
    result += localP.Dot( yV.Unit() ) * yA.Unit() * convUnits + A1;
    return result;
  };

  // ---------------------------------------------------------------------------
  TVector2 IcosProject( TVector3 pmtPos , int &segment ) {
    TVector3 pointOnSphere( pmtPos.X(), pmtPos.Y(), pmtPos.Z() );
    pointOnSphere = pointOnSphere.Unit();
    pointOnSphere.RotateX( -45.0 );
    // From http://www.rwgrayprojects.com/rbfnotes/polyhed/PolyhedraData/Icosahedralsahedron/Icosahedralsahedron.pdf
    const double t = ( 1.0 + sqrt( 5.0 ) ) / 2.0;
    const TVector3 V2 = TVector3( t * t, 0.0, t * t * t ).Unit();
    const TVector3 V6 = TVector3( -t * t, 0.0, t * t * t ).Unit();
    const TVector3 V12 = TVector3( 0.0, t * t * t, t * t ).Unit();
    const TVector3 V17 = TVector3( 0.0, -t * t * t, t * t ).Unit();
    const TVector3 V27 = TVector3( t * t * t, t * t, 0.0 ).Unit();
    const TVector3 V31 = TVector3( -t * t * t, t * t, 0.0 ).Unit();
    const TVector3 V33 = TVector3( -t * t * t, -t * t, 0.0 ).Unit();
    const TVector3 V37 = TVector3( t * t * t, -t * t, 0.0 ).Unit();
    const TVector3 V46 = TVector3( 0.0, t * t * t, -t * t ).Unit();
    const TVector3 V51 = TVector3( 0.0, -t * t * t, -t * t ).Unit();
    const TVector3 V54 = TVector3( t * t, 0.0, -t * t * t ).Unit();
    const TVector3 V58 = TVector3( -t * t, 0.0, -t * t * t ).Unit();

    // Faces {{ 2, 6, 17}, { 2, 12, 6}, { 2, 17, 37}, { 2, 37, 27}, { 2, 27, 12}, {37, 54, 27},
    // {27, 54, 46}, {27, 46, 12}, {12, 46, 31}, {12, 31, 6}, { 6, 31, 33}, { 6, 33, 17},
    // {17, 33, 51}, {17, 51, 37}, {37, 51, 54}, {58, 54, 51}, {58, 46, 54}, {58, 31, 46},
    // {58, 33, 31}, {58, 51, 33}}

    std::vector<TVector3> IcosahedralCentres;
    IcosahedralCentres.push_back( ( V2 + V6 + V17 ) * ( 1.0 / 3.0 ) );
    IcosahedralCentres.push_back( ( V2 + V12 + V6 ) * ( 1.0 / 3.0 ) );
    IcosahedralCentres.push_back( ( V2 + V17 + V37 ) * ( 1.0 / 3.0 ) );
    IcosahedralCentres.push_back( ( V2 + V37 + V27 ) * ( 1.0 / 3.0 ) );
    IcosahedralCentres.push_back( ( V2 + V27 + V12 ) * ( 1.0 / 3.0 ) );
    IcosahedralCentres.push_back( ( V37 + V54 + V27 ) * ( 1.0 / 3.0 ) );

    IcosahedralCentres.push_back( ( V27 + V54 + V46 ) * ( 1.0 / 3.0 ) );
    IcosahedralCentres.push_back( ( V27 + V46 + V12 ) * ( 1.0 / 3.0 ) );
    IcosahedralCentres.push_back( ( V12 + V46 + V31 ) * ( 1.0 / 3.0 ) );
    IcosahedralCentres.push_back( ( V12 + V31 + V6 ) * ( 1.0 / 3.0 ) );
    IcosahedralCentres.push_back( ( V6 + V31 + V33 ) * ( 1.0 / 3.0 ) );
    IcosahedralCentres.push_back( ( V6 + V33 + V17 ) * ( 1.0 / 3.0 ) );

    IcosahedralCentres.push_back( ( V17 + V33 + V51 ) * ( 1.0 / 3.0 ) );
    IcosahedralCentres.push_back( ( V17 + V51 + V37 ) * ( 1.0 / 3.0 ) );
    IcosahedralCentres.push_back( ( V37 + V51 + V54 ) * ( 1.0 / 3.0 ) );
    IcosahedralCentres.push_back( ( V58 + V54 + V51 ) * ( 1.0 / 3.0 ) );
    IcosahedralCentres.push_back( ( V58 + V46 + V54 ) * ( 1.0 / 3.0 ) );
    IcosahedralCentres.push_back( ( V58 + V31 + V46 ) * ( 1.0 / 3.0 ) );
    IcosahedralCentres.push_back( ( V58 + V33 + V31 ) * ( 1.0 / 3.0 ) );
    IcosahedralCentres.push_back( ( V58 + V51 + V33 ) * ( 1.0 / 3.0 ) );

    std::vector<double> distFromCentre;
    unsigned int uLoop;
    for( uLoop = 0; uLoop < IcosahedralCentres.size(); uLoop++ ){
      distFromCentre.push_back( ( IcosahedralCentres[uLoop] - pointOnSphere ).Mag() );
    }
    const int face = min_element( distFromCentre.begin(), distFromCentre.end() ) - distFromCentre.begin() + 1;
    TVector2 resultPosition;
    switch(face){
    case 1://{ 2, 6, 17}
      resultPosition = TransformCoord( V2, V6, V17, *fA2a, *fA6, *fA17a, pointOnSphere );
      break;
    case 2://{ 2, 12, 6}
      resultPosition = TransformCoord( V2, V12, V6, *fA2a, *fA12a, *fA6, pointOnSphere );
      break;
    case 3://{ 2, 17, 37}
      resultPosition = TransformCoord( V2, V17, V37, *fA2b, *fA17b, *fA37, pointOnSphere );
      break;
    case 4://{ 2, 37, 27}
      resultPosition = TransformCoord( V2, V37, V27, *fA2b, *fA37, *fA27, pointOnSphere );
      break;
    case 5://{ 2, 27, 12}
      resultPosition = TransformCoord( V2, V27, V12, *fA2b, *fA27, *fA12e, pointOnSphere );
      break;
    case 6://{37, 54, 27}
      resultPosition = TransformCoord( V37, V54, V27, *fA37, *fA54, *fA27, pointOnSphere );
      break;
    case 7://{27, 54, 46}
      resultPosition = TransformCoord( V27, V54, V46, *fA27, *fA54, *fA46, pointOnSphere );
      break;
    case 8://{27, 46, 12}
      resultPosition = TransformCoord( V27, V46, V12, *fA27, *fA46, *fA12d, pointOnSphere );
      break;
    case 9://{12, 46, 31}
      resultPosition = TransformCoord( V12, V46, V31, *fA12c, *fA46, *fA31, pointOnSphere );
      break;
    case 10://{12, 31, 6}
      resultPosition = TransformCoord( V12, V31, V6, *fA12b, *fA31, *fA6, pointOnSphere );
      break;
    case 11://{ 6, 31, 33}
      resultPosition = TransformCoord( V6, V31, V33, *fA6, *fA31, *fA33, pointOnSphere );
      break;
    case 12://{ 6, 33, 17}
      resultPosition = TransformCoord( V6, V33, V17, *fA6, *fA33, *fA17a, pointOnSphere );
      break;
    case 13://{17, 33, 51}
      resultPosition = TransformCoord( V17, V33, V51, *fA17a, *fA33, *fA51a, pointOnSphere );
      break;
    case 14://{17, 51, 37}
      resultPosition = TransformCoord( V17, V51, V37, *fA17b, *fA51e, *fA37, pointOnSphere );
      break;
    case 15://{37, 51, 54}
      resultPosition = TransformCoord( V37, V51, V54, *fA37, *fA51d, *fA54, pointOnSphere );
      break;
    case 16://{58, 54, 51}
      resultPosition = TransformCoord( V58, V54, V51, *fA58, *fA54, *fA51c, pointOnSphere );
      break;
    case 17://{58, 46, 54}
      resultPosition = TransformCoord( V58, V46, V54, *fA58, *fA46, *fA54, pointOnSphere );
      break;
    case 18://{58, 31, 46}
      resultPosition = TransformCoord( V58, V31, V46, *fA58, *fA31, *fA46, pointOnSphere );
      break;
    case 19://{58, 33, 31}
      resultPosition = TransformCoord( V58, V33, V31, *fA58, *fA33, *fA31, pointOnSphere );
      break;
    case 20://{58, 51, 33}
      resultPosition = TransformCoord( V58, V51, V33, *fA58, *fA51b, *fA33, pointOnSphere );
      break;
    }
    // output face, if needed
    segment = face;
    // 1 - x,y pos to project the same as node map
    return TVector2(1.0 - resultPosition.X(), 1.0 - 2.0 * resultPosition.Y() );
  };
  
  // ---------------------------------------------------------------------------
  void ReverseXAxis (TH1 *h) {
     // Remove the current axis
     h->GetXaxis()->SetLabelOffset(999);
     h->GetXaxis()->SetTickLength(0);
 
     // Redraw the new axis 
     gPad->Update();
     TGaxis *newaxis = new TGaxis(gPad->GetUxmax(), 
                                  gPad->GetUymin(),
                                  gPad->GetUxmin(),
                                  gPad->GetUymin(),
                                  h->GetXaxis()->GetXmin(),
                                  h->GetXaxis()->GetXmax(),
                                  510,"-");
     newaxis->SetLabelOffset(-0.03);
     newaxis->Draw();
  };

  // ---------------------------------------------------------------------------
  void ReverseYAxis (TH1 *h) {
     // Remove the current axis
     h->GetYaxis()->SetLabelOffset(999);
     h->GetYaxis()->SetTickLength(0);

     // Redraw the new axis 
     gPad->Update();
     TGaxis *newaxis = new TGaxis(gPad->GetUxmin(), 
                                  gPad->GetUymax(),
                                  gPad->GetUxmin()-0.001,
                                  gPad->GetUymin(),
                                  h->GetYaxis()->GetXmin(),
                                  h->GetYaxis()->GetXmax(),
                                  510,"+");
     newaxis->SetLabelOffset(-0.03);
     newaxis->Draw();
  };
  
}; // end of namespace func
