#include <TROOT.h>
#include <TMath.h>
#include <TVector2.h>
#include <TVector3.h>
#include <TGaxis.h>
#include <TGraph2D.h>
#include <TH1F.h>
#include "TF2.h"
#include "TPad.h"
using namespace std;

// Global constants
const double pi = TMath::Pi();
const TVector3 e1(1,0,0);
const TVector3 e2(0,1,0);
const TVector3 e3(0,0,1);

// Initialise functions
string printVector(const TVector3&);
void GetRotationAngles(const TVector3&, double&, double&);
void GetHotLimit(int*, int&);
void GetMaxColVal(const TVector3&, int*, int, int&, int&, const RAT::DU::PMTInfo&);
void FillHemisphere(const TVector3&, int*, int, TGraph**, TGraph2D*, int, int, const RAT::DU::PMTInfo&);
void DrawCircle(const TVector3&, double, TVector3**, int);
void FitLightSpot(TGraph2D*, double, double, double*);

// Display vector as string
string printVector(const TVector3& v) {
  string out = Form("(%8.2lf | %8.2lf | %8.2lf )", v.X(),  v.Y(), v.Z());
  return out.c_str();
}

// Rotate from given direction (z') to xyz-frame
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

// Get maximum allowed hits/PMT before qualifying as hot PMT
void GetHotLimit(int* pmthitcount, int NPMTS, int &maxnhit) { 
  const int NBINS = 100;
  int MAX_NHIT = *max_element(pmthitcount,pmthitcount+NPMTS);   // hottest PMT
  TH1F *hCount = new TH1F("hCount","",NBINS,0.,log10(MAX_NHIT+1.));
  for(int id=0; id<NPMTS; id++) {
    if (pmthitcount[id]==0) continue; 
    hCount->Fill(log10(pmthitcount[id]));
  }
  int i = hCount->GetMaximumBin(); // start at max. bin
  while (hCount->GetBinContent(i) > 1) i++;
  maxnhit = (int)round(pow(10.,hCount->GetBinLowEdge(i)));
  /*
  TCanvas *c1 = new TCanvas();
  c1->SetGrid();
  c1->SetLogy();
  hCount->Draw();
  c1->Print("hitcount.png");
  c1->Close();
  */
  if (hCount) delete hCount;
}

// Get maximum nhit values within both hemispheres, excluding hot PMTs
void GetMaxColVal(const TVector3& center, int* pmthitcount, int NPMTS, int &nearval, int &farval, const RAT::DU::PMTInfo& pmtinfo) {
  const int NBINS = 1e3;
  //int MAX_NHIT = *max_element(pmthitcount,pmthitcount+NPMTS);   // hottest PMT
  int hotlimit;
  GetHotLimit(pmthitcount, NPMTS, hotlimit);
  TH1F *hNear = new TH1F("hNear","",NBINS,0,hotlimit);
  TH1F *hFar  = new TH1F("hFar","",NBINS,0,hotlimit);
  TVector3 pmtpos, newpos;
  for(int id=0; id<NPMTS; id++) {
    pmtpos = pmtinfo.GetPosition(id);
    if (pmtpos.Mag()==0) continue;              // not a valid PMT position
    if (pmtinfo.GetType(id) != 1) continue;     // not a normal PMT
    if (pmthitcount[id] > hotlimit) continue;	// hot PMT
    if (center.Angle(pmtpos) <= pi/15.)         // narrow cone around central point
      hNear->Fill(pmthitcount[id]);
    else if (center.Angle(pmtpos) >= pi*14/15.) // same around opposite side
      hFar->Fill(pmthitcount[id]);
  }
  int j=1, k=1;
  while (hNear->Integral(1,j)/hNear->Integral(1,NBINS) < 0.99) j++;
  while (hFar->Integral(1,k)/hFar->Integral(1,NBINS) < 0.99) k++;
  //while (hNear->GetBinContent(j) > hNear->GetEntries()/1.e3) j++;
  //while (hFar->GetBinContent(k) > hFar->GetEntries()/1.e3) k++;
  nearval = hNear->GetXaxis()->GetBinCenter(j);
  farval = hFar->GetXaxis()->GetBinCenter(k);
  /*  
  TCanvas *c1 = new TCanvas();
  c1->Divide(2,1);
  c1->cd(1)->SetGrid();
  c1->cd(1)->SetLogy();
  hNear->Draw();
  c1->cd(2)->SetGrid();
  c1->cd(2)->SetLogy();
  hFar->Draw();
  c1->Print(Form("test_z=%d.png",(int)center.Z()));
  c1->Close();
  */
  if(hNear) delete hNear;
  if(hFar) delete hFar;
}

// Create detector view from above central point
void FillHemisphere(const TVector3& center, int* pmthitcount, int NPMTS, TGraph** dots, TGraph2D* graph, int NCOL, int MAXVAL, const RAT::DU::PMTInfo& pmtinfo) {
  
  // Get maximum value for color scale
  int HOTLIMIT;
  GetHotLimit(pmthitcount, NPMTS, HOTLIMIT);
  //GetMaxColVal(center, pmthitcount, NPMTS, MAXVAL, empty, pmtinfo);
   
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

    // Fill graph
    int step = (int)TMath::Ceil(pmthitcount[id]/(1.*MAXVAL/NCOL))+1;
    if (pmthitcount[id] >= MAXVAL) step = NCOL+1;   // cap color range
    if (pmthitcount[id] >= HOTLIMIT) step = 0;      // hot PMT
    //if (newpos.Perp()<1e3) printf("%s has step %d\n",printVector(newpos).c_str(),step);
    dotx[step][ndot[step]]=newpos.X()/1e3;
    doty[step][ndot[step]]=newpos.Y()/1e3;
    ndot[step]++;
    
    // Fill 2D graph (more effective?)
    if (pmtinfo.GetType(id) != 1) continue;     // not a normal PMT (remove OWLEs)
    if (pmthitcount[id]==0 || pmthitcount[id]>HOTLIMIT) continue;  // hot/off PMT
    if (newpos.Z() <= 0) continue;              // not in hemisphere (safety check)
    double xpt = newpos.X()/1e3;//*cos(pi/4)/cos(newpos.Theta());
    double ypt = newpos.Y()/1e3;//*cos(pi/4)/cos(newpos.Theta());
    //if(fabs(xpt)>10 || fabs(ypt)>10) continue;
    graph->SetPoint(counter, xpt, ypt, pmthitcount[id]);
    counter++;
    //for (int n=0; n<pmthitcount[id]; n++) 
    //  hist->Fill(newpos.X()/1e3*cos(pi/4)/cos(newpos.Theta()), newpos.Y()/1e3*cos(pi/4)/cos(newpos.Theta()));
  }

  for (int s=0; s<NCOL+2; s++) {
    int col;
    if      (s==0) col=1;  // hot PMTs
    else if (s==1) col=16; // off PMTs
    else           col=(int)(50.+(s-1)*(50./NCOL));
    if (col>100) printf("*** WARNING: col=%d\n",col);
    if (ndot[s]==0) {dots[s]=NULL; continue;}
    dots[s] = new TGraph(ndot[s],dotx[s],doty[s]);
    dots[s]->SetMarkerStyle(7);
    dots[s]->SetMarkerColor(col);
  }

}

void FitLightSpot(TGraph2D* graph, double radius, double angle, double* params) {
  int npts = graph->GetN();
  double *xpts = graph->GetX();
  double *ypts = graph->GetY();
  double *zpts = graph->GetZ();
  TGraph2D *graf = new TGraph2D();//graph->Clone();
  int pts=0;
  for (int n=0; n<npts; n++) {
    double rpt = sqrt(xpts[n]*xpts[n] + ypts[n]*ypts[n]);
    if (asin(rpt/radius) > angle/180*pi) continue;
    graf->SetPoint(pts,xpts[n],ypts[n],zpts[n]);
  }
  TF2 *fit = new TF2("gaus2d","[0]*TMath::Gaus(x,[1],[2])*TMath::Gaus(y,[3],[4])",-10,10,-10,10);
  fit->SetParameters(*max_element(zpts,zpts+npts),0.,2.,0.,2.);
  graf->Fit("gaus2d");
  for (int par=0; par<5; par++) params[par]=fit->GetParameter(par);
/*
  params[0] = fit->GetParameter(0); // Amplitude
  params[1] = fit->GetParameter(1); // x mean
  params[2] = fit->GetParameter(3); // y mean
  params[3] = fit->GetParameter(2)+fit->GetParameter(4))/2.; // 1 sigma deviation (x-y-averaged)
  params[4] = fit->GetParameter(4); // y deviation
*/
}

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

    vector<TVector3> IcosahedralCentres;
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

    vector<double> distFromCentre;
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
  
};
