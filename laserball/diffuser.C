// ---------------------------------------------------------
// Goal:          Simulate photon tracking through laserball diffuser flask
// Author:        Martti Nirkko, 24/05/2018
// Compile & run: clear && g++ -g -o diffuser.exe diffuser.C `root-config --cflags --libs` -I$RATROOT/include/libpq -I$RATROOT/include -L$RATROOT/lib -lRATEvent_Linux && ./diffuser.exe
// ---------------------------------------------------------

// Helper functions (includes everything else)
#include "../include/HelperFunc.C"

// Global constants
const int VERBOSE = 1; // verbosity (0-3)
const int NEVENTS = 1e6; // number of photons
const int NBINS = 100; // number of bins

// Laser
const double I0 = 1.; // laser intensity, normalised for now
const double NA = 0.2; // numerical aperture for fibre (from R. Ford's thesis)
const double lambda = 4.05e-5; // wavelength [cm]

// Diffuser ball - TODO REVIEW ALL PARAMETERS
const double flaskdiam = 106.; // diameter of the diffuser flask [mm]
const double R = flaskdiam/2.; // radius
const double na = 1.070; // refractive index (air) - TODO is actually 2:1 mixture of SO2:O2 at 1/3 atm!
const double ns = 1.404; // refractive index (SilGel)
const double ng = 1.530; // refractive index (glass)
// SO2 has n = 1.3047 (20 C) from http://www.microkat.gr/msds/Sulfur%20dioxide.html
// SO2 has n = 1.3396 (25 C) from https://pubchem.ncbi.nlm.nih.gov/compound/sulfur_dioxide
// O2 has n = 1.0003 (0 C) from http://www.kayelaby.npl.co.uk/general_physics/2_5/2_5_7.html

// Glass bubbles
const double r_bub = 2.e-3; // median radius of hollow glass microsphere [cm] - 50% percentile
const double d_bub = 8.8e-5; // approx. wall thickness of hollow sphere [cm] - specs from 3M
const double rho_bub = 2.25; // density of borosilicate glass [g/cm³] - specs from 3M
const double con_bub = 0.51e-3; // density of bubbles in silgel [g/cm³] - TODO use optimum value ~1.57 mg/mL

// Quartz rod
const double rodlen = 200.; // length [mm]
const double roddiam = 1.; // diameter [mm]
const double OFFSETS[] = {0.0, 2.5, 5.0, 7.5, 10.}; // injection point offsets
double OFFSET;

// Other
const double c0 = 299.792458; // speed of light in vacuum [mm/ns]
const double cs = c0/ns; // speed of light in SilGel [mm/ns]

// *****************************************************************************
// Function declarations
TVector3 diffuser(double& tracklen, double scatlen);
double distance_to_wall(TVector3&, TVector3&);
double reflection_prob(const double& n1, const double& n2, const double &impact);
TVector3 propagate(TVector3&, TVector3&, double& distance);
TVector3 scatter(TVector3&, TVector3&, double& impact);
TVector3 reflect_or_refract(TVector3&, TVector3&, double& impact, bool& exitflask);
double GetScatteringLength(double density);

// Other global objects (TODO - BAD PRACTICE)
TRandom3* gen = new TRandom3();
TF1* aperture = new TF1("aperture","cos(pi/2*x/[0])",0,NA);
TVector3 pos, dir, newpos, newdir, endpos, enddir;
TGraph trackS, trackT, trackU; // side view, top view
TH1D htrk("htrk","Single photon tracking (statistics);Distance between scatters [mm];Events",NBINS/2,0,5);
TFile outfile("diffuser.root","RECREATE");
int event;

// *****************************************************************************
// Main program
int main(int argc, char** argv) {

  // Set output verbosity (kPrint, kInfo, kWarning, kError, kBreak, kSysError, kFatal)
  gErrorIgnoreLevel = kWarning;

  // Scattering length of photons in diffuser medium
  double scatlen = 10*GetScatteringLength(con_bub); // mean free path [mm]
  
  // Initialisation
  gen->SetSeed(0);
  aperture->SetParameter(0,NA);
  TVector3 outdir;
  double length;
  gStyle->SetTitleOffset(1.1,"x");
  gStyle->SetTitleOffset(1.4,"y");

  // Loop over all offsets
  for (int it=0; it<sizeof(OFFSETS)/sizeof(double); it++) {
    OFFSET = OFFSETS[it];
    cout << endl << "SIMULATING DIFFUSER WITH INJECTION POINT OFFSET Z = " << OFFSET << " mm" << endl;

    string tlen = Form("hesc_z0=%.1f",OFFSET);
    string tphi = Form("hphi_z0=%.1f",OFFSET);
    string tcth = Form("hcth_z0=%.1f",OFFSET);
    string tang = Form("hang_z0=%.1f",OFFSET);
    TH1D hesc(tlen.c_str(),"Escape time from diffuser;t [ns];Events [#times 10^{3}]",10*NBINS,0,10.);
    TH1D hphi(tphi.c_str(),Form("Azimuthal distribution (z_{0} = %.1f mm);#phi [#pi];Events [#times 10^{3}]",OFFSET),NBINS,-1,1);
    TH1D hcth(tcth.c_str(),Form("Polar distribution (z_{0} = %.1f mm);cos(#theta) [ ];Events [#times 10^{3}]",OFFSET),NBINS,-1,1);
    TH2D hang(tang.c_str(),"Angular distribution of diffuser;#phi [#pi];cos(#theta) [ ]",NBINS/2,-1,1,NBINS/2,-1,1);
    
    // Generate events
    for (event=0; event<NEVENTS; event++) {
      printProgress(event,NEVENTS);
      outdir = diffuser(length, scatlen);
      hesc.Fill(length/cs); // assume constant speed in silicone, neglect time spend in bubbles
      hphi.Fill(outdir.Phi()/pi);
      hcth.Fill(cos(outdir.Theta()));
      hang.Fill(outdir.Phi()/pi,cos(outdir.Theta()));
    }
  
    // Plot distributions
    double scale = 1e-3;
    TCanvas c("c","",1200,1200);
    c.Divide(2,2);
    c.cd(1)->SetGrid();
    hphi.Scale(scale);
    hphi.SetMinimum(0);
    //hphi.SetAxisRange(0,1.25*NEVENTS/NBINS*scale,"Y");
    //hphi.GetYaxis()->SetTitleOffset(1.2);
    hphi.Draw();
    c.cd(2)->SetGrid();
    hcth.Scale(scale);
    hcth.SetMinimum(0);
    //hcth.SetAxisRange(0,1.25*NEVENTS/NBINS*scale,"Y");
    //hcth.GetYaxis()->SetTitleOffset(1.2);
    hcth.Draw();
    c.cd(3)->SetGrid();
    hang.Draw("colz");
    hang.SetMinimum(0);
    //hang.SetAxisRange(0,5*NEVENTS/NBINS/NBINS,"Z");
    //hang.GetYaxis()->SetTitleOffset(1.2);
    hang.Draw("colz same");
    c.cd(4)->SetGrid();
    hesc.Scale(scale);
    hesc.SetMinimum(0);
    //hesc.SetAxisRange(0,5*NEVENTS/NBINS*scale,"Y");
    hesc.Draw();
    c.Print(Form("diffuser_z0=%.1f.png",OFFSET));
    c.Print(Form("diffuser_z0=%.1f.pdf",OFFSET));
    c.Close();
 
    cout << "Time spread is " << getFWHM(&hesc) << " ns (FWHM)." << endl;

    outfile.Write();
  }
  outfile.Close();

  return 0;
}

// *****************************************************************************
// Track a single photon through the diffuser flask
// Inputs: none
// Output: escape direction of photon, seen from 6m sphere
TVector3 diffuser(double& tracklen, double scatlen) {
  
  // Photon position at injection point
  pos = e1;
  double A = pi*pow(roddiam/2.,2)*gen->Rndm(); // random within area
  pos.SetMag(sqrt(A/pi)); // radius [mm]
  pos.SetPhi(2*pi*gen->Rndm()); // angle [rad]
  pos += OFFSET*e3; // offset in z-direction

  // Photon direction at injection point
  dir = e3;
  dir.SetPhi(2*pi*gen->Rndm()); // azimuthal direction [0,2pi)
  double theta = aperture->GetRandom();
  dir.SetTheta(pi-theta); // downwards zenith direction [0,NA)
  
  // Intensity dropoff
  //TF1* inten = new TF1("intensity","[0]*exp(-x/[1])",0,R/10);
  //inten->SetParameter(0,I0);
  //inten->SetParameter(1,scatlen);
  
  // Photon tracking
  int step = 0;
  tracklen = 0;
  trackS.Set(0);
  trackT.Set(0);
  trackU.Set(0);
  // Injection point
  if (!event) {
    trackS.SetPoint(step,pos.X(),pos.Z());
    trackT.SetPoint(step,pos.X(),pos.Y());
    trackU.SetPoint(step,pos.Y(),pos.Z());
  }
  if (VERBOSE > 1) printf("ENTERING diffuser at (%8.3f %8.3f %8.3f), R=%6.3f\n",pos.X(),pos.Y(),pos.Z(),pos.Mag());
  bool exitflask = false;
  while (!exitflask) {
    if (VERBOSE > 2) printf("position (%8.3f %8.3f %8.3f), direction = (%8.3f %8.3f %8.3f)\n",pos.X(),pos.Y(),pos.Z(),dir.X(),dir.Y(),dir.Z());
    double y = gen->Rndm();
    double dsel = -scatlen*log(y);  // randomly selected distance in diffuser [mm]
    double dwall = distance_to_wall(pos,dir); // distance to flask wall [mm]
    double impact = sqrt(y);
    if (dsel < dwall) {
      // propagate selected distance
      newpos = propagate(pos,dir,dsel);
      tracklen += dsel;
      if (!event) htrk.Fill(dsel);
      // scatter in diffuser
      newdir = scatter(newpos,dir,impact);
      step++;
    } else {
      // propagate up to wall
      newpos = propagate(pos,dir,dwall);
      tracklen += dwall;
      if (!event) htrk.Fill(dwall);
      // exit or reflect internally
      newdir = reflect_or_refract(newpos,dir,impact,exitflask);
    }
    pos = newpos;
    dir = newdir;
    if(pos.Mag()>R+1e-6) printf("*** WARNING *** Position is (%8.3f %8.3f %8.3f), R=%6.3f\n",pos.X(),pos.Y(),pos.Z(),pos.Mag());
    if(fabs(dir.Mag()-1.)>1e-6) printf("*** WARNING *** Direction is (%8.3f %8.3f %8.3f), R=%6.3f\n",dir.X(),dir.Y(),dir.Z(),dir.Mag());
    if (!event) {
      trackS.SetPoint(step,pos.X(),pos.Z());
      trackT.SetPoint(step,pos.X(),pos.Y());
      trackU.SetPoint(step,pos.Y(),pos.Z());
    }
  }
  
  if (VERBOSE > 1) printf("EXITING diffuser after %d scatters at (%8.3f %8.3f %8.3f), R=%6.3f\n",step,pos.X(),pos.Y(),pos.Z(),pos.Mag());
  double RMAX = 6e3; // AV radius [mm]
  endpos = propagate(pos,dir,RMAX);
  if (VERBOSE > 1) printf("End point (%8.3f %8.3f %8.3f), R=%6.3f\n",endpos.X(),endpos.Y(),endpos.Z(),endpos.Mag());
  if (!event) {
    trackS.SetPoint(step+1,endpos.X(),endpos.Z());
    trackT.SetPoint(step+1,endpos.X(),endpos.Y());
    trackU.SetPoint(step+1,endpos.Y(),endpos.Z());
    
    double LIM=75;
    gStyle->SetOptStat(0);
    //gStyle->SetAxisLabelOffset?
    TCanvas c("c","",1200,1200);
    c.Divide(2,2);
    c.cd(1)->DrawFrame(-LIM,-LIM,LIM,LIM,"Single photon tracking (front view);X [mm];Z [mm]");
    // Flask
    TEllipse ball(0,0,R);
    ball.Draw("L same");
    // Rod
    TBox rod(-roddiam/2,OFFSET,roddiam/2,LIM); // must be within limits of frame
    rod.SetFillColor(0);
    rod.SetLineColor(4);
    rod.Draw("L same");
    // Photon track
    trackS.SetLineColor(2);
    trackS.Draw("L same");
    c.cd(3)->DrawFrame(-LIM,-LIM,LIM,LIM,"Single photon tracking (top view);X [mm];Y [mm]");
    // Flask
    ball.Draw("L same");
    // Rod
    TEllipse rodc(0,0,roddiam/2);
    rodc.SetFillColor(0);
    rodc.SetLineColor(4);
    rodc.Draw("L same");
    // Photon track
    trackT.SetLineColor(2);
    trackT.Draw("L same");
    c.cd(2)->DrawFrame(-LIM,-LIM,LIM,LIM,"Single photon tracking (side view);Y [mm];Z [mm]");
    // Flask
    ball.Draw("L same");
    // Rod
    rod.Draw("L same");
    // Photon track
    trackU.SetLineColor(2);
    trackU.Draw("L same");
    c.cd(4);
    htrk.GetXaxis()->SetTitleOffset(1.1);
    htrk.GetYaxis()->SetTitleOffset(1.4);
    htrk.Draw();
    c.Print(Form("photon_track_z0=%.1f.png",OFFSET));
    c.Print(Form("photon_track_z0=%.1f.pdf",OFFSET));
    c.Close();
  }

  return endpos.Unit();
}

// *****************************************************************************
// Calculate distance to inner flask boundary
double distance_to_wall(TVector3& pos, TVector3& dir) {
  double x = pos.X();
  double y = pos.Y();
  double z = pos.Z();
  double alpha = dir.X();
  double beta  = dir.Y();
  double gamma = dir.Z();
  double dist = -(x*alpha+y*beta+z*gamma)+sqrt(pow(x*alpha+y*beta+z*gamma,2)-(x*x+y*y+z*z-R*R));
  return dist;
}

// *****************************************************************************
// Reflection probability
double reflection_prob(const double& n1, const double& n2, const double& b) {
    // Consider transition between materials (n2 > n1)
    double ni = n1/n2;
    // Calculate Fresnel coefficients
    double rs = (sqrt(1.-b*b)-sqrt(ni*ni-b*b))/(sqrt(1.-b*b)+sqrt(ni*ni-b*b));
    double rp = (sqrt(ni*ni-b*b)-ni*ni*sqrt(1.-b*b))/(sqrt(ni*ni-b*b)+ni*ni*sqrt(1.-b*b));
    // Average over Fresnel coefficients:
    return 0.5*(rs*rs+rp*rp);
}

// *****************************************************************************
// Propagate through silicone
TVector3 propagate(TVector3& pos, TVector3& dir, double& dist) {
  double x = pos.X()+dir.X()*dist;
  double y = pos.Y()+dir.Y()*dist;
  double z = pos.Z()+dir.Z()*dist;
  return TVector3(x,y,z);
}

// *****************************************************************************
// Scatter off a glass bubble
TVector3 scatter(TVector3& pos, TVector3& dir, double& impact) {
  
  // Directional cosines
  double alpha = dir.X();
  double beta  = dir.Y();
  double gamma = dir.Z();
  
  // Impact parameter
  double b = impact;
  
  // Forward scattering angle depends on scattering modes A, B, C
  double costh;
  if (ng*b > 1) {
    // (A) forward scattering on glass bubble
    costh = 2.*b*b-1.;
  } else {
    // (B) scattering through glass bubble cavity
    costh = 2.*pow(ns*b*b+sqrt((1.-b*b)*(1.-ns*ns*b*b)),2)-1.;
    // (C) probability of reflecting despite nb <= 1
    double prob = reflection_prob(na,ng,b);
    double roll = gen->Rndm();
    if (roll < prob) costh = 2.*b*b-1.;
  }
  
  // Azimuthal angle chosen randomly [0,2*pi)
  double phi = 2.*pi*gen->Rndm();
  
  // New direction transformed back into global frame (R. Ford's thesis)
  double param = sqrt((1.-costh*costh)/(1.-gamma*gamma));
  double newalpha = alpha*costh + param*(alpha*gamma*cos(phi)-beta*sin(phi));
  double newbeta  = beta*costh + param*(beta*gamma*cos(phi)+alpha*sin(phi));
  double newgamma = gamma*costh - param*(1.-gamma*gamma)*cos(phi);  
  TVector3 scatdir(newalpha,newbeta,newgamma);
  return scatdir.Unit();
}

// *****************************************************************************
// Reflect at flask wall, or refract out of diffuser
TVector3 reflect_or_refract(TVector3& pos, TVector3& dir, double& impact, bool& exitflask) {

  // Photon position and directional cosines
  double x = pos.X();
  double y = pos.Y();
  double z = pos.Z();
  double alpha = dir.X();
  double beta  = dir.Y();
  double gamma = dir.Z();
  
  // Forward scattering angle depends on scattering modes A, B, C
  double prod = dir.Unit().Dot(pos.Unit());
  double sinOmega = ng*ng*(1.-fabs(prod));
  TVector3 scatdir;
  if (sinOmega > 1) {
    // (A) internal reflection
    scatdir = dir-2.*prod*pos;
    exitflask = false;
  } else {
    // (B) probability of reflecting despite sinOmega <= 1
    double prob = reflection_prob(ns,ng,impact);
    double roll = gen->Rndm();
    if (roll < prob) {
      scatdir = dir-2.*prod*pos;
      exitflask = false;
      if (VERBOSE > 1) printf("REFLECTING on diffuser wall at (%8.3f %8.3f %8.3f), R=%6.3f\n",x,y,z,pos.Mag());
    }
    // (C) refract out of the diffuser
    else {
      /*
      // Ford's thesis has an undefined 'm' in line (16), damn...
      double d = sqrt(1.-ng*ng*(1.-pow(fabs(prod),2)));
      double a = gamma*y - beta*z;
      double b = alpha*z - gamma*x;
      double f = (a*z*R*d)/(x*z*b-a*y*y);
      double g = (a*x*x+b*x*y*a*z*z)/(a*z*b-a*y*z);
      double c = (R*d+f*y)/x;
      double e = (z+g*y)/x;
      double gam1 = ((e*c+g*f) + sqrt(pow(e*c+g*f,2)-(e*e+g*g+1)*(c*c+f*f-1)))/(e*e+g*g+1);
      double gam2 = ((e*c+g*f) - sqrt(pow(e*c+g*f,2)-(e*e+g*g+1)*(c*c+f*f-1)))/(e*e+g*g+1);
      double m; // TODO - undefined variable!
      double bet1 = g*gam1 - m;
      double bet2 = g*gam2 - m;
      double alp1 = (R*d-bet1*y-gam1*z)/x;
      double alp2 = (R*d-bet2*y-gam2*z)/x;
      */
      
      // Change of plan: Find rotation angles for simple frame
      double rot_X, rot_Z;
      GetRotationAngles(pos.Unit(),rot_Z,rot_X);
      
      // Rotate original position to simple frame
      TVector3 rotpos = pos.Unit();
      rotpos.RotateZ(-rot_X);
      rotpos.RotateY(-rot_Z);
      
      // Rotate original direction to simple frame
      TVector3 rotdir = dir;
      rotdir.RotateZ(-rot_X);
      rotdir.RotateY(-rot_Z);
      
      // Rotate so that everything is in x-z plane (simple frame 2D)
      double angZ = atan(rotdir.Y()/rotdir.X());
      rotdir.RotateZ(-angZ);
      
      // Calculate direction after refraction
      double u = pos.Unit()*dir.Unit(); // cos(theta) in both frames)
      double alp = sqrt(1-pow(ns/ng,2)*(1-u*u));
      double bet = ns/ng*sqrt(1-u*u);
      if (bet*rotdir.X() < 0) bet*=-1; // change sign
      scatdir.SetXYZ(bet,0,alp);
      
      // Check that rotated vectors are all in x-z plane
      if (VERBOSE > 2) {
        printf("-----\n");
        printf("rotpos = (%8.3f %8.3f %8.3f), length=%8.3f\n",rotpos.X(),rotpos.Y(),rotpos.Z(),rotpos.Mag());
        printf("rotdir = (%8.3f %8.3f %8.3f), length=%8.3f\n",rotdir.X(),rotdir.Y(),rotdir.Z(),rotdir.Mag());
        printf("newdir = (%8.3f %8.3f %8.3f), length=%8.3f\n",scatdir.X(),scatdir.Y(),scatdir.Z(),scatdir.Mag());
      }
      
      // Rotate back last step
      rotdir.RotateZ(angZ);
      scatdir.RotateZ(angZ);
      
      // Rotate back to original frame
      rotpos.RotateUz(pos.Unit());
      rotdir.RotateUz(pos.Unit());
      scatdir.RotateUz(pos.Unit());
      
      // Check output vectors
      if (VERBOSE > 2) {
        printf("-----\n");
        printf("rotpos = (%8.3f %8.3f %8.3f), length=%8.3f\n",rotpos.X(),rotpos.Y(),rotpos.Z(),rotpos.Mag());
        printf("rotdir = (%8.3f %8.3f %8.3f), length=%8.3f\n",rotdir.X(),rotdir.Y(),rotdir.Z(),rotdir.Mag());
        printf("newdir = (%8.3f %8.3f %8.3f), length=%8.3f\n",scatdir.X(),scatdir.Y(),scatdir.Z(),scatdir.Mag());
        printf("-----\n");
      }
      
      // Verify output vectors are correct
      if (VERBOSE > 2) {
        printf("CHECK RESULTS\n");
        double d = sqrt(1.-pow(ns/ng,2)*(1.-pow(fabs(prod),2)));
        printf("(1) Snell's law: 0 = %e\n",scatdir*pos.Unit()-d);
        printf("(2) Unit length: 0 = %e\n",scatdir.Mag()-1);
        printf("(3) All planar:  0 = %e\n",scatdir.Dot(pos.Cross(dir)));
        printf("-----\n");
      }
      
      // Set break condition
      exitflask = true;
    }
  }
  
  return scatdir.Unit();
}

// *****************************************************************************
// Get scattering length from theory, given a bubble mass density
double GetScatteringLength(double density) {
  if (VERBOSE > 0) cout << "Provided mass density of glass bubbles: " << density << " g/cm^3" << endl;
  
  // Scattering cross section of glass microspheres
  // https://en.wikipedia.org/wiki/Anomalous_diffraction_theory
  // https://www.gfdl.noaa.gov/wp-content/uploads/files/user_files/pag/lecture2008/lecture3.pdf
  //double n = ns/na; // refractive index (outer vs inner)
  //double p = 4.*pi*r_bub*(n-1.)/lambda; // phase delay of wave through sphere
  //double Q = 2. - 4./p*sin(p) + 4./(p*p)*(1.-cos(p)); // scattering efficiency (same as extinction efficiency)
  double Q = 1; // artificial (ignore ADT)
  double sigma = pi*r_bub*r_bub*Q; // scattering cross section [cm²]
  if (VERBOSE > 0) cout << "Estimated scattering cross section: " << sigma << " cm^2" << endl;
  
  // Number density of glass microspheres (very approximate)
  double V = 4.*pi/3.*(pow(r_bub,3) - pow(r_bub-d_bub,3)); // volume of glass [cm³]
  double m = rho_bub*V; // mass of hollow sphere [g]
  double N = density/m; // number concentration [cm⁻³]
  if (VERBOSE > 0) cout << "Estimated number density of bubbles: " << N << "/cm^3" << endl;
  
  // Scattering parameters
  double eps = N*sigma; // scattering coefficient [cm⁻¹]
  double L = 1./eps; // scattering length [cm]
  if (VERBOSE > 0) cout << "Estimated scattering length: " << L << " cm" << endl;

  return L;
}
