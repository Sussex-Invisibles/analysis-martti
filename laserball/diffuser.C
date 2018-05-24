// ---------------------------------------------------------
// Goal:          Simulate photon tracking through laserball diffuser flask
// Author:        Martti Nirkko, 24/05/2018
// Compile & run: clear && g++ -g -o diffuser.exe diffuser.C `root-config --cflags --libs` -I$RATROOT/include/libpq -I$RATROOT/include -L$RATROOT/lib -lRATEvent_Linux && ./diffuser.exe
// ---------------------------------------------------------

// Helper functions (includes everything else)
#include "../HelperFunc.C"

// Global constants
const int VERBOSE = 0;

// Laser
const double I0 = 1.; // laser intensity, normalised for now
const double NA = 0.2; // numerical aperture for fibre (from R. Ford's thesis)

// Diffuser ball - TODO REVIEW ALL PARAMETERS
const double flaskdiam = 100.; // diameter of the diffuser flask [mm]
const double R = flaskdiam/2.; // radius
const double na = 1.0; // refractive index (air)
const double ns = 1.4; // refractive index (SilGel)
const double ng = 1.5; // refractive index (glass)
const double lambda = 1.; // optical property (depends on density)

// Quartz rod
const double rodlen = 200.; // length [mm]
const double roddiam = 1.; // diameter [mm]

// *****************************************************************************
// Function declarations
void diffuser();
double distance_to_wall(TVector3, TVector3);
double reflection_prob(double n1, double n2, double impact);
TVector3 propagate(TVector3, TVector3, double distance);
TVector3 scatter(TVector3, TVector3, double impact);
TVector3 reflect_or_refract(TVector3, TVector3, double impact, bool& exitflask);

// Random number generator
TRandom3* gen = new TRandom3();

// *****************************************************************************
// Main program
int main(int argc, char** argv) {
  diffuser();
  return 0;
}

// *****************************************************************************
// Track a single photon through the diffuser flask
// Inputs: laser intensity, quartz rod position, density of glass bubbles
// Output: intensity as function of solid angle
void diffuser() {
  gen->SetSeed(0);
  
  // Photon position at injection point
  TVector3 pos(1,0,0);
  double A = pi*pow(roddiam/2.,2)*gen->Rndm(); // random within area
  pos.SetMag(sqrt(A/pi)); // radius [mm]
  pos.SetPhi(2*pi*gen->Rndm()); // angle [rad]
  
  // Photon direction at injection point
  TVector3 dir(0,0,1);
  dir.SetPhi(2*pi*gen->Rndm()); // azimuthal direction [0,2pi)
  TF1* aperture = new TF1("aperture","cos(pi/2*x/[0])",0,NA);
  aperture->SetParameter(0,NA);
  double theta = aperture->GetRandom();
  dir.SetTheta(pi-theta); // downwards zenith direction [0,NA)
  
  // Intensity dropoff
  TF1* inten = new TF1("intensity","[0]*exp(-x/[1])",0,R/10);
  inten->SetParameter(0,I0);
  inten->SetParameter(1,lambda);
  
  /* If the selected distance is greater than d it is checked whether the photon
exits the diffuser or internally reflects. For a selected distance less than d a
scattering event will take place. */
  int step = 0;
  TVector3 newpos, newdir;
  TGraph track;
  track.SetPoint(step,pos.X(),pos.Z()); // initial point
  printf("ENTERING diffuser at (%8.3f %8.3f %8.3f), R=%6.3f\n",step,pos.X(),pos.Y(),pos.Z(),pos.Mag());
  bool exitflask = false;
  while (!exitflask) {
    if (VERBOSE) printf("position (%8.3f %8.3f %8.3f), direction = (%8.3f %8.3f %8.3f)\n",pos.X(),pos.Y(),pos.Z(),dir.X(),dir.Y(),dir.Z());
    double y = gen->Rndm();
    double dsel = -lambda*log(y);  // randomly selected distance in diffuser [mm]
    double dwall = distance_to_wall(pos,dir); // distance to flask wall [mm]
    if (dsel < dwall) {
      // propagate selected distance
      newpos = propagate(pos,dir,dsel);
      // scatter in diffuser
      newdir = scatter(newpos,dir,sqrt(y));
      step++;
    } else {
      // propagate up to wall
      newpos = propagate(pos,dir,dwall);
      // exit or reflect internally
      newdir = reflect_or_refract(newpos,dir,sqrt(y),exitflask);
    }
    pos = newpos;
    dir = newdir;
    track.SetPoint(step,pos.X(),pos.Z()); // ignore y for now
  }
  
  printf("EXITING diffuser after %d scatters at (%8.3f %8.3f %8.3f), R=%6.3f\n",step,pos.X(),pos.Y(),pos.Z(),pos.Mag());
  TVector3 endpos = propagate(pos,dir,6e3);
  track.SetPoint(step+1,endpos.X(),endpos.Z()); // ignore y for now
  printf("End point (%8.3f %8.3f %8.3f), R=%6.3f\n",endpos.X(),endpos.Y(),endpos.Z(),endpos.Mag());
  
  double LIM=75;
  TCanvas c("c","",800,800);
  c.DrawFrame(-LIM,-LIM,LIM,LIM,"Photon scattering in diffuser;X [mm];Z [mm]");
  // Flask
  TEllipse ball(0,0,R);
  ball.Draw("same");
  TBox rod(-roddiam/2,0,roddiam/2,LIM); // must be within limits of frame
  rod.SetFillColor(0);
  rod.SetLineColor(4);
  rod.Draw("L same");
  track.SetLineColor(2);
  track.Draw("L same");
  c.Print("diffuser.png");
  c.Close();
  
}

// *****************************************************************************
// Calculate distance to inner flask boundary
double distance_to_wall(TVector3 pos, TVector3 dir) {
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
double reflection_prob(double n1, double n2, double b) {
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
TVector3 propagate(TVector3 pos, TVector3 dir, double dist) {
  double x = pos.X()+dir.X()*dist;
  double y = pos.Y()+dir.Y()*dist;
  double z = pos.Z()+dir.Z()*dist;
  TVector3 newpos(x,y,z);
  return newpos;
}

// *****************************************************************************
// Scatter off a glass bubble
TVector3 scatter(TVector3 pos, TVector3 dir, double impact) {
  
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
  double newbeta  = beta*costh + param*(beta*gamma*cos(phi)-alpha*sin(phi));
  double newgamma = gamma*costh - param*(1.-gamma*gamma)*cos(phi);  
  TVector3 newdir(newalpha,newbeta,newgamma);
  return newdir.Unit();
}

// *****************************************************************************
// Reflect at flask wall, or refract out of diffuser
TVector3 reflect_or_refract(TVector3 pos, TVector3 dir, double impact, bool& exitflask) {

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
  TVector3 newdir;
  if (sinOmega > 1) {
    // (A) internal reflection
    newdir = dir-2.*prod*pos;
    exitflask = false;
  } else {
    // (B) probability of reflecting despite sinOmega <= 1
    double prob = reflection_prob(ns,ng,impact);
    double roll = gen->Rndm();
    if (roll < prob) {
      newdir = dir-2.*prod*pos;
      exitflask = false;
      printf("REFLECTING on diffuser wall at (%8.3f %8.3f %8.3f), R=%6.3f\n",x,y,z,pos.Mag());
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
      newdir.SetXYZ(bet,0,alp);
      
      // Check that rotated vectors are all in x-z plane
      if (VERBOSE) {
        printf("-----\n");
        printf("rotpos = (%8.3f %8.3f %8.3f), length=%8.3f\n",rotpos.X(),rotpos.Y(),rotpos.Z(),rotpos.Mag());
        printf("rotdir = (%8.3f %8.3f %8.3f), length=%8.3f\n",rotdir.X(),rotdir.Y(),rotdir.Z(),rotdir.Mag());
        printf("newdir = (%8.3f %8.3f %8.3f), length=%8.3f\n",newdir.X(),newdir.Y(),newdir.Z(),newdir.Mag());
      }
      
      // Rotate back last step
      rotdir.RotateZ(angZ);
      newdir.RotateZ(angZ);
      
      // Rotate back to original frame
      rotpos.RotateUz(pos.Unit());
      rotdir.RotateUz(pos.Unit());
      newdir.RotateUz(pos.Unit());
      
      // Check output vectors
      if (VERBOSE) {
        printf("-----\n");
        printf("rotpos = (%8.3f %8.3f %8.3f), length=%8.3f\n",rotpos.X(),rotpos.Y(),rotpos.Z(),rotpos.Mag());
        printf("rotdir = (%8.3f %8.3f %8.3f), length=%8.3f\n",rotdir.X(),rotdir.Y(),rotdir.Z(),rotdir.Mag());
        printf("newdir = (%8.3f %8.3f %8.3f), length=%8.3f\n",newdir.X(),newdir.Y(),newdir.Z(),newdir.Mag());
        printf("-----\n");
      }
      
      // Verify output vectors are correct
      if (VERBOSE) {
        printf("CHECK RESULTS\n");
        double d = sqrt(1.-pow(ns/ng,2)*(1.-pow(fabs(prod),2)));
        printf("(1) Snell's law: 0 = %e\n",newdir*pos.Unit()-d);
        printf("(2) Unit length: 0 = %e\n",newdir.Mag()-1);
        printf("(3) All planar:  0 = %e\n",newdir.Dot(pos.Cross(dir)));
        printf("-----\n");
      }
      
      // Set break condition
      exitflask = true;
    }
  }
  
  return newdir.Unit();
}

