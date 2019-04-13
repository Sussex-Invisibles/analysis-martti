// ---------------------------------------------------------
// Goal:            Plot TELLIE PIN readings vs subruns
// Author:          Martti Nirkko, 07/08/2017
// Compile & run:   clear && g++ -g -o trigger_check.exe trigger_check.C `root-config --cflags --libs` -I$RATROOT/include/libpq -I$RATROOT/include -L$RATROOT/lib -lRATEvent_Linux && ./trigger_check.exe
// ---------------------------------------------------------
#include "../include/HelperFunc.C"
#include <iomanip>
#include <ctime>
#include <ostream>

// Global constants
const int RUN_CLUSTER = 1;  // whether running on cluster (0=local)
const int VERBOSE = 0;      // verbosity flag
const int IS_MC = 0;        // Monte-Carlo flag 
const int NRUNS = 7;        // Number of runs
const int NSUBRUNS = 201;   // Number of subruns

// Initialise functions
using std::vector;
void get_run_data(int, bool, bool);

// Main program
using namespace std;
int main(int argc, char** argv) {
  
  // Test flag (0 = process all runs, else run number)
  const int TEST = (RUN_CLUSTER) ? 0 : 1;
  
  // Loop over all fibres in list
  //string input = "PIN_readings_2018-02-26.txt";
  string input = "PIN_readings_2018-03-12.txt";
  ifstream inp(input.c_str());
  if (!inp) { cerr<<"Failed to open "<<input<<endl; exit(1); }
  string line;
  int runnumbers[NRUNS]={0};
  int subrunhits[NRUNS][NSUBRUNS]={{0}};
  int subrunevs[NRUNS][NSUBRUNS]={{0}};
  int subrunpin[NRUNS][NSUBRUNS]={{0}};
  float subrunrms[NRUNS][NSUBRUNS]={{0}};
  float avgnhit[NRUNS][NSUBRUNS]={{0}};
  int run, subrun, pin;
  float rms;
  float mean[NRUNS]={0.}, err[NRUNS]={0.};
  int nruns=0, nsubruns[NRUNS]={0};
  for (int hdr=0; hdr<2; hdr++) {
    getline(inp,line);      // header
  }
  int lastrun = 0;
  TGraphErrors *gpins[NRUNS];
  TGraph *gmiss[NRUNS];
  TH1D *hdiff[NRUNS];
  for (int r=0; r<NRUNS; r++) {
    gpins[r] = new TGraphErrors();
    gmiss[r] = new TGraph();
    string htitle = Form("hdiff_%d",r);
    hdiff[r] = new TH1D(htitle.c_str(),"",200,0,4);
  }
  while (true) {
    inp >> run >> subrun >> pin >> rms;
    if (!inp.good()) break;
    if (TEST && lastrun!=0 && run!=lastrun) continue; // only want one run
    runnumbers[nruns] = run;
    if (run != lastrun) nruns++;
    if (nruns>NRUNS) break;
    subrunpin[nruns-1][subrun] = pin;
    subrunrms[nruns-1][subrun] = rms;
    if (nruns==1) gpins[nruns-1]->SetPoint(subrun,subrun,pin);
    else gpins[nruns-1]->SetPoint(subrun,subrun+0.5,pin);
    //gpins[nruns-1]->SetPointError(subrun,0,rms); // rms value
    //gpins[nruns-1]->SetPointError(subrun,0,rms/sqrt(NSUBRUNS)); // error on mean
  
    // Do the below only once per run
    if (run == lastrun) continue;
    lastrun = run;
    
    // Extract data from run file
    get_run_data(run, (bool)IS_MC, (bool)TEST);
 
    // Read data from text file
    string out = Form("./output/triggers_%d.out",run);
    ifstream in(out.c_str());
    string line;
    for (int hdr=0; hdr<2; hdr++) {
      getline(in,line);      // header
    }
    
    // Loop over input
    int subrun, gtid, nhit, lastsubrun=-1, lastgtid=-1;
    int trig, trig_tubii;
    long long time, diff, lasttime=0, lastdiff=0;
    bool bitflip=false;
    vector<int> badGTIDs;
    while (true) {
      in >> subrun >> gtid >> trig >> trig_tubii >> nhit >> time;
      if (!in.good()) break;
      if (!(trig & 0x8000)) { // EXTA trigger bit not set
        if(trig & 0x4000000) { // MISS trigger bit set
          cout << "WARNING: No EXTA trigger on GTID " << gtid << ", trigger word is " << TriggerToString(trig) << endl;
          badGTIDs.push_back(gtid);
        }
        continue;
      }
     
      subrunhits[nruns-1][subrun] = nhit;
      subrunevs[nruns-1][subrun]++;
  
      // Time difference between consecutive triggers
      bitflip = false;
      diff = time - lasttime;
      if (subrun!=lastsubrun) diff=0; // first event in subrun
      
      // Fix bit flips in 10 MHz clock, causing negative time differences
      if (lastdiff < 0) {
        bitflip = true;
        diff += lastdiff - 1e6; // fix time difference
        printf("INFO: Fixing bit flip in 10 MHz clock for GTID %7d: new diff is %.3f ms.\n",gtid-1,diff/1e6);
      }
      
      // Print a warning in case of unusual time differences    
      if (diff != 0 && (diff/1e6 > 1.4 || diff/1e6 < 0.9)) { // significant deviation from 1 ms
        printf("WARNING: Time difference between GTID %7d and %7d is %6.3f ms!",lastgtid,gtid,diff/1e6);
        if(gtid-lastgtid==1) { 
          if (bitflip) printf(" => Bit flip on GTID %7d.",gtid-1);
          else printf(" => Possible bit flip...");
        } else {
          if(gtid-lastgtid==2) printf(" => Missed GTID %7d!",gtid-1);
          else printf(" => Missed multiple GTIDs!!!");
        }
        printf("\n");
      }
      
      // Fill graphs
      hdiff[nruns-1]->Fill(diff/1e6);

      lastsubrun = subrun;
      lastgtid = gtid;
      lasttime = time;
      lastdiff = diff;
      
      //if (run==109799) printf("%d %d %lld %lld\n",lastsubrun,lastgtid,lasttime,lastdiff);
    }
    
    // Bad GTIDs
    cout << "Missed " << badGTIDs.size() << " GTIDs: ";
    for (int i=0; i<badGTIDs.size(); i++)
      cout << badGTIDs[i] << ' ';
    cout << endl;
    badGTIDs.clear(); 
  }
  
  // Subrun information
  if(!TEST) printf("*****\n");
  for (int r=0; r<NRUNS; r++) {
    for (int sr=0; sr<NSUBRUNS; sr++) {
      gmiss[r]->SetPoint(sr,sr,5000-subrunevs[r][sr]);
      if(subrunevs[r][sr]==0) {
        //printf("WARN: Run %6d subrun %3d has no events!\n", r, sr);
        continue;
      }
      avgnhit[r][sr] = (float)subrunhits[r][sr]/subrunevs[r][sr];
      if (VERBOSE) printf("Run %6d subrun %3d has %7d events with an average nhit of %5.2f.\n", run, sr, subrunevs[r][sr], avgnhit[r][sr]);
      mean[r] += subrunpin[r][sr];
      err[r] += pow(subrunrms[r][sr],2);
      nsubruns[r]++;
    }
    mean[r] /= NSUBRUNS;
    err[r] = sqrt(err[r]/NSUBRUNS);
    printf("Ran over %d subruns. Resulting PIN readout was %d +- %d.\n",nsubruns[r],(int)round(mean[r]),(int)round(err[r]));
  }
  if(!TEST) printf("*****\n");

  // ----------

  // Plotting section
  TCanvas *c = NULL;

  // Light stability during runs
  c = new TCanvas("c","",800,600);
  c->SetGrid();
  string title = "TELLIE internal readout stability (FT010A);Subrun number;PIN reading";
  c->DrawFrame(-0.025*NSUBRUNS,900,1.025*NSUBRUNS,1300,title.c_str())->GetYaxis()->SetTitleOffset(1.4);
  TLegend *leg = new TLegend(0.74,0.77,0.89,0.89);
  for(int r=2; r<NRUNS; r++) {
    gpins[r]->SetMarkerStyle(7);
    gpins[r]->SetMarkerColor(2*r-2);
    gpins[r]->SetLineColor(2*r-2);
    gpins[r]->Draw("P same");
    leg->AddEntry(gpins[r],std::to_string(runnumbers[r]).c_str(),"LP");
  }
  leg->Draw();
  c->Print("images/pin_vs_subrun.pdf");
  c->Print("images/pin_vs_subrun.png");
  c->Close();
  
  // Missing triggers during runs
  c = new TCanvas("c","",800,600);
  c->SetGrid();
  title = "TELLIE trigger stability (FT010A);Subrun number;Missing triggers";
  c->DrawFrame(-0.025*NSUBRUNS,-0.25,1.025*NSUBRUNS,5.25,title.c_str())->GetYaxis()->SetTitleOffset(1.4);
  leg = new TLegend(0.74,0.77,0.89,0.89);
  for(int r=2; r<NRUNS; r++) {
    gmiss[r]->SetMarkerStyle(8);
    gmiss[r]->SetMarkerColor(2*r-2);
    gmiss[r]->SetLineColor(2*r-2);
    gmiss[r]->Draw("P same");
    leg->AddEntry(gmiss[r],std::to_string(runnumbers[r]).c_str(),"LP");
  }
  leg->Draw();
  c->Print("images/trig_vs_subrun.pdf");
  c->Print("images/trig_vs_subrun.png");
  c->Close();

  // Difference in trigger times
  TText *txt = new TText();
  c = new TCanvas("c","",800,600);
  c->SetGrid();
  c->DrawFrame(-0.1,0.2,4.1,9e6,"Trigger time differences;#Delta t [ms];Counts");//->GetYaxis()->SetLog();
  c->SetLogy();
  leg = new TLegend(0.74,0.77,0.89,0.89);
  for(int r=2; r<NRUNS; r++) {
    hdiff[r]->SetLineWidth(2);
    hdiff[r]->SetLineColor(2*r-2);
    hdiff[r]->Draw("same");
    leg->AddEntry(hdiff[r],std::to_string(runnumbers[r]).c_str(),"LP");
  }
  leg->Draw();
  txt->DrawText(0,1.5*NSUBRUNS,"1st");
  txt->DrawText(0.90,1e4*NSUBRUNS,"EXTA");
  txt->DrawText(2.05,10,"MISS");
  c->Print("images/trig_difftime.pdf");
  c->Print("images/trig_difftime.png");
  c->Close();
  
  // End of program
  if (c) delete c;
  return 0;
}

// Get run data
void get_run_data(int run, bool isMC=false, bool TEST=false) {

  // ********************************************************************
  // Initialisation
  // ********************************************************************
  
  // Check files for given run
  if(!TEST) printf("*****\n");
  //printf("Checking files for run %d... ", run);
  //string fpath = (RUN_CLUSTER) ? "/lustre/scratch/epp/neutrino/snoplus/TELLIE_PCA_RUNS_PROCESSED" : "/home/nirkko/Software/SNOP/data";
  string fpath = (RUN_CLUSTER) ? "/lustre/scratch/epp/neutrino/snoplus/TELLIE_TEST_RUNS_PROCESSED" : "/home/nirkko/Software/SNOP/data";
  string fname = "";
  ifstream f;
  for (int pass=3;pass>=0;pass--) {
    fname = Form("%s/Analysis_r0000%d_s000_p00%d.root",fpath.c_str(),run,pass);
    f.open(fname.c_str());
    if (f.good()) break;
  }
  string out = Form("./output/triggers_%d.out",run);
  //ifstream f(fname.c_str());
  ifstream g(out.c_str());
  int scanned_file = 0;
  if(g.good()) {            // file extracted, but not processed
    printf("not processed! Generating plots for run %d.\n",run);
    scanned_file = 1;
  } else if(!f.good()) {    // file not downloaded
    printf("not downloaded! Skipping run %d.\n",run);
    return;
  } else {                  // file downloaded, but not processed
    printf("OK. Extracting data for run %d.\n",run);
  }
  
  // Output file
  if (scanned_file) return;
  FILE *output = fopen(out.c_str(),"w");
  fprintf(output,"Subrun GTID MTCA_trig TUBii_trig NHits Time\n");
  fprintf(output,"-------------------------------------------\n");

  // Aggro offline mode
  RAT::DB *db = RAT::DB::Get();
  db->SetAirplaneModeStatus(true);
  db->SetDefaultPlaneLockStatus(false);

  // Initialise RAT
  RAT::DU::DSReader dsreader(fname);

  // ********************************************************************
  // Sum PMT hit counts for entire run
  // ********************************************************************
  for(int iEntry=0; iEntry<dsreader.GetEntryCount(); iEntry++) {
    const RAT::DS::Entry& ds = dsreader.GetEntry(iEntry);
    printProgress(iEntry, dsreader.GetEntryCount());

    int subrun = ds.GetSubRunID();             // sub run this event belongs to
    
    // Loop over triggered events in each entry
    for(int iEv=0; iEv<ds.GetEVCount(); iEv++) { // mostly 1 event per entry
      const RAT::DS::EV& ev = ds.GetEV(iEv);
      
      // Trigger type
      int trig = ev.GetTrigType();        // MTCA+ trigger word
      int tubii = ev.GetTubiiTrig();      // TUBii trigger word
     
      // Event observables
      int gtid = ev.GetGTID();
      int nhits = ev.GetNhitsCleaned();   // calibrated PMT hits with removed crosstalk
      
      // Timing info
      RAT::DS::UniversalTime ut = ev.GetUniversalTime();
      unsigned int days = ut.GetDays();
      unsigned int secs = ut.GetSeconds();
      double nano = ut.GetNanoSeconds();
      long long fulltime = (86400*days + secs)*1e9 + nano;
      //cout << days << ", " << secs << ", " << nano << " = " << fulltime << " => ";
      //std::tm time = ut.GetTime();
      //cout << asctime(&time);
      
      // Write to file
      fprintf(output,"%3d %6d %10d %10d %4d %lld\n",subrun,gtid,trig,tubii,nhits,fulltime);
      
    } // event loop
    
  } // entry loop
  
  fclose(output);
  return;
}
