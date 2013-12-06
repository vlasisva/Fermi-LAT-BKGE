//Author: Vlasios Vasileiou <vlasisva@gmail.com>
#include "BackgroundEstimator/BKGE_Tools.h"
#include "TGraphAsymmErrors.h"


int TOOLS::Make_Burst_Plots(string DataClass, string FT1_FILE, string GRB_DIR, double RA_BURST, double DEC_BURST, double GRB_t0, 
      double Burst_Dur, TH1F* hROI, short int verbosity) {
  const double ENERGY_MIN=hROI->GetXaxis()->GetXmin();
  const double ENERGY_MAX=hROI->GetXaxis()->GetXmax();
  const int    ENERGY_BINS = hROI->GetNbinsX();

  char name[1000];
  double Burst_t1 = GRB_t0 + Burst_Dur;
  int MinCTBClassLevel  =GetCTBClassLevel(DataClass);
  int ConversionType    =GetConversionType(DataClass);
  string ConversionName = GetConversionName(DataClass);

  gROOT->cd("/");
  int iebin;
 
  int TrueEvents[ENERGY_BINS+2];
  memset(TrueEvents,0,sizeof(int)*(ENERGY_BINS+2));
  astro::SkyDir SCBurst = astro::SkyDir(RA_BURST, DEC_BURST ,astro::SkyDir::EQUATORIAL);
  astro::SkyDir SCEvent;

  float radius[ENERGY_BINS+2];
  for (int i=1;i<=ENERGY_BINS;i++) {
       radius[i] = hROI->GetBinContent(i)* DEG_TO_RAD;
  }

  TH1F * hCtsvsEnergy = (TH1F*)gROOT->Get("hCtsvsEnergy_Burst");
  if (hCtsvsEnergy) delete hCtsvsEnergy;  hCtsvsEnergy = new TH1F("hCtsvsEnergy_Burst","hCtsvsEnergy_Burst",ENERGY_BINS,ENERGY_MIN,ENERGY_MAX);
  hCtsvsEnergy->GetXaxis()->SetTitle("log_{10}(Energy/MeV)");
  hCtsvsEnergy->GetYaxis()->SetTitle("Counts/bin");

  fitsfile *fptr; 

  int EventClassMask=0;
  if (DataClass.find("P7")!=string::npos) {
     if (DataClass.find("TRANSIENT")!=string::npos) EventClassMask=int(pow(2.,0.));
  }

  
      int status = 0,anynul,hdutype;
      fits_open_file(&fptr, FT1_FILE.c_str(), READONLY, &status);
      if (status)  { fits_report_error(stderr, status);printf("%s: Can't open file %s\n",__FUNCTION__,FT1_FILE.c_str()); throw std::runtime_error("");}
      
      fits_movabs_hdu(fptr, 2, &hdutype, &status);if (status) fits_report_error(stderr, status);
      long nrows;int ncols;
      fits_get_num_rows(fptr, &nrows, &status); if (status) fits_report_error(stderr, status);
      fits_get_num_cols(fptr, &ncols, &status); if (status) fits_report_error(stderr, status);
      int format;
      if (DataClass.find("P7")!=string::npos) format=2;
      else if (ncols==22) format=1;
      else if (ncols==18) format=0;
      else {printf("%s: Unknown FITS file format file %s ncols=%d\n",__FUNCTION__,name,ncols); throw std::runtime_error("");}
      int col_time; fits_get_colnum(fptr, TRUE, "TIME", &col_time, &status);if (status) fits_report_error(stderr, status);
      int col_conv_type;fits_get_colnum(fptr, TRUE, "CONVERSION_TYPE", &col_conv_type, &status);if (status) fits_report_error(stderr, status);
      int col_energy;   fits_get_colnum(fptr, TRUE, "ENERGY", &col_energy, &status);if (status) fits_report_error(stderr, status);
      int col_evclass;   fits_get_colnum(fptr, TRUE, "EVENT_CLASS", &col_evclass, &status);if (status) fits_report_error(stderr, status);
      int col_ra;   fits_get_colnum(fptr, TRUE, "RA", &col_ra, &status); if (status) fits_report_error(stderr, status);
      int col_dec;   fits_get_colnum(fptr, TRUE, "DEC", &col_dec, &status); if (status) fits_report_error(stderr, status);

      /*
      int istart=-1;
      //perform a binary search to find starting event
      int min=0,max=nrows;

      
      while (1) {
            double PtTime;
            istart=(max+min)/2;
            //printf("%d %d istart=%d\n",min,max,istart);
            //printf("pttime=%f t0=%f\n",PtTime,GRB_t0);
            if (istart==0){
               istart=1;
               break;  //this is when the first event of the file is detected after our start time
            }  
            fits_read_col (fptr,TDOUBLE,col_time,istart, 1, 1, NULL,&PtTime, &anynul, &status); if (status) fits_report_error(stderr, status);
            if (PtTime>GRB_t0) max=istart; //go down
            else {
                if (istart==(nrows-1)) {printf("%s: istart==nrows? FT1 file seems to end before requested time interval\n",__FUNCTION__);throw std::runtime_error("");}
                fits_read_col (fptr,TDOUBLE,col_time,istart+1, 1, 1, NULL,&PtTime, &anynul, &status);if (status) fits_report_error(stderr, status);
                if (PtTime>GRB_t0) {istart++; break;} //found!
                else min=istart;//go up
            }
      }     
      */

      //a comment on the commented out code above and in the two lines below.. If the events were sorted in the FT1 file, then some parts would have been far faster
      //now we need to scan the whole file...
      
      for (long i=1;i<=nrows;i++) {
         double PtTime;
         status=0;
         fits_read_col (fptr,TDOUBLE,col_time,i, 1, 1, NULL,&PtTime, &anynul, &status); if (status) fits_report_error(stderr, status);

         //if (PtTime<GRB_t0) {printf("cont\n");continue;}
         //if (PtTime>Burst_t1) {printf("break\n");break;
         //CHECK CUTS
         if (PtTime<GRB_t0 || PtTime>Burst_t1) continue;
         
         if (format==2)  {
            int EventClass=0;
            fits_read_col (fptr,TINT,col_evclass,i, 1, 1, NULL,&EventClass, &anynul, &status); if (status) fits_report_error(stderr, status);
            //printf ("%d %d %d\n",EventClass,EventClassMask,EventClass%16);
            if (EventClass%16<EventClassMask) continue;
         }
         else {
            int CTBClassLevel;
            if      (format==0) fits_read_col (fptr,TINT,18,i, 1, 1, NULL,&CTBClassLevel, &anynul, &status); if (status) fits_report_error(stderr, status);
            else if (format==1) fits_read_col (fptr,TINT,15,i, 1, 1, NULL,&CTBClassLevel, &anynul, &status); if (status) fits_report_error(stderr, status);
            if (CTBClassLevel<MinCTBClassLevel) continue;
         }
         if (ConversionType!=-1) { //if !BOTH
            int aConversionType=0;
            fits_read_col (fptr,TINT,col_conv_type,i, 1, 1, NULL,&aConversionType, &anynul, &status); if (status) fits_report_error(stderr, status);
            if (aConversionType!=ConversionType) continue; //apply Conversion Type front/back
         }

         double FT1Energy;
         fits_read_col (fptr,TDOUBLE,col_energy,i, 1, 1, NULL,&FT1Energy, &anynul, &status); if (status) fits_report_error(stderr, status);
         FT1Energy=log10(FT1Energy);
         if (FT1Energy<ENERGY_MIN || FT1Energy>ENERGY_MAX) continue; //energy cut
  
         /////////////////////////////////////////////////////////////////////////

         double FT1RA,FT1DEC;
         fits_read_col (fptr,TDOUBLE,col_ra,i, 1, 1, NULL,&FT1RA, &anynul, &status); if (status) fits_report_error(stderr, status);
         fits_read_col (fptr,TDOUBLE,col_dec,i, 1, 1, NULL,&FT1DEC, &anynul, &status); if (status) fits_report_error(stderr, status);

         iebin=hROI->FindBin(FT1Energy);

         SCEvent = astro::SkyDir(FT1RA,FT1DEC,astro::SkyDir::EQUATORIAL);
         double EvAngDistance = SCBurst.difference(SCEvent);
         if (EvAngDistance>TMath::Pi()) EvAngDistance=TMath::Pi()-EvAngDistance;

         if (EvAngDistance>radius[iebin]) continue; 
         TrueEvents[iebin]++;
     }
     fits_close_file(fptr, &status);
   

  for (int ie=0;ie<=ENERGY_BINS+1;ie++) {
      hCtsvsEnergy->SetBinContent(ie,TrueEvents[ie]);
  }

  if (GRB_DIR!="") { //Save stuff 
     float SigError[2][ENERGY_BINS];
     memset(SigError,0,sizeof(SigError));
     float EnergyData[3][ENERGY_BINS];
     float TrueEvents_middlebins[ENERGY_BINS];
     for (int ie=0;ie<ENERGY_BINS+2;ie++) {
        if (ie>0 && ie<=ENERGY_BINS) { //no under/overflows
              TrueEvents_middlebins[ie-1]=TrueEvents[ie];
              EnergyData[0][ie-1]=hCtsvsEnergy->GetXaxis()->GetBinCenter(ie)-hCtsvsEnergy->GetXaxis()->GetBinLowEdge(ie);
              EnergyData[1][ie-1]=hCtsvsEnergy->GetXaxis()->GetBinCenter(ie);
              EnergyData[2][ie-1]=hCtsvsEnergy->GetXaxis()->GetBinUpEdge(ie)-hCtsvsEnergy->GetXaxis()->GetBinCenter(ie);
              SigError[1][ie-1]=PoissonErrorBar(1,TrueEvents[ie]);
              SigError[0][ie-1]=PoissonErrorBar(0,TrueEvents[ie]);
          }
     }

     TGraphAsymmErrors * gSignal = new TGraphAsymmErrors(ENERGY_BINS,EnergyData[1],TrueEvents_middlebins,
                      EnergyData[0],EnergyData[2],SigError[0],SigError[1]);
     gSignal->SetTitle("gSignal");

     char ResultsFilename[1000];
     sprintf(ResultsFilename,"%s/%s_sig_%.0f_%.0f.root",GRB_DIR.c_str(),DataClass.c_str(),pow(10,ENERGY_MIN),pow(10,ENERGY_MAX));
     TFile * fResults = new TFile(ResultsFilename,"RECREATE");

     gSignal->Write("gSignal");
     hROI->Write("hROI");
     hCtsvsEnergy->Write();

     sprintf(name,"RA/DEC %.3f %.3f",RA_BURST,DEC_BURST);
     TNamed Data = TNamed("Localization_Data",name);
     Data.Write();

     fResults->Close();
     
  }
  int SigEvents=0;
  for (int i=1;i<=ENERGY_BINS;i++) SigEvents+=TrueEvents[i];
  return SigEvents;
}

