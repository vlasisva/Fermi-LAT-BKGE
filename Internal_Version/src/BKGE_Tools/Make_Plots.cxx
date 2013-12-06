//Author: Vlasios Vasileiou <vlasisva@gmail.com>
#include "BackgroundEstimator/BKGE_Tools.h"
#include "TGraph.h"

void Interpolate(TH1F * hist, bool aRA=false, int verbosity=0);
void CalculateRockingAngle(TH1F* hRockingAnglevsTime, TH1F* hPtRazvsTime, TH1F* hPtDeczvsTime,  TH1F* hPtRaxvsTime, TH1F* hPtDecxvsTime, TH1F* hRAZenithvsTime, TH1F* hDecZenithvsTime);

int TOOLS::Make_Plots(double PreTime, double PostTime, double GRB_t0, string PlotsFile, string FT2_FILE, int verbosity, float TimeStep){

 if (PlotsFile=="") PlotsFile= GetS("OUTPUT_DIR")+"/Plots.root";
/*
 FILE * ftemp;
 FILE * ftemp = fopen(PlotsFile.c_str(),"r");
 if (ftemp) {
    fclose(ftemp); 
    printf("%s: File already exists %s, skipping..\n",__FUNCTION__,PlotsFile.c_str());
    return 0;
 }
*/
 double StartTime_orig,StopTime_orig;
 StartTime_orig=GRB_t0-PreTime;
 StopTime_orig =GRB_t0+PostTime;

 long nrows;
 int status=0,hdutype,anynul;
 fitsfile *fptr;

 double FT2_START,FT2_END;

 FILE* ftemp = fopen(FT2_FILE.c_str(),"r");
 if (ftemp) fclose (ftemp);
 else {
    printf("%s: Can't open FT2 file %s (file not found?)\n",__FUNCTION__,FT2_FILE.c_str()); 
    return -1;
 }


 //Open first file
 fits_open_file(&fptr, FT2_FILE.c_str(), READONLY, &status);
 fits_movabs_hdu(fptr, 2, &hdutype, &status);
 fits_get_num_rows(fptr, &nrows, &status);
 int col_start; fits_get_colnum(fptr, TRUE, "START", &col_start, &status);
 int col_stop; fits_get_colnum(fptr, TRUE, "STOP", &col_stop, &status);
 int col_razen;  fits_get_colnum(fptr, TRUE, "RA_ZENITH", &col_razen, &status);
 int col_deczen;  fits_get_colnum(fptr, TRUE, "DEC_ZENITH", &col_deczen, &status);
 int col_mcilwainl; fits_get_colnum(fptr, TRUE, "L_MCILWAIN", &col_mcilwainl, &status);
 int col_rascz; fits_get_colnum(fptr, TRUE, "RA_SCZ", &col_rascz, &status);
 int col_decscz; fits_get_colnum(fptr, TRUE, "DEC_SCZ", &col_decscz, &status);
 int col_rascx; fits_get_colnum(fptr, TRUE, "RA_SCX", &col_rascx, &status);
 int col_decscx; fits_get_colnum(fptr, TRUE, "DEC_SCX", &col_decscx, &status);

 fits_read_col (fptr,TDOUBLE,col_start,1, 1, 1, NULL,&FT2_START, &anynul, &status);
 fits_read_col (fptr,TDOUBLE,col_stop,nrows, 1, 1, NULL,&FT2_END, &anynul, &status);

 //decide if we have FT2 or FT2Seconds file
 double a1,a2;
 fits_read_col (fptr,TDOUBLE,col_start,1, 1, 1, NULL,&a1, &anynul, &status);
 fits_read_col (fptr,TDOUBLE,col_stop,1, 1, 1, NULL,&a2, &anynul, &status);
 int TimeExtension=0; //in sec
 bool FT2SECONDS=true;

 if (fabs((a2-a1)-1)>1.1) {
    TimeExtension=300; 
    FT2SECONDS=false;
    if (verbosity>1) printf("%s: Reading 30sec-step FT2 file %s\n",__FUNCTION__,FT2_FILE.c_str());
 }
 else if (verbosity>1) { printf("%s: Reading 1-sec step FT2 file %s\n",__FUNCTION__,FT2_FILE.c_str());}

 const double StartTime=StartTime_orig-TimeExtension;
 const double StopTime =StopTime_orig+TimeExtension;
 const int TimeBins      = (int)ceil((StopTime-StartTime)/TimeStep);
 const int TimeBins_orig = (int)ceil((StopTime_orig-StartTime_orig)/TimeStep);
 //this assumes we have 1sec time bins, will break otherwise

 TH1F hPtRazvsTime =  TH1F("hPtRazvsTime","PtRaz vs Time", TimeBins,StartTime,StopTime);
 TH1F hPtRaxvsTime = TH1F("hPtRaxvsTime","PtRax vs Time", TimeBins,StartTime,StopTime);
 TH1F hPtDeczvsTime = TH1F("hPtDeczvsTime","PtDecz vs Time", TimeBins,StartTime,StopTime);
 TH1F hPtDecxvsTime = TH1F("hPtDecxvsTime","PtDecx vs Time", TimeBins,StartTime,StopTime);
 TH1F hMcIlwainLvsTime = TH1F("hMcIlwainLvsTime","McIlwainLvsTime",TimeBins,StartTime,StopTime);
 TH1F hRAZenithvsTime = TH1F("hRAZenithvsTime","RAZenithvsTime",TimeBins,StartTime,StopTime);
 TH1F hDecZenithvsTime = TH1F("hDecZenithvsTime","DecZenithvsTime",TimeBins,StartTime,StopTime);

 
 if (FT2_START>hRAZenithvsTime.GetXaxis()->GetXmin()) {
  printf("%s: FT2 file starts at time %f and we need data from an earlier time %f to make the plots\n",
       __FUNCTION__,FT2_START,StartTime);
  throw std::runtime_error("");
 } 

 fits_read_col (fptr,TDOUBLE,col_stop,nrows, 1, 1, NULL,&FT2_END, &anynul, &status);
 if (FT2_END<hRAZenithvsTime.GetXaxis()->GetXmax()) {
      printf("%s: FT2 file stops at time %f and we need data up to a later time %f to make the plots\n",
           __FUNCTION__,FT2_END,StopTime);
     throw std::runtime_error("");
 } 

 int NBins=hRAZenithvsTime.GetNbinsX();
 double START,RA_ZENITH,DEC_ZENITH,PTRAZ,PTRAX,PTDECZ,PTDECX,MCILWAINL,STOP;
 
 for (long jj = 1; jj < nrows && !status; jj++) {
     int col_num;
     //if (jj%100==0){ printf("%5.2f\r",jj/float(nrows)); fflush(0);}

     fits_read_col (fptr,TDOUBLE,col_start,jj, 1, 1, NULL,&START, &anynul, &status);
     if (status) fits_report_error(stderr, status); 
     fits_read_col (fptr,TDOUBLE,col_stop,jj, 1, 1, NULL,&STOP, &anynul, &status);
     if (status) fits_report_error(stderr, status);
     long int ibin = hRAZenithvsTime.FindBin(START);
     long int ibinSTOP = hRAZenithvsTime.FindBin(STOP);

     if (ibinSTOP<1) continue;
     if (ibin>NBins) break;
    //Format of FT2 files seems to change all the time so let's try to be cautious...    

    fits_read_col (fptr,TDOUBLE,col_razen,jj, 1, 1, NULL,&RA_ZENITH, &anynul, &status);
    fits_read_col (fptr,TDOUBLE,col_deczen,jj, 1, 1, NULL,&DEC_ZENITH, &anynul, &status);
    fits_read_col (fptr,TDOUBLE,col_mcilwainl,jj, 1, 1, NULL,&MCILWAINL, &anynul, &status);
    fits_read_col (fptr,TDOUBLE,col_rascz,jj, 1, 1, NULL,&PTRAZ, &anynul, &status);
    fits_read_col (fptr,TDOUBLE,col_decscz,jj, 1, 1, NULL,&PTDECZ, &anynul, &status);
    fits_read_col (fptr,TDOUBLE,col_rascx,jj, 1, 1, NULL,&PTRAX, &anynul, &status);
    fits_read_col (fptr,TDOUBLE,col_decscx,jj, 1, 1, NULL,&PTDECX, &anynul, &status);    

    hPtRazvsTime.SetBinContent(ibin,PTRAZ);
    hPtRaxvsTime.SetBinContent(ibin,PTRAX);
    hPtDeczvsTime.SetBinContent(ibin,PTDECZ);
    hPtDecxvsTime.SetBinContent(ibin,PTDECX);
    hRAZenithvsTime.SetBinContent(ibin,RA_ZENITH);
    hDecZenithvsTime.SetBinContent(ibin,DEC_ZENITH);
    hMcIlwainLvsTime.SetBinContent(ibin,MCILWAINL);
    //Fill both START and STOP. If START of next event 
    //printf("ibin=%d ibinistop=%d\n",ibin,ibinSTOP);
    if (ibinSTOP!=ibin) {
      hPtRazvsTime.SetBinContent(ibinSTOP,PTRAZ);
      hPtRaxvsTime.SetBinContent(ibinSTOP,PTRAX);
      hPtDeczvsTime.SetBinContent(ibinSTOP,PTDECZ);
      hPtDecxvsTime.SetBinContent(ibinSTOP,PTDECX);
      hRAZenithvsTime.SetBinContent(ibinSTOP,RA_ZENITH);
      hDecZenithvsTime.SetBinContent(ibinSTOP,DEC_ZENITH);
      hMcIlwainLvsTime.SetBinContent(ibinSTOP,MCILWAINL);
    }
    if ((ibinSTOP-ibin)==2 && FT2SECONDS ) { //rounding errors fix
      hPtRazvsTime.SetBinContent(ibinSTOP-1,PTRAZ);
      hPtRaxvsTime.SetBinContent(ibinSTOP-1,PTRAX);
      hPtDeczvsTime.SetBinContent(ibinSTOP-1,PTDECZ);
      hPtDecxvsTime.SetBinContent(ibinSTOP-1,PTDECX);
      hRAZenithvsTime.SetBinContent(ibinSTOP-1,RA_ZENITH);
      hDecZenithvsTime.SetBinContent(ibinSTOP-1,DEC_ZENITH);
      hMcIlwainLvsTime.SetBinContent(ibinSTOP-1,MCILWAINL);
    }
    if (status) fits_report_error(stderr, status);
    //printf("%d %d %lf %lf %lf %lf j=%ld stat=%d %lf\n",ibin,ibinSTOP,START,STOP,hPtRazvsTime.GetXaxis()->GetBinLowEdge(3),hPtRazvsTime.GetXaxis()->GetBinUpEdge(3),jj,status,PTRAZ);

 }
 fits_close_file(fptr, &status);
 if (verbosity>2) printf("%s: Creating plots file %s\n",__FUNCTION__,PlotsFile.c_str());
 TFile * fout = new TFile(PlotsFile.c_str(),"RECREATE");

 if (!FT2SECONDS) {
   hPtRazvsTime.Write("hPtRazvsTime_uninterpolated");
   hPtRaxvsTime.Write("hPtRaxvsTime_uninterpolated");
   hPtDeczvsTime.Write("hPtDeczvsTime_uninterpolated");
   hPtDecxvsTime.Write("hPtDecxvsTime_uninterpolated");
   hRAZenithvsTime.Write("hRAZenithvsTime_uninterpolated");
   hDecZenithvsTime.Write("hDecZenithvsTime_uninterpolated");
   hMcIlwainLvsTime.Write("hMcIlwainLvsTime_uninterpolated");

   Interpolate(&hPtRazvsTime,true,verbosity);
   Interpolate(&hPtRaxvsTime,true,verbosity);
   Interpolate(&hPtDeczvsTime,false,verbosity);
   Interpolate(&hPtDecxvsTime,false,verbosity);
   Interpolate(&hRAZenithvsTime,true,verbosity);
   Interpolate(&hDecZenithvsTime,false,verbosity);
   Interpolate(&hMcIlwainLvsTime,false,verbosity);


   hPtRazvsTime.Write("hPtRazvsTime_interpolated");
   hPtRaxvsTime.Write("hPtRaxvsTime_interpolated");
   hPtDeczvsTime.Write("hPtDeczvsTime_interpolated");
   hPtDecxvsTime.Write("hPtDecxvsTime_interpolated");
   hRAZenithvsTime.Write("hRAZenithvsTime_interpolated");
   hDecZenithvsTime.Write("hDecZenithvsTime_interpolated");
   hMcIlwainLvsTime.Write("hMcIlwainLvsTime_interpolated");
 

   TH1F _hPtRazvsTime = TH1F("hPtRazvsTime","PtRaz vs Time", TimeBins_orig,StartTime_orig,StopTime_orig);
   TH1F _hPtRaxvsTime = TH1F("hPtRaxvsTime","PtRax vs Time", TimeBins_orig,StartTime_orig,StopTime_orig);
   TH1F _hPtDeczvsTime = TH1F("hPtDeczvsTime","PtDecz vs Time", TimeBins_orig,StartTime_orig,StopTime_orig);
   TH1F _hPtDecxvsTime = TH1F("hPtDecxvsTime","PtDecx vs Time", TimeBins_orig,StartTime_orig,StopTime_orig);
   TH1F _hMcIlwainLvsTime = TH1F("hMcIlwainLvsTime","McIlwainLvsTime",TimeBins_orig,StartTime_orig,StopTime_orig);
   TH1F _hRAZenithvsTime = TH1F("hRAZenithvsTime","RAZenithvsTime",TimeBins_orig,StartTime_orig,StopTime_orig);
   TH1F _hDecZenithvsTime = TH1F("hDecZenithvsTime","DecZenithvsTime",TimeBins_orig,StartTime_orig,StopTime_orig);
   TH1F _hRockingAnglevsTime  = TH1F("hRockingAnglevsTime","RockingAnglevsTime",TimeBins_orig,StartTime_orig,StopTime_orig);

   int ib;
   double atime;
   for (int i=1;i<=TimeBins_orig;i++) { //copy stuff from the extended filled plots to the short one
       atime=_hPtRazvsTime.GetBinCenter(i);
       ib   = hPtRazvsTime.FindBin(atime);
       _hPtRazvsTime.SetBinContent(i,hPtRazvsTime.GetBinContent(ib));
       _hPtRaxvsTime.SetBinContent(i,hPtRaxvsTime.GetBinContent(ib));
       _hPtDeczvsTime.SetBinContent(i,hPtDeczvsTime.GetBinContent(ib));
       _hPtDecxvsTime.SetBinContent(i,hPtDecxvsTime.GetBinContent(ib));
       _hRAZenithvsTime.SetBinContent(i,hRAZenithvsTime.GetBinContent(ib));
       _hDecZenithvsTime.SetBinContent(i,hDecZenithvsTime.GetBinContent(ib));
       _hMcIlwainLvsTime.SetBinContent(i,hMcIlwainLvsTime.GetBinContent(ib));
   }

   _hMcIlwainLvsTime.Write();
   _hDecZenithvsTime.Write();
   _hRAZenithvsTime.Write();
   _hPtDecxvsTime.Write();
   _hPtDeczvsTime.Write();
   _hPtRazvsTime.Write();
   _hPtRaxvsTime.Write();
   CalculateRockingAngle(&_hRockingAnglevsTime,&_hPtRazvsTime,&_hPtDeczvsTime, &_hPtRaxvsTime, &_hPtDecxvsTime, &_hRAZenithvsTime, &_hDecZenithvsTime);
   _hRockingAnglevsTime.Write();
  }
 else {
   TH1F hRockingAnglevsTime  = TH1F("hRockingAnglevsTime","RockingAnglevsTime",TimeBins,StartTime,StopTime);
   CalculateRockingAngle(&hRockingAnglevsTime,&hPtRazvsTime,&hPtDeczvsTime, &hPtRaxvsTime, &hPtDecxvsTime, &hRAZenithvsTime, &hDecZenithvsTime);
   hMcIlwainLvsTime.Write();
   hDecZenithvsTime.Write();
   hRAZenithvsTime.Write();
   hPtDecxvsTime.Write();
   hPtDeczvsTime.Write();
   hPtRazvsTime.Write();
   hPtRaxvsTime.Write();
   hRockingAnglevsTime.Write();
 }

 fout->Close();
 return 0;
}


//function that interpolates pointing plots
//takes care of 360.0 and 0.360 crossings of RAs
void Interpolate(TH1F* hist,bool aRA,int verbosity) {
 //if (verbosity>1) printf("Interpolate %s\n",hist->GetTitle());
 
 double X[200],Y[200];
 float aval;
 int ilast_filled=1;//last filled i
 short int nel=0;
 bool Process=false;
 TGraph g;

 for (int i=1;i<=hist->GetNbinsX();i++) {
    if (verbosity>1 && i%1000==0) {printf("%.2f \r",i/float(hist->GetNbinsX())); fflush(0);}
    if ((ilast_filled && (i-ilast_filled)>33 && nel<=1))  {  //there is a gap, reset table
        //  printf("long gap nel=%d i,ilast=%d %d %d time=%f %f\n",nel,i,ilast_filled,i-ilast_filled,hist->GetBinCenter(i),hist->GetBinCenter(ilast_filled));
         ilast_filled=0;
         nel=0;
    }
    if  (Process ||                              //There was a reversal or buffer full
        (ilast_filled && (i-ilast_filled)>33) || //gap found..interpolate
         (nel>1 && i==hist->GetNbinsX()) )  {      //last bin
          #ifdef DEBUG
          printf("processing at i=%d (%lf) nel=%d fill up to <%d %f Proc=%d ilast_filled=%d gap=%d\n",
                 i,hist->GetBinCenter(i),nel,int(floor(X[nel-1])),hist->GetBinCenter(int(floor(X[nel-1]))),Process,ilast_filled,i-ilast_filled);
          #endif
          g = TGraph(nel,X,Y);
          for (unsigned int ii=(unsigned int)floor(X[0])+1;ii<(unsigned int)floor(X[nel-1]);ii++) {
              if (aRA) {
                float x=g.Eval(ii,0,"S");
                if (x>360) x-=360;
                else if (x<0) x+=360;
                hist->SetBinContent(ii,x);
              }
              else hist->SetBinContent(ii,g.Eval(ii,0,"S"));
          } 

           
          /* //we don't care about the first bins now since we are working on an extended range
          if (first) {//extrapolate to fill the initial gap too
             for (int ii=1;ii<X[0];ii++) {
                if (aRA) {
                   float x=g.Eval(ii,0,"");
                   if (x>360) x-=360;
                   else if (x<0) x+=360;
                   hist->SetBinContent(ii,x);
                }
                else hist->SetBinContent(ii,g.Eval(ii,0,""));
              } 
              first=false;
          }
          */

          ilast_filled=0;
          Process=false;
          i=int(X[nel-1]);  //step back
          nel=0;
    }

    aval=hist->GetBinContent(i);
    if (aval) { //search for next entry
       if (aRA) { 
          if (nel && aval<40 && Y[nel-1]>250) { //reversal 
             aval+=360;
             Process=true;
          }
          else if (nel && aval>250 && Y[nel-1]<40) { //reversal 
             aval-=360;
             Process=true;
          }
       }

       X[nel]=i;
       Y[nel]=aval;
       nel++;
       ilast_filled=i;
       if (nel==200) {Process=true;}    
    }
 }

 /* //we don't care for the last bins since we are working on an extended range
 //fill last bins
 if (nel<=1) { //don't have enough data for extrapolation..find some
    nel=0;
    int ilast=0;
    for (int i=hist->GetNbinsX();i>0;i--) {    //search for last filled bin in the histo
      if (hist->GetBinContent(i)) {ilast=i; break;}
    }
    if (ilast==0) {printf("%s: ilast is still zero? \n",__FUNCTION__); hist->Write(); throw std::runtime_error("");}
    for (int i=ilast-1;i<=ilast;i++) { //just get two bins and extrapolate
        if ((Y[nel]=hist->GetBinContent(i))!=0) {
           X[nel]=i;
           nel++;
        }
    }
 }

 g = TGraph(nel,X,Y);
 for (int ii=floor(X[0])+1;ii<=hist->GetNbinsX();ii++) {
      if (aRA) {
        float x=g.Eval(ii,0,"");
        if (x>360) x-=360;
        else if (x<0) x+=360;
        hist->SetBinContent(ii,x);
      }
      else hist->SetBinContent(ii,g.Eval(ii,0,""));
 } 
 */

}



//this function calculates the rocking angle from raz,decz,razenith,deczenith and also checks if the histograms are good by calculating if PtZ and PtX are perpendicular
//as they should. If their angle is more than 5 degs it complains..
void CalculateRockingAngle(TH1F* hRockingAnglevsTime, TH1F* hPtRazvsTime, TH1F* hPtDeczvsTime, TH1F* hPtRaxvsTime, TH1F* hPtDecxvsTime, TH1F* hRAZenithvsTime, TH1F* hDecZenithvsTime) {

 //ROCKING ANGLE
 astro::SkyDir SCz,SCx,SCZenith;
 CLHEP::Hep3Vector LocalDir;
 astro::SkyDir local;

 CLHEP::HepRotation localToCelestial;
 localToCelestial.setTolerance(0.1);
 float PTRAZ,PTRAX,PTDECZ,PTDECX,RA_ZENITH,DEC_ZENITH;

 for (int i=1;i<=hRockingAnglevsTime->GetNbinsX();i++) {
      PTRAZ =  hPtRazvsTime->GetBinContent(i);
      if (PTRAZ==0) continue;
      PTRAX =  hPtRaxvsTime->GetBinContent(i);
      PTDECZ = hPtDeczvsTime->GetBinContent(i);
      PTDECX = hPtDecxvsTime->GetBinContent(i);
      SCz = astro::SkyDir(PTRAZ,PTDECZ,astro::SkyDir::EQUATORIAL);
      SCx = astro::SkyDir(PTRAX,PTDECX,astro::SkyDir::EQUATORIAL);

      RA_ZENITH = hRAZenithvsTime->GetBinContent(i);
      DEC_ZENITH= hDecZenithvsTime->GetBinContent(i);
      localToCelestial = CLHEP::HepRotation(SCx.dir(),(SCz.dir()).cross(SCx.dir()),SCz.dir());
     
      SCZenith = astro::SkyDir(RA_ZENITH,DEC_ZENITH,astro::SkyDir::EQUATORIAL);
      double tol=SCx.difference(SCz)/DEG_TO_RAD-90;
      if (fabs(tol)>2) { printf("%s: Accuracy of plots low.. %d %lf %f\n",__FUNCTION__,i,hPtRazvsTime->GetBinCenter(i),tol); }
      //if (fabs(tol)>3) throw std::runtime_error("");
      hRockingAnglevsTime->SetBinContent(i,SCZenith.difference(SCz)/DEG_TO_RAD);
  }
}


//ps. I hate 30-sec intervals..

