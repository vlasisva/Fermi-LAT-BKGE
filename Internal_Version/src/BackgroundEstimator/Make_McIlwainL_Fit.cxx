//Author: Vlasios Vasileiou <vlasisva@gmail.com>
// $Header: /nfs/slac/g/glast/ground/cvs/GRBAnalysis-scons/BackgroundEstimator/src/BackgroundEstimator/Make_McIlwainL_Fit.cxx,v 1.6 2013/11/28 13:20:18 vlasisva Exp $
#include "BackgroundEstimator/BackgroundEstimator_ext.h"

void BackgroundEstimator_ext::Make_McIlwainL_Fits(string FitsAllSkyFile){

  sprintf(name,"%s/Plots_%s.root",DataDir.c_str(),DataClass.c_str());
  TFile * fPlots = TFile::Open(name);
  TH1F* hPtRazvsTime     = (TH1F*)fPlots->Get("hPtRazvsTime;1");
  TH1F* hPtDeczvsTime    = (TH1F*)fPlots->Get("hPtDeczvsTime;1");
  TH1F* hMcIlwainLvsTime = (TH1F*)fPlots->Get("hMcIlwainLvsTime;1");
  TH1F* hRockingAnglevsTime = (TH1F*)fPlots->Get("hRockingAnglevsTime;1");
  if (!hPtRazvsTime || !hPtDeczvsTime || !hMcIlwainLvsTime) {printf("%s: Can't read plots from file %s\n",__FUNCTION__,name); throw std::runtime_error("");}

  
  TH1F * hTheta_away[Energy_Bins_user+1];
  sprintf(name,"%s/ThetaPhi_Fits_%s_%.1f.root",DataDir.c_str(),DataClass.c_str(),ThetaPhiFits_version);
  TFile * fThetaPhi_Fits = TFile::Open(name);
  for (int i=1;i<=Energy_Bins_user;i++) {
    sprintf(name,"hTheta_away_%d",i);
    hTheta_away[i]=(TH1F*)fThetaPhi_Fits->Get(name);
    if (!hTheta_away[i]) {printf("%s: no htheta_away_%d\n",__FUNCTION__,i); throw std::runtime_error("");}
  }
  
  //To avoid contamination from galactic-plane gamma-rays, only events when the LAT is pointing at a |GalacticLatitude|>B_Cut are considered.
  const float B_Cut=70;
  const float McIlwainLMin=0.983,McIlwainLMax=1.735;

  //the final rates are produced by upscaling the initially-calculated event rates (corresponding to FT1Theta<MaxTheta)
  //The upscaling factor is the inverse of the fraction of events having a FT1Theta<MaxTheta when looking away from the earth.
  //This factor comes from the ThetaPhi_Fits produced earlier.

  //I decreased maxtheta so that MaxTheta+RockAngle<=earth_limb_ztheta
  float MaxTheta=48;
  
  //We can use MaxTheta to also try to avoid the gamma rays from the galactic diffuse
  //A part of the theta>MaxTheta data (for the Phis that point towards the plane) have increased gamma ray contribution
  //By reducing MaxTheta we can reduce the gamma contamination of the McIlwainL fit

  //round up MaxTheta to the upper edge of the correspond bin in hTheta_away. 
  //This will avoid any accuracy problems in calculating the correction factor below.
  const int iMaxTheta = hTheta_away[1]->FindBin(MaxTheta);
  MaxTheta = hTheta_away[1]->GetXaxis()->GetBinUpEdge(iMaxTheta);
  const float MaxRockingAngle=FT1ZenithTheta_Cut-MaxTheta; //for ztheta_cut=100 and maxtheta=48, 3% of the observation time is rejected
  printf("%s: MaxTheta=%.2fdeg MaxRockingAngle=%f deg\n",__FUNCTION__,MaxTheta);

  unsigned short int McIlwainLBins[Energy_Bins_user+1];
  for (int i=0;i<=Energy_Bins_user;i++) McIlwainLBins[i]=10;
  //note if you make McIlwainLBins depend in i then it messes up hRateAll but this is OK

  TH1F *hEvents[Energy_Bins_user+1],*hEvents_Cut[Energy_Bins_user+1] ,*hRate[Energy_Bins_user+1],*hTime[Energy_Bins_user+1];
  
  TH1F*  hdt = new TH1F("hdt","Spectrum of time intervals between successive events",100,-6,5);
  hdt ->GetXaxis()->SetTitle("log_{10}(dt/sec)");
  hdt ->GetYaxis()->SetTitle("");

  for (int i=1;i<=Energy_Bins_user;i++) {
     sprintf(name,"hEvents_%d",i);
     hEvents[i] = new TH1F(name,"Number of events detected while at each McIlwainL bin",McIlwainLBins[i],McIlwainLMin,McIlwainLMax);
     hEvents[i]->GetXaxis()->SetTitle("McIlwainL");
     hEvents[i]->GetYaxis()->SetTitle("Events (events/bin)");

     sprintf(name,"hEvents_Cut_%d",i);
     hEvents_Cut[i] = new TH1F(name,"Number of events detected while at each McIlwainL bin with theta<ThetaCut",McIlwainLBins[i],McIlwainLMin,McIlwainLMax);
     hEvents_Cut[i]->GetXaxis()->SetTitle("McIlwainL");
     hEvents_Cut[i]->GetYaxis()->SetTitle("Events (events/bin)");

     sprintf(name,"hRate_%d",i);
     hRate[i] = new TH1F(name,"All sky event rate",McIlwainLBins[i],McIlwainLMin,McIlwainLMax);
     hRate[i]->GetXaxis()->SetTitle("McIlwainL");
     hRate[i]->GetYaxis()->SetTitle("Rate (Hz)");

     sprintf(name,"hTime_%d",i);
     hTime[i] = new TH1F(name,"Time spent at each McIlwainL bin",McIlwainLBins[i],McIlwainLMin,McIlwainLMax);
     hTime[i]->GetXaxis()->SetTitle("McIlwainL");
     hTime[i]->GetYaxis()->SetTitle("Time (sec/bin)");
  }

  TH1F * hRate_all = new TH1F("hRate_all","All sky event rate",McIlwainLBins[1],McIlwainLMin,McIlwainLMax);
  hRate_all->GetXaxis()->SetTitle("McIlwainL");
  hRate_all->GetYaxis()->SetTitle("Rate (Hz)");

  TH2F * Map = new TH2F("Map","Map",360,0,360,180,-90,90);

  double tlast=StartTime,dt;
  double PtL,PtB;
  short int iebin;
  int itimebin;

  fitsfile *fptr;

  FILE* ftemp = fopen(FitsAllSkyFile.c_str(),"r");
  int ifile=0,ibin;
  while (fscanf(ftemp,"%s",name)==1) {
     printf("%d\r",ifile); fflush(0); ifile++;
     int status=0,hdutype,anynul;
     fits_open_file(&fptr, name, READONLY, &status);
     fits_movabs_hdu(fptr, 2, &hdutype, &status);
     long nrows=0;int ncols=0;
     fits_get_num_rows(fptr, &nrows, &status);
     fits_get_num_cols(fptr, &ncols, &status);
     int format;
     if      (DataClass.find("P7")!=string::npos) format=DATA_FORMAT_P7;
     else if (DataClass.find("P6")!=string::npos){
        if      (ncols==22) format=DATA_FORMAT_P6_NEW;
        else if (ncols==18) format=DATA_FORMAT_P6_OLD; 
        else {printf("%s: unknown format file %s ncols=%d class=%s\n",__FUNCTION__,name,ncols,DataClass.c_str()); throw std::runtime_error("");}
     }
     else {printf("%s: Unknown fits file format file %s class=%s\n",__FUNCTION__,name,DataClass.c_str()); throw std::runtime_error("");}
 
     if (nrows<=0) {printf("%s: nrows=%ld\n", __FUNCTION__,nrows); throw std::runtime_error("");}
     for (long int i=1;i<=nrows;i++) {
         
         double PtTime;
         fits_read_col (fptr,TDOUBLE,10,i, 1, 1, NULL,&PtTime, &anynul, &status);
         dt=PtTime-tlast;
         tlast=PtTime;
         
         itimebin = hPtRazvsTime->FindBin(PtTime);
         if  (hPtRazvsTime->GetBinContent(itimebin)==0) {
            if      (hPtRazvsTime->GetBinContent(itimebin-1)!=0) ibin=itimebin-1;
            else if (hPtRazvsTime->GetBinContent(itimebin+1)!=0) ibin=itimebin+1;
            else    {printf("%s: there is a gap in the plots? %f %d\n",__FUNCTION__,PtTime,itimebin); continue;}
            itimebin=ibin;
         }
         if (hRockingAnglevsTime->GetBinContent(itimebin)>MaxRockingAngle) continue; //reject event if detected during a high rocking angle period
         
         TOOLS::Galactic((double)hPtRazvsTime->GetBinContent(itimebin),(double)hPtDeczvsTime->GetBinContent(itimebin),&PtL,&PtB);
         //exclude observations near the galactic plane
         if (fabs(PtB)<B_Cut || !PassesCuts(fptr,i,format)) continue;
         
         Map->Fill(PtL,PtB);
         hdt->Fill(log10(dt));

         if (dt<10) {
              float PtMcIlwainL=hMcIlwainLvsTime->GetBinContent(itimebin);
              for (int iee=1;iee<=Energy_Bins_user;iee++) hTime[iee]->Fill(PtMcIlwainL,dt);

              double FT1Energy;
              fits_read_col (fptr,TDOUBLE,1,i, 1, 1, NULL,&FT1Energy, &anynul, &status);
              iebin=Energy2Bin(FT1Energy);

              hEvents[iebin]->Fill(PtMcIlwainL);
              double FT1Theta;
              fits_read_col (fptr,TDOUBLE,6,i, 1, 1, NULL,&FT1Theta, &anynul, &status);
              if (FT1Theta<MaxTheta) hEvents_Cut[iebin]->Fill(PtMcIlwainL);
         }
     }
     fits_close_file(fptr, &status);
  }    
  fclose(ftemp);
     
  //CALCULATE RATES
  TCanvas * c1[Energy_Bins_user+1];
  sprintf(name,"%s/RateFit_%s_%.1f.root",DataDir.c_str(),DataClass.c_str(),RateFit_version);
  TFile * fout = new TFile(name,"RECREATE");
  for (int ie=1;ie<=Energy_Bins_user;ie++) {
     for (int i=1;i<=McIlwainLBins[ie];i++) {
         if (hTime[ie]->GetBinContent(i)) {
            hRate[ie]->SetBinContent(i,hEvents[ie]->GetBinContent(i)/hTime[ie]->GetBinContent(i));
            hRate[ie]->SetBinError(i,sqrt(hEvents[ie]->GetBinContent(i))/hTime[ie]->GetBinContent(i)); 
         }
     }
     sprintf(name,"hRate_alltheta_%d",ie);
     hRate[ie]->Write(name);

     for (int i=1;i<=McIlwainLBins[ie];i++) {
         if (hTime[ie]->GetBinContent(i)) {
            hRate[ie]->SetBinContent(i,hEvents_Cut[ie]->GetBinContent(i)/hTime[ie]->GetBinContent(i));
            hRate[ie]->SetBinError(i,sqrt(hEvents_Cut[ie]->GetBinContent(i))/hTime[ie]->GetBinContent(i)); //used for fitting
         }
     }
     sprintf(name,"hRate_withcut_uncorrected_%d",ie);
     hRate[ie]->Write(name);

     double FractionUnderThetaCut = hTheta_away[ie]->Integral(0,iMaxTheta)/hTheta_away[ie]->Integral();
     //printf("%d %f %f %f\n",ie,hTheta_away[ie]->Integral(0,MaxTheta),hTheta_away[ie]->Integral(),FractionUnderThetaCut);
     hRate[ie]->Scale(1./FractionUnderThetaCut);
     sprintf(name,"hRate_withcut_corrected_%d",ie);
     hRate[ie]->Write(name);


    //DRAW
    sprintf(name,"c_%d",ie);
    c1[ie] = new TCanvas(name,"R(McIlwainL) calculation");
    c1[ie]->Divide(2,2);
    c1[ie]->cd(1); hEvents[ie]->Draw();
    c1[ie]->cd(2); hTime[ie]->Draw();
    c1[ie]->cd(3); hdt->Draw();
    c1[ie]->cd(4);

    sprintf(name,"RateFit_%d",ie);
    TF1 * aRateFit;
    if (DataClass.find("DIFFUSE")!=string::npos) aRateFit = new TF1(name,"pol4");
    else aRateFit = new TF1(name,"pol5");
    hRate[ie]->Fit(name,"Q");

    c1[ie]->GetPad(3)->SetLogy();
    c1[ie]->Update();
    c1[ie]->Write();
    aRateFit->Write();
    aRateFit->Delete();
  }

  
  for (int ie=1;ie<=Energy_Bins_user;ie++) hRate_all->Add(hRate[ie]);
  hRate_all->Write();
  printf("       \r");


  //THE ALL-SKY RATE CALCULATED ABOVE IS FOR AN UNOBSTRUCTED BY A ZENITH-THETA CUT FOV
  //NOW THAT WE APPLY A ZENITH-THETA CUT (AND WE REJECT A PART OF THE SKY), WE HAVE TO CALCULATE THE FRACTIONAL 
  //DECREASE IN THE ALL-SKY RATE. WE CALCULATE THIS NUMBER FOR DIFFERENT ENERGY_BINS AND ROCKING ANGLES
  TH2F * hScaleFactor = new TH2F("hScaleFactor","hScaleFactor",180,0,180,Energy_Bins_user,1,Energy_Bins_user+1);
  const float FT1ZenithTheta_Cut_RAD=FT1ZenithTheta_Cut*DEG_TO_RAD;
  for (int i=1;i<=hScaleFactor->GetNbinsX();i++) {
     printf("%d\r",i); fflush(0);
     float rock_RAD=hScaleFactor->GetXaxis()->GetBinCenter(i)*DEG_TO_RAD;
     for (int j=1;j<=hScaleFactor->GetNbinsY();j++) {
       double eff=0,eff_sum=0;
       for (float theta_RAD=0;theta_RAD<80*DEG_TO_RAD;theta_RAD+=0.1*DEG_TO_RAD) { //theta integration
          //printf("%d %d %f\n",i,j,theta_RAD/DEG_TO_RAD);
          int thetabin=hTheta_away[j]->FindBin(theta_RAD/DEG_TO_RAD);
          float arate=hTheta_away[j]->GetBinContent(thetabin);
          for (float phi_RAD=0;phi_RAD<360*DEG_TO_RAD;phi_RAD+=1*DEG_TO_RAD) {
             double ZA=acos( sin(rock_RAD)*sin(theta_RAD)*cos(phi_RAD) + cos(rock_RAD)*cos(theta_RAD));
             if (ZA<=FT1ZenithTheta_Cut_RAD) eff+=arate;
             eff_sum+=arate;
          }//phi loop
       }//theta loop
       eff/=eff_sum; 
       hScaleFactor->SetBinContent(i,j,eff);
     } //energy loop
  }//rocking angle loop
  hScaleFactor->Write();  
  printf("       \r");


  sprintf(name,"%f",B_Cut);
  TNamed Data = TNamed("B_Cut",name);
  Data.Write();

  sprintf(name,"%f",FT1ZenithTheta_Cut);
  Data = TNamed("FT1ZenithThetaCut",name);
  Data.Write();

  sprintf(name,"%.2f",RateFit_version);
  Data = TNamed("version",name);
  Data.Write();

  Map->Write();

  fout->Close();
  fPlots->Close();
  fThetaPhi_Fits->Close();

}
