//Author: Vlasios Vasileiou <vlasisva@gmail.com>
//$Header: /nfs/slac/g/glast/ground/cvs/GRBAnalysis-scons/BackgroundEstimator/src/BackgroundEstimator/BackgroundEstimator.cxx,v 1.9 2013/11/09 08:40:12 vlasisva Exp $
#include "BackgroundEstimator/BackgroundEstimator.h"


#define DEBUG

ClassImp(BackgroundEstimator)

BackgroundEstimator::~BackgroundEstimator() 
{//do some cleanup
 if(fResidualOverExposure) fResidualOverExposure->Close();
 if(fRateFit)              fRateFit->Close();
 if(fThetaPhiFits)         fThetaPhiFits->Close();
 if(fCorrectionFactors)    fCorrectionFactors->Close();
}

BackgroundEstimator::BackgroundEstimator(string aClass, double EMin, double EMax, int EBins, bool initialize, bool ShowLogo):
Energy_Min_datafiles(0),Energy_Max_datafiles(0),Energy_Bins_datafiles(0),
Energy_Min_user(EMin),Energy_Max_user(EMax),Energy_Bins_user(EBins),UsingDefaultBinning(true),
DataClass(aClass),FT1ZenithTheta_Cut(100),EstimatorVersion(4),Residuals_version(2.0),RateFit_version(2.0),ThetaPhiFits_version(2.0),EastWest_version(2.0),TimeCorrectionFactors_version(3.0),
StartTime(0),EndTime(0),StopTime(0),TimeBins(0),BinSize(0.5),fResidualOverExposure(0),fRateFit(0),fThetaPhiFits(0),fCorrectionFactors(0)
{
 
  if (ShowLogo) {
   printf("*----------------------------------------------*\n");
   printf("|           Background Estimator (public)      |\n");
   printf("|              v%.0f September/2013               |\n",EstimatorVersion);
   printf("|                                              |\n");
   printf("| contact:     Vlasios Vasileiou               |\n");
   printf("|              vlasisva@gmail.com              |\n");
   printf("| http://arxiv.org/abs/1307.4284               |\n");
   printf("*----------------------------------------------*\n");
 }

  L_BINS = 720;
  B_BINS = 360;
  vector <string> VALID_CLASSES;
  VALID_CLASSES.push_back("P7TRANSIENT_V15");

  bool goodClass=false;
  for (unsigned int i=0;i<VALID_CLASSES.size();i++) {
     if (DataClass==VALID_CLASSES[i]) {goodClass=true; break;}
  }

  if (!goodClass) {
      printf("%s: Only \n",__FUNCTION__);
      for (unsigned int i=0;i<VALID_CLASSES.size();i++) printf(" %s",VALID_CLASSES[i].c_str());
      printf(" classes are supported\n");
      return;
  }

  EventClassMask=0;
  if (DataClass.find("P7")!=string::npos) {
     if (DataClass.find("TRANSIENT")!=string::npos) EventClassMask=int(pow(2.,0.));
  }

  if (ShowLogo) printf("%s: Using %s data class.\n",__FUNCTION__,DataClass.c_str());
  DataClassName_noConv=TOOLS::GetDataClassName_noConv(DataClass);
  ConversionName      =TOOLS::GetConversionName(DataClass);
  ConversionType      =TOOLS::GetConversionType(DataClass);
  DataClassVersion    =TOOLS::GetDataClassVersion(DataClass);
  MinCTBClassLevel    =TOOLS::GetCTBClassLevel(DataClass);
  DataDir             =TOOLS::GetS("BASEDIR")+"/gtgrb_data/Bkg_Estimator/";

  if (initialize) { //normal operation
      sprintf(name,"%sResidual_Over_Exposure_%s_%.1f.root",DataDir.c_str(),DataClass.c_str(),Residuals_version);
      fResidualOverExposure = TFile::Open(name);
      if (!fResidualOverExposure) {printf("%s: Data file %s cannot be read. Did you get the data files?\n",__FUNCTION__,name); exit(1);}

      if (sscanf(fResidualOverExposure->Get("Energy_Data")->GetTitle(),"%lf_%lf_%d",&Energy_Min_datafiles,&Energy_Max_datafiles,&Energy_Bins_datafiles)!=3)
          sscanf(fResidualOverExposure->Get("Energy_Data")->GetTitle(),"%lf-%lf-%d",&Energy_Min_datafiles,&Energy_Max_datafiles,&Energy_Bins_datafiles);
      

      if (Energy_Min_user<=0)  Energy_Min_user=Energy_Min_datafiles;
      else if (fabs(Energy_Min_user-Energy_Min_datafiles)>0.1) {
          if (ShowLogo) printf("%s: Energy_Min set to %f\n",__FUNCTION__,Energy_Min_user);
          UsingDefaultBinning=false;
      }
      if (Energy_Max_user<=0)  Energy_Max_user=Energy_Max_datafiles;
      else if (fabs(Energy_Max_user-Energy_Max_datafiles)>0.1) {
           if (ShowLogo) printf("%s: Energy_Max set to %f\n",__FUNCTION__,Energy_Max_user);
           UsingDefaultBinning=false;
      }
      if (Energy_Bins_user<=0) Energy_Bins_user=Energy_Bins_datafiles;
      else if (Energy_Bins_user!=Energy_Bins_datafiles) {
            if (ShowLogo) printf("%s: Energy_Bins set to %d\n",__FUNCTION__,Energy_Bins_user);
            UsingDefaultBinning=false;
      }
      
      if (ShowLogo) {
                                  printf("Data Files Energy: (%.1f-%.1f)MeV - %d bins\n",Energy_Min_datafiles,Energy_Max_datafiles,Energy_Bins_datafiles);
         if (UsingDefaultBinning) printf("User is using the default (shown above) energy configuration\n");
   	 else                     printf("User's Energy:     (%.1f-%.1f)MeV - %d bins\n",Energy_Min_user,Energy_Max_user,Energy_Bins_user);
      }
      
     
      sscanf(fResidualOverExposure->Get("FT1ZenithTheta_Cut")->GetTitle(),"%f",&FT1ZenithTheta_Cut);
      

      float aversion=atof(fResidualOverExposure->Get("version")->GetTitle());
      if (aversion<Residuals_version) {
           printf("%s: Warning! You are using a Residuals data file that has an older version than the latest (v%.2f) one.\n",__FUNCTION__,Residuals_version);
           Residuals_version = aversion;
      }

      sprintf(name,"%sRateFit_%s_%.1f.root",DataDir.c_str(),DataClass.c_str(),RateFit_version);
      fRateFit = TFile::Open(name);
      if (!fRateFit) {printf("%s: Data file %s cannot be read. Did you get the data files?\n",__FUNCTION__,name); exit(1);}
      
      aversion=atof(fRateFit->Get("version")->GetTitle());
      if (aversion<RateFit_version) {
           printf("%s: Warning! You are using a RateFit data file that has an older version than the latest (v%.1f) one.\n",__FUNCTION__,RateFit_version);
           RateFit_version = aversion;
      }     
      fRateFit->Close();
      
      sprintf(name,"%s/ThetaPhi_Fits_%s_%.1f.root",DataDir.c_str(),DataClass.c_str(),ThetaPhiFits_version);
      fThetaPhiFits = TFile::Open(name);
      if (!fThetaPhiFits) {printf("%s: Data file %s cannot be read. Did you get the data files?\n",__FUNCTION__,name); exit(1);}
      aversion=atof(fThetaPhiFits->Get("version")->GetTitle());

      if (aversion<RateFit_version) {
           printf("%s: Warning! You are using a ThetaPhiFits data file that has an older version than the latest (v%.1f) one.\n",__FUNCTION__,RateFit_version);
           ThetaPhiFits_version = aversion;
      }
      if (aversion<2.0) {
           printf("%s: You cannot use the bkge with such old data files. Please download the new ones.. \n",__FUNCTION__);
           exit(1);
      }
      
      fThetaPhiFits->Close();

      //Read Correction factors
      sprintf(name,"%s/TimeCorrectionFactor_%s_%.1f.root",DataDir.c_str(),DataClass.c_str(),TimeCorrectionFactors_version);
      fCorrectionFactors= TFile::Open(name);
      for (int iE=0;iE<Energy_Bins_datafiles;iE++) pRatiovsTime.push_back(0);
      
      if (!fCorrectionFactors) printf("%s: Can't open file %s\n",__FUNCTION__,name); 
      else {
         TCanvas * cCorr = (TCanvas*)fCorrectionFactors->Get("cRatiovsTime");
         if (cCorr) {
            for (int iE=0;iE<Energy_Bins_datafiles;iE++) {
               sprintf(name,"RatiovsTime_%d",iE);
               pRatiovsTime[iE]=(TProfile*)(cCorr->GetPad(iE+1)->FindObject(name));
               if (pRatiovsTime[iE]==NULL){ printf("%s: error reading TimeCorrectionFactors for %d\n",__FUNCTION__,iE); exit(1);}
            }      
         }
     } 
        
  }    
     
}


//#define DEBUG

void BackgroundEstimator::CreateDataFiles(string FitsAllSkyFilesList, string FT2_FILE, double StartTime_user, double EndTime_user, float ZenithThetaCut){

  FT1ZenithTheta_Cut=ZenithThetaCut;
  TimeStep=1.0; //don't change that - left for debugging purposes
  Energy_Min_datafiles=Energy_Min_user;
  Energy_Max_datafiles=Energy_Max_user;
  Energy_Bins_datafiles=Energy_Bins_user;

  printf("%s: ZenithTheta Cut=%.1f deg, BinSize=%.1f deg\n",__FUNCTION__,FT1ZenithTheta_Cut,BinSize);
  printf("%s: Energy Range = (%.1e, %.1e)MeV, Energy_Bins=%d\n",__FUNCTION__,Energy_Min_datafiles,Energy_Max_datafiles,Energy_Bins_datafiles);

  //Read starting and ending time of the data in the FT1 Fits Files
  fitsfile *fptr;
  int status = 0;	
  long nrows;
  int hdutype,anynul;
  FILE * ftemp;

  //1. PREPARE FILES  
  FILE * fOldFitsList = fopen(FitsAllSkyFilesList.c_str(),"r");
  if (!fOldFitsList) {printf("%s: file %s cannot be read.\n",__FUNCTION__,FitsAllSkyFilesList.c_str()); return;}

  //2.READ TIMES FROM GENERATED FILES    
  ftemp = fopen(FitsAllSkyFilesList.c_str(),"r");
  int res=fscanf(ftemp,"%s",name);
  status=0;
  fits_open_file(&fptr, name, READONLY, &status);
  if (status) fits_report_error(stderr, status);  status=0;
  fits_movabs_hdu(fptr, 2, &hdutype, &status);
  if (status) fits_report_error(stderr, status);  status=0;
  fits_read_col (fptr,TDOUBLE,10,1, 1, 1, NULL,&StartTime, &anynul, &status);
  if (status) fits_report_error(stderr, status);  status=0;
  fits_get_num_rows(fptr, &nrows, &status);
  if (status) fits_report_error(stderr, status);  status=0;
  fits_read_col (fptr,TDOUBLE,10,nrows, 1, 1, NULL,&EndTime, &anynul, &status);
  fits_close_file(fptr, &status);

  while (fscanf(ftemp,"%s",name)==1) {
     status=0;
     fits_open_file(&fptr, name, READONLY, &status);
     if (status) fits_report_error(stderr, status);
     fits_movabs_hdu(fptr, 2, &hdutype, &status);
     if (status) fits_report_error(stderr, status);  status=0;
     fits_get_num_rows(fptr, &nrows, &status);
     if (status) fits_report_error(stderr, status);  status=0;
     double anEndTime,aStartTime;
     fits_read_col (fptr,TDOUBLE,10,1, 1, 1, NULL,&aStartTime, &anynul, &status);
     if (status) fits_report_error(stderr, status);  status=0;
     fits_read_col (fptr,TDOUBLE,10,nrows, 1, 1, NULL,&anEndTime, &anynul, &status);
     if (status) fits_report_error(stderr, status);  status=0;
     fits_close_file(fptr, &status);
     if (status) fits_report_error(stderr, status);  status=0;
     printf("%s: %f (%f) to %f\n",name,aStartTime,aStartTime-EndTime,anEndTime);
     if (anEndTime>EndTime) EndTime=anEndTime;
  }
  
  fclose (ftemp);


  if (StartTime==EndTime) {printf("%s:StartTime==EndTime?? \n",__FUNCTION__); exit(1);} //this also includes when they are both zero
  if (StartTime_user!=0) {StartTime=StartTime_user; printf("%s: Override StartTime to %f\n",__FUNCTION__,StartTime);}
  if (EndTime_user!=0)   {EndTime=EndTime_user; printf("%s: Override EndTime to %f\n",__FUNCTION__,EndTime);}
  

  TimeBins = int((EndTime-StartTime)/TimeStep);
  StopTime = StartTime + TimeStep*TimeBins;
  printf("%s: Using MET:(%f-%f), dt=%e sec from all-sky file\n",__FUNCTION__,StartTime,StopTime,StopTime-StartTime);


  //3. MAKE PLOTS
  sprintf(name,"%s/Plots_%s.root",DataDir.c_str(),DataClass.c_str());
  ftemp = fopen(name,"r");
  if (!ftemp) {
       printf("%s: Making plots...\n",__FUNCTION__);
       double midtime=(StopTime+StartTime)/2;
       TOOLS::Make_Plots(midtime-StartTime+10,StopTime-midtime+10,midtime,name,FT2_FILE,99,10);
  }
  else  fclose(ftemp);
  
  //4. MAKE THETA/PHI DISTRIBUTIONS
  sprintf(name,"%s/ThetaPhi_Fits_%s_%.1f.root",DataDir.c_str(),DataClass.c_str(),ThetaPhiFits_version);
  ftemp = fopen(name,"r"); 
  if (!ftemp) { 
    printf("%s: Making Theta & Phi fits...\n",__FUNCTION__); 
    Make_ThetaPhi_Fits(FitsAllSkyFilesList); 
  } 
  else fclose(ftemp);

  //5. MAKE RATE FITS
  sprintf(name,"%s/RateFit_%s_%.1f.root",DataDir.c_str(),DataClass.c_str(),RateFit_version);
  ftemp = fopen(name,"r"); 
  if (!ftemp) { 
    printf("%s: Making Rate fit...\n",__FUNCTION__);
    Make_McIlwainL_Fits(FitsAllSkyFilesList); 
  } 
  else fclose(ftemp); 

      //3. RESIDUAL_OVER_EXPOSURE
       sprintf(name,"%s/Residual_Over_Exposure_%s_%.1f.root",DataDir.c_str(),DataClass.c_str(),Residuals_version);
       ftemp = fopen(name,"r"); 
       if (!ftemp) {
          //////////////////////////////////////////
          //5.1.Exposure

/*
             sprintf(name,"%s/ltCube_%s.fits",DataDir.c_str(),DataClass.c_str());
             ftemp=fopen(name,"r"); 
             if (!ftemp) {
                sprintf(name,"ltCube_%s.fits",DataDir.c_str(),DataClass.c_str());
                TOOLS::Run_gtltcube(Datadir, StartTime, EndTime, FT2_FILE, FT1ZenithTheta_Cut, 2, EventFile, name);
              } 
              else fclose(ftemp);
*/     
     
           //5.1. Calculate exposure maps
           sprintf(name,"%s/exposure_%s.fits",DataDir.c_str(),DataClass.c_str());
           ftemp = fopen(name,"r"); 
           if (!ftemp) { 
                char gtltcube_Filename[100];
                sprintf(gtltcube_Filename,"ltCube_%s.fits",DataClass.c_str());
                printf("%s: Creating exposure map %s...\n",__FUNCTION__,name); 
                TOOLS::Run_gtexpcube(DataDir, StartTime, EndTime, FT2_FILE, DataClass, FT1ZenithTheta_Cut, name, Energy_Min_datafiles,Energy_Max_datafiles, Energy_Bins_datafiles, 5, gtltcube_Filename);
                //sprintf(buffer,"gtexpcube infile=%s/ltCube_%s.fits evfile=@%s cmfile=NONE outfile=%s irfs=%s nxpix=1 nypix=1 pixscale=1 coordsys=GAL xref=0 yref=0 axisrot=0 proj=CAR emin=%f emax=%f enumbins=%d bincalc=CENTER clobber=yes 2>&1",DataDir.c_str(),DataClass.c_str(),FitsAllSkyFilesList.c_str(),name,DataClass.c_str(),Energy_Min_datafiles,Energy_Max_datafiles,Energy_Bins_datafiles);
          } 
          else fclose(ftemp); 
     
     
           //////////////////////////////////////////
           //5.2.Residual
           
           sprintf(name,"%s/Residual_%s.root",DataDir.c_str(),DataClass.c_str());
           ftemp = fopen(name,"r");
           if (!ftemp) { 
              printf("%s: Calculating Residuals...\n",__FUNCTION__); fflush(stdout); 
              CalcResiduals(FitsAllSkyFilesList);
           } 
           else fclose(ftemp); 
          
          ///////////////////////////////////////////
          //5.4. Calculate Residual Over Exposure
          TH2F * hExposureAllSky = new TH2F("hExposureAllSky","hExposureAllSky",L_BINS,-180,180,B_BINS,-90,90);
          sprintf(name,"%s/Residual_%s.root",DataDir.c_str(),DataClass.c_str()); 
          TFile * fResidual = TFile::Open(name);
          sprintf(name,"%s/Residual_Over_Exposure_%s_%.1f.root",DataDir.c_str(),DataClass.c_str(),Residuals_version); 
          TFile * fExposure = new TFile(name,"RECREATE");
     
          TH2D*  hSolidAngle = new TH2D("hSolidAngle","hSolidAngle",L_BINS,-180,180,B_BINS,-90,90); 
          double thetaphi=BinSize*DEG_TO_RAD; 
          double theta1,theta2,binarea; 
          for (int j=1;j<=B_BINS;j++) { 
             for (int i=1;i<=L_BINS;i++) { 
                theta1=(90+hSolidAngle->GetYaxis()->GetBinLowEdge(j))*DEG_TO_RAD; 
                theta2=(90+hSolidAngle->GetYaxis()->GetBinUpEdge(j))*DEG_TO_RAD; 
                binarea = (cos(theta1)-cos(theta2))*thetaphi; 
                hSolidAngle->SetBinContent(i,j,binarea); 
             } 
          }   
          hSolidAngle->Write(); 
      
          for (short unsigned int ie=1;ie<=Energy_Bins_datafiles;ie++) {
             sprintf(name,"hResidual_%d;1",ie); 
             TH2F * haResidual = (TH2F*)fResidual->Get(name);
             sprintf(name,"hResidual_Over_Exposure_%d",ie);
             TH2F *haResidualOverExposure = (TH2F*)haResidual->Clone(name);
             haResidualOverExposure->SetTitle("hResidual_Over_Exposure");
             haResidual->Delete();
     
             sprintf(name,"%s/exposure_%s.fits",DataDir.c_str(),DataClass.c_str());
             TOOLS::ReadExposureMap(name,hExposureAllSky,ie,1);
     
             sprintf(name,"hExposureAllSky_%d",ie);
             //hExposureAllSky->Write(name);
             haResidualOverExposure->Divide(hSolidAngle);
             haResidualOverExposure->Divide(hExposureAllSky);
     
             sprintf(name,"hResidual_Over_Exposure_%d",ie);
             haResidualOverExposure->Write(name);
             haResidualOverExposure->Delete();
     
             TOOLS::ProgressBar(ie-1,Energy_Bins_datafiles);
          }
          hExposureAllSky->Delete();
          fResidual->Close();

          sprintf(name,"%e_%e_%d",Energy_Min_datafiles,Energy_Max_datafiles,Energy_Bins_datafiles);
          TNamed Data = TNamed("Energy_Data",name);
          fExposure->cd(); Data.Write(); 

          sprintf(name,"%f",FT1ZenithTheta_Cut);
          Data = TNamed("FT1ZenithTheta_Cut",name);
          Data.Write(); 

          sprintf(name,"%.2f",Residuals_version);
          Data = TNamed("version",name);
          Data.Write();

          fExposure->Close();
       }
       else  fclose(ftemp);


};


unsigned short int BackgroundEstimator::Energy2Bin(float Energy){
  static const double dlEnergy = (log10(Energy_Max_datafiles)-log10(Energy_Min_datafiles))/Energy_Bins_datafiles;
  static const double lEnergy_Min = log10(Energy_Min_datafiles);

  if      (Energy<Energy_Min_datafiles) return 0;
  else if (Energy>=Energy_Max_datafiles) return Energy_Bins_datafiles+1;
  else return 1+(int)(floor((log10(Energy)-lEnergy_Min)/dlEnergy));
};


float BackgroundEstimator::Bin2Energy(unsigned short int bin){
  static const double dlEnergy = (log10(Energy_Max_datafiles)-log10(Energy_Min_datafiles))/Energy_Bins_datafiles;
  static const double lEnergy_Min = log10(Energy_Min_datafiles);
  if (bin==0) return Energy_Min_datafiles;
  else if (bin==Energy_Bins_datafiles+1) return Energy_Max_datafiles;
  else return pow(10,lEnergy_Min+(bin+0.5)*dlEnergy);
};


bool BackgroundEstimator::PassesCuts(fitsfile * fptr, long int i, int format) {
  int status=0,anynul;
  if (format==DATA_FORMAT_P7)  {
      static int EventClass=0;
      fits_read_col (fptr,TINT,15,i, 1, 1, NULL,&EventClass, &anynul, &status);
      //printf ("%d %d %d\n",EventClass,EventClassMask,EventClass%16);
      if (EventClass%16<EventClassMask) return false;
  }
  else {
    static int CTBClassLevel=0;
    if      (format==DATA_FORMAT_P6_OLD) fits_read_col (fptr,TINT,18,i, 1, 1, NULL,&CTBClassLevel, &anynul, &status);
    else if (format==DATA_FORMAT_P6_NEW) fits_read_col (fptr,TINT,15,i, 1, 1, NULL,&CTBClassLevel, &anynul, &status);
    else {printf("%s: what is going on?\n",__FUNCTION__); exit(1);}
    if   (CTBClassLevel<MinCTBClassLevel) return false; //apply class cut
  }



  if (ConversionType!=-1) { //if !BOTH
      static int aConversionType=0;
      fits_read_col (fptr,TINT,16,i, 1, 1, NULL,&aConversionType, &anynul, &status);
      //printf("conv %d %d\n",ConversionType,aConversionType);
      if (aConversionType!=ConversionType) return false; //apply Conversion Type front/back
  }

  static double FT1Energy=0;
  fits_read_col (fptr,TDOUBLE,1,i, 1, 1, NULL,&FT1Energy, &anynul, &status);
  //printf("energy %f %f %f\n",FT1Energy,Energy_Min_datafiles,Energy_Max_datafiles);
  if (FT1Energy<=Energy_Min_datafiles || FT1Energy>=Energy_Max_datafiles) return false; //energy cut

  static double FT1ZenithTheta=0;
  fits_read_col (fptr,TDOUBLE,8,i, 1, 1, NULL,&FT1ZenithTheta, &anynul, &status);
  //printf("zt %f %f\n",FT1ZenithTheta,FT1ZenithTheta_Cut);
  if (FT1ZenithTheta>=FT1ZenithTheta_Cut) return false; //ZTheta cut

  return true;
}

double BackgroundEstimator::GimmeCorrectionFactor(short int ie, double MET) {
  
  int idx=ie-1; //ie (the energy index) goes from 1 to =EnergyBins. idx goes from 0 to EnergyBins
  if (pRatiovsTime[idx]==0) return 0;

  if (ie<=0 || ie>Energy_Bins_datafiles) {printf("%s: energy bin that is out of bounds %d submitted\n",__FUNCTION__,ie); return 1;}
  int timebin=pRatiovsTime[idx]->FindBin(MET);
  if (timebin==0 || timebin>pRatiovsTime[idx]->GetNbinsX()) return 0;
  return pRatiovsTime[idx]->GetBinContent(timebin); //CorrectionCoeff = (est-act)/est 
  //printf("met=%f bin=%d cor=%f\n",MET,RatiovsTime[ie]->FindBin(MET),RatiovsTime[ie]->GetBinContent(RatiovsTime[ie]->FindBin(MET)));   
  
}

