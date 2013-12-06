//Author: Vlasios Vasileiou <vlasisva@gmail.com>
//$Header: /nfs/slac/g/glast/ground/cvs/GRBAnalysis-scons/BackgroundEstimator/src/BackgroundEstimator/BackgroundEstimator_ext.cxx,v 1.1 2013/11/28 13:20:18 vlasisva Exp $
#include "BackgroundEstimator/BackgroundEstimator_ext.h"

BackgroundEstimator_ext::~BackgroundEstimator_ext() 
{//do some cleanup
 if(fResidualOverExposure) fResidualOverExposure->Close();
 if(fRateFit)              fRateFit->Close();
 if(fThetaPhiFits)         fThetaPhiFits->Close();
 if(fCorrectionFactors)    fCorrectionFactors->Close();
}

//#define DEBUG

void BackgroundEstimator_ext::CreateDataFiles(string FitsAllSkyFilesList, string FT2_FILE, double StartTime_user, double EndTime_user, float ZenithThetaCut){

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


  if (StartTime==EndTime) {printf("%s:StartTime==EndTime?? \n",__FUNCTION__); throw std::runtime_error("");} //this also includes when they are both zero
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

