//Author: Vlasios Vasileiou <vlasisva@gmail.com>
// $Header: /nfs/slac/g/glast/ground/cvs/GRBAnalysis-scons/BackgroundEstimator/src/BackgroundEstimator/CalcResiduals.cxx,v 1.4 2013/11/13 07:56:16 vlasisva Exp $
#include "BackgroundEstimator/BackgroundEstimator.h"

#define DEBUG

ClassImp(BackgroundEstimator)

void BackgroundEstimator::CalcResiduals(string FitsAllSkyFile){

 const float hSimulatedSky_Earth_Map_Min_B=20;
 sprintf(name,"%.2f",Residuals_version);
 TNamed * File_Version = new TNamed("version",name);


 /////////////////////////////////////////////////////////////////

 TH2F * hSignal_True[Energy_Bins_datafiles+1],*hSignal_True_Earth[Energy_Bins_datafiles+1];
 sprintf(name,"%s/Residual_True_%s.root",DataDir.c_str(),DataClass.c_str());
 FILE * ftemp = fopen(name,"r");
 TFile * ftrue=NULL;
 if (ftemp) {
     ftrue = TFile::Open(name);
     for (int ie=1;ie<=Energy_Bins_datafiles;ie++) {
        sprintf(name,"hSignal_True_%d",ie);
        hSignal_True[ie] = (TH2F*)ftrue->Get(name);
        sprintf(name,"hSignal_True_Earth_%d",ie);
        hSignal_True_Earth[ie] = (TH2F*)ftrue->Get(name);
     }
     fclose(ftemp);
 }
 else {
    printf("%s: Filling Skymap\n",__FUNCTION__);
    ftrue = new TFile(name,"RECREATE");
    for (int ie=1;ie<=Energy_Bins_datafiles;ie++) {
        sprintf(name,"hSignal_True_%d",ie);
        hSignal_True[ie] = new TH2F(name,"hSignal_True",L_BINS,-180,180,B_BINS,-90,90);
        sprintf(name,"hSignal_True_Earth_%d",ie);
        hSignal_True_Earth[ie] = new TH2F(name,name,30,-110,110,30,-110,110);
        hSignal_True_Earth[ie]->SetName(name);
        hSignal_True_Earth[ie]->SetTitle(name);
        hSignal_True_Earth[ie]->SetContour(256);
    }
     
    fitsfile *fptr;
    int ifile=0;
    ftemp = fopen(FitsAllSkyFile.c_str(),"r");
    while (fscanf(ftemp,"%s",name)==1) {
       printf("%d\r",ifile); fflush(0);
       ifile++;
       long nrows;int ncols;
       int status=0,hdutype,anynul;
       fits_open_file(&fptr, name, READONLY, &status);
       if (status) {printf("error opening file %s\n",name); throw std::runtime_error("");}

       // Read event data
       fits_movabs_hdu(fptr, 2, &hdutype, &status);
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

        for (int i=1;i<=nrows;i++) {
 
             if (!PassesCuts(fptr,i,format)) continue;
    
             double FT1Energy;
             fits_read_col (fptr,TDOUBLE,1,i, 1, 1, NULL,&FT1Energy, &anynul, &status);
             short int iebin=Energy2Bin(FT1Energy);
             if (iebin==0 || iebin>Energy_Bins_datafiles) {
                printf("%s: iebin=%d e=%f skipped\n",__FUNCTION__,iebin,FT1Energy);
                continue;
             }
        
             double PtTime;
             fits_read_col (fptr,TDOUBLE,10,i, 1, 1, NULL,&PtTime, &anynul, &status);
             if (PtTime<StartTime) continue;
             if (PtTime>EndTime) break;
        
             double FT1L,FT1B;
             fits_read_col (fptr,TDOUBLE,4,i, 1, 1, NULL,&FT1L, &anynul, &status);
             fits_read_col (fptr,TDOUBLE,5,i, 1, 1, NULL,&FT1B, &anynul, &status);
             if (FT1L>180) FT1L-=360;

             hSignal_True[iebin]->Fill(FT1L,FT1B);

             //East-West stuff         
             if (fabs(FT1B)>hSimulatedSky_Earth_Map_Min_B) {
                double FT1EarthAzimuth,FT1ZenithTheta;
                //Do not include events with low galactic latitudes since we want to compare only CR signals
                //The simulated part should also have an identical cut.
                fits_read_col (fptr,TDOUBLE,8,i, 1, 1, NULL,&FT1ZenithTheta, &anynul, &status);
                fits_read_col (fptr,TDOUBLE,9,i, 1, 1, NULL,&FT1EarthAzimuth, &anynul, &status);
                double x=sin(FT1EarthAzimuth*TMath::Pi()/180.)*FT1ZenithTheta;
                double y=cos(FT1EarthAzimuth*TMath::Pi()/180.)*FT1ZenithTheta;
                hSignal_True_Earth[iebin]->Fill(x,y);
             }
        }

        fits_close_file(fptr, &status);
     }
  
     ftrue->cd();
     for (int ie=1;ie<=Energy_Bins_datafiles;ie++) {
        hSignal_True_Earth[ie]->Write(); 
        hSignal_True[ie]->Write();
     }
     File_Version->Write();
 }
 ////////////////////////////////////////////////////////////////


 TH2F * hSimulatedSky[Energy_Bins_datafiles+1]; 
 TH2F * hSimulatedSky_Earth[Energy_Bins_datafiles+1];
 sprintf(name,"%s/Residual_Sim_%s.root",DataDir.c_str(),DataClass.c_str());
 ftemp = fopen(name,"r");
 TFile * fsim=NULL;
 if (ftemp) {
    fsim= TFile::Open(name);
    for (int ie=1;ie<=Energy_Bins_datafiles;ie++) {
       sprintf(name,"hSimulatedSky_%d",ie);
       hSimulatedSky[ie]=(TH2F*)fsim->Get(name);
       sprintf(name,"hSimulatedSky_Earth_%d",ie);
       hSimulatedSky_Earth[ie]=(TH2F*)fsim->Get(name);
    }
    fclose(ftemp);
 }
 else {
     printf("%s: Simulating Sky\n",__FUNCTION__);
     printf("%s: --- Reading GTIs\n",__FUNCTION__);
     vector <double> GTI_Starts;
     vector <double> GTI_Ends;
     TOOLS::ReadGTI(GTI_Starts, GTI_Ends, FitsAllSkyFile, StartTime, EndTime);

     fsim = new TFile(name,"RECREATE");
     for (int ie=1;ie<=Energy_Bins_datafiles;ie++) {
        sprintf(name,"hSimulatedSky_%d",ie);
        hSimulatedSky[ie] = new TH2F(name,name,L_BINS,-180,180,B_BINS,-90,90);
        hSimulatedSky[ie]->SetContour(200);
        hSimulatedSky[ie]->GetXaxis()->SetTitle("Longitude (deg)");
        hSimulatedSky[ie]->GetYaxis()->SetTitle("Latitude (deg)");
        
        sprintf(name,"hSimulatedSky_Earth_%d",ie);
        hSimulatedSky_Earth[ie] = new TH2F(name,name,30,-110,110,30,-110,110);
        hSimulatedSky_Earth[ie]->SetContour(256);
     }
     
     sprintf(name,"%s/Plots_%s.root",DataDir.c_str(),DataClass.c_str()); 
     TFile * fPlots = TFile::Open(name);
     Plots_Struct myPlots_Struct;
     myPlots_Struct.hMcIlwainLvsTime    = (TH1F*)fPlots->Get("hMcIlwainLvsTime");  
     if (!myPlots_Struct.hMcIlwainLvsTime) {printf("%s: Can't read plots from file %s\n",__FUNCTION__,name); throw std::runtime_error("");}
     myPlots_Struct.hPtRazvsTime        = (TH1F*)fPlots->Get("hPtRazvsTime");
     myPlots_Struct.hPtDeczvsTime       = (TH1F*)fPlots->Get("hPtDeczvsTime");
     myPlots_Struct.hPtRaxvsTime        = (TH1F*)fPlots->Get("hPtRaxvsTime");
     myPlots_Struct.hPtDecxvsTime       = (TH1F*)fPlots->Get("hPtDecxvsTime");
     myPlots_Struct.hRAZenithvsTime     = (TH1F*)fPlots->Get("hRAZenithvsTime");
     myPlots_Struct.hDecZenithvsTime    = (TH1F*)fPlots->Get("hDecZenithvsTime");
     myPlots_Struct.hRockingAnglevsTime = (TH1F*)fPlots->Get("hRockingAnglevsTime");

     SimulateSky(myPlots_Struct, hSimulatedSky, GTI_Starts,GTI_Ends,Energy_Bins_datafiles, hSimulatedSky_Earth, 1000, hSimulatedSky_Earth_Map_Min_B);
     fPlots->Close();
     fsim->cd();
     for (int ie=1;ie<=Energy_Bins_datafiles;ie++) {
        if (hSimulatedSky[ie]->Integral()<=0) {printf("%s: Error simulating sky (result=%f)\n",__FUNCTION__,hSimulatedSky[ie]->Integral()); throw std::runtime_error("");}
        hSimulatedSky[ie]->Write();
        hSimulatedSky_Earth[ie]->Write();   
     }
     File_Version->Write();
 }

 ///////////////////////////////////////////////////////////////////////////////


 sprintf(name,"%s/Residual_%s.root",DataDir.c_str(),DataClass.c_str());
 ftemp = fopen(name,"r");
 if (ftemp) fclose (ftemp);
 else {
    TH2F * hResidual[Energy_Bins_datafiles];
    printf("%s: Calculating residual\n",__FUNCTION__);
    TFile * fResidual=new TFile(name,"RECREATE");
    for (int ie=1;ie<=Energy_Bins_datafiles;ie++) {
        sprintf(name,"hResidual_%d",ie);
        hResidual[ie]= new TH2F(name,name,L_BINS,-180,180,B_BINS,-90,90);
        hResidual[ie]->GetXaxis()->SetTitle("Longitude (deg)");
        hResidual[ie]->GetYaxis()->SetTitle("Latitude (deg)");
    
        hResidual[ie]->Add(hSignal_True[ie]);
        hResidual[ie]->Add(hSimulatedSky[ie],-1); //-1  = subtract
        hResidual[ie]->Write();  
    }
    File_Version->Write();
    fResidual->Close();
 }

 //////////////////////////////////////////////////////////////

 sprintf(name,"%s/EastWest_Correction_%s_%.1f.root",DataDir.c_str(),DataClass.c_str(),EastWest_version);
 ftemp = fopen(name,"r"); 
 if (ftemp) fclose(ftemp);
 else {
     printf("%s: Calculating EW Correction\n",__FUNCTION__);
     TFile* fCorrection_EW = new TFile(name,"RECREATE");
     for (int ie=1;ie<=Energy_Bins_datafiles;ie++) {
        //Scan the two maps for bins with low statistics. 
        //If any pair of bins has less then 50 events, set the bin of the true map equal to zero
        for (int i=1;i<=hSimulatedSky_Earth[ie]->GetNbinsX();i++) {
           for (int j=1;j<=hSimulatedSky_Earth[ie]->GetNbinsY();j++) {
               if (hSimulatedSky_Earth[ie]->GetBinContent(i,j)<100 || hSignal_True_Earth[ie]->GetBinContent(i,j)<100) 
                   hSignal_True_Earth[ie]->SetBinContent(i,j,0);
          }
       }
       hSignal_True_Earth[ie]->Divide(hSimulatedSky_Earth[ie]);
       sprintf(name,"hEW_Correction_TrueOverEst_%d",ie);
       hSignal_True_Earth[ie]->SetTitle(name);
       hSignal_True_Earth[ie]->SetName(name);
       hSignal_True_Earth[ie]->Write();
    }
    fCorrection_EW->Close();
 }
     

 ftrue->Close();
 fsim->Close();
}


