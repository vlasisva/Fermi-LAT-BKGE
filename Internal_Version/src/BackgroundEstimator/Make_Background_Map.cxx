//Author: Vlasios Vasileiou <vlasisva@gmail.com>
// $Header: /nfs/slac/g/glast/ground/cvs/GRBAnalysis-scons/BackgroundEstimator/src/BackgroundEstimator/Make_Background_Map.cxx,v 1.3 2011/10/03 12:05:14 vlasisva Exp $
#include "BackgroundEstimator/BackgroundEstimator.h"

ClassImp(BackgroundEstimator)

//Return codes 
//1 NO GTI
//2 NO EXPOSURE (i.e. burst out of FOV)
//0 ALL OK or file exists OK 
int BackgroundEstimator::Make_Background_Map(string FT1_FILE, string FT2_FILE, string GRB_DIR, double Burst_t0, double Burst_Dur,  int verbosity, bool Calc_Residual, bool Save_Earth_Coo_Map){

  const double Burst_t1 = Burst_t0+Burst_Dur;
           
  string ResultsFile = GRB_DIR+"/"+DataClass+"_BackgroundMaps.root";
  FILE * ftemp = fopen(ResultsFile.c_str(),"r");
  bool ZThetaChanged=false;
  if (ftemp) {
    fclose(ftemp); 
    TFile * fres = TFile::Open(ResultsFile.c_str());
    if  (fres->TestBit(TFile::kRecovered) || fres->GetNkeys()<=0) {
	printf("%s: File %s was not closed ok. Will overwrite\n",__FUNCTION__,ResultsFile.c_str()); 
    }
    else if (fres->Get("ERROR")) {
	printf("%s: ERROR code cound in file %s, skipping\n",__FUNCTION__,ResultsFile.c_str()); 
	fres->Close();
	return 1;
    }
    else if (fabs(atof(fres->Get("FT1ZenithTheta_Cut")->GetTitle())-FT1ZenithTheta_Cut)>.1) {
	printf("%s: Zenith Theta cuts different old=%s new=%.0f. Will resimulate sky\n",__FUNCTION__,fres->Get("FT1ZenithTheta_Cut")->GetTitle(),FT1ZenithTheta_Cut);
	ZThetaChanged=true;
    }
    else if (fres->GetNkeys()>1) {
	if (verbosity>1) {
	    printf("%s: File %s exists.. skipping\n",__FUNCTION__,ResultsFile.c_str());
	}
	fres->Close();
	return 0;
    }
    fres->Close();
  }

  TFile * fResults;
  int tries=0;
  ///Sometimes root does not manage to create the output file.. I keep looping until it manages?
  while (1) {
     string cmd="echo >"+ResultsFile;
     system(cmd.c_str());
     fResults = new TFile(ResultsFile.c_str(),"RECREATE");
     ftemp = fopen(ResultsFile.c_str(),"r");
     if (ftemp) {
        fclose(ftemp);
        break;
     }
     printf("%s: Could not create file %s. Retrying.. \n",__FUNCTION__,ResultsFile.c_str());
     tries++;
     if (tries>100){ printf("%s: Won't retry more.. exiting\n",__FUNCTION__); exit(1);}
  }
  
  
  TH2D* hSolidAngle = (TH2D*)fResidualOverExposure->Get("hSolidAngle");
  if (verbosity>1) printf("%s: Creating Background Map\n",__FUNCTION__);
  ///////////////////////////////////////////////////////

  //Check if data file provided has required data
  fitsfile *fptr;
  int status = 0;
  
  long nrows;
  int hdutype;
  double File_t0,File_t1;
  //Open first file
  char COMMENT[2000];
  fits_open_file(&fptr, FT1_FILE.c_str(), READONLY, &status);
  if (status) {
	printf("%s: Can't open fits file %s\n",__FUNCTION__,FT1_FILE.c_str());
        fits_report_error(stderr, status);
	exit(1);
  }
  fits_movabs_hdu(fptr, 2, &hdutype, &status);
  fits_get_num_rows(fptr, &nrows, &status);      
  //fits_read_col (fptr,TDOUBLE,10,1, 1, 1, NULL,&File_t0, &anynul, &status);
  status=0;
  fits_read_keyword(fptr, (char*)"TSTART",  name, COMMENT, &status);
  File_t0=atof(name);
  if (File_t0<1) sscanf(name,"'%lf'",&File_t0); //parsing failed so let's try with '' around the number

  status=0;
  fits_read_keyword(fptr, (char*)"TSTOP",  name, COMMENT, &status); 
  File_t1=atof(name);
  if (File_t1<1)   sscanf(name,"'%lf'",&File_t1); //parsing failed so let's try with '' around the number
  //printf("%d %s %f\n",status,name,File_t1);
  //fits_read_col (fptr,TDOUBLE,10,nrows, 1, 1, NULL,&File_t1, &anynul, &status);
  fits_close_file(fptr, &status);

  if (Burst_t0<File_t0) {printf("%s: data files start at time %f and you requested an estimation for an earlier time (%f).\n",__FUNCTION__,File_t0,Burst_t0); exit(1);}
  if (Burst_t1>File_t1) {
        printf("%s: data files end at time %f and you requested an estimation for a later time (%f).\n",__FUNCTION__,File_t1,Burst_t1);
        exit(1);
  }
  //////////////////////////////////////////////

  //Making plots
  string astring=GRB_DIR+"/Plots.root";
  ftemp = fopen(astring.c_str(),"r");
  if (ftemp && !ZThetaChanged) fclose (ftemp);
  else TOOLS::Make_Plots(0,Burst_Dur,Burst_t0,astring,FT2_FILE, verbosity);

  //Live-time cube
  if (verbosity>1) printf("%s: Creating live-time cube\n",__FUNCTION__);

  if (Calc_Residual) {
     sprintf(name,"%s/burst_ltCube.fits",GRB_DIR.c_str());
     ftemp = fopen(name,"r");
     if (ftemp && !ZThetaChanged) {
        printf("%s: File %s already exists.. skipping creation of it\n",__FUNCTION__,name);
        fclose(ftemp);
     }
     else TOOLS::Run_gtltcube(GRB_DIR, Burst_t0, Burst_t1, FT2_FILE, FT1ZenithTheta_Cut, verbosity);
  }
  //////////////////////////////////////////////////////

  //READ GTIs
  vector <double> GTI_Starts;
  vector <double> GTI_Ends;
  TOOLS::ReadGTI( GTI_Starts,GTI_Ends, FT1_FILE, Burst_t0, Burst_t1);
  
 
  if (GTI_Starts.size()==0) { 
       fResults->cd();
       printf("%s: No GTIs found\n",__FUNCTION__);
       TNamed Data = TNamed("ERROR","NO GTI");
       fResults->cd();
       Data.Write();
       fResults->Close();
       return 1;
  }

  if (verbosity>1) printf("%s: Read %d GTIs\n",__FUNCTION__,(int)GTI_Starts.size());
  /////////////////////////////////////////////////

  double SimEvents;
  TH2F * hSimulatedSky[Energy_Bins_datafiles+1];
  TH2F * hSimulatedSky_Earth[Energy_Bins_datafiles+1];
  hSimulatedSky_Earth[0]=NULL;
  for (int ie=1;ie<=Energy_Bins_datafiles;ie++) {

     if (Save_Earth_Coo_Map) {
        sprintf(name,"hSimulatedSky_Earth_%d",ie);
        hSimulatedSky_Earth[ie] = new TH2F(name,name,30,-110,110,30,-110,110);
        hSimulatedSky_Earth[ie]->SetContour(256);
     }
     else {
        hSimulatedSky_Earth[ie]=NULL;
     }

     sprintf(name,"hSimulatedSky_%d",ie);
     hSimulatedSky[ie] = new TH2F(name,name,L_BINS,-180,180,B_BINS,-90,90);
     hSimulatedSky[ie]->SetContour(200);
  }
  TH2F* hExposureBurst=NULL,*hExposureBurst_1deg=NULL,*hFinalBackground=NULL,*hFinalResidual=NULL,*hFinalSimulatedSky=NULL;
  if (Calc_Residual) {
     hExposureBurst = new TH2F("hExposureBurst","hExposureBurst",L_BINS,-180,180,B_BINS,-90,90);
     hExposureBurst_1deg = new TH2F("hExposureBurst_1deg","hExposureBurst_1deg",360,-180,180,180,-90,90);
     hFinalBackground = new TH2F("hFinalBackground","hFinalBackground",L_BINS,-180,180,B_BINS,-90,90);
     hFinalBackground->SetContour(200);
     hFinalResidual = new TH2F("hFinalResidual","hFinalResidual",L_BINS,-180,180,B_BINS,-90,90);
     hFinalResidual->SetContour(200);
     hFinalSimulatedSky = new TH2F("hFinalSimulatedSky","hFinalSimulatedSky",L_BINS,-180,180,B_BINS,-90,90);
     hFinalSimulatedSky->SetContour(200);
  }


  //EXPOSURE
  if (verbosity) printf("%s: Calculating exposure\n",__FUNCTION__);
  if (Calc_Residual) {
     sprintf(name,"%s/%s_burst_exposure.fits",GRB_DIR.c_str(),DataClass.c_str());
     ftemp = fopen(name,"r");
     if (ftemp && !ZThetaChanged) fclose(ftemp);
     else TOOLS::Run_gtexpcube(GRB_DIR, Burst_t0, Burst_t1, FT2_FILE, DataClass, FT1ZenithTheta_Cut, name, Energy_Min_datafiles, Energy_Max_datafiles, Energy_Bins_datafiles,verbosity);
  }

  sprintf(name,"%s/Plots.root",GRB_DIR.c_str());
  TFile * fPlots = TFile::Open(name);
  Plots_Struct myPlots_Struct;
  myPlots_Struct.hMcIlwainLvsTime    = (TH1F*)fPlots->Get("hMcIlwainLvsTime");
  if (!myPlots_Struct.hMcIlwainLvsTime) {printf("%s: Can't read plots from file %s\n",__FUNCTION__,name); exit(1);}
  myPlots_Struct.hPtRazvsTime        = (TH1F*)fPlots->Get("hPtRazvsTime");
  myPlots_Struct.hPtDeczvsTime       = (TH1F*)fPlots->Get("hPtDeczvsTime");
  myPlots_Struct.hPtRaxvsTime        = (TH1F*)fPlots->Get("hPtRaxvsTime");
  myPlots_Struct.hPtDecxvsTime       = (TH1F*)fPlots->Get("hPtDecxvsTime");
  myPlots_Struct.hRAZenithvsTime     = (TH1F*)fPlots->Get("hRAZenithvsTime");
  myPlots_Struct.hDecZenithvsTime    = (TH1F*)fPlots->Get("hDecZenithvsTime");
  myPlots_Struct.hRockingAnglevsTime = (TH1F*)fPlots->Get("hRockingAnglevsTime");

  if (verbosity>1) printf("%s: Simulating CR background\n",__FUNCTION__);
  SimulateSky(myPlots_Struct, hSimulatedSky, GTI_Starts,GTI_Ends,Energy_Bins_datafiles, hSimulatedSky_Earth, 0, 0);
  
  for (int ie=1;ie<=Energy_Bins_datafiles;ie++) { 
     if (verbosity>0) TOOLS::ProgressBar(ie-1,Energy_Bins_datafiles);
     
     //EXPOSURE
     if (Calc_Residual) {
        sprintf(name,"%s/%s_burst_exposure.fits",GRB_DIR.c_str(),DataClass.c_str());
        TOOLS::ReadExposureMap(name,hExposureBurst,ie,verbosity);
        TOOLS::ReadExposureMap(name,hExposureBurst_1deg,ie,verbosity);
        if (hExposureBurst->Integral()==0) {
    	    printf("%s: Integral of exposure is zero? bailing out.. Check or delete file %s\n",__FUNCTION__,name);
    	    return 2;
    	}
     }


     fResults->cd();
     //sprintf(name,"hSimulatedSky[ie]_nosolid_%d",ie);
     //hSimulatedSky[ie]->Write(name);

     SimEvents = hSimulatedSky[ie]->GetEntries();
     if (SimEvents<=0) {printf("%s: Problem simulating sky (return %f)\n",__FUNCTION__,hSimulatedSky[ie]->Integral()); exit(1);}
 
     /////////////////////////////////////////////////////
     hSimulatedSky[ie]->Divide(hSolidAngle);
     fResults->cd();
     //RESIDUAL
     if (Calc_Residual) {
        sprintf(name,"hResidual_Over_Exposure_%d;1",ie);
        TH2F* htemp = (TH2F*)fResidualOverExposure->Get(name);
        sprintf(name,"hResidualBurst_%d",ie);
        TH2F* hResidualBurst = (TH2F*)htemp->Clone(name);
        htemp->Delete();
        hResidualBurst->SetTitle(name);
        hResidualBurst->Multiply(hExposureBurst);
        
	hFinalResidual->Add(hResidualBurst);
    	hFinalSimulatedSky->Add(hSimulatedSky[ie]);
	sprintf(name,"hExposure_Burst_%d",ie);
	hExposureBurst->Write(name);
	sprintf(name,"hExposure_Burst_1deg_%d",ie);
        hExposureBurst_1deg->Write(name);
        sprintf(name,"hResidualBurst_%d",ie);
        hResidualBurst->Write(name);
        hResidualBurst->Delete();
    }

    sprintf(name,"hSimulatedSky_%d",ie);
    hSimulatedSky[ie]->SetTitle(name);
    hSimulatedSky[ie]->Write(name);
    if (hSimulatedSky_Earth[ie]) hSimulatedSky_Earth[ie]->Write(name);

  }
  if (Calc_Residual) {
     hSolidAngle->Write();
     hFinalBackground->Add(hFinalSimulatedSky);
     hFinalBackground->Add(hFinalResidual);
     hFinalBackground->Write();   hFinalBackground->Delete();
     hFinalResidual->Write();     hFinalResidual->Delete();
     hFinalSimulatedSky->Write(); hFinalSimulatedSky->Delete();
  }
  printf("\n");
  //////////////////////

  sprintf(name,"%e-%e-%d",Energy_Min_datafiles,Energy_Max_datafiles,Energy_Bins_datafiles);
  TNamed Data = TNamed("Energy_Data",name);
  Data.Write();

  sprintf(name,"%.2f-%.2f",Burst_t0,Burst_Dur);
  Data = TNamed("Time_Data",name);
  Data.Write();

  sprintf(name,"%.1f",FT1ZenithTheta_Cut);
  Data = TNamed("FT1ZenithTheta_Cut",name);
  Data.Write();

  sprintf(name,"%.2f",EstimatorVersion);
  Data = TNamed("Estimator_Version",name);
  Data.Write();

  sprintf(name,"%.1f %.1f %.1f %.1f %.1f",Residuals_version,RateFit_version,ThetaPhiFits_version,TimeCorrectionFactors_version,EastWest_version);
  Data = TNamed("DataFiles_Version",name);
  Data.Write();

  Data = TNamed("GRB_NAME",(TOOLS::GetS("GRB_NAME")).c_str());  Data.Write();
  Data = TNamed("FT1_FILE",FT1_FILE.c_str());  Data.Write();
  Data = TNamed("FT2_FILE",FT2_FILE.c_str());  Data.Write();
  Data = TNamed("DataClass",DataClass.c_str());  Data.Write();

  if (Calc_Residual) {
    hExposureBurst->Delete();
    hExposureBurst_1deg->Delete();
  }
  for (int ie=1;ie<=Energy_Bins_datafiles;ie++) {
     hSimulatedSky[ie]->Delete();
     if (hSimulatedSky_Earth[ie]) hSimulatedSky_Earth[ie]->Delete();
  }

  fResults->Close();
  fPlots->Close();
  return 0;
}

