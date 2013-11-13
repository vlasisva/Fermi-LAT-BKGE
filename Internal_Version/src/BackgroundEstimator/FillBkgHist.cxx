//Author: Vlasios Vasileiou <vlasisva@gmail.com>
#include "BackgroundEstimator/BackgroundEstimator.h"
#include "TGraph.h"


//return codes 
//-1 error 
//0 all good

//type 1=Galactic Gammas, 2=CR+ExtraG, 3=both
int BackgroundEstimator::FillBackgroundHist(string GRB_DIR, TH1F * hROI, double RA_BURST, double DEC_BURST, short int type, int verbosity) {
  string ResultsFile = GRB_DIR+"/"+DataClass+"_BackgroundMaps.root";
  TFile * fResults = TFile::Open(ResultsFile.c_str());
  if (!fResults || fResults->GetNkeys()==0) {
        printf("%s: Problem opening file %s\n",__FUNCTION__,ResultsFile.c_str()); 
        return -1;
  }
  double MET; sscanf(fResults->Get("Time_Data")->GetTitle(),"%lf-%*f",&MET);
  double FT1ZenithTheta_Cut_file;  sscanf(fResults->Get("FT1ZenithTheta_Cut")->GetTitle(),"%lf",&FT1ZenithTheta_Cut_file);

  if (fabs(FT1ZenithTheta_Cut_file-FT1ZenithTheta_Cut)>0.1) {
       printf("%s: Different ZTheta cuts! configured=%f, in bkg file=%f\n",__FUNCTION__,FT1ZenithTheta_Cut,FT1ZenithTheta_Cut_file);
       return -1;
  }

  if (fabs(pow(10,hROI->GetXaxis()->GetXmin())-Energy_Min_user)>0.1 || fabs(pow(10,hROI->GetXaxis()->GetXmax())-Energy_Max_user)>0.1 || hROI->GetNbinsX()!=Energy_Bins_user) {
     printf("%s: Different energy configuration between hROI (%f/%f/%d) and user-passed values (%f/%f/%d) \n",__FUNCTION__,
         pow(10,hROI->GetXaxis()->GetXmin()),pow(10,hROI->GetXaxis()->GetXmax()),hROI->GetNbinsX(),Energy_Min_user,Energy_Max_user,Energy_Bins_user);
    return -1;
  }


  TH1F hBkg = TH1F("hCtsvsEnergy_Est","Background Estimate",hROI->GetNbinsX(),hROI->GetXaxis()->GetXmin(),hROI->GetXaxis()->GetXmax());
  hBkg.GetXaxis()->SetTitle("log_{10}	(Energy/MeV)");
  hBkg.GetYaxis()->SetTitle("Events/bin");

  double L_BURST,B_BURST;
  TOOLS::Galactic(RA_BURST,DEC_BURST,&L_BURST,&B_BURST);

  if (L_BURST>180) L_BURST-=360;

  TH2F * hMap[Energy_Bins_datafiles+2];

  for (int ie=1;ie<=Energy_Bins_datafiles;ie++) {
      if (type==1) {
          sprintf(name,"hResidualBurst_%d;1",ie);
          hMap[ie]=(TH2F*)fResults->Get(name);
          if (!hMap[ie]){fResults->Close(); printf ("%s: no %s!\n",__FUNCTION__,name); return -1;}
      }
      else if (type==2) {
          sprintf(name,"hSimulatedSky_%d;1",ie);
          hMap[ie]=(TH2F*)fResults->Get(name);
          if (!hMap[ie]){fResults->Close(); printf ("%s: no %s!\n",__FUNCTION__,name);return -1;}
      }
      else if (type==3) {
          sprintf(name,"hSimulatedSky_%d;1",ie);
          TH2F* htemp=(TH2F*)fResults->Get(name);
          if (!htemp){fResults->Close(); printf ("%s: no %s!\n",__FUNCTION__,name);return -1;}

          sprintf(name,"hResidualBurst_%d;1",ie);
          hMap[ie]=(TH2F*)fResults->Get(name);
          if (!hMap[ie]){fResults->Close(); printf ("%s: no %s!\n",__FUNCTION__,name);return -1;}
    
          hMap[ie]->Add(htemp);
          htemp->Delete();
      }
      else {printf("%s: Unknown particle-species type %d\n",__FUNCTION__,type); throw std::runtime_error("");}
  }

  //CALCULATE SPECTRAL TEMPLATE
  bool OutOfFOV=false;
  TGraph gSpectralTemplate;
  TH1F hBkg_old = TH1F("hbkg_Old","hbkg_old",Energy_Bins_datafiles,log10(Energy_Min_datafiles),log10(Energy_Max_datafiles));
  if (verbosity>=5) printf("%s: Default binning: %d\n",__FUNCTION__,UsingDefaultBinning);
  if (!UsingDefaultBinning) {
        double BKG_old[Energy_Bins_datafiles],E_old[Energy_Bins_datafiles];

        for (int ibin_old = 1;ibin_old<=Energy_Bins_datafiles;ibin_old++) {
             E_old[ibin_old-1]  = hBkg_old.GetBinCenter(ibin_old);
             float CorrectionFactor=GimmeCorrectionFactor(ibin_old,MET);
             if (fabs(CorrectionFactor)>0.3) printf("%s: weird correction factor.. %f\n",__FUNCTION__,CorrectionFactor);

             BKG_old[ibin_old-1]= TOOLS::Integrate(hMap[ibin_old], L_BURST,B_BURST, 10)*(1-CorrectionFactor);
             if (type!=1 && BKG_old[ibin_old-1]<=0) {
                    printf("%s: weird background %f at bin %d (type=%d) correction factor=%f\n",__FUNCTION__,BKG_old[ibin_old-1],ibin_old,type,CorrectionFactor);
                    //if (BKG_old[ibin_old-1]) {
                      // sprintf(name,"hResidualBurst_%d;1",ibin_old);
                      // TH2F * htemp=(TH2F*)fResults->Get(name);
                      // printf("gamma=%f\n",TOOLS::Integrate(htemp, L_BURST,B_BURST, 10));
                      // sprintf(name,"hSimulatedSky_%d;1",ibin_old);
                      // htemp = (TH2F*)fResults->Get(name);
                      // printf("CR=%f\n",TOOLS::Integrate(htemp, L_BURST,B_BURST, 10));
                      // delete htemp;
                    //}
             }
             /*
             if (type!=0 && BKG_old[ibin_old-1]==0 && ibin_old==1) {
                   printf("%s: weird background %f at bin %d (type=%d) -- Assuming the source is outside the FOV\n",__FUNCTION__,BKG_old[ibin_old-1],ibin_old,type);
                   for (int ibin_old_ = 1;ibin_old_<=Energy_Bins_datafiles;ibin_old_++) {
                      E_old[ibin_old_-1]  = hBkg_old.GetBinCenter(ibin_old_);
                      BKG_old[ibin_old_-1]= 0;
                   }
                   OutOfFOV=true;
                   break;
             }
             */
        }

        gSpectralTemplate = TGraph(Energy_Bins_datafiles,E_old,BKG_old);
        TH1F hSpectralTemplate = TH1F("hSpectralTemplate","hSpectralTemplate",200,log10(Energy_Min_datafiles),log10(Energy_Max_datafiles));
        for (int i=1;i<=200;i++) {
            double interpolated_value;
            double le=hSpectralTemplate.GetBinCenter(i);
            if (le>log10(1.9)) interpolated_value= gSpectralTemplate.Eval(le,0,"S");
            else interpolated_value= gSpectralTemplate.Eval(le); //linear extrapolation
            hSpectralTemplate.SetBinContent(i,interpolated_value);
        }
  }
  if (OutOfFOV) {
     for (int i_new=1;i_new<=hBkg.GetNbinsX();i_new++)
         hBkg.SetBinContent(i_new,0);
  }
  else if (UsingDefaultBinning) {
      for (int i_new=1;i_new<=hBkg.GetNbinsX();i_new++) {  //Loop over energy bins
          double BKG;
          float CorrectionFactor=GimmeCorrectionFactor(i_new,MET);
          if (fabs(CorrectionFactor)>0.3) {printf("%s: weird correction factor.. %f\n",__FUNCTION__,CorrectionFactor); throw std::runtime_error("");}
          BKG = TOOLS::Integrate(hMap[i_new], L_BURST, B_BURST, hROI->GetBinContent(i_new))*(1-CorrectionFactor);
          if (BKG<=0 && type!=1) {
              printf("%s: A bkg bin (%d) was negative or zero. Setting it to the value of the previous bin.\n",__FUNCTION__,i_new);
              BKG=hBkg.GetBinContent(i_new-1);
          }
          hBkg.SetBinContent(i_new,BKG);  
          if (verbosity) TOOLS::ProgressBar(i_new-1,hBkg.GetNbinsX());
      }
  }
  else {
    for (int i_new=1;i_new<=hBkg.GetNbinsX();i_new++) {  //Loop over new energy bins
         double bin_lE_min = hBkg.GetXaxis()->GetBinLowEdge(i_new)
            ,bin_lE_max = hBkg.GetXaxis()->GetBinUpEdge(i_new);
   
         int bin_min_Old = hBkg_old.FindBin(bin_lE_min),
             bin_max_Old = hBkg_old.FindBin(bin_lE_max); 
         if (bin_min_Old==0)          bin_min_Old=1;
         if (bin_max_Old>Energy_Bins_datafiles) bin_max_Old=Energy_Bins_datafiles;
         double BKG=0;

         for (int ibin_old = bin_min_Old;ibin_old<=bin_max_Old;ibin_old++) {
               //printf("%d %f %f\n",ibin_old,pow(10,bin_lE_min),pow(10,bin_lE_max));
               double BKG_WHOLE_OLD_BIN=TOOLS::Integrate(hMap[ibin_old], L_BURST, B_BURST, hROI->GetBinContent(i_new));
               BKG_WHOLE_OLD_BIN*=(1-GimmeCorrectionFactor(ibin_old,MET));
               double bin_lE_min_old = hBkg_old.GetXaxis()->GetBinLowEdge(ibin_old)
                     ,bin_lE_max_old = hBkg_old.GetXaxis()->GetBinUpEdge(ibin_old);
               //printf("ibin_old %d %f %f new :%d %f %f BKG_WHOLE=%e ROI=%.1e\n",ibin_old,bin_lE_min_old,bin_lE_max_old,
               //                                          i_new,bin_lE_min, bin_lE_max,BKG_WHOLE_OLD_BIN,hROI->GetBinContent(i_new));
            //see if we have a subset or the whole bin
               if (bin_lE_min_old>=bin_lE_min && bin_lE_max_old<=bin_lE_max) {
                   //printf("whole bin\n");
                   BKG+=BKG_WHOLE_OLD_BIN; //old bin is full contained 
               }
               else {
                   double Intersection_lE_min,Intersection_lE_max;
                   Intersection_lE_min=std::max(bin_lE_min_old,bin_lE_min);
                   Intersection_lE_max=std::min(bin_lE_max_old,bin_lE_max);
                   const double tolerance=0.0001;
                   if (fabs(Intersection_lE_min-Intersection_lE_max)<tolerance) continue; //didn't find overlap between old and new bins
                   //printf("intersection at %f %f\n",Intersection_lE_min, Intersection_lE_max);
                   double sum=0,sum_all=0;
                   const int npoints=200;
                   double lE_int   =Intersection_lE_min,
                   Delta_lE_int    =(Intersection_lE_max-Intersection_lE_min),
                   dlE_int         =Delta_lE_int/npoints;
                   //Here I calculate the "ScalingFactor" to scale the contents of one of the original bins to the new bin
                   for (int i=0;i<npoints;i++) {
                       double interpolated_value; 
                       if (lE_int>log10(1.9)) interpolated_value=gSpectralTemplate.Eval(lE_int,0,"S"); 
                       else                   interpolated_value=gSpectralTemplate.Eval(lE_int); 
                       sum+=interpolated_value*dlE_int;
                       //printf("---%d le=%f gs=%f dlE_in=%f sum=%f\n",i,lE_int,gSpectralTemplate.Eval(lE_int,0,"S"),dlE_int,sum);
                       lE_int+=dlE_int;
                   } 

                   lE_int =bin_lE_min_old;
                   dlE_int=hBkg_old.GetXaxis()->GetBinWidth(ibin_old)/npoints;
                   for (int i=0;i<npoints;i++) {
                       sum_all+=gSpectralTemplate.Eval(lE_int,0,"S")*dlE_int;
                       //printf("+++ %d le=%f gs=%f dlE_in=%f sum=%f\n",i,lE_int,gSpectralTemplate.Eval(lE_int,0,"S"),dlE_int,sum_all);
                       lE_int+=dlE_int;
                   } 

                   double ScalingFactor = sum/sum_all;
                   //printf("sum=%e %e scalingfactor=%f\n",sum,sum_all,ScalingFactor);
                   BKG+= BKG_WHOLE_OLD_BIN*ScalingFactor;
               }
          }
          //printf("fill %e at %d\n",BKG,i_new);
          if (BKG<0 && type!=1) { 
              printf("%s: A bkg bin (%d) was negative. Setting it to the value of the previous bin.\n",__FUNCTION__,i_new);
              BKG=hBkg.GetBinContent(i_new-1);
          }
          hBkg.SetBinContent(i_new,BKG);

          if (verbosity) TOOLS::ProgressBar(i_new-1,hBkg.GetNbinsX());
     }
  }
   
  if (GRB_DIR!="") {

     string bkgtype;
     if      (type==1) bkgtype="_GALGAMMAS";
     else if (type==2) bkgtype="_CR_EGAL";
     else if (type==3) bkgtype="";
     char OutputFilename[1000];
     sprintf(OutputFilename,"%s/%s_bkg_%.0f_%.0f%s.root",GRB_DIR.c_str(),DataClass.c_str(),pow(10,hROI->GetXaxis()->GetXmin()),pow(10,hROI->GetXaxis()->GetXmax()),bkgtype.c_str());
   
     TFile * fBkg = new TFile(OutputFilename,"RECREATE");
     TH1F * hExposure = (TH1F*)hROI->Clone("hExposure");
     hExposure->GetXaxis()->SetTitle("log_{10}(Energy/MeV)");
     hExposure->GetYaxis()->SetTitle("Exposure (cm^{2} sec)");
     hExposure->SetTitle("Exposure");
     TOOLS::CalcExposure(fResults, L_BURST,B_BURST, FT1ZenithTheta_Cut, hExposure, GRB_DIR, verbosity);
     hExposure->Write();
     hExposure->Delete();
     hBkg.Write();
     hROI->Write("hROI");
     gSpectralTemplate.Write("gSpectralTemplate");

     sprintf(name,"RA/DEC %.3f %.3f",RA_BURST,DEC_BURST);
     TNamed Data = TNamed("Localization_Data",name);
     Data.Write();

     ((TNamed*)fResults->Get("Estimator_Version"))->Write();
     ((TNamed*)fResults->Get("DataFiles_Version"))->Write();
     ((TNamed*)fResults->Get("FT1ZenithTheta_Cut"))->Write();
     ((TNamed*)fResults->Get("Time_Data"))->Write();
     ((TNamed*)fResults->Get("Energy_Data"))->Write();

     fBkg->Close();
  }
     fResults->Close();
/*
  //Some diagnostic stuff
  for (int i=1;i<=Energy_Bins_datafiles;i++) {
    sprintf(name,"%s/integratedmap_%d.root",GRB_DIR.c_str(),i);
    double bkg= Integrate(hMap[i], L_BURST, B_BURST, hROI->GetBinContent(i),string(name));
    hBkg->SetBinContent(i,bkg);
    printf("old %d %f\n",i,bkg);
  }
  sprintf(name,"%s/bkg_orig.root",GRB_DIR.c_str());
  TFile * fbkg_orig = new TFile(name,"RECREATE");
  hBkg.Write();
  hROI->Write();
  fbkg_orig->Close();  
*/

 return 0;

}

