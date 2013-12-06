//Author: Vlasios Vasileiou <vlasisva@gmail.com>
//$Header: /nfs/slac/g/glast/ground/cvs/GRBAnalysis-scons/BackgroundEstimator/src/BackgroundEstimator/Calc_TimeCorrectionFactors.cxx,v 1.4 2013/11/28 13:20:18 vlasisva Exp $

#include "BackgroundEstimator/BackgroundEstimator_ext.h"
#include <algorithm>
#include "TProfile.h"
#include "TLine.h"
#include "TGraphErrors.h"

void BKGE_NS_EXT::Calc_TimeCorrectionFactors(vector<string> GRB_folders, vector <double> METs, string Dataclass, double MinE, double MaxE, int NBins){
 gStyle->SetLineScalePS(0.01);
 gStyle->SetPadTopMargin(0.07);
 gStyle->SetPadLeftMargin(0.13);
 gStyle->SetPadRightMargin(0.011);
 gStyle->SetPadBottomMargin(0.15);

 const int DURATION=10000;
 char name[1000],name2[1000];

 const double MinMET=*min_element(METs.begin(),METs.end())+DURATION/2;
 const double MaxMET=*max_element(METs.begin(),METs.end())+DURATION/2;
 const double lMinE = log10(MinE);
 const double lMaxE = log10(MaxE);
 const double dlMinE = (lMaxE-lMinE)/NBins;

 TProfile * hRatiovsTime[NBins];
 TGraph * gRatiovsTime[NBins];
 for (int i=0;i<NBins;i++) {
    hRatiovsTime[i]= NULL;
 }

 TCanvas * cRatiovsTime = new TCanvas("cRatiovsTime","cRatiovsTime",1400,400); 
 cRatiovsTime->Divide(5,4,1e-3,1e-3);
 double Sig[NBins][METs.size()],Back[NBins][METs.size()];
 double ratio[NBins][METs.size()],ratio_err[NBins][METs.size()];

 for (unsigned int iSig=0;iSig<METs.size();iSig++) {
 
      sprintf(name,"%s/%s_bkg_%.0f_%.0f.root",GRB_folders[iSig].c_str(),Dataclass.c_str(),MinE,MaxE);
      TFile * fEst = TFile::Open(name);
      if (!fEst) { printf("no %s\n",name);  continue;    }
      TH1F * hEst = (TH1F*)fEst->Get("hCtsvsEnergy_Est");
      if (!hEst){printf("no hBurst in %s\n",name);  fEst->Close(); continue;}

      sprintf(name,"%s/%s_sig_%.0f_%.0f.root",GRB_folders[iSig].c_str(),Dataclass.c_str(),MinE,MaxE);
      TFile * fBurst = TFile::Open(name);
      if (!fBurst) { printf("no %s\n",name);     continue;    }
      TH1F * hBurst = (TH1F*)fBurst->Get("hCtsvsEnergy_Burst");
      if (!hBurst){printf("no hBurst in %s\n",name);  fBurst->Close(); continue;}
      
      for (int iE=0;iE<NBins;iE++) {
      /*
          if (iE==19) {
             Back[iE][iSig]=hEst->Integral(1,10);
             Sig[iE][iSig] =hBurst->Integral(1,10);
             ratio[iE][iSig]=ratio_err[iE][iSig]=0;
             if (Back[iE][iSig]>40 && Sig[iE][iSig]!=0) {
                ratio[iE][iSig]    = Back[iE][iSig]/Sig[iE][iSig];
             }
             continue;
          }
      */
          Back[iE][iSig]=hEst->GetBinContent(iE+1);
          Sig[iE][iSig] =hBurst->GetBinContent(iE+1);
          ratio[iE][iSig]=ratio_err[iE][iSig]=0;
          if (Back[iE][iSig]>20 && Sig[iE][iSig]!=0) {
//             ratio[iE][iSig]    = Back[iE][iSig]/Sig[iE][iSig];
             ratio[iE][iSig]    = (Back[iE][iSig]-Sig[iE][iSig])/Back[iE][iSig];
             //ratio[iE][iSig]    = (Back[iE][iSig]-Sig[iE][iSig])/Back[iE][iSig];
             //ratio_err[iE][iSig]= ratio[iE][iSig]/sqrt(Back[iE][iSig]);
          } 
          //else printf("skip %d %f\n",iE,Back[iE][iSig]);
      }
      fBurst->Close();
      fEst->Close();
 }
 
 
  for (int iE=0;iE<NBins;iE++) {
     //Make hRatiovsTime Profiles
     int bins=50;
     sprintf(name,"RatiovsTime_%d",iE);
     sprintf(name2,"Ratio vs Time, %.0f<E<%.0f MeV",pow(10,lMinE+(iE)*dlMinE),pow(10,lMinE+(iE+1)*dlMinE));

     while (bins>=5) {
        if (hRatiovsTime[iE]) hRatiovsTime[iE]->Delete();
        hRatiovsTime[iE]= new TProfile(name,name2,bins,MinMET,MaxMET);
        for (unsigned int iSig=0;iSig<METs.size();iSig++) {
             if (ratio[iE][iSig]!=0) hRatiovsTime[iE]->Fill(METs[iSig]+DURATION/2,ratio[iE][iSig]);
        }         
        bool hist_ok=true;
        for (int ibin=1;ibin<=bins && hist_ok;ibin++) {
            if (hRatiovsTime[iE]->GetBinError(ibin)>0.05) {
               //printf("ie=%d bins=%d bin=%d err=%f\n",iE,ibin,bins,hRatiovsTime[iE]->GetBinError(ibin));
               hist_ok=false;
            }   
        }
        if (hist_ok==true)  break;
        bins--;   
     }
     //for (int i=1;i<=hRatiovsTime[iE]->GetNbinsX();i++) printf("%d %f\n",i,hRatiovsTime[iE]->GetBinContent(i));
    
     cRatiovsTime->cd(iE+1);
//     TH1F * hh = cRatiovsTime->DrawFrame(MinMET*0.99,0.75,MaxMET*1.01,1.25);
     //TH1F * hh = cRatiovsTime->DrawFrame(MinMET*0.99,-0.25,MaxMET*1.01,0.25);
     
     //TH1F * hh = cRatiovsTime->DrawFrame(MinMET,-0.2,MaxMET,0.2);
     TLine * l1 = new TLine(MinMET,0,MaxMET,0);
     hRatiovsTime[iE]->SetLineWidth(3);

     gRatiovsTime[iE] = new TGraph(); 
    
     
     int ip=0;
     for (int is=0;is<METs.size();is++) {
        if (ratio[iE][is]!=0) {
           gRatiovsTime[iE]->SetPoint(ip,METs[is]+DURATION/2,ratio[iE][is]);
           //gRatiovsTime[iE]->SetPointError(ip,0,ratio_err[iE][is]);
           ip++;
        }   
     }
     gRatiovsTime[iE]->GetXaxis()->SetTimeDisplay(1);
     gStyle->SetTimeOffset(978336000);
     gRatiovsTime[iE]->GetXaxis()->SetTimeFormat("%d/%m/%y");
     gRatiovsTime[iE]->GetXaxis()->SetTitle("Date");
     gRatiovsTime[iE]->GetYaxis()->SetTitle("(N_{est}-N_{act})/N_{est}");
     sprintf(name,"fCorrection_%d",iE);
     //gRatiovsTime[iE]->SetLineColor(2);
     gRatiovsTime[iE]->SetMarkerStyle(20);
     gRatiovsTime[iE]->SetMarkerColor(kGray);
     gRatiovsTime[iE]->SetMarkerSize(0.4);
     
     //gRatiovsTime[iE]->Fit(name,"0WQ");
     
     gRatiovsTime[iE]->Draw("AP");
     hRatiovsTime[iE]->SetLineColor(2); hRatiovsTime[iE]->SetMarkerColor(1);
     hRatiovsTime[iE]->Draw("CPSAME");
     
     l1->Draw("SAME");
     
  }
    //Save Correction Factors
    sprintf(name,"TimeCorrectionFactors_%s_v2.0.root",Dataclass.c_str());
    TFile *fout = new TFile(name,"RECREATE");
    cRatiovsTime->Write();
    fout->Close();
    printf("Correction factors saved in file %s\n",name);
}


