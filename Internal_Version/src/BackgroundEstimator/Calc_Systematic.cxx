//Author: Vlasios Vasileiou <vlasisva@gmail.com>

#include "BackgroundEstimator/BackgroundEstimator.h"
#include <algorithm>
#include "TProfile.h"
#include "TLine.h"
#include "TGraph.h"
#include <vector>
#include "TGraphErrors.h"
#include "TStyle.h"

void BKGE_NS::Calc_Systematic(vector<string> GRB_folders, vector <double> METs, string Dataclass, double MinE, double MaxE, int NBins, int CalcSystematic){
 gStyle->SetOptTitle(1);
 gStyle->SetOptStat(1);
 char name[1000];

 const double lMinE = log10(MinE);
 const double lMaxE = log10(MaxE);
 const double dlE = (lMaxE-lMinE)/NBins;


 vector <double> Sig[NBins];//double because of TMath::QUantiles
 vector <float> Back[NBins];
 vector <float> ratio[NBins];
 

 for (unsigned int iTest=0;iTest<METs.size();iTest++) {
 
      sprintf(name,"%s/%s_bkg_%.0f_%.0f.root",GRB_folders[iTest].c_str(),Dataclass.c_str(),MinE,MaxE);
      TFile * fEst = TFile::Open(name);
      if (!fEst) { printf("no %s\n",name);  continue;    }
      TH1F * hEst = (TH1F*)fEst->Get("hCtsvsEnergy_Est");
      if (!hEst){printf("no hBurst in %s\n",name);  fEst->Close(); continue;}

      sprintf(name,"%s/%s_sig_%.0f_%.0f.root",GRB_folders[iTest].c_str(),Dataclass.c_str(),MinE,MaxE);
      TFile * fBurst = TFile::Open(name);
      if (!fBurst) { printf("no %s\n",name);     continue;    }
      TH1F * hBurst = (TH1F*)fBurst->Get("hCtsvsEnergy_Burst");
      if (!hBurst){printf("no hBurst in %s\n",name);  fBurst->Close(); continue;}
     
      for (int iE=0;iE<NBins;iE++) {
/*
          if (iE==14) {
              double aback=hEst->Integral();
              double asig=hBurst->Integral(); 
                    
              if (aback<10 || asig==0) continue; //don't fill faint cases
              int ip=Back[iE].size();
              Back[iE].push_back(aback);
              Sig[iE].push_back(asig);
          
          }
  */
          double aback=hEst->GetBinContent(iE+1);
          double asig=hBurst->GetBinContent(iE+1);
                    
          //if (aback<30 || asig==0) continue; //don't fill faint cases
          if (aback<10 || asig==0) continue; //don't fill faint cases
          //int ip=Back[iE].size();
          Back[iE].push_back(aback);
          Sig[iE].push_back(asig);
      }
      fBurst->Close();
      fEst->Close();
 
 }

 TGraph  * gErrorvsE = new TGraph();
 TGraphErrors *gBiasvsE = new TGraphErrors();
 TGraphErrors *gRMSvsE = new TGraphErrors();

 for (int iE=0;iE<NBins;iE++) {
// for (int iE=14;iE<=14;iE++) {
/*
   for (int iE=0;iE<NBins;iE++) {
     if (ratio[iE].size()==0) continue;
     float med=TMath::Median(ratio[iE].size(),&ratio[iE].front());
     
       //printf("%d unbias back - divide by %f\n",iE,AveR[iE]);
       //unbias back and sig
     for (int is=0;is<Back[iE].size();is++) {
          //Back[iE][is]/=AveR[iE];
          Back[iE][is]/=med;
     }
   }  
 */


     int ip=gErrorvsE->GetN();
     float anerror=0;
     float anE=pow(float(10),lMinE+(ip+0.5)*dlE);
     const int nEntries=Back[iE].size();   
 
     if (nEntries<20){ printf("not enough entries %d\n",nEntries); continue;}

     anerror=DecomposeError(nEntries,&Back[iE].front(),&Sig[iE].front(),true);
     gErrorvsE->SetPoint(ip,anE,anerror); //syst. error
   
     //float mean_ratio=TMath::Mean(nEntries,&ratio[iE].front())-1;
     //float rms       =TMath::RMS(nEntries,&ratio[iE].front());
     //float mean_ratio_err=rms/sqrt(ratio[iE].size());
     //float median_sig=TMath::Median(nEntries,&Sig[iE].front());

     //double sig_q,q=0.50;   
     //TMath::Quantiles(nEntries,1,&Sig[iE].front(),&sig_q,&q,false);
   
   //printf("%d rms=%f stat_exp=%f aee=%f (med_sig=%f)\n",ip,rms,1/sqrt(sig_q),aee,sig_q);
   
     //gRMSvsE->SetPoint(ip, anE,rms);  //width of ratio dist   
     //gBiasvsE->SetPoint(ip, anE,mean_ratio); //mean of ratio dist
     //gBiasvsE->SetPointError(ip, 0,mean_ratio_err); //mean of ratio dist
     
   }
/*
   if (CalcSystematic>1) {
 
   float Back_all_E[METs.size()];
   double Sig_all_E[METs.size()];
   vector <float> ratio_all_E;
   for (int is=0;is<METs.size();is++) {
      Back_all_E[is]=Sig_all_E[is]=0;
      for (int iE=0;iE<NBins;iE++) {
         Back_all_E[is]+=Back[iE][is];
         Sig_all_E[is] +=Sig[iE][is];
      }
      if (Back_all_E[is]>400) ratio_all_E.push_back(Back_all_E[is]/Sig_all_E[is]);
   }  
    
   int nn=ratio_all_E.size();
   //printf("%f  %f %f\n",TMath::Mean(nn,&ratio_all_E.front()),TMath::RMS(nn,&ratio_all_E.front()),TMath::Median(nn,&ratio_all_E.front()));
   DecomposeError(METs.size(),Back_all_E,Sig_all_E,true);
 }
 */
 

 
   gErrorvsE->GetYaxis()->SetTitle("Syst. Error (\%)");
   gErrorvsE->GetXaxis()->SetTitle("Energy (MeV)");
    
   gBiasvsE->GetYaxis()->SetTitle("Bias (\%)");
   gBiasvsE->GetXaxis()->SetTitle("Energy (MeV)");
   gRMSvsE->GetYaxis()->SetTitle("RMS of ratio distribution");
   gRMSvsE->GetXaxis()->SetTitle("Energy (MeV)");

  
   TCanvas * cFinal = new TCanvas("cFinal","cFinal");
   cFinal->Divide(2,2);
 
   cFinal->cd(1); 
   cFinal->GetPad(1)->SetGridy();
   cFinal->GetPad(1)->SetLogx();
   gErrorvsE->GetYaxis()->SetRangeUser(0,0.4);
   gErrorvsE->Draw("A*");
     

   cFinal->cd(2); 
   cFinal->GetPad(2)->SetGridy();
   cFinal->GetPad(2)->SetLogx();
   gRMSvsE->GetYaxis()->SetRangeUser(0,0.4);
   gRMSvsE->Draw("A*");
    
   cFinal->cd(3); 
   cFinal->GetPad(3)->SetGridy();
   cFinal->GetPad(3)->SetLogx();
   gBiasvsE->GetYaxis()->SetRangeUser(-0.1,0.1);
   gBiasvsE->Draw("A*");
 


}


