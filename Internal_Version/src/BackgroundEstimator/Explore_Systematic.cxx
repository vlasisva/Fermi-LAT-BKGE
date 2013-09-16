//Author: Vlasios Vasileiou <vlasisva@gmail.com>
//$Header: Exp $

#include "BackgroundEstimator/BackgroundEstimator.h"
#include <algorithm>
#include "TProfile.h"
#include "TLine.h"
#include "TGraph.h"
#include <vector>
#include "TGraphErrors.h"
#include "TStyle.h"
#include "TMarker.h"
#include "TColor.h"

void BKGE_NS::Explore_Systematic(vector<string> GRB_folders, vector <double> METs, vector <double> GRB_L, vector <double> GRB_B, vector <double> GRB_Theta, vector <double> GRB_ZTheta, vector <double> GRB_ZPhi, vector <double> GRB_McIlwainL, vector <double> GRB_ROCK, string Dataclass, double MinE, double MaxE, int NBins){
 gStyle->SetOptTitle(0);
 gStyle->SetOptStat(1);
  const int prof_bins=10;
 char name[1000],name2[1000];

 const double MinMET=*min_element(METs.begin(),METs.end());
 const double MaxMET=*max_element(METs.begin(),METs.end());

 const double lMinE = log10(MinE);
 const double lMaxE = log10(MaxE);
 const double dlE = (lMaxE-lMinE)/NBins;


 vector <double> Sig[NBins];//double because of TMath::QUantiles
 vector <double> Back[NBins];
 vector <float> Factor[NBins];
 vector <int> idx[NBins];
 

 const UInt_t Number = 3;
 Double_t Red[Number]    = { 1.00, 0.00, 0.00};
 Double_t Green[Number]  = { 0.00, 1.00, 0.00};
 Double_t Blue[Number]   = { 0.00, 0.00, 1.00};
 Double_t Length[Number] = { 0.00, 0.50, 1.00 };
 const Int_t palette_colors=256;
 const int palette_idx=TColor::CreateGradientColorTable(Number,Length,Blue,Green,Red,palette_colors);
 

 TH1F * hFactor[NBins];
 for (int iE=0;iE<15;iE++) {
    sprintf(name,"hFactor_%d",iE);
    sprintf(name2,"Bin #%d %.0f<E<%.0f MeV",iE,pow(10,lMinE+(iE)*dlE),pow(10,lMinE+(iE+1)*dlE));
    hFactor[iE] = new TH1F(name,name2,50,-2,2);
    hFactor[iE]->GetXaxis()->SetTitle("(N_{est}-N_{act})/N_{est}");
 }

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
      if (GRB_L[iTest]>180) GRB_L[iTest]-=360;

      for (int iE=0;iE<NBins;iE++) {
         if (iE==14) {
              float aback=hEst->Integral();
              float asig=hBurst->Integral(); 
              float aFactor=(aback-asig)/aback;
              if (fabs(aFactor>2) || aback<=0){
                    printf("too large aggregate Factor ra=%f b=%e s=%e\n",aFactor,aback,asig);
                    printf("%s\n",name);
                    continue;
              }      
              //if (aback<1000) continue;
              //if (asig==0) continue; //don't fill faint cases
              //printf("%f %f\n",aback,asig);
            
              Back[iE].push_back(aback);              Sig[iE].push_back(asig);
              idx[iE].push_back(iTest);
              Factor[iE].push_back(aFactor);
              hFactor[iE]->Fill(aFactor);
              //if ((aback>1800 && (aback-asig)/aback<-0.1)) printf("%s\n",name);
              continue;
         }
          float aback=hEst->GetBinContent(iE+1);
          float asig=hBurst->GetBinContent(iE+1);
          if (aback<5) continue;        

          float aFactor=(aback-asig)/aback;
          //if ((aback+asig)<10 || asig==0) continue; //don't fill faint cases
          //if ((aback+asig)<10 || asig==0) continue; //don't fill faint cases
          //if (iE==0 && aback/asig<1.4) break;
          //printf("%s\n",name);
          Back[iE].push_back(aback);
          Sig[iE].push_back(asig);
          idx[iE].push_back(iTest);
          Factor[iE].push_back(aFactor);
          hFactor[iE]->Fill(aFactor);
          //Factor[iE].push_back(aback/asig);
      }
      fBurst->Close();
      fEst->Close();
 }
 
 //Done reading data
//////////////////////////////////////////////////////////////////////////////////////////////
 gStyle->SetLineScalePS(0.01);
  gStyle->SetPadTopMargin(0.07);
  gStyle->SetPadLeftMargin(0.13);
  gStyle->SetPadRightMargin(0.011);
  gStyle->SetPadBottomMargin(0.15);
 
  TCanvas * cFactor = new TCanvas("cFactor","Factor",1400,400);
  cFactor->Divide(5,3,0.001,0.001);
 
  TCanvas * cFactorvsTheta = new TCanvas("cFactorVsTheta","Factor vs Theta",1400,400);
  cFactorvsTheta->Divide(5,3,0.001,0.001);

  TCanvas * cFactorvsROCK = new TCanvas("cFactorVsROCK","Factor vs RockingAngle",1400,400);
  cFactorvsROCK->Divide(5,3,0.001,0.001);

  
  TCanvas * cFactorvsZTheta = new TCanvas("cFactorVsZTheta","Factor vs ZTheta",1400,400);
  cFactorvsZTheta->Divide(5,3,0.001,0.001);
  
  TCanvas * cFactorvsZPhi = new TCanvas("cFactorVsZPhi","Factor vs ZPhi",1400,400);
  cFactorvsZPhi->Divide(5,3,0.0001,0.0001);

  TCanvas * cFactorvsMcIlwainL = new TCanvas("cFactorVsMcIlwainL"," Factor vs McIlwainL",1400,400);
  cFactorvsMcIlwainL->Divide(5,3,0.0001,0.0001);

  TCanvas * cDiffcsAct = new TCanvas("cDiffcsAct","Factor vs N_est",1400,400);
  cDiffcsAct->Divide(5,3,0.0001,0.0001);

  TCanvas * cDiffcsAct2 = new TCanvas("cDiffcsAct2","Factor vs L/B",1400,400);
  cDiffcsAct2->Divide(5,3,0.0001,0.0001);

  TCanvas * cErrorvsEst = new TCanvas("cErrorvsEst","Factor vs N_est",1400,400);
  cErrorvsEst->Divide(5,3,0.0001,0.0001);

  TCanvas * cSigvsBack = new TCanvas("cSigvsBack","Signal vs Background",1400,400);
  cSigvsBack->Divide(5,3,0.0001,0.0001);

  TCanvas * cFactorvsTimeMap = new TCanvas("cFactorvsTimeMap","Factor vs Time",1400,400);
  cFactorvsTimeMap->Divide(5,3,0.0001,0.0001);

  TCanvas * cFactorvsBkgAndTheta = new TCanvas("cFactorvsBkgAndTheta","Factor vs BkgAndTheta",1400,400);
  cFactorvsBkgAndTheta ->Divide(5,3,0.0001,0.0001);


  TLine * l1;
  
  TGraph *   gFactorvsEst[NBins];



 TH2F * hFactor2vsE[NBins],*hFactorvsSig[NBins],*hFactor2vsE_Norm[NBins],*hSigvsBack[NBins];
 TH2F * hFactorvsTimeMap[NBins];
 TGraph * gFactorvsTheta[NBins],*gFactorvsZTheta[NBins],*gFactorvsZPhi[NBins],*gFactorvsMcIlwainL[NBins];
 TGraph * gFactorvsROCK[NBins]; 
  
 for (int iE=0;iE<15;iE++) {

   cFactor->cd(iE+1); hFactor[iE]->Draw();

    sprintf(name,"FactorvsTheta_%d",iE);
    sprintf(name2,"Bin #%d %.0f<E<%.0f MeV",iE,pow(10,lMinE+(iE)*dlE),pow(10,lMinE+(iE+1)*dlE));
    if (iE==14) sprintf(name2,"All energy bins");
    gFactorvsTheta[iE] = new TGraph();    gFactorvsTheta[iE]->SetTitle(name2);

    sprintf(name,"FactorvsMcIlwainL_%d",iE);
    gFactorvsMcIlwainL[iE] = new TGraph();    gFactorvsMcIlwainL[iE]->SetTitle(name2);

    sprintf(name,"FactorvsRock_%d",iE);
    gFactorvsROCK[iE] = new TGraph();    gFactorvsROCK[iE]->SetTitle(name2);


    sprintf(name,"FactorvsZTheta_%d",iE);
    gFactorvsZTheta[iE] = new TGraph();  gFactorvsZTheta[iE]->SetTitle(name2);

    sprintf(name,"FactorvsZPhi_%d",iE);
    gFactorvsZPhi[iE] = new TGraph();  gFactorvsZPhi[iE]->SetTitle(name2);

    sprintf(name,"Factor2vsE_%d",iE);
    hFactor2vsE[iE]= new TH2F(name,name2,20,-180,180,20,-90,90);
    hFactor2vsE[iE]->GetXaxis()->SetTitle("L (deg)");
    hFactor2vsE[iE]->GetYaxis()->SetTitle("B (deg)");
    hFactor2vsE[iE]->SetContour(200);

    sprintf(name,"Factor2vsE_Norm__%d",iE);
    hFactor2vsE_Norm[iE]= (TH2F*)hFactor2vsE[iE]->Clone(name);

    sprintf(name,"FactorvsBack_%d",iE);
    hFactorvsSig[iE]= new TH2F(name,name2,50,0,log10(5000),30,-50,50);
    hFactorvsSig[iE]->GetXaxis()->SetTitle("log_{10}(Signal)");
    hFactorvsSig[iE]->GetYaxis()->SetTitle("Factor Est/Real");
    hFactorvsSig[iE]->SetContour(200);

    sprintf(name,"SigvsBack_%d",iE);
    hSigvsBack[iE]= new TH2F(name,name2,30,0,log10(5000),30,0,log10(5000));
    hSigvsBack[iE]->GetXaxis()->SetTitle("log_{10}(Measured Signal)");
    hSigvsBack[iE]->GetYaxis()->SetTitle("log_{10}(Estimated Background)");
    hSigvsBack[iE]->SetContour(200);

    sprintf(name,"FactorvsTimeMap_%d",iE);
    hFactorvsTimeMap[iE]= new TH2F(name,name2,45,MinMET,MaxMET,50,-0.5,0.5);

    sprintf(name,"DiffDist_%d",iE);
    gFactorvsEst[iE]= new TGraph();
    gFactorvsEst[iE]->SetName(name);
    gFactorvsEst[iE]->SetTitle(name2);

 }


   
  //Fill some data   

  for (int iE=0;iE<15;iE++) {
     const unsigned int nEntries=Back[iE].size();
     if (!nEntries) {printf("cont %d\n",iE); continue;}
     cFactorvsBkgAndTheta->cd(iE+1);
     TH1F * h = (TH1F*)cFactorvsBkgAndTheta->DrawFrame(0,0,TMath::MaxElement(nEntries,&Back[iE].front()),100);
     h->GetXaxis()->SetTitle("Estimated background");
     h->GetYaxis()->SetTitle("LAT Theta (deg)");
     
     for (unsigned int iPass=0;iPass<nEntries;iPass++) {
          
          //if (iE==6 && Factor[iE][iPass]>2) printf("---------%s\n",name);
          hFactor2vsE[iE]->Fill(GRB_L[idx[iE][iPass]],GRB_B[idx[iE][iPass]],Factor[iE][iPass]);
          hFactor2vsE_Norm[iE]->Fill(GRB_L[idx[iE][iPass]],GRB_B[idx[iE][iPass]]);
         
          hFactorvsSig[iE]->Fill(log10(Back[iE][iPass]),Factor[iE][iPass]);
          hSigvsBack[iE]  ->Fill(log10(Sig[iE][iPass]),log10(Back[iE][iPass]));
          hFactorvsTimeMap[iE]->Fill(METs[idx[iE][iPass]],Factor[iE][iPass]);
          gFactorvsTheta[iE]->SetPoint(iPass,GRB_Theta[idx[iE][iPass]],Factor[iE][iPass]);
          gFactorvsMcIlwainL[iE]->SetPoint(iPass,GRB_McIlwainL[idx[iE][iPass]],Factor[iE][iPass]);
          gFactorvsZTheta[iE]->SetPoint(iPass,GRB_ZTheta[idx[iE][iPass]],Factor[iE][iPass]);
          gFactorvsZPhi[iE]->SetPoint(iPass,GRB_ZPhi[idx[iE][iPass]],Factor[iE][iPass]);      
          gFactorvsEst[iE]->SetPoint(iPass,Back[iE][iPass],Factor[iE][iPass]);
          gFactorvsROCK[iE]->SetPoint(iPass,GRB_ROCK[idx[iE][iPass]],Factor[iE][iPass]);
          //gFactorvsEst[iE]->SetPoint(iPass,Sig[iE][iPass],Factor[iE][iPass]);
          
          
          TMarker * mkr = new TMarker(Back[iE][iPass],GRB_Theta[idx[iE][iPass]],20);
          int color_idx=int(palette_colors*(Factor[iE][iPass]+1)/2.0); //-1.0--1.0
          if (color_idx<0) color_idx=0;
          else if (color_idx>=palette_colors) color_idx=palette_colors-1;
          mkr->SetMarkerSize(0.5);
          mkr->SetMarkerColor(palette_idx+color_idx);
          mkr->Draw();
     }
  
  
       double cl[2]={0.02,0.98};
       double q[2];
       //TMath::Quantiles(nEntries,2,&Sig[iE].front(),q,cl,false);
       TMath::Quantiles(nEntries,2,&Back[iE].front(),q,cl,false);
       TProfile * hp = new TProfile(name,name2,prof_bins,q[0],q[1]);
       
       cDiffcsAct->cd(iE+1);
   
       for (int id=0;id<nEntries;id++)  hp->Fill(Back[iE][id],Factor[iE][id]);
       //for (int id=0;id<nEntries;id++)  hp->Fill(Sig[iE][id],Factor[iE][id]);
       TGraphErrors * gp = new TGraphErrors(); gp->SetMarkerColor((2)); gp->SetLineColor(2);
       for (int ib=0;ib<prof_bins;ib++) {
          float amean=hp->GetBinContent(ib+1);   
          gp->SetPoint(ib,hp->GetXaxis()->GetBinCenter(ib+1),amean);
          float aerror=hp->GetBinError(ib+1);
          gp->SetPointError(ib,hp->GetXaxis()->GetBinWidth(ib+1)/2.,aerror);
          //printf("%d %f %f\n",ib,amean,aerror);
       }
       
       gFactorvsEst[iE]->SetMarkerColor(kGray);
       gFactorvsEst[iE]->Draw("AP");   gFactorvsEst[iE]->GetYaxis()->SetRangeUser(-0.5,0.5);
       gp->Draw("PSAME");
       gFactorvsEst[iE]->GetYaxis()->SetTitle("(N_{est}-N_{act})/N_{est}");
       gFactorvsEst[iE]->GetXaxis()->SetTitle("N_{est}");
       
     
     cDiffcsAct2->cd(iE+1);
     hFactor2vsE[iE]->Divide(hFactor2vsE_Norm[iE]);
     hFactor2vsE[iE]->GetZaxis()->SetRangeUser(-1,1);
     hFactor2vsE[iE]->Draw("COLZ");

     
     cSigvsBack->cd(iE+1);
     l1 = new TLine(hSigvsBack[iE]->GetXaxis()->GetXmin(),hSigvsBack[iE]->GetYaxis()->GetXmin(),
                           hSigvsBack[iE]->GetXaxis()->GetXmax(),hSigvsBack[iE]->GetYaxis()->GetXmax());
     hSigvsBack[iE]->Draw("COLZ");
     l1->Draw("SAME");

     cErrorvsEst->cd(iE+1);

     l1 = new TLine(hFactorvsSig[iE]->GetXaxis()->GetXmin(),0,hFactorvsSig[iE]->GetXaxis()->GetXmax(),0);     
     hFactorvsSig[iE]->Draw("COLZ");
     l1->Draw("SAME");


     cFactorvsTimeMap->cd(iE+1);
     hFactorvsTimeMap[iE]->SetContour(64);
     hFactorvsTimeMap[iE]->Draw("COLZ");
 
     cFactorvsTheta->cd(iE+1); 
     sprintf(name,"hp_Factorvstheta_%d",iE);

     hp = new TProfile(name,name,prof_bins,0,80);
     gp = new TGraphErrors();gp->SetMarkerColor((2)); gp->SetLineColor(2);
     
     sprintf(name,"hm_Factorvstheta_%d",iE);
     
     for (int ii=0;ii<nEntries;ii++) hp->Fill(GRB_Theta[idx[iE][ii]],Factor[iE][ii]);     
     
     for (int ib=0;ib<prof_bins;ib++) {
        float amean=hp->GetBinContent(ib+1);   
        gp->SetPoint(ib,hp->GetXaxis()->GetBinCenter(ib+1),amean);
        float aerror=hp->GetBinError(ib+1);
        gp->SetPointError(ib,hp->GetXaxis()->GetBinWidth(ib+1)/2.,aerror);
        //gp->SetPointError(ib,hp->GetXaxis()->GetBinWidth(ib+1)/2.,aerror);
     }
     
        
     gFactorvsTheta[iE]->SetMarkerStyle(6);
     gFactorvsTheta[iE]->SetMarkerColor(kGray);
     gFactorvsTheta[iE]->Draw("AP"); gFactorvsTheta[iE]->GetYaxis()->SetRangeUser(-0.5,0.5);
     gFactorvsTheta[iE]->GetXaxis()->SetRangeUser(0,80); 
     gFactorvsTheta[iE]->GetYaxis()->SetTitle("(N_{est}-N_{act})/N_{est}");
     gFactorvsTheta[iE]->GetXaxis()->SetTitle("#theta_{LAT} (deg)");
     
     TLine * ltheta = new TLine(0,0,80,0); ltheta->Draw("SAME");
     gp->Draw("PSAME");


//////////
     cFactorvsMcIlwainL->cd(iE+1); 
     sprintf(name,"hp_FactorvsMcIlwainL_%d",iE);

     hp = new TProfile(name,name,prof_bins,0.97,1.8);
     gp = new TGraphErrors();gp->SetMarkerColor((2)); gp->SetLineColor(2);
     
     sprintf(name,"hm_FactorvsMcIlwainL_%d",iE);
     
     for (int ii=0;ii<nEntries;ii++) hp->Fill(GRB_McIlwainL[idx[iE][ii]],Factor[iE][ii]);     
     
     for (int ib=0;ib<prof_bins;ib++) {
        float amean=hp->GetBinContent(ib+1);   
        gp->SetPoint(ib,hp->GetXaxis()->GetBinCenter(ib+1),amean);
        float aerror=hp->GetBinError(ib+1);
        gp->SetPointError(ib,hp->GetXaxis()->GetBinWidth(ib+1)/2.,aerror);
        //gp->SetPointError(ib,hp->GetXaxis()->GetBinWidth(ib+1)/2.,aerror);
     }

     gFactorvsMcIlwainL[iE]->SetMarkerColor(kGray);
     gFactorvsMcIlwainL[iE]->SetMarkerStyle(6);
     gFactorvsMcIlwainL[iE]->Draw("AP"); gFactorvsMcIlwainL[iE]->GetYaxis()->SetRangeUser(-0.5,0.5);
     gFactorvsMcIlwainL[iE]->GetXaxis()->SetRangeUser(0.97,1.8);
     gFactorvsMcIlwainL[iE]->GetXaxis()->SetRangeUser(0.9,1.8); 
     gFactorvsMcIlwainL[iE]->GetYaxis()->SetTitle("(N_{est}-N_{act})/N_{est}");
     gFactorvsMcIlwainL[iE]->GetXaxis()->SetTitle("McIlwainL");
     
     TLine * lMcIlwainL = new TLine(0.9,0,1.80,0); lMcIlwainL->Draw("SAME");
     gp->Draw("PSAME");



/////


     
     
     cFactorvsZTheta->cd(iE+1); 
     sprintf(name,"hp_Factorvsztheta_%d",iE);
     hp = new TProfile(name,name,prof_bins,0,110);
     gp = new TGraphErrors();gp->SetMarkerColor((2)); gp->SetLineColor(2);
     for (unsigned int i=0;i<nEntries;i++)hp->Fill(GRB_ZTheta[idx[iE][i]],Factor[iE][i]);
     
     
     for (int ib=0;ib<prof_bins;ib++) {
        float amean=hp->GetBinContent(ib+1);   
        gp->SetPoint(ib,hp->GetXaxis()->GetBinCenter(ib+1),amean);
        float aerror=hp->GetBinError(ib+1);
        gp->SetPointError(ib,hp->GetXaxis()->GetBinWidth(ib+1)/2.,aerror);
     }     
     
     gFactorvsZTheta[iE]->SetMarkerColor(kGray);
     gFactorvsZTheta[iE]->SetMarkerStyle(6);
     gFactorvsZTheta[iE]->Draw("AP");     gFactorvsZTheta[iE]->GetYaxis()->SetRangeUser(-0.5,0.5);
     gFactorvsZTheta[iE]->GetXaxis()->SetRangeUser(0,110);
     gFactorvsZTheta[iE]->GetYaxis()->SetTitle("(N_{est}-N_{act})/N_{est}");
     gFactorvsZTheta[iE]->GetXaxis()->SetTitle("#theta_{#oplus} (deg)");
     gp->Draw("PSAME");
     TLine * lztheta = new TLine(0,0,110,0); lztheta->Draw("SAME");



     cFactorvsZPhi->cd(iE+1); 
     sprintf(name,"hp_FactorvsZPhi_%d",iE);
     hp = new TProfile(name,name,prof_bins,0,360);
     gp = new TGraphErrors();gp->SetMarkerColor(2); gp->SetLineColor(2);
     for (int ii=0;ii<nEntries;ii++)hp->Fill(GRB_ZPhi[idx[iE][ii]],Factor[iE][ii]);
      
     
     for (int ib=0;ib<prof_bins;ib++) {
        float amean=hp->GetBinContent(ib+1);   
        gp->SetPoint(ib,hp->GetXaxis()->GetBinCenter(ib+1),amean);
        float aerror=hp->GetBinError(ib+1);
        gp->SetPointError(ib,hp->GetXaxis()->GetBinWidth(ib+1)/2.,aerror);
     }     

     gFactorvsZPhi[iE]->SetMarkerColor(kGray);
     gFactorvsZPhi[iE]->SetMarkerStyle(6);
     gFactorvsZPhi[iE]->Draw("AP");   gFactorvsZPhi[iE]->GetYaxis()->SetRangeUser(-0.5,0.5);
      gFactorvsZPhi[iE]->GetXaxis()->SetRangeUser(0,360);
     gFactorvsZPhi[iE]->GetYaxis()->SetTitle("(N_{est}-N_{act})/N_{est}");
     gFactorvsZPhi[iE]->GetXaxis()->SetTitle("#phi_{#oplus} (deg)");
     gp->Draw("PSAME");
     TLine * lZPhi = new TLine(0,0,360,0); lZPhi->Draw("SAME");
////////////////////


     cFactorvsROCK->cd(iE+1); 
     sprintf(name,"hp_FactorvsROCK_%d",iE);
     hp = new TProfile(name,name,prof_bins,0,180);
     gp = new TGraphErrors();gp->SetMarkerColor(2); gp->SetLineColor(2);
     for (int ii=0;ii<nEntries;ii++)hp->Fill(GRB_ROCK[idx[iE][ii]],Factor[iE][ii]);
      
     for (int ib=0;ib<prof_bins;ib++) {
        float amean=hp->GetBinContent(ib+1);   
        gp->SetPoint(ib,hp->GetXaxis()->GetBinCenter(ib+1),amean);
        float aerror=hp->GetBinError(ib+1);
        gp->SetPointError(ib,hp->GetXaxis()->GetBinWidth(ib+1)/2.,aerror);
     }     


     gFactorvsROCK[iE]->SetMarkerColor(kGray);
     gFactorvsROCK[iE]->SetMarkerStyle(6);
     gFactorvsROCK[iE]->Draw("AP");   gFactorvsROCK[iE]->GetYaxis()->SetRangeUser(-0.5,0.5);
      gFactorvsROCK[iE]->GetXaxis()->SetRangeUser(0,110);
     gFactorvsROCK[iE]->GetYaxis()->SetTitle("(N_{est}-N_{act})/N_{est}");
     gFactorvsROCK[iE]->GetXaxis()->SetTitle("#theta_{rock} (deg)");
     gp->Draw("PSAME");
     TLine * lROCK = new TLine(0,0,110,0); lROCK->Draw("SAME");

  }


  cFactorvsTheta->Update();cFactorvsTheta->GetPad(15)->SaveAs("Factor_vs_thetalat_100_15deg.eps");
  cFactorvsZTheta->Update();cFactorvsZTheta->GetPad(15)->SaveAs("Factor_vs_ztheta_100_15deg.eps");  
  cFactorvsZPhi->Update();cFactorvsZPhi->GetPad(15)->SaveAs("Factor_vs_zphi_100_15deg.eps");  
   cFactorvsMcIlwainL->Update(); cFactorvsMcIlwainL->GetPad(15)->SaveAs("Factor_vs_mcilwainl_100_15deg.eps");  
   cFactorvsROCK->Update(); cFactorvsROCK->GetPad(15)->SaveAs("Factor_vs_ROCK_100_15deg.eps");  
   cDiffcsAct->Update(); 
   cDiffcsAct2->Update();
   cErrorvsEst->Update();
   cSigvsBack->Update(); 
   cFactorvsTimeMap->Update(); 
   cFactorvsBkgAndTheta->Update();


}


