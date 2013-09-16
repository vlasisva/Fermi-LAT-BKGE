// Author: Vlasios Vasileiou <vlasisva@gmail.com>
// $Header: /nfs/slac/g/glast/ground/cvs/GRBAnalysis-scons/BackgroundEstimator/src/BKGE_Tools/FeldmanCousins.cxx,v 1.1.1.1 2011/06/02 19:41:04 jrb Exp $
/* This macro calculates upper limits according to Feldman-Cousins 1998 */
//Author Vlasios Vasileiou 2009

#include "BackgroundEstimator/BKGE_Tools.h"
#include "TGraph.h"
#include "Math/ProbFuncMathCore.h"
#include "Math/PdfFuncMathCore.h"

typedef struct {
  double R;
  int n,mu_best;
  double Pnu_mu,Pnu_mu_best;
} adatum;

bool myfunction (adatum i,adatum j);

double TOOLS::FeldmanCousins(float CL, int SIG, double BKG, int MuMin,int MuMax,float dmu,bool Draw) {

 char name[1000];

 //Intervals for the confidence belt
 TH1F * h0 = new TH1F("h0","Limits",MuMax-MuMin,MuMin,MuMax);
 TH1F * h1 = new TH1F("h1","Limits",MuMax-MuMin,MuMin,MuMax);
 h1->GetXaxis()->SetTitle("Measured number of events");
 h1->GetYaxis()->SetTitle("Signal mean #mu");

 TCanvas * can = NULL;
 TFile * f = NULL;
 TNamed * DataNamed = NULL;
 if (Draw) {
   can = new TCanvas("cFC","cFC");
   can->Divide(1,2);
   DataNamed = new TNamed("Results","");
   f = new TFile("/tmp/UL.root","RECREATE");
 }

 const int i_max = (int)((MuMax-MuMin)/dmu);
 TGraph g[i_max];

 float mu=MuMin;

 //create the confidence belt
 for (int i=0;i<i_max;i++) {
    //printf("n\tPnu_mu\tmu_best\tP_nu_mubest\tR\n");
    vector <adatum> Data;
    adatum datum;
    for (datum.n=0;datum.n<i_max;datum.n++) {
      datum.Pnu_mu = ROOT::Math::poisson_pdf(datum.n,BKG+mu);
      datum.mu_best = int(datum.n-BKG);
      if ( datum.mu_best<0) datum.mu_best=0;
      datum.Pnu_mu_best = ROOT::Math::poisson_pdf(datum.n,BKG+datum.mu_best);
      datum.R=datum.Pnu_mu/datum.Pnu_mu_best;
      Data.push_back(datum);
      //printf("%d\t%.3f\t%d\t%.3f\t\t%.2e\n",datum.n,datum.Pnu_mu,datum.mu_best,datum.Pnu_mu_best,datum.R);
    }
    std::sort(Data.begin(),Data.end(),myfunction);

    double CLsum=0;
    double x[2]={9999,0},y[2]={mu,mu};
    for (int j=0;j<i_max;j++) {
       CLsum+=Data[j].Pnu_mu;
       if (x[1]<Data[j].n) x[1]=Data[j].n;
       if (x[0]>Data[j].n) x[0]=Data[j].n;
       if (CLsum>CL) break;
    }
    //printf("%f\n",x[1]);
    g[i] = TGraph(2,x,y);
    mu+=dmu;
 }
 //////////////////////////////////////////////////

 for (int n=0;n<=h0->GetNbinsX();n++) {
    double amu,x0,x1;
    for (int i=0;i<i_max;i++) {
       g[i].GetPoint(1,x1,amu);
       if (x1>n) break;
    }
    h0->SetBinContent(n+1,amu);
    for (int i=0;i<i_max;i++) {
       g[i].GetPoint(0,x0,amu);
       //printf("%d %d %f %f \n",n,i,x0,amu);
       if (x0>n) break;
    }
    h1->SetBinContent(n+1,amu);
 }


 if (Draw) {
   can->cd(2);
   h1->Draw();
   h0->Draw("SAME");

   //Belt
   double x[2]={MuMin,MuMax},y[2]={MuMin,MuMax};
   can->cd(1);

   for (int i=0;i<i_max;i++) {
      if (i==0) {    
        g[i].SetTitle("Confidence Belt");
        g[i].GetXaxis()->SetTitle("Measured number of events");
        g[i].GetYaxis()->SetTitle("Signal #mu");
     //  g[i].GetXaxis()->SetLimits(MuMin,MuMax);
        g[i].GetYaxis()->SetRangeUser(MuMin,MuMax);
        g[i].Draw("AC");
     }
     else g[i].Draw("CSAME");
   } 
 

   x[0]=x[1]=SIG;
   y[0]=MuMin;y[1]=MuMax;
   TGraph gMeas = TGraph(2,x,y);
   gMeas.SetLineColor(2);
   gMeas.SetLineWidth(2);
   gMeas.Draw("CSAME");
   can->cd(2); gMeas.Draw("CSAME");

   can->Write();
 
   //printf("%s\n",name);
   DataNamed->SetTitle(name);
   DataNamed->Write();
   delete can;
   f->Close();
 }

 double UL=h1->GetBinContent(SIG+1);
 delete h0;
 delete h1;
 return UL;
}


bool myfunction (adatum i,adatum j) {
 return (i.R>j.R);
}



