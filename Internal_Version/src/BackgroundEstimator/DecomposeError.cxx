#include "BackgroundEstimator/BackgroundEstimator.h"
#include "TLegend.h"
#include "TGraph.h"
#include "Math/ProbFuncMathCore.h"
#include "Math/PdfFuncMathCore.h"
#include "Math/QuantFuncMathCore.h"
#include <algorithm>
#include "TLine.h"

double WeightedP(int nev, float bkg, float dbkg);


float  BKGE_NS::DecomposeError(int ent0, float Back0[], double MeasSig0[], bool Plot){

 
 const int bins=100;
 const int binss=100;

 vector <double> MeasSig;
 vector <double> Back;

 //select for analysis only cases that the expected bkg >5 events
 for (int i=0;i<ent0;i++) {
    if (Back0[i]>20){
       MeasSig.push_back(MeasSig0[i]);
       Back.push_back(Back0[i]); 
    }
 }
 unsigned int nData=Back.size();

 if (nData<10) {printf ("%s: not enough data points.. returning\n",__FUNCTION__); return 0;}

 const int MaxSig=TMath::MaxElement(nData,&MeasSig.front());
 const int MinSig=TMath::MinElement(nData,&MeasSig.front());
 printf("nData=%d MinSig/MaxSig= %d/%d ,",nData,MinSig,MaxSig);

 TCanvas * c = 0;
 TH1F * hBack=0, * hSig=0, * hSigma=0, *hProb=0;
 TGraph *gQ_fit=0;
 TF1 * fitgaus= 0;
 TGraph * gQ = 0;
 TGraph * gCum=0;
 TGraph * gProbvsBack=0;
 if (Plot) {
   c= new TCanvas("cDecompose","cDecompose",1024,768);
   c->Divide(2,3);

   c->cd(2); 
   hProb = new TH1F("hProb","Probability distribution",binss,0,1);
   hProb->GetXaxis()->SetTitle("Probability"); 
   hProb->SetMinimum(0);
   hProb->Draw();
  
   c->cd(1); 
   hBack = new TH1F("hBack","Background Estimates",bins,MinSig,MaxSig);
   hBack->GetXaxis()->SetTitle("Estimated number of bkg events");
   hSig = new TH1F("hSig","Signal distributions",bins,MinSig,MaxSig);
   hSig->GetXaxis()->SetTitle("Number of signal events");
   for (unsigned int i=0;i<nData;i++) {
       hSig->Fill(MeasSig[i]);
       hBack->Fill(Back[i]);
   }  
  
   c->GetPad(1)->SetLogy();
   hSig->Draw(); 
   hBack->SetLineColor(3);
   hBack->Draw("SAME"); 

   TLegend * l = new TLegend(0.7,0.8,0.9,0.9);
   l->AddEntry(hSig,"Measured signal","l");
   l->AddEntry(hBack,"Background","l");
   l->Draw();
   
   c->cd(4);   
   hSigma = new TH1F("hSigma","Significance distribution",30,-5,5);
   hSigma->GetXaxis()->SetTitle("Significance (sigma)");
   hSigma->GetYaxis()->SetTitle("dN/dP");
      
   c->cd(3);
      
   fitgaus=new TF1("fit","gaus");
   fitgaus->SetParameters(1,0,1);
   fitgaus->FixParameter(1,0);
   fitgaus->FixParameter(2,1);

    c->cd(5);
    gQ = new TGraph();
    gQ->SetTitle("Kolmogorov-test result");
    gQ->GetXaxis()->SetTitle("Test systematic error on the bkg estimate");
    gQ->GetYaxis()->SetTitle("Kolmogorov-test result");
    gQ->Draw("A*");      
    gQ->SetMinimum(0);
 /*
    gQ_fit = new TGraph();
    gQ_fit->SetTitle("Kolmogorov-test result");
    gQ_fit->GetXaxis()->SetTitle("Test systematic error on the bkg estimate");
    gQ_fit->GetYaxis()->SetTitle("Gaussian Fit result");
    gQ_fit->SetMarkerColor(2);
    gQ_fit->Draw("*SAME");      
*/
 }

 short nerrors=0;
 float BestQ=0,BestProb=0;
 const int nData_null=1000;
 vector <double> CorrectedPValues_null(nData_null);
 for (int i=0;i<nData_null;i++) CorrectedPValues_null[i]=gRandom->Uniform();
 sort(CorrectedPValues_null.begin(),CorrectedPValues_null.end()); 
 
 for (float anerror=0;anerror<0.2;anerror+=0.01) { //loop over systematic error on the bkg
   
    if (Plot) { hProb->Reset();hSigma->Reset();}
   
    vector <double> CorrectedPValues(nData);
   
    for (unsigned int i=0;i<nData;i++) { //process entries
        float aBkgError=anerror*Back[i];
       
      
        /*
        int iback=floor(Back[i]+0.5);
        if (MeasSig[i]<=iback) for (int isig=MeasSig[i];isig<iback;isig++)       CorrectedIntP+=WeightedP(isig,Back[i],aBkgError);
        else                   for (int isig=iback+1;isig<=MeasSig[i];isig++)    CorrectedIntP+=WeightedP(isig,Back[i],aBkgError);
        CorrectedIntP*=2;
        CorrectedIntP+=WeightedP(iback,iback,aBkgError);
        */

        if (anerror<1e-4) {
           CorrectedPValues[i]=ROOT::Math::poisson_cdf_c(MeasSig[i],Back[i]);//ROOT::Math::poisson_pdf(MeasSig[i],Back[i])/2.;
          //CorrectedPValues[i]=ROOT::Math::poisson_cdf_c(MeasSig[i],Back[i]);//-ROOT::Math::poisson_pdf(MeasSig[i],Back[i]);
        }   
        else {
          for (int asig=0;asig<MeasSig[i];asig++) CorrectedPValues[i]+=WeightedP(asig,Back[i],aBkgError);
          //CorrectedPValues[i]+=WeightedP(MeasSig[i],Back[i],aBkgError)/2.;    
          //CorrectedPValues[i]+=WeightedP(MeasSig[i],Back[i],aBkgError);
        
          if (CorrectedPValues[i]>1.000001) {
            printf("correctedP>1? p=%f sig=%f bk=%f\n",CorrectedPValues[i],MeasSig[i],Back[i]);
            CorrectedPValues[i]=0.9999;
          }    
          if (CorrectedPValues[i]<=0) {
             printf("correctedP<0?\n");
             //printf("corp=%e MeasSig[i]=%.1f  back=%f iback=%d\n",CorrectedPValues[i],MeasSig[i],Back[i],iback); //check to see if denominator is ~1 as it should be
            //CorrectedPValues[i]=1e-1;     
           }
            
          //we calculated integral from 0 to MeasSig (i.e. the cdf) and the probability we want to use is the cdf_c
          CorrectedPValues[i]=1-CorrectedPValues[i]; 
       }      
       //printf("%f \n",CorrectedPValues[i]/(ROOT::Math::poisson_cdf_c(MeasSig[i],Back[i])+ROOT::Math::poisson_pdf(MeasSig[i],Back[i])/2.)); 
       //CorrectedPValues[i]=ROOT::Math::poisson_cdf_c(MeasSig[i],Back[i])+ROOT::Math::poisson_pdf(MeasSig[i],Back[i])/2.; 
        
        //CorrectedPValues[i]=TOOLS::WeightedP(MeasSig[i],Back[i],aBkgError);
        //printf("%d: sig=%d back=%f bkgerror=%f corrP=%e\n",i,MeasSig[i],Back[i],aBkgError,CorrectedPValues[i]);

        
        if (Plot) {
           hProb->Fill(CorrectedPValues[i]);  
           //if (MeasSig[i]<Back[i]) CorrectedPValues[i]*=-1;
           hSigma->Fill(ROOT::Math::gaussian_quantile_c(CorrectedPValues[i],1));
        }

     }
     
      sort(CorrectedPValues.begin(),CorrectedPValues.end());
      float Q=TMath::KolmogorovTest(nData,&CorrectedPValues.front(),nData_null,&CorrectedPValues_null.front(),"");
      if (Q>BestQ) {BestQ=Q;BestProb=anerror;}
       
  
      if (Plot) {
         hProb->Scale(1./nData,"width");
         c->cd(3); 
         if (gCum) delete gCum;
         gCum = new TGraph();
         for (int ip=1;ip<=hProb->GetNbinsX();ip++) gCum->SetPoint(ip-1,hProb->GetXaxis()->GetBinUpEdge(ip),hProb->Integral(1,ip)/hProb->Integral());
         gCum->SetMarkerStyle(4);
         gCum->SetMarkerSize(0.2);
         
         gCum->GetXaxis()->SetTitle("Probability P'");
         gCum->GetYaxis()->SetTitle("Number of events with P<P'");
         gCum->Draw("AP"); 
         TLine * ll = new TLine(0,0,1,1);
         ll->SetLineColor(2); ll->Draw("SAME");
         
         c->cd(4);
         //c->GetPad(4)->SetLogy();
         hSigma->Fit("fit","BLQ");
        
         gQ->SetPoint(nerrors,anerror,Q);
         //gQ_fit->SetPoint(nerrors,anerror,fitgaus->GetProb());
         for (int ip=1;ip<=6;ip++) c->GetPad(ip)->Modified() ;
         
         if (gProbvsBack) delete gProbvsBack;
         //gProbvsBack=new TGraph(nData,&CorrectedPValues.front(),&Back.front());
         //c->cd(6);  gProbvsBack->Draw("A*");
         
         c->Update();
       
            
      }  
    
      nerrors++;

   }


 if (Plot) {
    getchar();
   delete hBack;
   delete hSig;
   delete hSigma;
   
   delete c;
   delete fitgaus;
   delete hProb;  
   delete gQ; delete gQ_fit;
 }
 printf ("e=%f\n", BestProb);

 return BestProb;

}   

//(differential) Poisson probability of detecting nev, while expecting bkg+=dbkg
double WeightedP(int nev, float bkg, float dbkg) { 

    double PSum=0,WSum=0;
    float BkgMin=TMath::Max(0.,bkg-3*dbkg),
          BkgMax=bkg+3*dbkg,
          daback=dbkg/10.;
    
    
    for (float aback=BkgMin;aback<=BkgMax;aback+=daback) {
       double weight = ROOT::Math::gaussian_pdf(aback,dbkg,bkg);
       double ap     = ROOT::Math::poisson_pdf(nev,aback);
       PSum+=weight*ap*daback;
       WSum+=weight*daback;
    }
    if (WSum<0.9 || WSum>1) {printf("%s: wsum=%e %d %e %e \n",__FUNCTION__,WSum,nev,bkg,dbkg); exit(1);}
    return PSum/WSum;
}
