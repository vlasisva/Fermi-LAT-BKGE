//Author: Vlasios Vasileiou <vlasisva@gmail.com>
#include "BackgroundEstimator/BKGE_Tools.h"

//Integrate a background map over an ROI
//hMap is a map of the background per solid angle
double TOOLS::Integrate(TH2F * hMap, float L_BURST, float B_BURST, float ROI_RADIUS, string MAPNAME) {

  if (ROI_RADIUS>180) {
     printf("%s: ROI_RADIUS (%f)>180? \n",__FUNCTION__,ROI_RADIUS); 
     exit(1);
  }

  if (!hMap) {printf("%s: hMAP is null!!\n",__FUNCTION__); exit(1);}
  TH2F* hIntegrantMap=NULL;
  TFile * fout = NULL;
  if (MAPNAME!="") {
     fout = new TFile(MAPNAME.c_str(),"RECREATE");
     hIntegrantMap = new TH2F("hIntegrantMap","hIntegrantMap",hMap->GetNbinsX(),-180,180,hMap->GetNbinsY(),-90,90);       
  }

  B_BURST+=90;  //move from latitude theta range (-90,90) to CLHEP theta range (0,180) 
  L_BURST*=DEG_TO_RAD;
  B_BURST*=DEG_TO_RAD;
   
  CLHEP::Hep3Vector gBurst;
  gBurst.setRThetaPhi(1,B_BURST,L_BURST);
  const CLHEP::Hep3Vector gBurst_Perp = gBurst.orthogonal();

  float Integrator_dtheta=0.25;
  if (ROI_RADIUS<0.5) Integrator_dtheta=0.01;

  ROI_RADIUS*=DEG_TO_RAD;
  double sum=0;

  for (double theta_local=0;theta_local<ROI_RADIUS;theta_local+=Integrator_dtheta*DEG_TO_RAD) {
      //dphi step changes with theta so that we don't sample as finely near the pole (theta=0)
      //this can be optimized more
      double dphi = (1 - 0.49*(1.0-pow(cos(theta_local),2))) *DEG_TO_RAD;
      //printf("%f %f\n",theta_local/DEG_TO_RAD,dphi/DEG_TO_RAD);
      for (double phi_local=-TMath::Pi();phi_local<TMath::Pi();phi_local+=dphi) {
          double val,plotTheta,plotPhi;
          int bin;
          CLHEP::Hep3Vector gPoint=gBurst;
          gPoint.rotate(theta_local,gBurst_Perp);
          gPoint.rotate(phi_local,gBurst);

          plotTheta=gPoint.theta()/DEG_TO_RAD;
          plotPhi=gPoint.phi()/DEG_TO_RAD;
          bin=hMap->FindBin(plotPhi,plotTheta-90);
          val=hMap->GetBinContent(bin);
          //printf("%f %f %f %f val=%f\n",theta_local/DEG_TO_RAD,phi_local/DEG_TO_RAD,gPoint.theta()/DEG_TO_RAD,gPoint.phi()/DEG_TO_RAD,val);
          if (MAPNAME!="") {
              hIntegrantMap->Fill(plotPhi,plotTheta-90,(Integrator_dtheta*DEG_TO_RAD*dphi)*sin(theta_local)*val);
          }
          sum+=(Integrator_dtheta*DEG_TO_RAD*dphi)*sin(theta_local)*val;
      }
 }


 if (MAPNAME!="") {
   hIntegrantMap->SetContour(200);
   hIntegrantMap->Write();
   printf("%f %f\n",sum,hIntegrantMap->Integral());
   hIntegrantMap->Delete();
   fout->Close();

 }
 return sum;
 
}

