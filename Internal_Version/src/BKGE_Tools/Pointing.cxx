//Author: Vlasios Vasileiou <vlasisva@gmail.com>
//$Header: /nfs/slac/g/glast/ground/cvs/GRBAnalysis-scons/BackgroundEstimator/src/BKGE_Tools/Pointing.cxx,v 1.2 2011/06/13 16:59:55 vlasisva Exp $
#include "BackgroundEstimator/BKGE_Tools.h"
#include "TLine.h"

TCanvas * TOOLS::MakePointingPlots(double PreTime, double PostTime, double GRB_t0, double RA, double DEC, string FT2_FILE) {

 if (RA<-999 || DEC<-999) {
    RA =  Get("GRB_RA");
    DEC = Get("GRB_DEC");
 }
 char PlotsFile[] = "/tmp/plots.root";
 if (Make_Plots(PreTime,PostTime,GRB_t0,PlotsFile,FT2_FILE)) return NULL;
 
 const int NbinsX=int(PostTime+PreTime);
 const double T_MIN=-PreTime;
 const double T_MAX=+PostTime;

 CLHEP::HepRotation inverse;
 inverse.setTolerance(0.1);
 gROOT->cd("/");
 TH1F * hTheta = (TH1F*)gROOT->Get("hTheta");if (hTheta) delete hTheta;
 TH1F * hPhi = (TH1F*)gROOT->Get("hPhi");  if (hPhi) delete hPhi;
 TH1F * hZTheta = (TH1F*)gROOT->Get("hZTheta");  if (hZTheta) delete hZTheta;
 TH1F * hRock = (TH1F*)gROOT->Get("hRockingAngle");  if (hRock) delete hRock;
 TH1F * hMcIlwainL = (TH1F*)gROOT->Get("hMcIlwainL");  if (hMcIlwainL) delete hMcIlwainL;
 TH1F * hMcIlwainB = (TH1F*)gROOT->Get("hMcIlwainB");  if (hMcIlwainB) delete hMcIlwainB;
 hTheta = new TH1F("hTheta","Theta",NbinsX,T_MIN,T_MAX);
 hPhi = new TH1F("hPhi","Phi",NbinsX,T_MIN,T_MAX);
 hZTheta = new TH1F("hZTheta","Zenith Theta",NbinsX,T_MIN,T_MAX);
 hRock = new TH1F("hRockingAngle","Rocking Angle",NbinsX,T_MIN,T_MAX);
 hMcIlwainL = new TH1F("hMcIlwainL","McIlwain L",NbinsX,T_MIN,T_MAX);
 hMcIlwainB = new TH1F("hMcIlwainB","McIlwain B",NbinsX,T_MIN,T_MAX);
 hTheta->GetYaxis()->SetTitle("Theta (deg)");
 hZTheta->GetYaxis()->SetTitle("Zenith Theta (deg)");
 hRock->GetYaxis()->SetTitle("Rocking Angle (deg)");
 hMcIlwainL->GetYaxis()->SetTitle("McIlwain L");
 hMcIlwainB->GetYaxis()->SetTitle("McIlwain B");

 hTheta->GetXaxis()->SetTitle("t-t_{burst} (sec)");
 hPhi->GetXaxis()->SetTitle("t-t_{burst} (sec)");
 hZTheta->GetXaxis()->SetTitle("t-t_{burst} (sec)");
 hRock->GetXaxis()->SetTitle("t-t_{burst} (sec)");
 hMcIlwainL->GetXaxis()->SetTitle("t-t_{burst} (sec)");
 hMcIlwainB->GetXaxis()->SetTitle("t-t_{burst} (sec)");

 astro::SkyDir DirBurst = astro::SkyDir(RA,DEC,astro::SkyDir::EQUATORIAL);
 CLHEP::Hep3Vector dir = DirBurst();

 TFile * fBurst = TFile::Open(PlotsFile,"r");
 TH1F * hPtRazvsTime = (TH1F*)fBurst->Get("hPtRazvsTime");
 if (!hPtRazvsTime) {printf("%s: Something went wrong with making the pointing plots.. \n",__FUNCTION__); return NULL;}
 TH1F * hPtRaxvsTime = (TH1F*)fBurst->Get("hPtRaxvsTime");
 TH1F * hPtDeczvsTime = (TH1F*)fBurst->Get("hPtDeczvsTime");
 TH1F * hPtDecxvsTime = (TH1F*)fBurst->Get("hPtDecxvsTime");
 TH1F * hRAZenithvsTime = (TH1F*)fBurst->Get("hRAZenithvsTime");
 TH1F * hDecZenithvsTime = (TH1F*)fBurst->Get("hDecZenithvsTime");
 TH1F * hRock_=(TH1F*)fBurst->Get("hRockingAnglevsTime");
 TH1F * hMcIlwainL_=(TH1F*)fBurst->Get("hMcIlwainLvsTime");
 TH1F * hMcIlwainB_=(TH1F*)fBurst->Get("hMcIlwainBvsTime");

 using namespace astro;
 for (int i=1;i<=NbinsX;i++) {
     hRock->SetBinContent(i,hRock_->GetBinContent(i));
     hMcIlwainL->SetBinContent(i,hMcIlwainL_->GetBinContent(i));
     hMcIlwainB->SetBinContent(i,hMcIlwainB_->GetBinContent(i));
     if (hPtRazvsTime->GetBinContent(i)==0) continue; //not in GTI
     astro::SkyDir SCz(hPtRazvsTime->GetBinContent(i),hPtDeczvsTime->GetBinContent(i),astro::SkyDir::EQUATORIAL);
     astro::SkyDir SCx(hPtRaxvsTime->GetBinContent(i),hPtDecxvsTime->GetBinContent(i),astro::SkyDir::EQUATORIAL);
     astro::SkyDir SCEarth(hRAZenithvsTime->GetBinContent(i),hDecZenithvsTime->GetBinContent(i),astro::SkyDir::EQUATORIAL);
     //printf("%f %f %f %f\n",hPtRazvsTime->GetBinContent(i),hPtRaxvsTime->GetBinContent(i),hPtDeczvsTime->GetBinContent(i),h
     PointingTransform p = PointingTransform(SCz,SCx);
     inverse = CLHEP::inverseOf(p.localToCelestial()); 
     SkyDir LocalDir = SkyDir(inverse*DirBurst());
     double y=LocalDir().y();
     double x=LocalDir().x();
     double phi=atan2(y,x)*RAD_TO_DEG;
     //double theta=acos(LocalDir().z())*RAD_TO_DEG; 

     hTheta->SetBinContent(i,SCz.difference(DirBurst)/DEG_TO_RAD);
     hZTheta->SetBinContent(i,SCEarth.difference(DirBurst)/DEG_TO_RAD);
     hPhi->SetBinContent(i,phi);

 }
 fBurst->Close();
 TCanvas * cPoint  = new TCanvas("cPoint","Pointing",500,800); cPoint->Divide(1,3);
 TCanvas * cPoint2 = new TCanvas("cPoint2","Pointing2",500,800); cPoint2->Divide(1,3);

 cPoint->cd(1); hTheta->SetLineColor(4); hTheta->SetLineWidth(2); hTheta->GetYaxis()->SetRangeUser(0,180);  hTheta->Draw();
 cPoint->cd(2); hZTheta->SetLineColor(4);hZTheta->SetLineWidth(2); hZTheta->GetYaxis()->SetRangeUser(0,180); hZTheta->Draw();
 cPoint->cd(3); hPhi->SetLineColor(4);hPhi->SetLineWidth(2); hPhi->Draw();

 cPoint2->cd(1); hRock->SetLineColor(4); hRock->SetLineWidth(2); hRock->GetYaxis()->SetRangeUser(0,100); hRock->Draw();
 cPoint2->cd(2); hMcIlwainL->SetLineColor(4); hMcIlwainL->SetLineWidth(2); hMcIlwainL->GetYaxis()->SetRangeUser(0.9,1.8); hMcIlwainL->Draw();
 cPoint2->cd(3); hMcIlwainB->SetLineColor(4); hMcIlwainB->SetLineWidth(2); hMcIlwainB->GetYaxis()->SetRangeUser(0.8,6.2); hMcIlwainB->Draw();
 for (int i=0;i<3;i++) cPoint->GetPad(i+1)->SetGridy();
 for (int i=0;i<3;i++) cPoint2->GetPad(i+1)->SetGridy();

 TLine * l = new TLine(0,hTheta->GetMinimum(),0,hTheta->GetMaximum());
 l->SetLineColor(3);  cPoint->cd(1); l->Draw();
 if (hTheta->GetMaximum()>70) {
    l = new TLine(-PreTime,70,PostTime,70);
    l->SetLineColor(2);
    l->SetLineStyle(2);
    cPoint->cd(1);l->Draw();
 }

 if (hZTheta->GetMaximum()>105) {
    l = new TLine(-PreTime,105,PostTime,105);
    l->SetLineColor(2);
    l->SetLineStyle(2);
    cPoint->cd(2);l->Draw();
 }

 l = new TLine(0,hZTheta->GetMinimum(),0,hZTheta->GetMaximum());
 l->SetLineColor(3);  cPoint->cd(2); l->Draw();

 l = new TLine(0,hRock->GetMinimum(),0,hRock->GetMaximum());
 l->SetLineColor(3);  cPoint2->cd(1); l->Draw();

 //remove(PlotsFile);

 printf("%s: At t0: Theta=%.1f, ZTheta=%.1f, Rocking Angle=%.1f deg\n",__FUNCTION__,
     hTheta->GetBinContent(NbinsX/2), hZTheta->GetBinContent(NbinsX/2), hRock->GetBinContent(NbinsX/2));
 /*
 char name[1000];
 string GRB_NAME=TOOLS::GetS("GRB_NAME");
 sprintf(name,"Results/%s/%s_Pointing1.png",GRB_NAME.c_str(),GRB_NAME.c_str());
 cPoint->SaveAs(name);
 sprintf(name,"Results/%s/%s_Pointing2.png",GRB_NAME.c_str(),GRB_NAME.c_str());
 cPoint2->SaveAs(name);
 */
 return cPoint;
}


void TOOLS::GetThetaPhi(double &theta, double &phi, double &ZTheta, double MET, string FT2_FILE, float GRB_RA, float GRB_DEC) {
  if (GRB_RA<-400)  GRB_RA=TOOLS::Get("GRB_RA");
  if (GRB_DEC<-400) GRB_DEC=TOOLS::Get("GRB_DEC");

  char PlotsFile[1000];
  sprintf(PlotsFile,"/tmp/plots_%d.root",rand());
  if (Make_Plots(1,1,MET,PlotsFile,FT2_FILE)){
     theta=phi=ZTheta=-999; 
     printf("%s: Problem making plots.. returning\n",__FUNCTION__);
     return;
  }
  TFile * fPlots = TFile::Open(PlotsFile);
  TH1F* hPtRazvsTime = (TH1F*)fPlots->Get("hPtRazvsTime");
  if (!hPtRazvsTime) {
      printf("%s: Can't read histograms from file %s\n",__FUNCTION__,PlotsFile);
      theta=phi=ZTheta=-9999;
      remove(PlotsFile);
      fPlots->Close();
      return;
  }
  float PtRaz = hPtRazvsTime->GetBinContent(1);
  float PtRax = ((TH1F*)fPlots->Get("hPtRaxvsTime"))->GetBinContent(1);
  float PtDecx = ((TH1F*)fPlots->Get("hPtDecxvsTime"))->GetBinContent(1);
  float PtDecz = ((TH1F*)fPlots->Get("hPtDeczvsTime"))->GetBinContent(1);
  float RAZenith = ((TH1F*)fPlots->Get("hRAZenithvsTime"))->GetBinContent(1);
  float DecZenith = ((TH1F*)fPlots->Get("hDecZenithvsTime"))->GetBinContent(1);

  if (PtRaz==0 || PtRax==0 || PtDecx==0 || PtDecz==0) {
      printf("%s: Error reading pointing data (time=%lf)\n",__FUNCTION__,MET);
      theta=phi=ZTheta=-9999;
      remove(PlotsFile);
      fPlots->Close();
      return;
  }
  fPlots->Close();

  double x,y;
  astro::SkyDir SCz(PtRaz,PtDecz,astro::SkyDir::EQUATORIAL);
  astro::SkyDir SCx(PtRax,PtDecx,astro::SkyDir::EQUATORIAL);
  astro::PointingTransform p = astro::PointingTransform(SCz,SCx);
  astro::SkyDir ZenithDir(RAZenith,DecZenith,astro::SkyDir::EQUATORIAL);

  astro::SkyDir DirBurst = astro::SkyDir(GRB_RA,GRB_DEC,astro::SkyDir::EQUATORIAL);
  CLHEP::Hep3Vector dir = DirBurst();
  CLHEP::HepRotation inverse = CLHEP::inverseOf(p.localToCelestial());
  astro::SkyDir LocalDir = astro::SkyDir(inverse*DirBurst(),astro::SkyDir::EQUATORIAL);
  y=LocalDir().y();
  x=LocalDir().x();
  phi=atan2(y,x)*RAD_TO_DEG;
  if (phi<0) phi+=360.;
  theta=acos(LocalDir().z())*RAD_TO_DEG;
  ZTheta=ZenithDir.difference(DirBurst)*RAD_TO_DEG;

  remove(PlotsFile);
  fPlots->Close();
}



float TOOLS::GimmeEarthAzimuth(float C1, float C2, float C1Zenith, float C2Zenith, astro::SkyDir::CoordSystem coordinate) {
  //C1/C2 are RA/DEC or L/B

  astro::SkyDir zenith = astro::SkyDir(C1Zenith,C2Zenith,coordinate);
  Hep3Vector north_pole(0,0,1);
  Hep3Vector east_dir(north_pole.cross(zenith()).unit() ); // east is perp to north_pole and zenith
  Hep3Vector north_dir(zenith().cross(east_dir));
  astro::SkyDir Dir = astro::SkyDir(C1,C2,coordinate);
  float azimuth=atan2( Dir().dot(east_dir), Dir().dot(north_dir) ); // z is north, heading.
  if( azimuth <0) azimuth += 2*M_PI; // to 0-360 deg.
  return azimuth/DEG_TO_RAD;
}
