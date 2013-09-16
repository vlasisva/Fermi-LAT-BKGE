#include "BackgroundEstimator/BKGE_Tools.h"
#include "astro/JulianDate.h"
#include "astro/SolarSystem.h"


void TOOLS::GBM_AppendPointingPlots(string PlotsFile) {

  TFile * fPlots = TFile::Open(PlotsFile.c_str(),"update");
  if (fPlots->Get("hGBMSolarAngle_0vsTime")) { //plots file already has GBM data in
     fPlots->Close();
     return;
  }
  
  TH1F * hPtRazvsTime=(TH1F*)fPlots->Get("hPtRazvsTime");
  TH1F * hPtDeczvsTime=(TH1F*)fPlots->Get("hPtDeczvsTime");
  TH1F * hPtRaxvsTime=(TH1F*)fPlots->Get("hPtRaxvsTime");
  TH1F * hPtDecxvsTime=(TH1F*)fPlots->Get("hPtDecxvsTime");
  TH1F * hPtRaZenithvsTime=(TH1F*)fPlots->Get("hRAZenithvsTime");
  TH1F * hPtDecZenithvsTime=(TH1F*)fPlots->Get("hDecZenithvsTime");
  
  char name[1000];
  TH1F * hSolarAngle[14],*hPtGal_L[14],*hPtGal_B[14],*hZenithAngle[14], *hAzimuthAngle[14];

  const int bins=hPtRazvsTime->GetNbinsX();
  const float xmin=hPtRazvsTime->GetXaxis()->GetXmin();
  const float xmax=hPtRazvsTime->GetXaxis()->GetXmax();
  
  for (int idet=0;idet<14;idet++) { //loop over 14 GBM dets

      sprintf(name,"hGBMSolarAngle_%dvsTime",idet);
      hSolarAngle[idet]= new TH1F(name,name,bins,xmin,xmax);

      sprintf(name,"hGBMPtGal_L_%dvsTime",idet);
      hPtGal_L[idet]= new TH1F(name,name,bins,xmin,xmax);

      sprintf(name,"hGBMPtGal_B_%dvsTime",idet);
      hPtGal_B[idet]= new TH1F(name,name,bins,xmin,xmax);

      sprintf(name,"hGBMEarthZenithAngle_%dvsTime",idet);
      hZenithAngle[idet]= new TH1F(name,name,bins,xmin,xmax);

      sprintf(name,"hGBMEarthAzimuthAngle_%dvsTime",idet);
      hAzimuthAngle[idet]= new TH1F(name,name,bins,xmin,xmax);
      
      for (int ib=1;ib<=bins;ib++) {
         float PtRaz=hPtRazvsTime->GetBinContent(ib);
         if (PtRaz==0) continue;
         float PtRax=hPtRaxvsTime->GetBinContent(ib);
         float PtDecz=hPtDeczvsTime->GetBinContent(ib);
         float PtDecx=hPtDecxvsTime->GetBinContent(ib);
         double time_met=hPtRazvsTime->GetBinCenter(ib);
         float sangle=GetGBMDet_SolarAngle(time_met, idet, PtRaz, PtDecz, PtRax, PtDecx);
         hSolarAngle[idet]->SetBinContent(ib,sangle);
         
         float PtGBM_L,PtGBM_B;
         GetGBMDet_GalPt(idet, PtGBM_L, PtGBM_B, PtRaz, PtDecz, PtRax, PtDecx);
         hPtGal_L[idet]->SetBinContent(ib,PtGBM_L);
         hPtGal_B[idet]->SetBinContent(ib,PtGBM_B);
         
         float PtRaZenith=hPtRaZenithvsTime->GetBinContent(ib);
         float PtDecZenith=hPtDecZenithvsTime->GetBinContent(ib);
        
         float EarthZenith,EarthAzimuth;
         GetGBMDet_EarthCoordinates(idet, PtRaZenith, PtDecZenith, PtRaz, PtDecz,PtRax,  PtDecx, EarthZenith, EarthAzimuth);
         hZenithAngle[idet]->SetBinContent(ib,EarthZenith);
         hAzimuthAngle[idet]->SetBinContent(ib,EarthAzimuth);
      }
      hSolarAngle[idet]->Write();
      hPtGal_L[idet]->Write();
      hPtGal_B[idet]->Write();
      hZenithAngle[idet]->Write();
      hAzimuthAngle[idet]->Write();
  } 
  
  fPlots->Close();


}
///////////////////////////////////////////////////////////////////


///////////////////////////////////////////////////////////////////////////////////
int TOOLS::PlotGBMDet_Angle(int nDet, float RA, float DEC, TH1F* h, string Plots_File){

  TFile * fPlots = TFile::Open(Plots_File.c_str());
  TH1F * hPtRazvsTime=(TH1F*)fPlots->Get("hPtRazvsTime");
  TH1F * hPtDeczvsTime=(TH1F*)fPlots->Get("hPtDeczvsTime");
  TH1F * hPtRaxvsTime=(TH1F*)fPlots->Get("hPtRaxvsTime");
  TH1F * hPtDecxvsTime=(TH1F*)fPlots->Get("hPtDecxvsTime");
  printf("%d bins\n",h->GetNbinsX());
  for (int i=1;i<=h->GetNbinsX();i++) {
       int itimebin=hPtRazvsTime->GetXaxis()->FindBin(h->GetXaxis()->GetBinCenter(i));
       float PtRaz=hPtRazvsTime->GetBinContent(itimebin);
       printf("ptraz=%f\n",PtRaz);
       if (PtRaz==0) continue;
       float PtDecz=hPtDeczvsTime->GetBinContent(itimebin);
       float PtRax=hPtRaxvsTime->GetBinContent(itimebin);
       float PtDecx=hPtDecxvsTime->GetBinContent(itimebin);
       float angle=GetGBMDet_Angle(nDet, RA, DEC, PtRaz, PtDecz, PtRax, PtDecx);
       printf("i=%d itimebin=%d %f %f %f %f %f %f angle=%f\n",i,itimebin,RA,DEC,PtRaz, PtDecz, PtRax, PtDecx,angle);
       h->SetBinContent(i,angle);
  }
  fPlots->Close();
  return 0;  
}

/////////////////////////////////////////////////////////////////////////////
int TOOLS::PlotGBMDet_ZenithAngle(int nDet, TH1F* h, string Plots_File){

  TFile * fPlots = TFile::Open(Plots_File.c_str());
  TH1F * hPtRazvsTime=(TH1F*)fPlots->Get("hPtRazvsTime");
  TH1F * hPtDeczvsTime=(TH1F*)fPlots->Get("hPtDeczvsTime");
  TH1F * hPtRaxvsTime=(TH1F*)fPlots->Get("hPtRaxvsTime");
  TH1F * hPtDecxvsTime=(TH1F*)fPlots->Get("hPtDecxvsTime");
  TH1F * hPtRaZenithvsTime=(TH1F*)fPlots->Get("hRAZenithvsTime");
  TH1F * hPtDecZenithvsTime=(TH1F*)fPlots->Get("hDecZenithvsTime");

  for (int i=1;i<=h->GetNbinsX();i++) {
       int itimebin=hPtRazvsTime->GetXaxis()->FindBin(h->GetXaxis()->GetBinCenter(i));
       float PtRaz=hPtRazvsTime->GetBinContent(itimebin);
       if (PtRaz==0) continue;
       float PtDecz=hPtDeczvsTime->GetBinContent(itimebin);
       float PtRax=hPtRaxvsTime->GetBinContent(itimebin);
       float PtDecx=hPtDecxvsTime->GetBinContent(itimebin);
       float PtRaZenith=hPtRaZenithvsTime->GetBinContent(itimebin);
       float PtDecZenith=hPtDecZenithvsTime->GetBinContent(itimebin);
       float angle=GetGBMDet_Angle(nDet, PtRaZenith, PtDecZenith, PtRaz, PtDecz, PtRax, PtDecx);
       h->SetBinContent(i,angle);
  }
  fPlots->Close();
  return 0;  
}

//////////////////////////////////////////////////////////////////////////////////////////////////
int TOOLS::PlotGBMDet_GalCenterAngle(int nDet, TH1F* h, string Plots_File){

  TFile * fPlots = TFile::Open(Plots_File.c_str());
  TH1F * hPtRazvsTime=(TH1F*)fPlots->Get("hPtRazvsTime");
  TH1F * hPtDeczvsTime=(TH1F*)fPlots->Get("hPtDeczvsTime");
  TH1F * hPtRaxvsTime=(TH1F*)fPlots->Get("hPtRaxvsTime");
  TH1F * hPtDecxvsTime=(TH1F*)fPlots->Get("hPtDecxvsTime");

  for (int i=1;i<=h->GetNbinsX();i++) {
       int itimebin=hPtRazvsTime->GetXaxis()->FindBin(h->GetXaxis()->GetBinCenter(i));
       float PtRaz=hPtRazvsTime->GetBinContent(itimebin);
       if (PtRaz==0) continue;
       float PtDecz=hPtDeczvsTime->GetBinContent(itimebin);
       float PtRax=hPtRaxvsTime->GetBinContent(itimebin);
       float PtDecx=hPtDecxvsTime->GetBinContent(itimebin);
       h->SetBinContent(i,GetGBMDet_Angle(nDet, 266.36445, -28.994293, PtRaz, PtDecz, PtRax, PtDecx));
  }
  fPlots->Close();
  return 0;  
}
//////////////////////////////////////////////////////////////////

float TOOLS::GetGBMDet_SolarAngle(double time_met, int nDet, float PtRaz, float PtDecz, float PtRax, float PtDecx){

    astro::JulianDate julian(astro::JulianDate::missionStart()+time_met/astro::JulianDate::secondsPerDay);
    //    astro::JulianDate julian2(time);
    //std::cout<<" test "<<julian<<" "<<julian2<<" "<<std::endl;
    //astro::SkyDir moondir=moon.direction(julian);
    //astro::SkyDir earthdir=earth.direction(julian);
    astro::SolarSystem sun(astro::SolarSystem::SUN);
    astro::SkyDir sundir=sun.direction(julian);
    return GetGBMDet_Angle(nDet, sundir.ra(), sundir.dec(), PtRaz, PtDecz, PtRax,       PtDecx);
}

/////////////////////////////////////////////////////////////////
int TOOLS::PlotSolar_Angle(int nDet,  TH1F* h, string Plots_File){

  TFile * fPlots = TFile::Open(Plots_File.c_str());
  TH1F * hPtRazvsTime=(TH1F*)fPlots->Get("hPtRazvsTime");
  TH1F * hPtDeczvsTime=(TH1F*)fPlots->Get("hPtDeczvsTime");
  TH1F * hPtRaxvsTime=(TH1F*)fPlots->Get("hPtRaxvsTime");
  TH1F * hPtDecxvsTime=(TH1F*)fPlots->Get("hPtDecxvsTime");
  for (int i=1;i<=h->GetNbinsX();i++) {
       int itimebin=hPtRazvsTime->GetXaxis()->FindBin(h->GetXaxis()->GetBinCenter(i));
       float PtRaz=hPtRazvsTime->GetBinContent(itimebin);
       if (PtRaz==0) continue;
       float PtDecz=hPtDeczvsTime->GetBinContent(itimebin);
       float PtRax=hPtRaxvsTime->GetBinContent(itimebin);
       float PtDecx=hPtDecxvsTime->GetBinContent(itimebin);
       h->SetBinContent(i,GetGBMDet_SolarAngle(h->GetXaxis()->GetBinCenter(i), nDet, PtRaz,  PtDecz, PtRax, PtDecx));
  }
  fPlots->Close();
  return 0;  
}

//////////////////////////////////////////////////////////////////////////
float TOOLS::GetGBMDet_Angle(int nDet, float RA, float DEC, float PtRaz, float PtDecz, float PtRax, float PtDecx){

 astro::SkyDir SCz = astro::SkyDir(PtRaz,PtDecz,astro::SkyDir::EQUATORIAL);
 astro::SkyDir SCx = astro::SkyDir(PtRax,PtDecx,astro::SkyDir::EQUATORIAL);
 float SCz_SCx_diff=SCz.difference(SCx)/DEG_TO_RAD;
 if (fabs(SCz_SCx_diff-90)>5) {
     printf("%s: SCz and SCz not normal; angle=%f\n",__FUNCTION__,SCz_SCx_diff);
     
 }
 CLHEP::Hep3Vector GBMDir_local;
 GBMDir_local.setRThetaPhi(1,GBMDet_Theta_Azimuth[nDet][0]*DEG_TO_RAD,GBMDet_Theta_Azimuth[nDet][1]*DEG_TO_RAD);
 CLHEP::HepRotation localToCelestial = HepRotation(SCx.dir(),(SCz.dir()).cross(SCx.dir()),SCz.dir());
 localToCelestial.setTolerance(0.2);
 astro::SkyDir GBMDir_cel = astro::SkyDir(localToCelestial*GBMDir_local);

 astro::SkyDir Source_cel = astro::SkyDir(RA,DEC,astro::SkyDir::EQUATORIAL);
 
/*
 astro::PointingTransform p = astro::PointingTransform(SCz,SCx);
 CLHEP::HepRotation inverse = inverseOf(p.localToCelestial());
 astro::SkyDir SDSource_Local = astro::SkyDir(inverse*SDSource());
*/
 return GBMDir_cel.difference(Source_cel)/DEG_TO_RAD;

}

////////////////////////////////////////////////////////////////////////////////
void TOOLS::GetGBMDet_EarthCoordinates(int nDet, float RAZenith, float DecZenith, float PtRaz, float PtDecz, float PtRax, float PtDecx, float & EarthZenith, float &EarthAzimuth){

  astro::SkyDir SCz = astro::SkyDir(PtRaz,PtDecz,astro::SkyDir::EQUATORIAL);
  astro::SkyDir SCx = astro::SkyDir(PtRax,PtDecx,astro::SkyDir::EQUATORIAL);

  CLHEP::Hep3Vector GBMDir_local;
  GBMDir_local.setRThetaPhi(1,GBMDet_Theta_Azimuth[nDet][0]*DEG_TO_RAD,GBMDet_Theta_Azimuth[nDet][1]*DEG_TO_RAD);
  CLHEP::HepRotation localToCelestial = HepRotation(SCx.dir(),(SCz.dir()).cross(SCx.dir()),SCz.dir());
  localToCelestial.setTolerance(0.2);
  astro::SkyDir GBMDir_cel = astro::SkyDir(localToCelestial*GBMDir_local);

  astro::SkyDir zenith = astro::SkyDir(RAZenith, DecZenith, astro::SkyDir::EQUATORIAL);
  Hep3Vector north_pole(0,0,1);
  Hep3Vector east_dir(north_pole.cross(zenith()).unit() ); // east is perp to north_pole and zenith
  Hep3Vector north_dir(zenith().cross(east_dir));
  
  EarthAzimuth=atan2( GBMDir_cel().dot(east_dir), GBMDir_cel().dot(north_dir) ); // z is north, heading.
  if( EarthAzimuth <0) EarthAzimuth += 2*M_PI; // to 0-360 deg.
  EarthAzimuth/=DEG_TO_RAD;

  EarthZenith=GBMDir_cel.difference(zenith)/DEG_TO_RAD;

}


/////////////////////////////////////////////////////////////////////////////

void TOOLS::GetGBMDet_GalPt(int nDet, float & PtGBM_L, float & PtGBM_B, float PtRaz, float PtDecz, float PtRax, float PtDecx){

 astro::SkyDir SCz = astro::SkyDir(PtRaz,PtDecz,astro::SkyDir::EQUATORIAL);
 astro::SkyDir SCx = astro::SkyDir(PtRax,PtDecx,astro::SkyDir::EQUATORIAL);
 float SCz_SCx_diff=SCz.difference(SCx)/DEG_TO_RAD;
 if (fabs(SCz_SCx_diff-90)>5) {
     printf("%s: SCz and SCz not normal; angle=%f\n",__FUNCTION__,SCz_SCx_diff);
     
 }
 CLHEP::Hep3Vector GBMDir_local;
 GBMDir_local.setRThetaPhi(1,GBMDet_Theta_Azimuth[nDet][0]*DEG_TO_RAD,GBMDet_Theta_Azimuth[nDet][1]*DEG_TO_RAD);
 CLHEP::HepRotation localToCelestial = HepRotation(SCx.dir(),(SCz.dir()).cross(SCx.dir()),SCz.dir());
 localToCelestial.setTolerance(0.2);
 astro::SkyDir GBMDir_cel = astro::SkyDir(localToCelestial*GBMDir_local);
 PtGBM_L=GBMDir_cel.l();
 PtGBM_B=GBMDir_cel.b();
 
}

/////////////////////////////////////////////////////////////
int TOOLS::PlotGBMGal_Coordinates(int nDet, TH1F* hL, TH1F* hB, string Plots_File){

  if (hL->GetNbinsX()!=hB->GetNbinsX() || hL->GetXaxis()->GetXmin()!=hB->GetXaxis()->GetXmin()|| hL->GetXaxis()->GetXmax()!=hB->GetXaxis()->GetXmax()) return -1;
  
  TFile * fPlots = TFile::Open(Plots_File.c_str());
  TH1F * hPtRazvsTime=(TH1F*)fPlots->Get("hPtRazvsTime");
  TH1F * hPtDeczvsTime=(TH1F*)fPlots->Get("hPtDeczvsTime");
  TH1F * hPtRaxvsTime=(TH1F*)fPlots->Get("hPtRaxvsTime");
  TH1F * hPtDecxvsTime=(TH1F*)fPlots->Get("hPtDecxvsTime");
  for (int i=1;i<=hL->GetNbinsX();i++) {
       int itimebin=hPtRazvsTime->GetXaxis()->FindBin(hL->GetXaxis()->GetBinCenter(i));
       float PtRaz=hPtRazvsTime->GetBinContent(itimebin);
       if (PtRaz==0) continue;
       float PtDecz=hPtDeczvsTime->GetBinContent(itimebin);
       float PtRax=hPtRaxvsTime->GetBinContent(itimebin);
       float PtDecx=hPtDecxvsTime->GetBinContent(itimebin);
       float PtGBM_L,PtGBM_B;
       TOOLS::GetGBMDet_GalPt(nDet,PtGBM_L, PtGBM_B, PtRaz, PtDecz, PtRax, PtDecx);
       if (PtGBM_B>=360 || PtGBM_L>=360) printf("too high %f %f %f %f %f %f\n",PtGBM_B,PtGBM_L,PtRaz,PtDecz,PtRax,PtDecx);
       //printf("put %f %f to %d\n",PtGBM_L,PtGBM_B,i);
       hL->SetBinContent(i,PtGBM_L);
       hB->SetBinContent(i,PtGBM_B);
  }
  fPlots->Close();
  return 0;  
}


