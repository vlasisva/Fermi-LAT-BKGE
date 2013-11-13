// Author: Vlasios Vasileiou <vlasisva@gmail.com>
#include "BackgroundEstimator/BKGE_Tools.h"
#include "rootIrfLoader/rootIrfLoader/Psf.h"
#include "rootIrfLoader/rootIrfLoader/Aeff.h"
#include <algorithm>

void TOOLS::CalculatePSF(TH1F * hROI, double MET, string FT2_FILE, string DATACLASS, float Containment, float MaxRadius) {
  double theta,phi,ztheta;

  GetThetaPhi(theta,phi,ztheta,MET,FT2_FILE);
  if (theta>70) {
    printf("%s: Theta too high (%f). Calculating PSF for theta=70deg\n",__FUNCTION__,theta);
    theta=70;
  }
  else if (theta<0){
    printf("%s: Cannot read pointing data. Calculating PSF for theta=40deg\n",__FUNCTION__);
    theta=40;
  }
  CalculatePSF_ThetaPhi(hROI, theta, phi, DATACLASS, Containment, MaxRadius);
}

void TOOLS::CalculatePSF_ThetaPhi(TH1F* hROI, float theta, float phi, string DataClass, float Containment, float MaxRadius) {

  if (Containment<0) Containment=Get("ROI_CONTAINMENT");
  if (MaxRadius<0)   MaxRadius  =Get("ROI_MAX_RADIUS");
  double LocalizationError      =Get("ROI_LOCALIZATION_ERROR");

  char name[1000];
  if (theta>80) {
     printf("%s: Theta is >80, can't calculate ROI.\n",__FUNCTION__);
     throw std::runtime_error("");
  }

  rootIrfLoader::Psf PSF_FRONT  =  rootIrfLoader::Psf(DataClass+string("::FRONT"));
  rootIrfLoader::Aeff EA_FRONT  =  rootIrfLoader::Aeff(DataClass+string("::FRONT"));
  rootIrfLoader::Psf PSF_BACK   =  rootIrfLoader::Psf(DataClass+string("::FRONT"));
  rootIrfLoader::Aeff EA_BACK   =  rootIrfLoader::Aeff(DataClass+string("::FRONT"));



  double act;
  sprintf(name,"%.0f%% Containment PSF",Containment*100);
  double E;
  for (int i=1;i<=hROI->GetNbinsX();i++) {
     E = pow(10,hROI->GetBinCenter(i));
     float radius=0;
     act=0;
        double AEFront = EA_FRONT(E,theta,phi);
        double AEBack  = EA_BACK(E,theta,phi);
        if (AEFront<=0 || AEBack<=0) {printf("%s: Effective area weird! AEfron=%f aeback=%f theta=%f phi=%f E=%f theta=%f\n",__FUNCTION__,AEFront,AEBack,theta,phi,E,theta); throw std::runtime_error("");}

        do { 
            float FrontInt=PSF_FRONT.angularIntegral(E,theta,phi,radius);
            float BackInt =PSF_BACK.angularIntegral(E,theta,phi,radius);
            act = (AEBack*BackInt + AEFront*FrontInt)/(AEBack+AEFront);
            radius+=0.1;
        } while (Containment>act);
        radius-=2*0.1;
        do {
            float FrontInt=PSF_FRONT.angularIntegral(E,theta,phi,radius);
            float BackInt =PSF_BACK.angularIntegral(E,theta,phi,radius);
            act = (AEBack*BackInt + AEFront*FrontInt)/(AEBack+AEFront);
            radius+=0.01;
        } while (Containment>act);
        radius-=2*0.01;
        do {
            float FrontInt=PSF_FRONT.angularIntegral(E,theta,phi,radius);
            float BackInt =PSF_BACK.angularIntegral(E,theta,phi,radius);
            act = (AEBack*BackInt + AEFront*FrontInt)/(AEBack+AEFront);
            radius+=0.001;
        } while (Containment>act);

     //printf("%d Energy=%f theta=%f phi=%f radius=%f ct=%f\n",i,E,theta,phi,radius,act);
     if (!(radius>0)) {printf("%s: ROI Radius invalid!! RADIUS=%f E=%f theta=%f phi=%f\n",__FUNCTION__,radius,E,theta,phi); throw std::runtime_error("");}
     hROI->SetBinContent(i,std::min(sqrt(pow(radius,2)+pow(LocalizationError,2)),(double)MaxRadius));
  }
}



int TOOLS::ReadROI_File(TH1F * hROI, string filename) {
   printf("%s: Using ROI Radius from file %s.\n",__FUNCTION__,filename.c_str());
   FILE * fROI = fopen(filename.c_str(),"r");
   if (!fROI) {printf("%s: Can't open ROI-data file %s\n",__FUNCTION__,filename.c_str()); return -1;}
   float aroi;
   for (int i=1;i<=hROI->GetNbinsX();i++) {
      if (fscanf(fROI,"%f",&aroi)!=1) {
          printf("%s: Error reading file %s\nDid you provide exactly %d values?\n",__FUNCTION__,filename.c_str(),hROI->GetNbinsX());
          return -1;
      }
      hROI->SetBinContent(i,aroi);
   }
   if (fscanf(fROI,"%f",&aroi)==1) {
       printf("%s: Warning: File %s contains more valus than the available ENERGY_BINS (%d)\n",__FUNCTION__,filename.c_str(),hROI->GetNbinsX());
       return -1;
   }
   fclose (fROI);
   return 0;
}

