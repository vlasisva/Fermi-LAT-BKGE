//author vlasisva@gmail.com
#include "BackgroundEstimator/BackgroundEstimator.h"

int BKGE_NS::MakeGtLikeTemplate(float gtlike_ROI, string GRB_DIR, string DATACLASS) {

 const float error = TOOLS::Get("BKG_ESTIMATE_ERROR");

 char name[1000];
 const char BKG_COMPONENT[3][30]={"GALGAMMAS","CR_EGAL","TOTAL"};
 BackgroundEstimator * Est = new BackgroundEstimator(DATACLASS,-1,-1,-1,false);
 
 const double ENERGY_MIN_DEFAULT =Est->Energy_Min_datafiles;
 const double ENERGY_MAX_DEFAULT =Est->Energy_Max_datafiles;
 const int    ENERGY_BINS_DEFAULT=Est->Energy_Bins_datafiles;

 TH1F hFlatROI = TH1F("hFlatROI","hFlatROI",ENERGY_BINS_DEFAULT,log10(ENERGY_MIN_DEFAULT),log10(ENERGY_MAX_DEFAULT));
 for (int ib=0;ib<=ENERGY_BINS_DEFAULT;ib++) hFlatROI.SetBinContent(ib,gtlike_ROI);
 
 TH1F * hEstParticles[3];
 TFile * fEstParticles[3];
 
 //Create and read estimated background component files
 if (Est->FillBackgroundHist(GRB_DIR, &hFlatROI,TOOLS::Get("GRB_RA"),TOOLS::Get("GRB_DEC"),1)) return -1; //only gamma
 if (Est->FillBackgroundHist(GRB_DIR, &hFlatROI,TOOLS::Get("GRB_RA"),TOOLS::Get("GRB_DEC"),2)) return -1; //CR + extragalactic
 for (int iParType=0;iParType<2;iParType++) {//here open only the two components (do not open the total root file since it is for a 95% ROI instead of flat gtlike one)
      sprintf(name,"%s/%s_bkg_%.0f_%.0f_%s.root",GRB_DIR.c_str(),Est->DataClass.c_str(),ENERGY_MIN_DEFAULT,ENERGY_MAX_DEFAULT,BKG_COMPONENT[iParType]);
      fEstParticles[iParType]= TFile::Open(name);
      if (!fEstParticles[iParType]) {
          printf ("%s: Cannot open file %s\n",__FUNCTION__,name);
          return -1;
      }
      hEstParticles[iParType] = (TH1F*)fEstParticles[iParType]->Get("hCtsvsEnergy_Est");
 }
 //for iParType==2 (total), manually add the two components. I do not open the total bkg file since that is not for the glike ROI but some 95% ROI
 hEstParticles[2]=(TH1F*)hEstParticles[0]->Clone();
 hEstParticles[2]->Add(hEstParticles[1]);

 TH1F * hExposure = (TH1F*)fEstParticles[0]->Get("hExposure");
 if (!hExposure) {
    sprintf(name,"%s/%s_bkg_%.0f_%.0f_%s.root",GRB_DIR.c_str(),Est->DataClass.c_str(),ENERGY_MIN_DEFAULT,ENERGY_MAX_DEFAULT,BKG_COMPONENT[0]);
    printf("%s: Can't find hExposure in file %s\n",__FUNCTION__,name);
    return -1;
 }
 for (int iParType=0;iParType<3;iParType++) {
     sprintf(name,"%s/%s_gtlike_%s_%s.txt",GRB_DIR.c_str(),Est->DataClass.c_str(),(TOOLS::GetS("GRB_NAME")).c_str(),BKG_COMPONENT[iParType]);
     FILE *fgtlike = fopen(name,"w");
     if (hEstParticles[iParType]->GetNbinsX()!=ENERGY_BINS_DEFAULT) {
          printf("hEstParticles bins %d != default bins %d\n",hEstParticles[iParType]->GetNbinsX(),ENERGY_BINS_DEFAULT);
	  throw std::runtime_error("");
     }
     for (int ie=1;ie<=ENERGY_BINS_DEFAULT;ie++) {
        double Flux = hEstParticles[iParType]->GetBinContent(ie); //counts per bin
        double dE = pow(10,hExposure->GetXaxis()->GetBinUpEdge(ie))-pow(10,hExposure->GetXaxis()->GetBinLowEdge(ie));
        Flux/=dE; //counts/MeV
        Flux/=hExposure->GetBinContent(ie); //counts/mev/sec/cm2
        double ROI_AREA = 2*TMath::Pi()*(1-cos(gtlike_ROI*DEG_TO_RAD)); //area of the ROI in sr
	Flux/=ROI_AREA;
    	fprintf(fgtlike,"%f %e %e\n",pow(10,hExposure->GetXaxis()->GetBinCenter(ie)),Flux,Flux*error);
    	//printf("%d %f %f %f %f\n",ie,hEstParticles[iParType]->GetBinContent(ie),dE,hExposure->GetBinContent(ie),ROI_AREA);
    	//printf("%f %e %e\n",pow(10,hExposure->GetXaxis()->GetBinLowEdge(ie)),Flux,Flux*error);
     }
     fclose(fgtlike);

 }
 for (int iParType=0;iParType<2;iParType++) fEstParticles[iParType]->Close();
 return 0;
}
