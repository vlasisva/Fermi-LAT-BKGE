// Author: Vlasios Vasileiou <vlasisva@gmail.com>
// $Header: /nfs/slac/g/glast/ground/cvs/GRBAnalysis-scons/BackgroundEstimator/src/BKGE_Tools/CalcExposure.cxx,v 1.1.1.1 2011/06/02 19:41:04 jrb Exp $
#include "BackgroundEstimator/BKGE_Tools.h"
#include "rootIrfLoader/rootIrfLoader/Aeff.h"

void TOOLS::CalcExposure(TFile * fResults, float L_BURST, float B_BURST, float FT1ZenithTheta_Cut, TH1F * hExposure, string GRB_DIR, int verbosity){

  char name[1000];
  //This code averages the contents of the exposure map by circularly scanning it over the ROI area. 
  TH2F * hExposureMap = NULL;

  CLHEP::Hep3Vector gBurst;
  gBurst.setRThetaPhi(1,(B_BURST+90)*DEG_TO_RAD,L_BURST*DEG_TO_RAD);
  const CLHEP::Hep3Vector gBurst_Perp = gBurst.orthogonal();
  const double lEnergy_Min_user=hExposure->GetXaxis()->GetXmin();
  const double lEnergy_Max_user=hExposure->GetXaxis()->GetXmax();
  const double Energy_Min_user=pow(10,lEnergy_Min_user);
  const double Energy_Max_user=pow(10,lEnergy_Max_user);
  const int Energy_Bins_user=hExposure->GetNbinsX();

  int Energy_Bins;
  double Energy_Max,Energy_Min;

  sscanf(fResults->Get("Energy_Data")->GetTitle(),"%lf-%lf-%d",&Energy_Min,&Energy_Max,&Energy_Bins);
  int ie_max=std::max(Energy_Bins,Energy_Bins_user);
  bool Default_Binning=true;

  string GRB_NAME,FT1_FILE,FT2_FILE,DataClass;
  double Burst_t0,Burst_dur;
  if (Energy_Bins!=Energy_Bins_user || fabs(log10(Energy_Min)/lEnergy_Min_user-1)>0.01 || fabs(log10(Energy_Max)/lEnergy_Max_user-1)>0.01 ) {
     printf("%s: You are using a non-default energy binning. Have to recalculate the exposures.. \n",__FUNCTION__);
     printf("%s: Default energy binning: %d bins from %.0f to %.0f MeV.. \n",__FUNCTION__,Energy_Bins,Energy_Min,Energy_Max);
     printf("%s: Your energy binning   : %d bins from %.0f to %.0f MeV.. \n",__FUNCTION__,Energy_Bins_user,Energy_Min_user,Energy_Max_user); 
     Default_Binning=false;
     if (fResults->Get("GRB_NAME")==0) {
         printf("%s: You are using a custom binning and a background file generated with a previous version of the bkg estimator. \n",__FUNCTION__);
         printf("%s: Please use the default binning or regenerate your background files and try again \n",__FUNCTION__);
         exit(1);
     }
     GRB_NAME     = (((TNamed*)fResults->Get("GRB_NAME"))->GetTitle());
     FT1_FILE     = (((TNamed*)fResults->Get("FT1_FILE"))->GetTitle());
     FT2_FILE     = (((TNamed*)fResults->Get("FT2_FILE"))->GetTitle());
     DataClass    = (((TNamed*)fResults->Get("DataClass"))->GetTitle());
     sscanf((((TNamed*)fResults->Get("Time_Data"))->GetTitle()),"%lf-%lf",&Burst_t0,&Burst_dur);
     hExposureMap = new TH2F("hExposureMap","hExposureMap",360,-180,180,180,-90,90); //1deg
  }

  for (int ie=1;ie<=ie_max;ie++) {
     ProgressBar(ie-1,ie_max);
     if (!Default_Binning) {
        sprintf(name,"%s/Burst_Exposure_user.fits",GRB_DIR.c_str());
        if (ie==1) TOOLS::Run_gtexpcube(GRB_DIR, Burst_t0, Burst_t0+Burst_dur, FT2_FILE, DataClass, FT1ZenithTheta_Cut, name, Energy_Min_user, Energy_Max_user, ie_max, verbosity);

        ReadExposureMap(name,hExposureMap,ie,verbosity);
     }
     else { //default binning
        sprintf(name,"hExposure_Burst_1deg_%d",ie);
        hExposureMap=(TH2F*)fResults->Get(name);
     }
     float ROI_RADIUS = hExposure->GetBinContent(ie)*DEG_TO_RAD; //hexposure comes filled with the roi values

     int nadd=0;
     double sum=0;

     for (float theta_local=0;theta_local<ROI_RADIUS;theta_local+=ROI_RADIUS/20) {
         //dphi step changes with theta so that we don't sample as finely near the pole (theta=0)
         //this can be optimized more
         int phisteps= (int)abs(200*sin(theta_local)/sin(12*TMath::Pi()/180.));
         //if (phisteps==0 && theta_local) {printf("%s: phisteps=0? theta_local=%f\n",__FUNCTION__,theta_local); exit(1);}
         float dphi = 2*TMath::Pi()/phisteps;
         //printf("%f %f\n",theta_local/DEG_TO_RAD,dphi/DEG_TO_RAD);
         for (float phi_local=0;phi_local<2*TMath::Pi();phi_local+=dphi) {
             CLHEP::Hep3Vector gPoint=gBurst;
             gPoint.rotate(theta_local,gBurst_Perp);
             gPoint.rotate(phi_local,gBurst);
        
             float plotTheta=gPoint.theta()/DEG_TO_RAD;
             float plotPhi=gPoint.phi()/DEG_TO_RAD;
             float val=hExposureMap->GetBinContent(hExposureMap->FindBin(plotPhi,plotTheta-90));
             //printf("%f %f %f %f val=%f sum=%f\n",theta_local/DEG_TO_RAD,phi_local/DEG_TO_RAD,gPoint.theta()/DEG_TO_RAD,gPoint.phi()/DEG_TO_RAD,val,sum);
             sum+=val;
             nadd++;
         }
     }
     //printf ("divide by %d\n",nadd);
     sum/=nadd;

     //int bin=hExposureMap->FindBin(L_BURST,B_BURST);
     hExposure->SetBinContent(ie,sum);
  }
  delete hExposureMap;

}


double TOOLS::CalcSpectrallyWeightedExposure(TH1F * hExposure, double a) {
    double AEFF_SpectWgt_nom=0,denom=0;
    for (int i=1;i<=hExposure->GetNbinsX();i++) {
       double EbinMax=pow(10,hExposure->GetXaxis()->GetBinUpEdge(i));
       double EbinMin=pow(10,hExposure->GetXaxis()->GetBinLowEdge(i));
       double weight;
       if (a!=-1) weight=(powl(EbinMax,a+1)-powl(EbinMin,a+1))/(a+1);
       else       weight=log(EbinMax)-log(EbinMin);
       float expo=hExposure->GetBinContent(i);
       if (expo<0) {
	    printf("%s: Something is wrong with the exposure value of bin %d (exp=%f) \n",__FUNCTION__,i,expo); 
	    for (int ii=1;ii<=hExposure->GetNbinsX();ii++) printf("%d %e\n",i,hExposure->GetBinContent(i));
	    exit (1);
       }
       AEFF_SpectWgt_nom+=expo*weight;
       denom+=weight;
    }
    if (AEFF_SpectWgt_nom<0 || denom<=0) {printf("%s: Something's wrong nom=%f denom=%f\n",__FUNCTION__,AEFF_SpectWgt_nom,denom); exit(1);}
    AEFF_SpectWgt_nom/=denom;
    return AEFF_SpectWgt_nom;
}


double TOOLS::CalcMeanEnergy(double MinEnergy, double MaxEnergy, double a) {
    double nom=0,denom=0;
    if (a!=-2) nom=(pow(MaxEnergy,a+2)-pow(MinEnergy,a+2))/(a+2);
    else       nom=log(MaxEnergy)-log(MinEnergy);

    if (a!=-1) denom=(pow(MaxEnergy,a+1)-pow(MinEnergy,a+1))/(a+1);
    else       denom=log(MaxEnergy)-log(MinEnergy);
    return nom/denom;
}


double TOOLS::CalcMeanEnergy(double MinEnergy,double PeakEnergy, double MaxEnergy, double a,double b) {

    if (MinEnergy>MaxEnergy || PeakEnergy<MinEnergy || PeakEnergy>MaxEnergy || a<=b ) {
       printf("%s: Check your input\n",__FUNCTION__);
       return -1;
    }
    double nom=0,denom=0;
    const int steps=200;
    const double dlE = (log10(MaxEnergy)-log10(MinEnergy))/steps;
    const double lBreak = log10((a-b)*PeakEnergy);
    for (int iE=0;iE<steps;iE++) {
        double lE=log10(MinEnergy)+(iE+0.5)*dlE;
        double E=pow(10,lE); //in keV
        double Band_Ne;
        if (lE<=lBreak) {
            Band_Ne=pow(E/100,a)*exp(-E/PeakEnergy);
        }
        else {
            Band_Ne=pow((a-b)*PeakEnergy/100,a-b)*exp(b-a)*pow(E/100,b);
        }
        double dE=pow(10,lE+0.5*dlE)-pow(10,lE-0.5*dlE);
        nom+=E*Band_Ne*dE;
        denom+=Band_Ne*dE;
        //printf("%d de=%f MeanE=%f band=%e dlE=%e\n",iE,dE,MeanE,Band_Ne,dlE);
    }
    nom/=denom;
    return nom;
}



#include "TGraph.h"


void TOOLS::PlotBand(float MinEnergy, float PeakEnergy, float MaxEnergy, float a,float b,double IntFlux) {

    const double keV2erg=1.60217646e-9;
    //const double erg2keV =1/keV2erg;

    double MeanEnergy_keV=CalcMeanEnergy(MinEnergy,PeakEnergy,MaxEnergy,a,b);

    if (MinEnergy>MaxEnergy || a<=b ) {
       printf("%s: Check your input\n",__FUNCTION__);
       return ;
    } 
    printf("%s: Integrated flux from %.1f to %.1f keV is %e erg/cm2/sec\n",__FUNCTION__,MinEnergy,MaxEnergy,IntFlux);
    double integral=0;
    const int steps=500;
    const double dlE = (log10(MaxEnergy)-log10(MinEnergy))/steps;
    const double lBreak = log10((a-b)*PeakEnergy);
    double lE_point[steps],nufnu_point[steps],nufnu_point_erg[steps];
    for (int iE=0;iE<steps;iE++) {
        double lE=log10(MinEnergy)+(iE+0.5)*dlE;
        double E=pow(10,lE); //in keV
        double E_erg=E*keV2erg;
        lE_point[iE]=E;
        double Band_Ne;
        if (lE<=lBreak && 0) {
            Band_Ne=pow(float(E/100.),float(a))*exp(-E/PeakEnergy);
        }
        else {
            Band_Ne=pow(float((a-b)*PeakEnergy/100.),float(a-b))*exp(b-a)*pow(float(E/100.),float(b));
        }
        double dE=pow(10,lE+0.5*dlE)-pow(10,lE-0.5*dlE);
        integral+=Band_Ne*dE;//in photon/cm2/sec
        nufnu_point[iE]=Band_Ne*E*dE;
        nufnu_point_erg[iE]=Band_Ne*E*E_erg;
        //printf("%d de=%f E=%e band=%e dlE=%e %e-%e\n",iE,dE,E,Band_Ne,dlE,pow(10,lE-0.5*dlE),pow(10,lE+0.5*dlE));
    }
    integral*=MeanEnergy_keV;
    double Norm_erg=IntFlux/(integral);
    double Norm    =IntFlux/integral;
    for (int iE=0;iE<steps;iE++) {
          nufnu_point_erg[iE]*=Norm_erg;
          nufnu_point[iE]*=Norm;
    }

    TCanvas * c = new TCanvas("cBand","cBand");
    c->Divide(1,2);
    c->GetPad(1)->SetLogx();
    c->GetPad(2)->SetLogx();
//    c->GetPad(1)->SetLogy();
//    c->GetPad(2)->SetLogy();

    TGraph * gBand_keV = new TGraph(steps,lE_point,nufnu_point);
    TGraph * gBand_erg = new TGraph(steps,lE_point,nufnu_point_erg);
    char name[1000];
    sprintf(name,"a=%.2f, b=%.2f, Peak Energy = %.0f keV Integral=%.2e erg/cm2/sec",a,b,PeakEnergy,IntFlux);
    gBand_keV->GetXaxis()->SetTitle("Energy (KeV)");
    gBand_erg->GetXaxis()->SetTitle("Energy (KeV)");
    gBand_keV->GetYaxis()->SetTitle("#nu F_{#nu} (keV cm^{-2} s^{-1}");
    gBand_erg->GetYaxis()->SetTitle("#nu F_{#nu} (erg cm^{-2} s^{-1}");
    gBand_keV->SetTitle(name);
    gBand_erg->SetTitle(name);
    c->cd(1); gBand_keV->Draw("AP");
    c->cd(2); gBand_erg->Draw("AP");


}

