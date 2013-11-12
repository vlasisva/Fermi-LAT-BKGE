//Author: Vlasios Vasileiou <vlasisva@gmail.com>
#include "BackgroundEstimator/BackgroundEstimator.h"
#include <cmath>
#include <ctime>

//#define NO_EAST_WEST
//#define SAVE_DEBUG_MAPS

void BackgroundEstimator::SimulateSky(Plots_Struct myPlots_Struct, TH2F * hSimulatedSky[], vector <double> GTI_Start, vector <double> GTI_End, const int nEnergy, TH2F* hSimulatedSky_Earth[], float TimeStep_user, float hSimulatedSky_Earth_Map_Min_B) {
  //hSimulatedSky_Earth does not have a [0] element

  const float L_BinWidth=360./L_BINS;
  const float B_BinWidth=180./B_BINS;
  
  //Make a temp buffer skymap
  TH2F *htemp[nEnergy+1];
  TH2F *htemp_earth[nEnergy+1];
  for (int iEnergy=1;iEnergy<=nEnergy;iEnergy++) {
      sprintf(name,"htemp_%d",iEnergy);
      htemp[iEnergy] = new TH2F(name,name,L_BINS,-180,180,B_BINS,-90,90);

      if (hSimulatedSky_Earth[1]) {
         sprintf(name,"htemp_earth_%d",iEnergy);
         htemp_earth[iEnergy] = (TH2F*) hSimulatedSky_Earth[1]->Clone();
         htemp_earth[iEnergy]->SetTitle(name);
         htemp_earth[iEnergy]->SetName(name);
      }
  }

  //Get East-West Corrections
  static bool HaveEastWest=false;
  TH2F * hEastWest[nEnergy+1];
  TFile * fEastWest=NULL;
  sprintf(name,"%s/EastWest_Correction_%s_%.1f.root",DataDir.c_str(),DataClass.c_str(),EastWest_version);
  fEastWest = TFile::Open(name);
  if (fEastWest) {
     HaveEastWest=true;
     for (int iEnergy=1;iEnergy<=nEnergy;iEnergy++) {
        sprintf(name,"hEW_Correction_TrueOverEst_%d",iEnergy);
        hEastWest[iEnergy] = (TH2F*)fEastWest->Get(name);
        if  (!hEastWest[iEnergy]) { printf("%s: Can't find %s. Will not apply East West corrections\n",__FUNCTION__,name); HaveEastWest=false; break;}         
     }
  }
  

  //Get Ratefits
  sprintf(name,"%s/RateFit_%s_%.1f.root",DataDir.c_str(),DataClass.c_str(),RateFit_version);
  TFile * fRates = TFile::Open(name);

  TF1 * RateFit[nEnergy+1];
  for (int iEnergy=1;iEnergy<=nEnergy;iEnergy++) {
     sprintf(name,"RateFit_%d;1",iEnergy);
     RateFit[iEnergy] = (TF1*)fRates->Get(name);
  }
  TH1F * hScaleFactor =(TH1F*)fRates->Get("hScaleFactor");

  sprintf(name,"%s/ThetaPhi_Fits_%s_%.1f.root",DataDir.c_str(),DataClass.c_str(),ThetaPhiFits_version);
  TFile * fThetaPhi_Fits = TFile::Open(name);
  TH1F* hThetaPhi_rescaled[nEnergy+1][5+1];
  for (int iEnergy=1;iEnergy<=nEnergy;iEnergy++) {
     for (int iPhi=1;iPhi<=5;iPhi++) {
        char name2[100];
        sprintf(name2,"cc_%d",iEnergy);
        sprintf(name,"hProjection_Coarse_rescaled_%d_%d",iPhi,iEnergy);
        hThetaPhi_rescaled[iEnergy][iPhi] = (TH1F*)(((TCanvas*)fThetaPhi_Fits->Get(name2))->GetPad(iPhi)->FindObject(name));
        if (!hThetaPhi_rescaled[iEnergy][iPhi]) {printf("%s: Can't find %s in %s\n",__FUNCTION__,name,name2); exit(1);}
     }
  }
  
    
  CLHEP::HepRotation inverse;
  CLHEP::Hep3Vector LocalDir;
  astro::SkyDir local;
  
 
  HepRotation localToCelestial;
  localToCelestial.setTolerance(0.1);
  
  
  unsigned int igti=0;
  double TIME_START= myPlots_Struct.hMcIlwainLvsTime->GetXaxis()->GetXmin();
  double TIME_0  = TIME_START;
  double TIME_END= myPlots_Struct.hMcIlwainLvsTime->GetXaxis()->GetXmax();

  float TimeStep;
  if (TimeStep_user) TimeStep=TimeStep_user;
  else if ((TIME_END-TIME_START)<200) TimeStep=5;
  else TimeStep =30;
  
  #ifdef SAVE_DEBUG_MAPS
  TFile * fjunk = new TFile("/tmp/junk.root","RECREATE");
  #endif 

  
  //time_t start,end;
  //double dif;      
  //time (&start);  
  bool showed_rocking_angle_warning=false;
  
  while (1) {
      double fraction_done=(TIME_0-TIME_START)/(TIME_END-TIME_START);
      TOOLS::ProgressBar(int(fraction_done*40),40);
      /* 
      time (&end);

      dif = difftime (end,start);
      double remaining = dif/fraction_done-dif;
      printf("%f dif=%f  %f hrs remaining\r",fraction_done,dif,remaining/60./60.); fflush(0);
      */
      
      //Check if TIME_0 is in GTI
      if (fabs(TIME_0-TIME_END)<0.01) {
          //printf("TIME_0 reached TIME_END, done!\n");
          break;
      }

      if (TIME_0<GTI_Start[igti]) {
         TIME_0=GTI_Start[igti];
         //printf("advanced TIME_0 to the start of gti %d %f\n",igti,TIME_0);
      }
      bool BailOut=false;
      while (TIME_0>=GTI_End[igti]) { //TIME_0 not in gti
         if (igti==GTI_Start.size()-1) {
             //printf("TIME_0=%f out of GTIs. Current GTI %d %f-%f\n",TIME_0,igti,GTI_Start[igti],GTI_End[igti]);
             BailOut=true;
             break;
         }
         else {
            igti++;
            TIME_0=GTI_Start[igti];
            //printf("advanced GTI to %d TIME_0=%f GTI=%f-%f\n",igti,TIME_0,GTI_Start[igti],GTI_End[igti]);            
         }
      }
      if (BailOut) break;
      
      //find a TIME_1
      float aTimeStep=TimeStep;
      double TIME_1=TIME_0+aTimeStep;
      //printf("tentative TIME_1=%f timestep=%f\n",TIME_1,aTimeStep);
      if (TIME_1>GTI_End[igti]) {
         aTimeStep=GTI_End[igti]-TIME_0;
         TIME_1=TIME_0+aTimeStep;
         //printf("Adjusted TIME_1 to %f and timestep to %f (GTI_End=%f TIME_END=%f)\n",TIME_1,aTimeStep,GTI_End[igti],TIME_END);
      }
      
      if (aTimeStep<0.001) {
         TIME_0=TIME_1;
         continue;
      }
      
      int i_0=myPlots_Struct.hMcIlwainLvsTime->GetXaxis()->FindBin(TIME_0);
      int i_1=myPlots_Struct.hMcIlwainLvsTime->GetXaxis()->FindBin(TIME_1);
      int imid=(i_0+i_1)/2;
      
      //if (i_0%100==0) {printf("%.3f   %f - %f   \r",i_0/float(TimeBins),TIME_1,GTI_End[GTI_End.size()-2]);fflush(0);}
      //2.CALCULATE RATE    
      float RockingAngle = myPlots_Struct.hRockingAnglevsTime->GetBinContent(imid);
      if (RockingAngle<=0) {
             printf("%s: rocking angle weird rock=%f i_0=%d time=%f\n",__FUNCTION__,RockingAngle, i_0,TIME_0); 
             printf("%s: GTI data start\end: %f %f\n",__FUNCTION__,GTI_Start[igti],GTI_End[igti]);
             exit(1);
      }
      if (RockingAngle>70 && !showed_rocking_angle_warning) {
        printf("%s: WARNING - the rocking angle of the spacecraft exceeded 70degrees during your observation (ARR?). \n",__FUNCTION__);
        printf("For high rocking-angle observations, the contamination from the Earth limb becomes higher. \n");
        printf("Since the BKGE does not fully take into account such contamination, its final estimates may be erroneously lower\n");
        printf("See discussion in associated publication http://arxiv.org/abs/1307.4284\n");
        printf("To see for how long the rocking angle was at high values see the file Plots.root in the output directory\n");
        showed_rocking_angle_warning=true;
      }
    
      float McIlwainL=myPlots_Struct.hMcIlwainLvsTime->GetBinContent(imid);
      if (!McIlwainL) {
          //printf("%s: MCIlwainL=0? bin=%d time=%f\n",__FUNCTION__,i_0,TIME_0);
          exit(1);
      }
    
      float AllSkyRate[nEnergy+1];
      for (int iEnergy=1;iEnergy<=nEnergy;iEnergy++) {
          float ScaleFactor = hScaleFactor->GetBinContent(hScaleFactor->GetXaxis()->FindBin(RockingAngle),iEnergy);
          //printf("ie=%d rock=%f fac=%f\n",iEnergy,RockingAngle,ScaleFactor);
          if (ScaleFactor<=0 && RockingAngle<170) {printf("%s: Scalefactor=%f, rocking angle=%f TIME_0=%f i_0=%d gti=%d %f/%f\n",__FUNCTION__,ScaleFactor,RockingAngle,TIME_0,i_0,igti,GTI_Start[igti],GTI_End[igti]); exit(1);}
          //printf("%d ratefit=%f mcilwainl=%f scalefactor=%f rockingangle=%f\n",iEnergy,RateFit[iEnergy]->Eval(McIlwainL),McIlwainL,ScaleFactor,RockingAngle);
          AllSkyRate[iEnergy] = RateFit[iEnergy]->Eval(McIlwainL)*ScaleFactor;
      }      

      
      //printf("%d %d %f %f\n",i_0,TimeBins,ScaleFactor,AllSkyRate);
      //3.PREPARE VECTORS
      float PtRaz =  myPlots_Struct.hPtRazvsTime->GetBinContent(imid);
      float PtRax =  myPlots_Struct.hPtRaxvsTime->GetBinContent(imid);
      float PtDecz = myPlots_Struct.hPtDeczvsTime->GetBinContent(imid);
      float PtDecx = myPlots_Struct.hPtDecxvsTime->GetBinContent(imid);
      if (PtRaz==0 || PtRax==0 || PtDecz==0 || PtDecx==0){ printf("%s: there is a gap? %d %f %f %f %f \n",__FUNCTION__,imid,PtRaz,PtRax,PtDecz,PtDecx); exit(1);}
    
      astro::SkyDir SCz = astro::SkyDir(PtRaz,PtDecz,astro::SkyDir::EQUATORIAL);
      astro::SkyDir SCx = astro::SkyDir(PtRax,PtDecx,astro::SkyDir::EQUATORIAL);
      float PtBz_rad=SCz.b()*DEG_TO_RAD;          
      localToCelestial = HepRotation(SCx.dir(),(SCz.dir()).cross(SCx.dir()),SCz.dir());
          
      float RAZenith = myPlots_Struct.hRAZenithvsTime->GetBinContent(imid);
      float DecZenith= myPlots_Struct.hDecZenithvsTime->GetBinContent(imid);
      astro::SkyDir SCZenith = astro::SkyDir(RAZenith,DecZenith,astro::SkyDir::EQUATORIAL);
    
      //float RockingAngle=SCZenith.difference(SCz)/DEG_TO_RAD;

      if (RAZenith==0 || DecZenith==0) {printf("%s: FT2 has a gap! %d time=%f %f %f\n",__FUNCTION__,imid,myPlots_Struct.hRAZenithvsTime->GetBinCenter(imid),RAZenith,DecZenith); exit(1);}
      
      const float ThetaMax=80*DEG_TO_RAD;
      
      //Because we fill the east west maps using a FT1B>20degs cut, we need to keep a separate sum of the whole sky in this EW_Map_Integral variable
      //For the galactic coordinates map, there are no FT1B/L cuts, so we normalize them using their integral
      double EW_Map_Integral_low_lat[nEnergy];
      if (hSimulatedSky_Earth[1]) {
         for (int iEnergy=1;iEnergy<=nEnergy;iEnergy++) EW_Map_Integral_low_lat[iEnergy]=0;
      }
      
      for (int iB=1;iB<=B_BINS;iB++) {  
         float Bin_B=-90+(iB-0.5)*B_BinWidth;
         //htemp[1]->GetYaxis()->GetBinCenter(iB);
         
         if (fabs(Bin_B-SCz.b())>80) continue;
         int iL_Min,iL_Max;
             
         float Bin_B_rad=Bin_B*DEG_TO_RAD;
         
         double kernel=sqrt((pow(sin(ThetaMax/2.),2)-pow(sin((Bin_B_rad-PtBz_rad)/2.),2))/(cos(Bin_B_rad)*cos(PtBz_rad)));
         //printf("Pt %f/%f binB %f k=%f dl=%f\n",SCz.b(),SCz.l(),Bin_B,kernel,2*asin(kernel)/DEG_TO_RAD);
         if (isnan(kernel)) continue; //out of FOV
         else if (kernel>1) {        //the whole RA range is covered
              iL_Min=1;iL_Max=L_BINS;
         }
         else {   //in FOV
             float dd=2*asin(kernel)/DEG_TO_RAD;
             if (dd>180) {printf("dd>180 %f\n",dd); exit(1);}
             float L_Min=SCz.l()-dd;
             float L_Max=SCz.l()+dd;

             
             if (L_Min<-180)L_Min+=360;
             if (L_Min>180) L_Min-=360;
             if (L_Max>180) L_Max-=360;
             if (L_Max<-180)L_Max+=360;
             //if (L_Min>L_Max) {float temp=L_Min; L_Min=L_Max; L_Max=temp;}
             iL_Min=1+int((L_Min+180)/L_BinWidth);
             //=htemp[1]->GetXaxis()->FindBin(L_Min);
             iL_Max=1+int((L_Max+180)/L_BinWidth);
             //=htemp[1]->GetXaxis()->FindBin(L_Max);
             if (iL_Min>iL_Max) iL_Max+=L_BINS;
             //printf("%f %f %d %d\n",L_Min,L_Max,iL_Min,iL_Max);
         }
         float Bin_SolidAngle=cos(Bin_B*DEG_TO_RAD);
         for (int aiL=iL_Min;aiL<=iL_Max;aiL++) {
              //this trick with iL and aiL is for cases where the goodrange passes the L=+180degs border.
              int iL=aiL;
              if (iL>L_BINS) iL-=L_BINS;
                  
              float Bin_L=-180+(iL-0.5)*L_BinWidth;              
              //=htemp[1]->GetXaxis()->GetBinCenter(iL);

              astro::SkyDir SCBin = astro::SkyDir(Bin_L,Bin_B,astro::SkyDir::GALACTIC);
                  
              float BinZTheta_rad=SCBin.difference(SCZenith);
              float BinZTheta=BinZTheta_rad/DEG_TO_RAD;
              
              if (BinZTheta>FT1ZenithTheta_Cut) continue;
                            
              
              astro::PointingTransform p = astro::PointingTransform(SCz,SCx);
              //intf("%f %f\n",FT1Theta*DEG_TO_RAD,FT1Phi*DEG_TO_RAD);

              CLHEP::Hep3Vector dir = SCBin();
              CLHEP::HepRotation inverse = inverseOf(p.localToCelestial());
              astro::SkyDir LocalDir = astro::SkyDir(inverse*SCBin());

              float yy=LocalDir().y();
              float xx=LocalDir().x();
              float BinPhi=atan2(yy,xx)*RAD_TO_DEG;
              if (BinPhi<0) BinPhi+=360;
              float BinTheta=acos(LocalDir().z())/DEG_TO_RAD;
              if (BinTheta>84 || BinPhi<0 || BinPhi>360){  
                 printf("bintheta=%f BinPhi=%f L=%f/%d PtL/B=%f/%f \n",BinTheta,BinPhi,Bin_L,iL,SCz.l(),SCz.b()); 
                 BinTheta=80;
              }
              /////////////////////////////////////////////////
                                          
              //float BinTheta_rad=SCBin.difference(SCz);
              
              int iPhi=1+int(((int(BinPhi/45)%2)?45-fmod(BinPhi,45):fmod(BinPhi,45))/9);
              if (iPhi==6) iPhi=5;
              if (iPhi<=0 || iPhi>6) {printf("%f %d\n",BinPhi,iPhi); iPhi=1;}
                           
              float x=0,y=0;
              int EW_Bin=0;//This is the bin of the east west map

              if (HaveEastWest || hSimulatedSky_Earth[1]) {
                  float EarthAzimuth=DEG_TO_RAD*TOOLS::GimmeEarthAzimuth(Bin_L,Bin_B,SCZenith.l(),SCZenith.b(), astro::SkyDir::GALACTIC);
                  x=sin(EarthAzimuth)*BinZTheta;
                  y=cos(EarthAzimuth)*BinZTheta;
                  if (HaveEastWest) EW_Bin = hEastWest[1]->FindBin(x,y);
              }
              //The above assumes that hEastWest and hSimulatedSky_Earth have all the same bins in all energies
            
              
              for (int iEnergy=1;iEnergy<=nEnergy;iEnergy++) {
                       int itheta_phi_bin=hThetaPhi_rescaled[iEnergy][iPhi]->FindBin(BinTheta);
                       //rate = relative rate 
                       float rate=hThetaPhi_rescaled[iEnergy][iPhi]->GetBinContent(itheta_phi_bin)*Bin_SolidAngle;
                       //if (iEnergy==1) printf("L/B=%f/%f rate=%f BinTheta/Phi= %f/%f thetafit=%f phifit=%f area_fac=%f iphibin=%d\n",Bin_L,Bin_B,rate,BinTheta,BinPhi,ThetaFit[iphibin][iEnergy]->Eval(BinTheta),PhiFit[iEnergy]->Eval(BinPhi),cos(Bin_B*DEG_TO_RAD),iphibin);
                       if (rate<0) continue;
                       
                       #ifndef NO_EAST_WEST
            	       if (HaveEastWest) {      
                           float EWCorrection = hEastWest[iEnergy]->GetBinContent(EW_Bin);
                           if (EWCorrection!=0) rate*=EWCorrection;
                       }
                       #endif
                       
                       //The simulated earth skymap is EW-Corrected (if applicable)
                       if (hSimulatedSky_Earth[1]) { //fill EW simulated map
                              if (fabs(Bin_B)> hSimulatedSky_Earth_Map_Min_B) {
                                  htemp_earth[iEnergy]->Fill(x,y,rate);
                              }
                              else EW_Map_Integral_low_lat[iEnergy]+=rate;
                       }
                                  
                       htemp[iEnergy]->SetBinContent(iL,iB,rate);
                       //htemp[iEnergy]->SetBinContent(iL,iB,iPhi);
                     
                       //printf("%d %f %d %d\n",iEnergy,BinTheta,iPhi,itheta_phi_bin);
             }
          }
      }

      for (int iEnergy=1;iEnergy<=nEnergy;iEnergy++) {
          //printf("allskyrate=%f atimestep=%f integral=%f\n",AllSkyRate[iEnergy],aTimeStep,htemp[iEnergy]->Integral());
          #ifdef SAVE_DEBUG_MAPS
          sprintf(name,"htemp_%d_%d",iEnergy,i_0); htemp[iEnergy]->Write(name);
          if (hSimulatedSky_Earth[1]) {sprintf(name,"htemp_earth_%d_%d",iEnergy,i_0); htemp_earth[iEnergy]->Write(name);}
          #endif
          htemp[iEnergy]->Scale(AllSkyRate[iEnergy]*aTimeStep/htemp[iEnergy]->Integral());          
          if (hSimulatedSky_Earth[1]){ htemp_earth[iEnergy]->Scale(AllSkyRate[iEnergy]*aTimeStep/(EW_Map_Integral_low_lat[iEnergy]+htemp_earth[iEnergy]->Integral()));}

          
          #ifdef SAVE_DEBUG_MAPS
          sprintf(name,"htemp_scaled_%d_%d",iEnergy,i_0); htemp[iEnergy]->Write(name);
          if (hSimulatedSky_Earth[1]) {sprintf(name,"htemp_earth_scaled_%d_%d",iEnergy,i_0); htemp_earth[iEnergy]->Write(name);}
          #endif
                    
          hSimulatedSky[iEnergy]->Add(htemp[iEnergy]);
          if (hSimulatedSky_Earth[1]) hSimulatedSky_Earth[iEnergy]->Add(htemp_earth[iEnergy]);
          
          #ifdef SAVE_DEBUG_MAPS
          sprintf(name,"hsim_%d_%d",iEnergy,i_0); hSimulatedSky[iEnergy]->Write(name);
          if (hSimulatedSky_Earth[1]) {sprintf(name,"hsim_earth_%d_%d",iEnergy,i_0); hSimulatedSky_Earth[iEnergy]->Write(name);}
          #endif
          //printf("hsim integral=%f\n",hSimulatedSky->Integral());
        
          htemp[iEnergy]->Reset();
          if (hSimulatedSky_Earth[1]) htemp_earth[iEnergy]->Reset();
      }
      TIME_0=TIME_1;
  } //for timebins
  
  for (int iEnergy=1;iEnergy<=nEnergy;iEnergy++) {
     htemp[iEnergy]->Delete();
     if (hSimulatedSky_Earth[1]) htemp_earth[iEnergy]->Delete();
  }
  fRates->Close();
  fThetaPhi_Fits->Close();
  if (fEastWest) fEastWest->Close();
  #ifdef SAVE_DEBUG_MAPS
  fjunk->Close();
  #endif
}



