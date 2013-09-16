//Author: Vlasios Vasileiou <vlasisva@gmail.com>
// $Header: /nfs/slac/g/glast/ground/cvs/GRBAnalysis-scons/BackgroundEstimator/src/BackgroundEstimator/Make_ThetaPhi_Fits.cxx,v 1.3 2011/10/03 12:05:14 vlasisva Exp $
#include "BackgroundEstimator/BackgroundEstimator.h"

ClassImp(BackgroundEstimator)

void MakeFits(string DataDir, string DataDirs, int Energy_Bins_user, int nPhi, TFile * fout);

void BackgroundEstimator::Make_ThetaPhi_Fits(string FitsAllSkyFile){

  sprintf(name,"%s/ThetaPhi_Fits_%s_%.1f.root",DataDir.c_str(),DataClass.c_str(),ThetaPhiFits_version);
  FILE * ftemp = fopen(name,"r");
  TFile * fout;
  if (ftemp) {
     fout  = TFile::Open(name);
     fclose (ftemp);
     MakeFits(DataDir,DataClass,Energy_Bins_user,5,fout);
     fout->Close();
     return;
  }
  fout  = new TFile(name,"RECREATE");
     

  sprintf(name,"%s/Plots_%s.root",DataDir.c_str(),DataClass.c_str());
  TFile * fPlots = TFile::Open(name);
  TH1F* hPtRazvsTime = (TH1F*)fPlots->Get("hPtRazvsTime");
  TH1F* hPtDeczvsTime = (TH1F*)fPlots->Get("hPtDeczvsTime");
  TH1F* hPtRaxvsTime = (TH1F*)fPlots->Get("hPtRaxvsTime");
  TH1F* hPtDecxvsTime = (TH1F*)fPlots->Get("hPtDecxvsTime");
  TH1F* hRAZenithvsTime = (TH1F*)fPlots->Get("hRAZenithvsTime");
  TH1F* hDecZenithvsTime = (TH1F*)fPlots->Get("hDecZenithvsTime");

  TH2F * hThetaPhiOctant[Energy_Bins_user+1];
  TH1F * hTheta_away[Energy_Bins_user+1],*hTheta_towards[Energy_Bins_user+1];
  fout->cd();
  const float MinB=70;

  TH2F* hThetavsPhi_restricted = new TH2F("hThetavsPhi_restricted","ThetavsPhi with a ztheta cut",360,0,360,40,0,80);
  TH2F* hThetavsPhi_unrestricted = new TH2F("hThetavsPhi_ur","ThetavsPhi with no zt cut and looking away from the earth",180,0,360,40,0,80);
  TH1F* hZTheta_away = new TH1F("hZTheta_away","ZTheta_looking away from the earth",120,0,120);
  TH1F* hZTheta_towards = new TH1F("hZTheta_towards","hZTheta looking towards the earth",120,0,120);
  TH1F* hPhi_away = new TH1F("hPhi_away","Phi looking away from the earth",360,0,360);
  TH1F* hPhi_towards = new TH1F("hPhi_towards","Phi looking towards the earth",360,0,360);

  for (int i=1;i<=Energy_Bins_user;i++) {
     sprintf(name,"hThetaPhiOctant_%d",i);
     hThetaPhiOctant[i]=new TH2F(name,"hThetaPhiOctant",5,0,45,80,0,80);
     hThetaPhiOctant[i]->SetContour(256);
     
     sprintf(name,"hTheta_away_%d",i);
     hTheta_away[i] = new TH1F(name,"Theta looking away from the earth",90,0,90);//don't make the number of bins variable -- it will mess up the rate correction-factor calculation

     sprintf(name,"hTheta_towards_%d",i);
     hTheta_towards[i] = new TH1F(name,"Theta looking towards the earth",90,0,90);

  }
 
  fitsfile *fptr;
 
  Hep3Vector localdir;
  const double maxtheta=80*DEG_TO_RAD;
  const double maxZtheta=FT1ZenithTheta_Cut*DEG_TO_RAD;
  HepRotation localToCelestial;
  localToCelestial.setTolerance(0.1);
  astro::SkyDir SCZenith,SCz,SCx,locGAL;    

  ftemp = fopen(FitsAllSkyFile.c_str(),"r");
  int ifile=0;
  int ibin;
  while (fscanf(ftemp,"%s",name)==1) {
    printf("%d %s\r",ifile,name); fflush(0);
    ifile++;
    long nrows;int ncols;
    int status=0,hdutype,anynul;
    fits_open_file(&fptr, name, READONLY, &status);
    if (status) {printf("%s: error opening file %s\n",__FUNCTION__,name); exit(1);}
    fits_movabs_hdu(fptr, 2, &hdutype, &status);
    fits_get_num_rows(fptr, &nrows, &status);
    fits_get_num_cols(fptr, &ncols, &status);
    int format;
    if      (DataClass.find("P7")!=string::npos) format=DATA_FORMAT_P7;
    else if (DataClass.find("P6")!=string::npos){
        if      (ncols==22) format=DATA_FORMAT_P6_NEW;
        else if (ncols==18) format=DATA_FORMAT_P6_OLD;
        else {printf("%s: unknown format file %s ncols=%d class=%s\n",__FUNCTION__,name,ncols,DataClass.c_str()); exit(1);}
    }
    else {printf("%s: Unknown fits file format file %s class=%s\n",__FUNCTION__,name,DataClass.c_str()); exit(1);}

    for (long i=1;i<=nrows;i++) {
        if (i%100000==0) {printf("%ld \ %ld   \r",i,nrows); fflush(0);} 
    
        double PtTime;
        fits_read_col (fptr,TDOUBLE,10,i, 1, 1, NULL,&PtTime, &anynul, &status);
        
        int itimebin = hPtRazvsTime->FindBin(PtTime);

        if  (hPtRazvsTime->GetBinContent(itimebin)==0) {
           if      (hPtRazvsTime->GetBinContent(itimebin-1)!=0) ibin=itimebin-1;
           else if (hPtRazvsTime->GetBinContent(itimebin+1)!=0) ibin=itimebin+1;
           else    {printf("%s: there is a gap in the plots? %f %d\n",__FUNCTION__,PtTime,itimebin); continue;}
           itimebin=ibin;
        }
        double PtRaz  = hPtRazvsTime->GetBinContent(itimebin);
        double PtDecz = hPtDeczvsTime->GetBinContent(itimebin);
        double PtL,PtB;
        TOOLS::Galactic(PtRaz,PtDecz,&PtL,&PtB);
        
        if (fabs(PtB)>MinB && PassesCuts(fptr,i,format)) {
        
            double FT1ZenithTheta,FT1Theta,FT1Phi;
            fits_read_col (fptr,TDOUBLE,8,i, 1, 1, NULL,&FT1ZenithTheta, &anynul, &status);
            fits_read_col (fptr,TDOUBLE,6,i, 1, 1, NULL,&FT1Theta, &anynul, &status);
            fits_read_col (fptr,TDOUBLE,7,i, 1, 1, NULL,&FT1Phi, &anynul, &status);

            double FT1Energy;
            fits_read_col (fptr,TDOUBLE,1,i, 1, 1, NULL,&FT1Energy, &anynul, &status);
            short ebin=Energy2Bin(FT1Energy);
            if (ebin==0) {
               printf("ebin=0? Ft1energy=%f i=%d \n",FT1Energy,i);
               exit(1);
            } 
            float RAZENITH  = hRAZenithvsTime->GetBinContent(itimebin);
            float DECZENITH = hDecZenithvsTime->GetBinContent(itimebin);
            float PtRax     = hPtRaxvsTime->GetBinContent(itimebin);
            float PtDecx    = hPtDecxvsTime->GetBinContent(itimebin);
            SCZenith = astro::SkyDir(RAZENITH,DECZENITH,astro::SkyDir::EQUATORIAL);

            SCz = astro::SkyDir(PtRaz,PtDecz,astro::SkyDir::EQUATORIAL);
            SCx = astro::SkyDir(PtRax,PtDecx,astro::SkyDir::EQUATORIAL);
            localToCelestial = HepRotation(SCx.dir(),(SCz.dir()).cross(SCx.dir()),SCz.dir());

            localdir.setRThetaPhi(1,maxtheta,FT1Phi*DEG_TO_RAD);
            locGAL=astro::SkyDir(localToCelestial*localdir,astro::SkyDir::EQUATORIAL);

            //if an event with the same phi but with theta=maxtheta had a ztheta>maxZtheta then this event is towards the earth 
            if (SCZenith.difference(locGAL)>maxZtheta) {
                hThetavsPhi_restricted->Fill(FT1Phi,FT1Theta); 
                hZTheta_towards->Fill(FT1ZenithTheta);
                hPhi_towards->Fill(FT1Phi);
                hTheta_towards[ebin]->Fill(FT1Theta);
            }
            else {
               float PhiOffset=(int(FT1Phi/45)%2)?45-fmod(FT1Phi,45):fmod(FT1Phi,45);
               if (PhiOffset<0 || PhiOffset>45) printf("%f %f\n",PhiOffset,FT1Phi);
               hThetaPhiOctant[ebin]->Fill(PhiOffset,FT1Theta);
               hThetavsPhi_unrestricted->Fill(FT1Phi,FT1Theta); 
               hZTheta_away->Fill(FT1ZenithTheta);
               hTheta_away[ebin]->Fill(FT1Theta);
               hPhi_away->Fill(FT1Phi);
            } 
        }
       //if (i%1000==0) { printf("%.0f \r",100*i/(float)imax);  fflush(stdout); }
   }
   fits_close_file(fptr, &status);
  }    
  fclose (ftemp);
  printf("done reading files\n");
  for (int i=0;i<=hThetavsPhi_unrestricted->GetNbinsX();i++) {
     double sum1=0,sum2=0;
     for (int j=0;j<=hThetavsPhi_unrestricted->GetNbinsY();j++) {
       sum1+=hThetavsPhi_unrestricted->GetBinContent(i,j);
       sum2+=hThetavsPhi_restricted->GetBinContent(i,j);
     }

     for (int j=0;j<=hThetavsPhi_unrestricted->GetNbinsY();j++) {
         if (sum1) hThetavsPhi_unrestricted->SetBinContent(i,j,hThetavsPhi_unrestricted->GetBinContent(i,j)/sum1);
         if (sum2) hThetavsPhi_restricted->SetBinContent(i,j,hThetavsPhi_restricted->GetBinContent(i,j)/sum2);
     } 
  }
  
  printf("writing plots\n");
  for (int ie=1;ie<=Energy_Bins_user;ie++) {
     hThetaPhiOctant[ie]->Write();
     hTheta_away[ie]->Write();
     hTheta_towards[ie]->Write();
  }
  hPhi_towards->Write();
  hPhi_away->Write();
  hThetavsPhi_unrestricted->Write();
  hThetavsPhi_restricted->Write();
  hZTheta_away->Write();
  hZTheta_towards->Write();
  
  fPlots->Close();

  sprintf(name,"%.2f",ThetaPhiFits_version);
  TNamed * Version_TNamed = new TNamed("version",name);
  Version_TNamed->Write();
  
  MakeFits(DataDir,DataClass,Energy_Bins_user,5,fout);
  fout->Close();
}



void MakeFits(string DataDir, string DataClass, int Energy_Bins_user, int nPhi, TFile * fout){
  char name[1000];
   
  TCanvas * cc[Energy_Bins_user+1];
  
  fout->cd();
  
  TCanvas * c_report[Energy_Bins_user+1];
  for (int ie=1;ie<=Energy_Bins_user;ie++) {
  
     sprintf(name,"creport_%d",ie);
     c_report[ie]=new TCanvas(name,name);
     c_report[ie]->Divide(3,2);
      
     //get fine map 
     sprintf(name,"hThetaPhiOctant_%d",ie); 
     TH2F * hThetaPhi_Octant = (TH2F*)fout->Get(name); 
       
     c_report[ie]->cd(1);
     hThetaPhi_Octant->SetContour(256);
     hThetaPhi_Octant->Draw("COLZ");

     //make fine projections and plot them
     TH1F * hP_fine[nPhi+1];
     for (int iPhi=1;iPhi<=nPhi;iPhi++) {
         sprintf(name,"hProjection_fine_%d_%d",iPhi,ie);
         hP_fine[iPhi] = (TH1F*)hThetaPhi_Octant->ProjectionY(name,iPhi,iPhi);
         for (int itheta=1;itheta<=nPhi;itheta++) hP_fine[iPhi]->SetBinError(itheta,sqrt(hP_fine[iPhi]->GetBinContent(itheta)));
         hP_fine[iPhi]->SetLineColor(iPhi);
         hP_fine[iPhi]->SetMaximum(hThetaPhi_Octant->GetMaximum());
         c_report[ie]->cd(1+3);
         if (iPhi==1) hP_fine[iPhi]->Draw();
         else hP_fine[iPhi]->Draw("SAME");
     }
     
     //rebin fine map to increase bin statistics
     sprintf(name,"cc_%d",ie);
     cc[ie]=new TCanvas(name,name,1024,768);
     cc[ie]->Divide(3,2);
     
     TH1F * hP_coarse[nPhi+1],*hP_coarse_rescaled[nPhi+1];
     for (int iPhi=1;iPhi<=nPhi;iPhi++) {
        vector <double> xx,yy;
        xx.push_back(hP_fine[iPhi]->GetXaxis()->GetXmin());
        
        double sum=0;
        float theta_start=0;
        for (int iTheta=1;iTheta<=hThetaPhi_Octant->GetNbinsY();iTheta++) {
            sum+=hP_fine[iPhi]->GetBinContent(iTheta);
            float theta_end=hP_fine[iPhi]->GetXaxis()->GetBinUpEdge(iTheta);
            printf("ie=%d iTheta=%d %f-%f sum=%f \n",ie,iTheta,theta_start,theta_end,sum);
            float dtheta=theta_end-theta_start;
            if (sum && dtheta>2 &&
                ((1/sqrt(sum))<0.05 || 
                (theta_end>=50 && dtheta>10)||
                (theta_end<50 && dtheta>20))) {
                
                xx.push_back(hP_fine[iPhi]->GetXaxis()->GetBinUpEdge(iTheta));

                yy.push_back(sum);
                //printf("push xx=%f yy=%f\n",xx.back(),yy.back());
                sum=0;
                theta_start=hP_fine[iPhi]->GetXaxis()->GetBinUpEdge(iTheta);
            }
            
         }
         if (sum) {
            xx.push_back(hP_fine[iPhi]->GetXaxis()->GetXmax());
            yy.push_back(sum); 
            //printf("push xx=%f yy=%f\n",xx.back(),yy.back());
         }
         
         hP_coarse[iPhi]=new TH1F();
         hP_coarse[iPhi]->SetBins(xx.size()-1,&xx.front());
         hP_coarse_rescaled[iPhi]=new TH1F();
         hP_coarse_rescaled[iPhi]->SetBins(xx.size()-1,&xx.front());

         for (unsigned int i=1;i<xx.size();i++) {
            hP_coarse[iPhi]->SetBinContent(i,yy[i-1]);
            hP_coarse[iPhi]->SetBinError(i,sqrt(yy[i-1]));
            float theta_rad_0=hP_coarse[iPhi]->GetXaxis()->GetBinLowEdge(i)*DEG_TO_RAD;         
            float theta_rad_1=hP_coarse[iPhi]->GetXaxis()->GetBinUpEdge(i)*DEG_TO_RAD;         
            double dc=cos(theta_rad_0)-cos(theta_rad_1);
            //printf("dc=%f thetas=%f %f i=%d\n",dc,theta_rad_0/DEG_TO_RAD,theta_rad_1/DEG_TO_RAD,i);
            hP_coarse_rescaled[iPhi]->SetBinContent(i,yy[i-1]/dc);
            hP_coarse_rescaled[iPhi]->SetBinError(i,sqrt(yy[i-1])/dc);
         }
         hP_coarse_rescaled[iPhi]->Smooth(1);
         
         c_report[ie]->cd(2+3);
         sprintf(name,"hProjection_Coarse_%d_%d",iPhi,ie);
         hP_coarse[iPhi]->SetTitle(name);
         hP_coarse[iPhi]->SetName(name);
         hP_coarse[iPhi]->SetLineColor(iPhi);
         hP_coarse[iPhi]->SetMinimum(0);
         
         if (iPhi==1) hP_coarse[iPhi]->Draw();
         else         hP_coarse[iPhi]->Draw("SAME");
         
         cc[ie]->cd(iPhi);
         sprintf(name,"hProjection_Coarse_rescaled_%d_%d",iPhi,ie);
         hP_coarse_rescaled[iPhi]->SetTitle(name);
         hP_coarse_rescaled[iPhi]->SetName(name);
         hP_coarse_rescaled[iPhi]->SetLineColor(iPhi);
         hP_coarse_rescaled[iPhi]->SetMinimum(0);
         hP_coarse_rescaled[iPhi]->Draw();
         
         c_report[ie]->cd(3+3);
         if (iPhi==1) hP_coarse_rescaled[iPhi]->Draw();
         else         hP_coarse_rescaled[iPhi]->Draw("SAME");
         
     }
    
      c_report[ie]->Write();
      cc[ie]->Write();
    
  }
  
  

}

