//Author: Vlasios Vasileiou <vlasisva@gmail.com>
//$Header: /nfs/slac/g/glast/ground/cvs/GRBAnalysis-scons/BackgroundEstimator/BackgroundEstimator/BackgroundEstimator.h,v 1.5 2011/10/05 14:36:17 vlasisva Exp $
#ifndef _BackgroundEstimator_H
#define _BackgroundEstimator_H
#include "TProfile.h"
#include "include/BKGE_Tools.h"

//this enum helps with monitoring the different data formats
enum {DATA_FORMAT_P6_OLD,DATA_FORMAT_P6_NEW,DATA_FORMAT_P7};

/// LAT Background estimation for transient events

namespace BKGE_NS{
    ///Calculate the background map
    int CalculateBackground(string Interval_name, double MET, double DURATION, string FT1_FILE, string FT2_FILE, string DATACLASS, double Energy_Min_user, double Energy_Max_user, int Energy_Bins_user, int verbosity=1,bool Calc_Residual=true);
    ///Integrate the background map over the ROI
    string PlotBackground(string Interval_name, double MET, double DURATION, string FT1_FILE, string FT2_FILE, string DATACLASS, double Energy_Min_user, double Energy_Max_user, int Energy_Bins_user, bool OverwritePlots=true, int verbosity=1, double MET_FOR_THETA=-1);
    int MakeGtLikeTemplate(float gtlike_ROI, string GRB_DIR, string DATACLASS);
};


class BackgroundEstimator{
  struct Plots_Struct{
     TH1F* hPtRazvsTime,* hPtDeczvsTime,* hPtRaxvsTime,* hPtDecxvsTime,* hRAZenithvsTime,* hDecZenithvsTime;
     TH1F* hRockingAnglevsTime,* hMcIlwainLvsTime;
  } ;

  public:
    BackgroundEstimator(string aClass, double EMin=-1, double EMax=-1, int EBins=-1, bool ShowLogo=true);
    ~BackgroundEstimator();
   
    ///Calculate a background skymap. This is the first part of the bkg estimation
    int Make_Background_Map(string FT1_FILE, string FT2File, string GRB_DIR, double Burst_t0, double Burst_Dur,int verbosity=1, bool Calc_Residual=true, bool Save_Earth_Coo_Map=false); 

    ///Integrate a background map over the ROI to produce the final bkg estimate
    int FillBackgroundHist(string GRB_DIR, TH1F * hROI_Max, double RA_BURST, double DEC_BURST, short int type, int verbosity=0);

    ///Min and max energy in MeV of the datafiles
    double Energy_Min_datafiles, Energy_Max_datafiles;
    ///Number of logarithmically-spaced bins for the bkg estimate in the datafiles
    int Energy_Bins_datafiles;
    ///Min and max energy in MeV of the user
    double Energy_Min_user, Energy_Max_user;
    ///Number of logarithmically-spaced bins for the bkg estimate of the user
    int Energy_Bins_user;
    ///Zenith Theta cut -- can be used to override the default cut

    ///Flag that shows if the user is using the energy binning of the data files (set by the bkge)
    bool UsingDefaultBinning;

    ///Dataclass (e.g. P6_V3_TRANSIENT::FRONT)
    string DataClass;
    string DataClassName_noConv;
    string ConversionName;
    string DataClassVersion;
    int ConversionType;

    unsigned short int Energy2Bin(float Energy);
    float Bin2Energy(unsigned short int bin);

 private:
    void SimulateSky(Plots_Struct myPlots_Struct, TH2F * hSimulatedSky[], vector <double> GTI_Start, vector <double> GTI_End, const int nEnergy, TH2F* hSimulatedSky_Earth[]=0, float TimeStep_user=0,float hSimulatedSky_Earth_Map_Min_B=20);
    double GimmeCorrectionFactor(short int ie, double MET);

    bool PassesCuts(fitsfile * fptr, long int i, int format);

    ///Correction factors
    vector <TProfile*> pRatiovsTime;

    float EstimatorVersion,Residuals_version,RateFit_version,ThetaPhiFits_version,EastWest_version,TimeCorrectionFactors_version;
    double StartTime, EndTime, StopTime;
    long int TimeBins;
    float FT1ZenithTheta_Cut;
    int MinCTBClassLevel; //for P6
    int EventClassMask;   //for P7
    float BinSize;
    double TimeStep;
    int L_BINS, B_BINS; ///Number of longitude and latitude map bins 
    TFile * fResidualOverExposure,*fRateFit,*fThetaPhiFits,*fCorrectionFactors;

    char name[2000];
    string DataDir;
    ClassDef(BackgroundEstimator,1)

};

#endif

