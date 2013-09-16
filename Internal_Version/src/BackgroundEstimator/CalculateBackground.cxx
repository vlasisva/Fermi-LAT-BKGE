//Author: Vlasios Vasileiou <vlasisva@gmail.com>
// $Header: /nfs/slac/g/glast/ground/cvs/GRBAnalysis-scons/BackgroundEstimator/src/BackgroundEstimator/CalculateBackground.cxx,v 1.1 2011/10/05 14:36:17 vlasisva Exp $
#include "BackgroundEstimator/BackgroundEstimator.h"
#include <sys/stat.h>
#include <sys/types.h>

//if Energy_Min_user,Energy_Max_user,Energy_Bins_user<=0 then their default values are used 
int BKGE_NS::CalculateBackground(string Interval_name, double MET, double DURATION, string FT1_FILE, string FT2_FILE, string DATACLASS, double Energy_Min_user, double Energy_Max_user, int Energy_Bins_user, int verbosity, bool Calc_Residual){

 static bool ShowLogo=true;
 if (verbosity==0) ShowLogo=false;

 string path = TOOLS::GetS("OUTPUT_DIR")+"/Bkg_Estimates/";
 mkdir(path.c_str(),S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
 string DataClassName=DATACLASS.substr(0,DATACLASS.find("::"));
 BackgroundEstimator * Est = new BackgroundEstimator(DATACLASS, Energy_Min_user, Energy_Max_user, Energy_Bins_user, true,ShowLogo);
  

 string GRB_DIR;
 GRB_DIR=TOOLS::GetS("OUTPUT_DIR")+"/Bkg_Estimates/"+Interval_name+"/";
 mkdir(GRB_DIR.c_str(),S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
 if (verbosity>1) printf("%s (MET/DUR)=(%.1f %.1f) \n",Interval_name.c_str(),MET,DURATION);
 int result=Est->Make_Background_Map(FT1_FILE, FT2_FILE, GRB_DIR, MET, DURATION, verbosity,Calc_Residual);
 
 delete Est;
 ShowLogo=false; 
 return result;
}
