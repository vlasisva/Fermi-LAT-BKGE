//Author: Vlasios Vasileiou <vlasisva@gmail.com>
#include "include/BackgroundEstimator.h"
#include <sys/stat.h>
#include <sys/types.h>

//if Energy_Min_user,Energy_Max_user,Energy_Bins_user<=0 then their default values are used 
int BKGE_NS::CalculateBackground(string Interval_name, double MET, double DURATION, string FT1_FILE, string FT2_FILE, string DATACLASS, double Energy_Min_user, double Energy_Max_user, int Energy_Bins_user, int verbosity, bool Calc_Residual){

 static bool ShowLogo=true;
 if (verbosity==0) ShowLogo=false;

 string path = TOOLS::GetS("OUTPUT_DIR");
 mkdir(path.c_str(),S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
 string DataClassName=DATACLASS.substr(0,DATACLASS.find("::"));
 BackgroundEstimator * Est = new BackgroundEstimator(DATACLASS, Energy_Min_user, Energy_Max_user, Energy_Bins_user, ShowLogo);
  
 string GRB_DIR=TOOLS::GetS("OUTPUT_DIR")+ "/" + TOOLS::GetS("GRB_NAME"); mkdir(GRB_DIR.c_str(),S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
 GRB_DIR= GRB_DIR + "/"+ Interval_name; mkdir(GRB_DIR.c_str(),S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);

 if (verbosity>1) printf("%s (MET/DUR)=(%.1f %.1f) \n",Interval_name.c_str(),MET,DURATION);
 int result=Est->Make_Background_Map(FT1_FILE, FT2_FILE, GRB_DIR, MET, DURATION, verbosity,Calc_Residual);
 
 delete Est;
 ShowLogo=false; 
 return result;
}
