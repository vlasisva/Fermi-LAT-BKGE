//Author: Vlasios Vasileiou <vlasisva@gmail.com>
// $Header: /nfs/slac/g/glast/ground/cvs/GRBAnalysis-scons/BackgroundEstimator/src/BKGE_Tools/tools.cxx,v 1.2 2011/09/14 14:57:11 vlasisva Exp $
#include "BackgroundEstimator/BKGE_Tools.h"

void TOOLS::ProgressBar(const short int current, const short int max) {
  string out="[";
  for (int ii=0;ii<=current;ii++)    out+='*';
  for (int ii=current;ii<max-1;ii++) out+=' ';
  out+=']';
  printf("%s\r",out.c_str());
  fflush(stdout);
}

int TOOLS::GetClassMask(string DataClass) { //this is for the bitfield operations in P7
  if       (DataClass.find("TRANSIENT")!=string::npos) return 4;
  else if  (DataClass.find("SOURCE")  !=string::npos)  return 8;
  else     {printf("%s:Unknown class %s\n",__FUNCTION__,DataClass.c_str()); exit(1);}
}

int TOOLS::GetCTBClassLevel(string DataClass) {
  if       (DataClass.find("TRANSIENT")!=string::npos) return 1;
  else if  (DataClass.find("DIFFUSE")  !=string::npos) return 3;
  else     {printf("%s:Unknown class %s\n",__FUNCTION__,DataClass.c_str()); exit(1);}
}

string TOOLS::GetDataClassName_noConv(string DataClass) {
  if       (DataClass.find("TRANSIENT")!=string::npos) return "TRANSIENT";
  else if  (DataClass.find("DIFFUSE")!=string::npos)   return "DIFFUSE";
  else if  (DataClass.find("SOURCE")!=string::npos)    return "SOURCE";
  else {printf("%s:Unknown class %s\n",__FUNCTION__,DataClass.c_str()); exit(1);}
}

int TOOLS::GetConversionType(string DataClass) {
  if      (DataClass.find("FRONT")!=string::npos) return 0;
  else if (DataClass.find("BACK")!=string::npos)  return 1;
  else                                            return -1;
}

string TOOLS::GetConversionName(string DataClass) {
  if      (DataClass.find("FRONT")!=string::npos) return "FRONT";
  else if (DataClass.find("BACK")!=string::npos)  return "BACK";
  else                                            return "BOTH";
}

string TOOLS::GetDataClassVersion(string DataClass) {
    //see if we have P6 or P7 first
    if (DataClass.find("P7")!=string::npos) { //P7
       return DataClass.substr(DataClass.find("_")+1,2);
    }
    else { //P6
       return DataClass.substr(0,DataClass.find("_",3)+1);
    }
}





double TOOLS::GimmeLumDistance(double z) {

  double Wm=0.27,Wl=0.73,Wr=0,//,Wr=0,
         H0=71*3.240779289e-20, //H0 in 1/sec
         Wk=1-Wm-Wl-Wr;
  double c=3e8;//m/sec

  double ZInt=0;
  double dz=z/10000.;
  for (double zz=0;zz<z;zz+=dz) {
     ZInt+=dz/sqrt(Wr*pow(1+zz,4)+Wm*pow(1+zz,3)+Wk*pow(1+zz,2)+Wl);
  }
  double LumDist=c/H0*(1+z)*ZInt;
 return LumDist;
}


