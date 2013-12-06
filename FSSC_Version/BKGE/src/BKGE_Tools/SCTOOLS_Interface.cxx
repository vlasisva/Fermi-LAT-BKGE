#include "BackgroundEstimator/BKGE_Tools.h"

void TOOLS::Run_gtltcube(string GRB_DIR, double TMin, double TMax, string FT2_FILE, float FT1ZenithTheta_Cut, int verbosity, string gtltcube_Filename) {

  //gtltcube takes parameter binsize for SCTOOLS v<=r14 and binsz for v>=15
  //Have to choose which parameter to use -- I don't want to rely on parsing the version number of science tools
  //so I just check if there exists a parameter binsz
  char buffer[1000];char name[1000];
  buffer[0]='\0'; //empty char array
  sprintf(name,"plist gtltcube |grep binsize");
  FILE * pipe = popen(name,"r");
  char *cha =fgets(buffer,sizeof(buffer),pipe);
  pclose(pipe);
  //printf("--%s--\n",buffer);
  if (gtltcube_Filename=="") gtltcube_Filename="burst_ltCube.fits";
  if (strcmp(buffer,"")) sprintf(buffer,"gtltcube scfile=%s outfile=%s/%s dcostheta=0.025 binsize=1 zmax=%f chatter=4 phibins=2 tmin=%f tmax=%f evfile='' ",
     FT2_FILE.c_str(),GRB_DIR.c_str(),gtltcube_Filename.c_str(),FT1ZenithTheta_Cut,TMin,TMax);
  else                  sprintf(buffer,"gtltcube scfile=%s outfile=%s/%s dcostheta=0.025 binsz=1 zmax=%f chatter=4 phibins=2 tmin=%f tmax=%f evfile='' ",
     FT2_FILE.c_str(),GRB_DIR.c_str(),gtltcube_Filename.c_str(),FT1ZenithTheta_Cut,TMin,TMax);
 
 if (verbosity>3 ) { 
     pipe = popen(buffer, "r");
     printf("%s: Executing command: %s\n",__FUNCTION__,buffer);
     while (fgets(buffer,sizeof(buffer),pipe)) printf("|%s",buffer);
     pclose(pipe);
 }
  else {
     sprintf(buffer,"%s >/dev/null 2>&1",buffer);
     system(buffer);
 }

}

void TOOLS::Run_gtexpcube(string GRB_DIR,  double TMin, double TMax, string FT2_FILE, string DATACLASS, float FT1ZenithTheta_Cut, char * Outfile, float Energy_Min, float Energy_Max, int Energy_Bins, int verbosity, string gtltcube_Filename) {
 //first check if the ltcube exists
  char name[1000];
  if (gtltcube_Filename=="") gtltcube_Filename="burst_ltCube.fits";
  sprintf(name,"%s/%s",GRB_DIR.c_str(),gtltcube_Filename.c_str());
  FILE * ftemp  = fopen(name,"r");
  if (ftemp) fclose(ftemp);
  else Run_gtltcube(GRB_DIR, TMin, TMax, FT2_FILE, FT1ZenithTheta_Cut, verbosity, gtltcube_Filename);

  char buffer[1000];
  buffer[0]='\0'; //empty char array
     
  sprintf(buffer,"gtexpcube2 infile=%s/%s cmap=none outfile=%s irfs=%s nxpix=360 nypix=180 binsz=1 coordsys=GAL xref=0 yref=0 axisrot=0 proj=CAR emin=%f emax=%f enumbins=%d bincalc=CENTER chatter=4 ignorephi=yes thmax=180 thmin=0",
    	GRB_DIR.c_str(),gtltcube_Filename.c_str(),Outfile,DATACLASS.c_str(),Energy_Min,Energy_Max,Energy_Bins);

  if (verbosity>3) {
        FILE * pipe = popen(buffer, "r");
        printf("%s: Executing command: %s\n",__FUNCTION__,buffer);
        while (fgets(buffer,sizeof(buffer),pipe))   printf("|%s",buffer);
        pclose(pipe);
  }
  else  {
        sprintf(buffer,"%s >/dev/null 2>&1",buffer);
        system(buffer);
  }
}
