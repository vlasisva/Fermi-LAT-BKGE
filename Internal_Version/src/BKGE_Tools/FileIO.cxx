//Author: Vlasios Vasileiou <vlasisva@gmail.com>
// $Header: /nfs/slac/g/glast/ground/cvs/GRBAnalysis-scons/BackgroundEstimator/src/BKGE_Tools/FileIO.cxx,v 1.1.1.1 2011/06/02 19:41:04 jrb Exp $
#include "BackgroundEstimator/BKGE_Tools.h"

//this script creates a file (named OFILENAME) that contains a list of fits files
//found in the FITSDIR that include the time range MET, MET+DURATION
vector<string> TOOLS::MakeFitsList(double MET, double DURATION, string FITSDIR, int verbosity) {

  vector <string> myList;
  if (verbosity) printf("%s: Making fitslist file for MET=%lf DUR=%lf, FITSDIR.c_str()=%s \n",__FUNCTION__,MET,DURATION,FITSDIR.c_str());
  string astring;
  char name[1000];
  vector <double> RUNS;
  sprintf(name,"ls -1 %s/r0*_ft1.fit",FITSDIR.c_str());
  FILE *pipe = popen(name, "r");
  while (fgets(name,sizeof(name),pipe)) {
     astring.assign(name);
     astring = astring.substr(1+astring.rfind("/"));
     double arun;
     sscanf(astring.c_str(),"r%10lf%*s",&arun);
     RUNS.push_back(arun);
  }

  unsigned int lasti=0;
  int nfound=0;
  for (unsigned int i=0;i<RUNS.size()-1;i++) {
      if (RUNS[i+1]<MET || RUNS[i]>(MET+DURATION)) continue;
      sprintf(name,"%s/r0%.0lf_ft1.fit",FITSDIR.c_str(),RUNS[i]);
      myList.push_back((string)name);
      lasti=i;
      nfound++;
      //printf("i=%d r0%lf_ft1.fit (next=%.0lf)\n",i,RUNS[i],RUNS[i+1]);
  }
  if (nfound==0 || lasti>RUNS.size()-1) {printf("Could not find enough fits files for the requested period\n"); exit(1);}
  sprintf(name,"%s/r0%.0lf_ft1.fit",FITSDIR.c_str(),RUNS[lasti+1]);
  myList.push_back((string)name);
  return myList;
}

