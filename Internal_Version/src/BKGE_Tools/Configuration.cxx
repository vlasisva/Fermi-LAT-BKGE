// Author: Vlasios Vasileiou <vlasisva@gmail.com>
// $Header: /nfs/slac/g/glast/ground/cvs/GRBAnalysis-scons/BackgroundEstimator/src/BKGE_Tools/Configuration.cxx,v 1.2 2011/09/23 16:30:10 vlasisva Exp $
#include "BackgroundEstimator/BKGE_Tools.h"

void TOOLS::PrintConfig() {
  TDirectory * DataDir = gROOT->GetDirectory("/DataDir");
  if (!DataDir) {printf("%s: Can't find DataDir!\n",__FUNCTION__); return;}
  DataDir->ls();
}

double TOOLS::Get(string name) {
  TDirectory * DataDir = gROOT->GetDirectory("/DataDir");
  if (!DataDir) {printf("Can't find DataDir!\n"); exit(1);}
  TNamed * t=((TNamed*)DataDir->FindObject(name.c_str()));
  if (!t) {printf("%s:Can't find TNamed %s \n",__FUNCTION__,name.c_str()); exit(1);}
  return atof(t->GetTitle());
}


void TOOLS::Set(string name,double val) {
  TDirectory * DataDir = gROOT->GetDirectory("/DataDir");
  if (!DataDir) DataDir= new TDirectory("DataDir","DataDir");
  TNamed * t=((TNamed*)DataDir->FindObject(name.c_str()));
  char vals[100];
  sprintf(vals,"%f",val);
  if (!t) {
     t = new TNamed(name,vals);
     DataDir->Add(t);
  }
  else t->SetTitle(vals);
}


void TOOLS::Set(string name,string vals) {
  TDirectory * DataDir = gROOT->GetDirectory("/DataDir");
  if (!DataDir) DataDir= new TDirectory("DataDir","DataDir");
  TNamed * t=((TNamed*)DataDir->FindObject(name.c_str()));
  if (!t) {
     t = new TNamed(name,vals);
     DataDir->Add(t);
  }
  else t->SetTitle(vals.c_str());
  if (name=="DATA_CLASS") {
       //For P7 and !transient continue
       if (!(vals.find("P7")!=string::npos and vals.find("TRANSIENT")==string::npos )) {
          Set("_DataClassName_noConv",GetDataClassName_noConv(vals));
          Set("_ConversionName",      GetConversionName(vals));
          Set("_DataClassVersion",    GetDataClassVersion(vals));
          Set("_CTBClassLevel",    GetCTBClassLevel(vals));       
       }
  }     

}


string TOOLS::GetS(string name) {
  TDirectory * DataDir = gROOT->GetDirectory("/DataDir");
  if (!DataDir) {printf("Can't find DataDir!\n"); exit(1);}
  TNamed * t=((TNamed*)DataDir->FindObject(name.c_str()));
  if (!t) {printf("%s:Can't find TNamed %s \n",__FUNCTION__,name.c_str()); exit(1);}
  return t->GetTitle();
}

void TOOLS::WriteConfig(string ConfigFile) {
    TDirectory * DataDir = gROOT->GetDirectory("/DataDir");
    TFile * fout = new TFile(ConfigFile.c_str(),"RECREATE");
    DataDir->Write();
    fout->Close();
}

void TOOLS::LoadConfig(string ConfigFile, bool ls) {
 TDirectory * DataDir = gROOT->GetDirectory("DataDir");
 if (!DataDir) DataDir= new TDirectory("DataDir","DataDir");

 char name[1000];
 sprintf(name,"%s/ConfigFiles/%s",(TOOLS::GetS("GANGSTER_DIR")).c_str(),ConfigFile.c_str());
 FILE * ftemp = fopen(name,"r");
 if (!ftemp) {printf("%s: Can't open file %s for parsing\n",__FUNCTION__,name); return;}
 
 char varname[1000],val[1000],line[1000];
 TNamed * named;
 while (fgets(line,sizeof(line),ftemp)) {
    if (line[0]=='#') continue;
    if (sscanf(line,"%s %s",varname,val)!=2) continue;
    if (!strcmp(varname,"") || !strcmp(val,"")) continue;
    if (!strcmp(varname,"INCLUDE")) {LoadConfig(val,false); continue;}

    named = (TNamed*)DataDir->FindObject(varname);
    if (named) {
       if (strcmp(named->GetTitle(),val)) {
          printf("%s:Resetting %s to %s from %s\n",__FUNCTION__,varname,named->GetTitle(),val);
          named->SetTitle(val);
       }
    }
    else {
       named=new TNamed(varname,val);
       DataDir->Add(named);
    }
    if (!strcmp(varname,"DATA_CLASS")) {
       Set("_DataClassName_noConv",GetDataClassName_noConv(val));
       Set("_ConversionName",      GetConversionName(val));
       Set("_DataClassVersion",    GetDataClassVersion(val));
       Set("_CTBClassLevel",    GetCTBClassLevel(val));
    }
 }
 if (ls) DataDir->ls();
}


