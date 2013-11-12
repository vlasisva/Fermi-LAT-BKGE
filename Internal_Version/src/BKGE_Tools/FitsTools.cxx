#include "BackgroundEstimator/BKGE_Tools.h"

void TOOLS::ReadGTI(vector <double>& GTI_Starts, vector <double>& GTI_Ends, string FitsAllSkyFile, double StartTime, double EndTime) {


 FILE* ftemp = fopen(FitsAllSkyFile.c_str(),"r");
 if (!ftemp){ printf("%s: file %s does not exist\n",__FUNCTION__,FitsAllSkyFile.c_str()); exit(1);}
 
 fitsfile *fptr;
 vector <string> FitsFiles;
 if (FitsAllSkyFile.find(".fit")!=string::npos) FitsFiles.push_back(FitsAllSkyFile); //this is a fits file.. 
 else {
    char name[1000];
    while (fscanf(ftemp,"%s",name)==1) FitsFiles.push_back(name);
 }
 fclose(ftemp);
 
 for (int i=0;i<FitsFiles.size();i++) {    
    int status=0,hdutype,anynul;
    fits_open_file(&fptr, FitsFiles[i].c_str(), READONLY, &status);
    if (status) {printf("%s: error opening file %s\n",__FUNCTION__,FitsFiles[i].c_str()); exit(1);}

    //Read GTIs 
    long int GTIs;
    fits_movabs_hdu(fptr, 3, &hdutype, &status);
    fits_get_num_rows(fptr, &GTIs, &status);
    double agti0,agti1;
    for (long int i=1;i<=GTIs;i++) {
       fits_read_col (fptr,TDOUBLE,1,i, 1, 1, NULL,&agti0, &anynul, &status);
       fits_read_col (fptr,TDOUBLE,2,i, 1, 1, NULL,&agti1, &anynul, &status);
       if (agti1<StartTime || agti0>EndTime) continue; //skip external GTIs
       if (fabs(agti1-agti0)<=0.01) {printf("%s: Skipping short GTI %f-%f dt=%f\n",__FUNCTION__,agti0,agti1,agti1-agti0); continue;}
       if (StartTime && agti0<StartTime) agti0=StartTime; //advance start of GTI to match starting time
       if (EndTime && agti1>EndTime)     agti1=EndTime;
       GTI_Starts.push_back(agti0);
       GTI_Ends.push_back(agti1);
    }
    fits_close_file(fptr, &status);
 }

}
