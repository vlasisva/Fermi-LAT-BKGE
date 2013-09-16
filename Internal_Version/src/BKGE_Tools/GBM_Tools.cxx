#include "BackgroundEstimator/BKGE_Tools.h"
#include "astro/JulianDate.h"
#include "astro/SolarSystem.h"
#include <climits>

//returns 0 if ok and -1 if failure
int TOOLS::ReadCTIME(string CTIME_File, short MinChannel, short MaxChannel, vector <unsigned short>& Counts, vector<double>& Time, vector <double>& Exposure, double &GBM_TSTART, double &GBM_TSTOP, int verbosity ) {
 
  FILE * ftemp = fopen(CTIME_File.c_str(),"r");
  if (ftemp) fclose(ftemp);
  else {
     printf("%s: Can't find %s\n",__FUNCTION__,CTIME_File.c_str()); //probably there are multiple version numbers (?)
     return -1;
  }

  fitsfile *fptr;
  int status = 0;
  int anynul,hdutype;
  long nrows;

  if (fits_open_file(&fptr,CTIME_File.c_str(), READONLY, &status)) {
      printf("%s: Error opening file %s\n",__FUNCTION__,CTIME_File.c_str());
      return -1;
  }

  fits_movabs_hdu(fptr, 3, &hdutype, &status);
  status=0;
  if (fits_get_num_rows(fptr, &nrows, &status)) fits_report_error(stderr, status);
  if (nrows<=0) {
      printf("%s:bad rows %ld in file %s\n",__FUNCTION__,nrows,CTIME_File.c_str());
      fits_close_file(fptr, &status);
      return -1;
  }

  char COMMENT[1000];
  char name[1000];
  fits_read_keyword(fptr, (char*)"TSTART",  name, COMMENT, &status); GBM_TSTART=atof(name);
  fits_read_keyword(fptr, (char*)"TSTOP",  name, COMMENT, &status); GBM_TSTOP=atof(name);
 

  //printf("%s: nrows=%ld\n",__FUNCTION__,nrows);
  Time.reserve(nrows);
  Counts.reserve(nrows);
  Exposure.reserve(nrows);

  for (long int jj = 1; jj <= nrows; jj++) {
        //if (jj%10000==0) {printf("%.2f\r",float(jj)/nrows); fflush(0);}

        float aexposure;
        if (fits_read_col (fptr,TFLOAT,2,jj, 1, 1, NULL,&aexposure, &anynul, &status)) { fits_report_error(stderr, status); return -1;}

        short quality;
        if (fits_read_col (fptr,TSHORT,3,jj, 1, 1, NULL,&quality, &anynul, &status)) { fits_report_error(stderr, status); return -1;}
        if (quality) {
            if (verbosity) printf("%s: skip bad quality %d bin =%ld\n",__FUNCTION__,quality,jj);
            continue;
        }

        short COUNTS[8];
        if (fits_read_col (fptr,TSHORT,1,jj, 1, 8, NULL,COUNTS, &anynul, &status)) { fits_report_error(stderr, status); return -1;}
        int sum=0;
        for (int iel=MinChannel;iel<=MaxChannel;iel++) sum+=COUNTS[iel];
        if (sum<=0) {
            //printf("%s: skip bad sum %d bin =%ld\n",__FUNCTION__,sum,jj);
            continue;
        }
        Exposure.push_back(aexposure);
        if (sum>=USHRT_MAX) {
            printf("%s: Counts overflow!! counts=%d max=%d\n",__FUNCTION__,sum,SHRT_MAX);
            sum=SHRT_MAX-1;
        }
        Counts.push_back(sum);

        double tstart;//,tstop;
        if (fits_read_col (fptr,TDOUBLE,4,jj, 1, 1, NULL,&tstart, &anynul, &status)) { fits_report_error(stderr, status); return -1;}
        Time.push_back(tstart);
  }
  fits_close_file(fptr, &status);

  if (Time.size()==0) {
      printf("%s: did not manage to read any time intervals (nrows=%ld)\n",__FUNCTION__,nrows);
      return -1;
  }

  return 0;

}

////////////////////////////////////////////////////////////////////////////////

void TOOLS::GetGBM_StartStop(string file, double &GBM_TSTART, double &GBM_TSTOP) {

  fitsfile *fptr;
  int status = 0;
  int hdutype;

  if (fits_open_file(&fptr,file.c_str(), READONLY, &status)) {
      printf("%s: Error opening file %s\n",__FUNCTION__,file.c_str());
      return;
  }

  fits_movabs_hdu(fptr, 3, &hdutype, &status);
  status=0;
  char COMMENT[1000];
  char name[1000];
  fits_read_keyword(fptr, (char*)"TSTART",  name, COMMENT, &status); GBM_TSTART=atof(name);
  fits_read_keyword(fptr, (char*)"TSTOP",  name, COMMENT, &status); GBM_TSTOP=atof(name);

  fits_close_file(fptr, &status);


}

