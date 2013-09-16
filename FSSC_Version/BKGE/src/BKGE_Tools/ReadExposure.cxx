//Author: Vlasios Vasileiou <vlasisva@gmail.com>
#include "fitsio.h" 
#include "include/BKGE_Tools.h"
//#define DEBUG

#define buffsize 362
void printerror( int status);

void TOOLS::ReadExposureMap(string ExposureFilename, TH2F * hExposure, int ie, const short int verbosity){

    hExposure->Reset();
    fitsfile *fptr;       /* pointer to the FITS file, defined in fitsio.h */
    int status=0,  nfound, anynull;
    long  fpixel, npixels, naxes[2];

    double nullval, buffer[buffsize];

    TH2D * htemp = new TH2D("htemp","htemp",360,0,360,180,-90,90);
#ifdef DEBUG
    TH2D * hArea = (TH2D*) hExposure->Clone("hArea");
    hArea->SetTitle("Area per bin");
#endif

    if ( fits_open_file(&fptr,ExposureFilename.c_str() , READONLY, &status) ) printerror( status );

    /* read the NAXIS1 and NAXIS2 keyword to get image size */
    char name[]="NAXIS";
    if ( fits_read_keys_lng(fptr,name, 1, 2, naxes, &nfound, &status) )  printerror( status );
    if (naxes[0]>=buffsize) {printf("%s: Please increase buffer size\n",__FUNCTION__);exit(1);}
    npixels  = naxes[0] * naxes[1];         /* number of pixels in the image */
    fpixel   = 1+(ie-1)*npixels;
    nullval  = 0;                /* don't check for null values in the image */
    for (int icol=0;icol<naxes[1];icol++) {
      if ( fits_read_img(fptr, TDOUBLE, fpixel, naxes[0], &nullval, buffer, &anynull, &status) ) printerror( status );
      for (int ii = 0; ii < naxes[0]; ii++)  {
         //printf("%d %d %f\n",ii,icol,buffer[ii]);
         //hExposure->SetBinContent(-ii+naxes[0],1+icol,buffer[ii]);
         htemp->SetBinContent(-ii+naxes[0],1+icol,buffer[ii]);
      }

      fpixel  += naxes[0];    /* next pixel to be read in image */
    }
    if ( fits_close_file(fptr, &status) )   printerror( status );
    for (int i=1;i<=hExposure->GetNbinsX();i++) {
       float L=hExposure->GetXaxis()->GetBinCenter(i)+180;
       while (L<0) L+=360;
       for (int j=1;j<=hExposure->GetNbinsY();j++) {
          float B=hExposure->GetYaxis()->GetBinCenter(j);
          int bin=htemp->FindBin(L,B);
          hExposure->SetBinContent(i,j,htemp->GetBinContent(bin));
       }
    }

    delete htemp;
}


void printerror( int status){

    if (status)
    {
       fits_report_error(stderr, status); /* print error report */
       exit( status );    /* terminate the program, returning error status */
    }
    return;
}

