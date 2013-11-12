//Author: Vlasios Vasileiou <vlasisva@gmail.com>
//$Header: /nfs/slac/g/glast/ground/cvs/GRBAnalysis-scons/BackgroundEstimator/BackgroundEstimator/BKGE_Tools.h,v 1.6 2013/10/25 12:51:12 vlasisva Exp $

#ifndef _BKGE_Tools_H
#define _BKGE_Tools_H

#include "TTree.h"
#include "TH1F.h"
#include "TH1I.h"
#include "TFile.h"
#include <cmath>
#include <cstdlib>
#include "TF1.h"
#include "TH2F.h"
#include "astro/PointingTransform.h"
#include "astro/SkyDir.h"
#include "TMath.h"
#include "TRandom.h"
#include "TStyle.h"
#include "TROOT.h"
#include "TCanvas.h"
#include <cstring>
#include "fitsio.h"
#include <vector>
#include "TCanvas.h"
const double RAD_TO_DEG=57.2957795130;
const double DEG_TO_RAD=0.01745329255;


/*! \brief Miscellaneous support tools.
*/
using namespace std;
using namespace CLHEP;

namespace TOOLS{

    /*! \brief Fill hROI with sqrt(PSF^2+LocalizationError^2) -- calculates theta & phi automatically
       * It calculates the 'Containment' radius of the PSF given by the DATA_CLASS
       * configuration parameter evaluated at the (theta,phi) of the source at time MET.\n\n
       * The Localization error is added in quadrature. The maximum radius is MaxRadius.
       * The results are filled in the hROI histogram.\n\n
       * If values for MET, Containment and MaxRadius are not given, then the values from the currently loaded configuration are used.
    */
    void CalculatePSF(TH1F * hROI ,double MET, string FT2_FILE, string DATACLASS, float Containment=-1, float MaxRadius=-1);

    /*! \brief Fill hROI with sqrt(PSF^2+LocalizationError^2) -- for given theta & phi
       * It calculates the 'Containment' radius of the PSF given by the DATA_CLASS
       * configuration parameter evaluated at the (theta,phi) of the source at time MET.\n\n
       * The Localization error is added in quadrature. The maximum radius is MaxRadius.
       * The results are filled in the hROI histogram.\n\n
       * If values for MET, Containment and MaxRadius are not given, then the values from the currently loaded configuration are used.
    */
    void CalculatePSF_ThetaPhi(TH1F * hROI, float theta, float phi, string DATACLASS, float Containment=-1, float MaxRadius=-1);

    /*! \brief Fill hROI with ROI data from a file
       * The file should have the same number of lines as the number of bins of hROI.\n\n
       * Format: one ROI size per line, and in degrees.
    */
    int ReadROI_File(TH1F * hROI, string filename);

    ///Show a prgress bar filled with a fraction current/max
    void ProgressBar(const short int max, const short int current);

    //Upper Limits helpers
    ///Return the spectrally weighted exposure from hExposure assuming an index a (dN/dE=E0*E^-a)
    double CalcSpectrallyWeightedExposure(TH1F * hExposure, double a);
    ///Return the average energy for a spectrum dN/dE=E0*E^-a
    double CalcMeanEnergy(double MinEnergy, double MaxEnergy, double a);
    ///Return the average energy for a Band spectrum. Give energies in keV
    double CalcMeanEnergy(double MinEnergy, double PeakEnergy, double MaxEnergy, double a, double b);

    //Map Manipulations/Skymaps/etc
    ///Read the exposure map from a fits file and save it in TH2F hExposure
    void ReadExposureMap(string ExposureFilename, TH2F * hExposure, int ie, const short int verbosity);
    ///Average the exposure over an ROI and save it in hExposure
    void CalcExposure(TFile * fResults, float L_BURST, float B_BURST, float FT1ZenithTheta_Cut, TH1F * hExposure, string GRB_DIR, int verbosity);
    int Make_Burst_Plots(string DataClass, string FT1_FILE, string GRB_DIR, double RA_BURST, double DEC_BURST, double GRB_t0, double Burst_Dur, TH1F* hROI, short int verbosity=1);

    ///Integrate a skymap that shows ev/sr over an ROI to produce events in the ROI
    double Integrate(TH2F * hMap, float L_BURST, float B_BURST, float ROI_RADIUS, string MAPNAME="");

    ///Make Theta,Phi,ZTheta,RockingAngle vs Time plots
    int Make_Plots(double PreTime, double PostTime, double Burst_t0, string Filename, string FT2_FILE, int verbosity=1, float TimeStep=1);

    ///Get theta phi and zenith theta coordinates
    void GetThetaPhi(double &theta, double &phi, double &ZTheta, double MET, string FT2_FILE, float RA=-999, float DEC=-999);
    ///Generates plots of theta and phi vs time for the source
    TCanvas* MakePointingPlots(double PreTime, double PostTime, double Burst_t0, double RA, double DEC, string FT2_FILE);
    ///Return EarthAzimuth angle of coordinate RA/DEC
    float GimmeEarthAzimuth(float C1, float C2, float C1Zenith, float C2Zenith, astro::SkyDir::CoordSystem coordinate=astro::SkyDir::EQUATORIAL);

    //Coordinate Transformations
    ///Convert from Equatorial to Galactic coordinates
    void Galactic(float ra, float dec, float *glon, float *glat);
    ///Convert from Equatorial to Galactic coordinates
    void Galactic(double ra, double dec, double *glon, double *glat);
    ///Convert from Galactic to Equatorial coordinates
    void unGalactic(float lon, float lat, float* ra, float* dec);
    ///Convert from Galactic to Equatorial coordinates
    void unGalactic(double lon, double lat, double* ra, double* dec);

    //Data-Class Functions
    int GetCTBClassLevel(string DataClass);
    int GetClassMask(string DataClass);
    string GetDataClassName_noConv(string DataClass);
    int  GetConversionType(string DataClass);
    string GetConversionName(string DataClass);
    string GetDataClassVersion(string DataClass);

    ///Return a Poisson probability for the case of an error in the bkg estimate (a la GRB080825c Bayesian method)
    double WeightedP(int nev, double bkg, double dbkg);
    ///Return an asymmetric "1sigma" Poisson error. Parameter "type" sets the kind of error (upper/lower)
    double PoissonErrorBar(short int type, int events);

    ///Get the value of a numeric parameter
    double Get(string name);
    ///Get the value of a string parameter
    string GetS(string name);
    ///Set the value of a numeric parameter
    void Set(string name,double val);
    ///Set the value of a string parameter
    void Set(string name,string vals);
    ///Load a configuration file
    void LoadConfig(string ConfigFile,bool ls=true);
    ///Print the current configuration
    void PrintConfig();
    void WriteConfig(string ConfigFile);

    void Run_gtltcube(string GRB_DIR,  double TMin, double TMax, string FT2_FILE, float FT1ZenithTheta_Cut, int verbosity, string gtltcube_Filename="");
    void Run_gtexpcube(string GRB_DIR, double TMin, double TMax, string FT2_FILE, string DATACLASS, float FT1ZenithTheta_Cut, char * Outfile, float Energy_Min,float Energy_Max, int Energy_Bins, int verbosity, string gtltcube_Filename="");

    void ReadGTI(vector <double>& GTI_Starts, vector <double>& GTI_Ends, string FitsAllSkyFile, double StartTime, double EndTime);

};

#endif

