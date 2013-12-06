//Author: Vlasios Vasileiou <vlasisva@gmail.com>
//$Header: /nfs/slac/g/glast/ground/cvs/GRBAnalysis-scons/BackgroundEstimator/BackgroundEstimator/BackgroundEstimator_ext.h,v 1.1 2013/11/28 13:20:18 vlasisva Exp $
#ifndef _BackgroundEstimator_ext_H
#define _BackgroundEstimator_ext_H

#include "BackgroundEstimator.h"

namespace BKGE_NS_EXT{
    void Calc_TimeCorrectionFactors(vector<string> GRB_folders, vector  <double> METs, string Dataclass, double MinE, double MaxE, int NBins);
};


class BackgroundEstimator_ext : public BackgroundEstimator{
  public:
    BackgroundEstimator_ext(string aClass, double EMin=-1, double EMax=-1, int EBins=-1, bool initialize=true, bool ShowLogo=true):
    BackgroundEstimator(aClass,EMin,EMax,EBins,initialize,ShowLogo){};
    ~BackgroundEstimator_ext();
    
    ///Create the data files used for the background estimation. Normal users don't need to run that
    void CreateDataFiles(string FitsAllSkyFilesList, string FT2_FILE, double StartTime=0, double EndTime=0, float ZenithThetaCut=100); 

  private:
   ///Data Files
    void CalcResiduals(string FitsAllSkyFile);
    void Make_McIlwainL_Fits(string FitsAllSkyFile);
    void Make_ThetaPhi_Fits(string FitsAllSkyFile);
};

#endif

