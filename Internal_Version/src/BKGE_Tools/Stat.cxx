//Author: Vlasios Vasileiou <vlasisva@gmail.com>
// $Header: /nfs/slac/g/glast/ground/cvs/GRBAnalysis-scons/BackgroundEstimator/src/BKGE_Tools/Stat.cxx,v 1.3 2013/11/13 07:56:16 vlasisva Exp $
#include "BackgroundEstimator/BKGE_Tools.h"
#include "Math/ProbFuncMathCore.h"
#include "Math/PdfFuncMathCore.h"


double TOOLS::PoissonErrorBar(short int type, int events) {
 if (type!=0 && type!=1) {printf("%s: unknown error type %d\n",__FUNCTION__,type); throw std::runtime_error("");}
 //type==0 is bottom error bar, type==1 is top
 static bool first=true;
 static double ParError[2][3]; //0 is low, 1 is high
 if (first) {
    ParError[0][0]=-8.646473e-01;
    ParError[0][1]=7.718411e-01;
    ParError[0][2]=2.894903e-03;

    ParError[1][0]=1.421087e+00;
    ParError[1][1]=1.224625e+00;
    ParError[1][2]=-2.816522e-03;
    first=false;
 }
 if (events>30)      return sqrt((float)events);
 else if (events==0) return 0;
 double error =0;
 for (short int ipar=0;ipar<3;ipar++) error+=ParError[type][ipar]*pow((float)events,ipar);
 error = fabs(error-events);
 if (type==0 && error>events) error=events;
 return error;

}


//This gives a Poisson probability with a nexp that has a finite uncertainty
//the returned probability is the integral from ndet to infinity
double TOOLS::WeightedP(int nev, double bkg, double dbkg) {
    if (dbkg==0) return ROOT::Math::poisson_cdf_c(nev,bkg)+ROOT::Math::poisson_pdf(nev,bkg);

    double PSum=0,WSum=0;
    double Start=bkg-3*dbkg;

    const double daback=dbkg/100.;
    if (Start<0) Start=daback;
    
    for (double aback=Start;aback<=bkg+3*dbkg;aback+=daback) {
       double weight = ROOT::Math::gaussian_pdf(aback,dbkg,bkg);
       double ap = ROOT::Math::poisson_cdf_c(nev,aback)+ROOT::Math::poisson_pdf(nev,aback);
       PSum+=weight*ap*daback;
       WSum+=weight*daback;
    }
    PSum/=WSum;
    if (WSum<0.95 || WSum>1) {printf("wsum=%e %e %e %e \n",WSum,PSum,bkg,dbkg);}
    return PSum;
}



