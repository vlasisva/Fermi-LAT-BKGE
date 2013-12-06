// Author: Vlasios Vasileiou <vlasisva@gmail.com>
// $Header: /nfs/slac/g/glast/ground/cvs/GRBAnalysis-scons/BackgroundEstimator/src/BKGE_Tools/CoordinateTransformations.cxx,v 1.1.1.1 2011/06/02 19:41:04 jrb Exp $
#include "BackgroundEstimator/BKGE_Tools.h"



//ok the next two functions could have been merged using a template but
//I couldn't figure out how to define the template in the root dictionary
void TOOLS::Galactic(float ra, float dec, float* glon, float* glat){
  double aglon,aglat;
  TOOLS::Galactic((double)ra,double(dec),&aglon,&aglat);
  *glon=aglon;
  *glat=aglat;
}

//This is taken (borrowed) fron Milagro analysis code
void TOOLS::Galactic(double ra, double dec, double* glon, double* glat){

    int i;
    double rar,decr,x,y,z,rxy2,rxy;
    double vect[3],vectgal[3];
    static double jgal[3][3] = {{-0.054875539726,-0.873437108010,-0.483834985808},
            {0.494109453312,-0.444829589425, 0.746982251810},
            {-0.867666135858,-0.198076386122, 0.455983795705}};

    if (ra>180.) ra-=360.;
    decr= dec*DEG_TO_RAD;
    rar = ra*DEG_TO_RAD;

    // convert spherical to Cartesian coord
    vect[0] = cos (rar) * cos (decr);
    vect[1] = sin (rar) * cos (decr);
    vect[2] = sin (decr);

    for (i = 0; i < 3; i++) {
       vectgal[i] = vect[0]*jgal[i][0] + vect[1]*jgal[i][1] + vect[2]*jgal[i][2];
    }

    x = vectgal[0];
    y = vectgal[1];
    z = vectgal[2];
    *glon = atan2 (y, x);
    if (*glon < 0.) *glon = *glon + 6.283185307179586;

    rxy2 = x*x + y*y;
    rxy = sqrt (rxy2);
    *glat = atan2 (z, rxy);

    *glon = *glon/DEG_TO_RAD;
    *glat = *glat/DEG_TO_RAD;
//     if (*glon>180.) *glon -=360.;
//     if (*glon<=-180.) *glon +=360.;

}


//ok the next two functions could have been merged using a template but
//I couldn't figure out how to define the template in the root dictionary
void TOOLS::unGalactic(float lon, float lat, float* ra, float* dec){
  double ara,adec;
  TOOLS::unGalactic((double)lon,double(lat),&ara,&adec);
  *ra=ara;
  *dec=adec;
}


void TOOLS::unGalactic(double lon, double lat, double* ra, double* dec){

    double sind,y,x;
    const double gPoleRA = 3.3660333;
    const double gPoleDEC = 0.4734773;
    const double gpAscNode = 0.575959;  //Ascending node of GP on equator (rad)

    lon *= DEG_TO_RAD;
    lat *= DEG_TO_RAD;


    sind = cos(lat)*cos(gPoleDEC)*sin(lon-gpAscNode) + sin(lat)*sin(gPoleDEC);
    y = cos(lat)*cos(lon - gpAscNode);
    x = sin(lat)*cos(gPoleDEC) - cos(lat)*sin(gPoleDEC)*sin(lon-gpAscNode);

    *ra = atan(y/x);
    *dec = asin(sind);

    //convert back to degrees
    *dec/= DEG_TO_RAD;
    *ra /= DEG_TO_RAD;

    //Require longitude >=0 and <360
    while(*ra < 0)   *ra+=360;
    while(*ra>=360)  *ra-=360;
    //correct inverse tan calc for quadrant ambiguity
    if(x<0. && y<0.){
        //3rd quadrant, require longitude between 180 and 270 deg
        while(*ra < 180)         *ra +=180;
        while(*ra > 270)         *ra -=180;
        if(*ra >270 || *ra <180) printf("LON out of bounds.  Q3\n");
        }
    else if(x<0. && y>=0.){
        //2nd quadrant, require longitude between 90 and 180 deg
        while(*ra > 180)         *ra -=180;
        while(*ra < 90)          *ra +=180;
        }
    else if(x>=0. && y<0.){
        //4th quadrant, require longitude between 270 and 360 deg
        while(*ra < 270)         *ra +=180;
        if(*ra >360 || *ra <270) printf("LON out of bounds.  Q4\n");
        }
    else{
        //1st quadrant, require longitude between 0 and 90 deg
        while(*ra > 90)          *ra -=180;
        if(*ra >90 || *ra <0) printf("LON out of bounds.  Q1\n");
        }

    *ra += gPoleRA/DEG_TO_RAD;
    if(*ra>=360)  *ra-=360.;

}
