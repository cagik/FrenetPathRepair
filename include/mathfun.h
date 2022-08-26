#ifndef MATHFUN_H
#define MATHFUN_H

#include <cmath>
#include <vector>

#include "data_sturct.h"


using namespace std;

namespace Tage {

    float getGPSDistance(double lat1,double lng1,double lat2,double lng2);
    float getAngleByTwoGps(double lat1, double lng1, double lat2, double lng2);
    void transGpsIntoInertia(const GPSInfo &my_gps, const GPSInfo &target_point,float *inertia_x, float *inertia_y);
    GPSInfo transXYToGps(double lat,double lng,float head,double x,double y);
}
#endif // MATHFUN_H
