#include "mathfun.h"

namespace Tage {
    float getGPSDistance(double lat1, double lng1, double lat2, double lng2)
    {
        float disn;
        double _x,_y;
        _x=fabs(lng1-lng2)*111700*cos(lat1/180*M_PI);
        _y=fabs(lat1-lat2)*111700;

        disn=sqrt(_x*_x+_y*_y);
        return disn;
    }

    float getAngleByTwoGps(double lat1, double lng1, double lat2, double lng2)
    {
        float result_angle;
        float dx = (lng2 - lng1) * 111700 *cos(lat1/180 * M_PI);
        float dy = (lat2 - lat1) * 111700;

        result_angle = atan2(dx,dy)/M_PI*180;//计算与y轴正方向的夹角(-180~180)
        result_angle = (result_angle <0) ? 360+result_angle:result_angle;

        return result_angle;
    }

    void transGpsIntoInertia(const GPSInfo &my_gps, const GPSInfo &target_point,
                                             float *inertia_x, float *inertia_y)
    {

        float dis = getGPSDistance(my_gps.iLatitude, my_gps.iLongitude, target_point.iLatitude, target_point.iLongitude);//两点之间绝对位置
        float azimuth = getAngleByTwoGps(my_gps.iLatitude,my_gps.iLongitude,target_point.iLatitude,target_point.iLongitude);//由车辆的gps点到目标点的航向,正北为0度,0~360顺时针
        float d_angle = azimuth - my_gps.iHead;
        if(d_angle >=180)
            d_angle -=360;
        if(d_angle <-180)
            d_angle +=360;
        *inertia_x = dis*sin(d_angle*M_PI/180);
        *inertia_y = dis*cos(d_angle*M_PI/180);
    }

    GPSInfo transXYToGps(double lat, double lng, float head, double x, double y)
    {
        GPSInfo gps_object;
        double angle,temp_head;
        double dis = sqrt(x*x + y*y);
        if(y == 0)
        {
            if(x == 0)
            {
               temp_head = head;
            }
            else if(x > 0)
            {
               temp_head = 90+head;
            }
            else if(x < 0)
            {
               temp_head = head-90;
            }
        }
        else if(y > 0)
        {
            if(x == 0)
                temp_head = head;
            else
            {
                angle = atan(x/y)*180/M_PI;
                temp_head = angle+head;
            }
        }
        else if(y <0)
        {
            if(x == 0)
                temp_head = head+180;
            else
            {
                angle = atan(x/y)*180/M_PI;
                temp_head = 180+angle+head;
            }
        }

        double dx = dis*cos((90-temp_head)*M_PI/180);
        double dy = dis*sin((90-temp_head)*M_PI/180);
        gps_object.iLatitude = lat+dy/111700;
        gps_object.iLongitude = lng+dx/111700/cos(lat*M_PI/180);
        return gps_object;
    }

    }

