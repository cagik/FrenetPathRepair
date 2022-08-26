#ifndef DATA_STURCT_H
#define DATA_STURCT_H

#include <vector>

typedef struct
{
    unsigned char	road_type:4;
    unsigned char	right_wall:1;
    unsigned char	left_wall:1;
    unsigned char	road_direction:2;
    unsigned char	reserve:8;
} Property;
typedef union
{
    Property data;
    char byte[sizeof(Property)];
}PropertyUnion;


struct RRTResult{
    std::vector<double> x;
    std::vector<double> y;
    std::vector<double> phi;
};

typedef struct _GPSInfo
{
  unsigned short int road_id;
  float iHead;//方位角
  double iLatitude;//纬度
  double iLongitude;//经度
  float altitude;//高度
  float velocity;
  float left_width;
  float right_width;
  PropertyUnion property;
  float iDistance;//距离
  float iCurvature;//曲率
  float dist_origin;//距离原点的距离
  double theta;     //与东偏向
  double xEN;
  double yEN;
}GPSInfo;

struct VehicleParameter{
    float vehicle_length;
    float vehicle_width;
};

#endif // DATA_STURCT_H
