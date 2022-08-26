#include <iostream>
#include <jsoncpp/json/json.h>
#include <fstream>
#include <string>
#include <sstream>

#include <iomanip>

#include <vector>

#include "pathrepairplanner.h"

using namespace std;

int main(){

    vector<pair<double,double>> leftBoundStream;
    vector<pair<double,double>> rightBoundStream;
    vector<pair<double,double>> refPathStream;
    vector<pair<double,double>> refinePathStream;
    vector<pair<double,double>> errorPoint_Mid_Stream;
    vector<pair<double,double>> errorPoint_StartAndEnd_Stream;

    Json::Value root;
    Json::Reader reader;
    ifstream jsonFile;
    jsonFile.open("/home/cavata/workspace/pathRefine/bin/new1.json");

    if(!jsonFile.good()){
        cout << "error" << endl;
    }

    if(reader.parse(jsonFile, root)){
        cout << "succ" << endl;
    }
    else{
        cout << "failed" << endl;
    }

    Json::Value msg = root["msg"];

    Json::Value dis = msg["distanceSetting"];

    Json::Value leftBoundJson = msg["leftBoundInfo"];
    Json::Value rightBoundJson = msg["rightBoundInfo"];
    Json::Value mainPathJson = msg["mainPathInfo"];

    std::vector<GPSInfo> mainPath;
    std::vector<GPSInfo> leftBound;
    std::vector<GPSInfo> rightBound;


    if(rightBoundJson.isArray())
    {   
        rightBoundStream.resize(rightBoundJson.size());
        for(int i = 0; i < rightBoundJson.size(); i++){
            GPSInfo tmp;
            tmp.iLatitude = rightBoundJson[i]["latitude"].asDouble();
            tmp.iLongitude = rightBoundJson[i]["longitude"].asDouble();
            leftBound.push_back(tmp);
            rightBoundStream[i].first = rightBoundJson[i]["latitude"].asDouble();
            rightBoundStream[i].second = rightBoundJson[i]["longitude"].asDouble();
        }
    }

    if(leftBoundJson.isArray())
    {
        leftBoundStream.resize(leftBoundJson.size());
        for(int i = 0; i < leftBoundJson.size(); i++){
            GPSInfo tmp;
            tmp.iLatitude = leftBoundJson[i]["latitude"].asDouble();
            tmp.iLongitude = leftBoundJson[i]["longitude"].asDouble();
            rightBound.push_back(tmp);
            leftBoundStream[i].first = leftBoundJson[i]["latitude"].asDouble();
            leftBoundStream[i].second = leftBoundJson[i]["longitude"].asDouble();
        }
    }

    if(mainPathJson.isArray())
    {
        refPathStream.resize(mainPathJson.size());
        for(int i = 0; i < mainPathJson.size(); i++){
            GPSInfo tmp;
            tmp.iLatitude = mainPathJson[i]["latitude"].asDouble();
            tmp.iLongitude = mainPathJson[i]["longitude"].asDouble();
            tmp.iHead = mainPathJson[i]["head"].asDouble();
            mainPath.push_back(tmp);
            refPathStream[i].first = mainPathJson[i]["latitude"].asDouble();
            refPathStream[i].second = mainPathJson[i]["longitude"].asDouble();
        }
    }


    Json::Value veh_json = msg["vehicleParam"];

    VehicleParameter vehicleParam;

    vehicleParam.vehicle_length = veh_json["length"].asDouble();
    vehicleParam.vehicle_width = veh_json["width"].asDouble();
    double distanceSet = msg["distanceSetting"].asDouble();


    vector<GPSInfo> refine_result;
    vector<GPSInfo> left_edge_trajectory;
    vector<GPSInfo> right_edge_trajectory;

    int point_index_select = 0;
    int point_end_index_select  = mainPath.size();
    std::vector<std::pair<double,double>> edge_bound;


    int result_flag_ = 0;

    vector<GPSInfo> error_point_se;
    vector<GPSInfo> error_point_mid;

    PathRepairPlanner *planner;
    if(!planner->plan(mainPath, leftBound, rightBound, point_index_select, point_end_index_select, vehicleParam,
                      &refine_result, &left_edge_trajectory, &right_edge_trajectory, edge_bound, distanceSet, result_flag_,
                      error_point_se, error_point_mid))
    {
        std::cout << "plan failed" << std::endl;
    }
    else{
        std::cout << "plane success" << std::endl;
    }
    
    cout << "flag" << result_flag_ << endl;

    refinePathStream.resize(refine_result.size());
    for(size_t i = 0; i < refine_result.size(); i++){
        refinePathStream[i].first = refine_result[i].iLatitude;
        refinePathStream[i].second = refine_result[i].iLongitude;
    }

    errorPoint_StartAndEnd_Stream.resize(error_point_se.size());
    for(size_t i = 0; i < error_point_se.size(); i++){
        errorPoint_StartAndEnd_Stream[i].first = error_point_se[i].iLatitude;
        errorPoint_StartAndEnd_Stream[i].second = error_point_se[i].iLongitude;
    }


    errorPoint_Mid_Stream.resize(error_point_mid.size());
    for(size_t i = 0; i < error_point_mid.size(); i++){
        errorPoint_Mid_Stream[i].first = error_point_mid[i].iLatitude;
        errorPoint_Mid_Stream[i].second = error_point_mid[i].iLongitude;
    }

    DataToTxt("./data/leftBound.txt", leftBoundStream);
    DataToTxt("./data/rightBound.txt", rightBoundStream);
    DataToTxt("./data/refPath.txt", refPathStream);
    DataToTxt("./data/refinePath.txt", refinePathStream);
    DataToTxt("./data/errorPoint_Mid.txt", errorPoint_Mid_Stream);
    DataToTxt("./data/errorPoint_StartAndEnd.txt", errorPoint_StartAndEnd_Stream);

    jsonFile.close();


    return 0;
}