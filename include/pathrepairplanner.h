#ifndef PATHREPAIRPLANNER_H
#define PATHREPAIRPLANNER_H
#include <array>
#include <iostream>
#include <eigen3/Eigen/Core>
#include <eigen3/Eigen/Dense>
#include <set>
#include <vector>
#include <osqp/osqp.h>


#include <iomanip>

#include "mathfun.h"
#include "data_sturct.h"
#include "KDTree.hpp"
#include "TxtStream.hpp"

using namespace std;


typedef struct _Coffient
{
  double s0;
  double x_a0;
  double x_a1;
  double x_a2;
  double x_a3;
  double y_a0;
  double y_a1;
  double y_a2;
  double y_a3;
} Coffient;

template <typename T>
T* CopyData(const std::vector<T>& vec) {
  T* data = new T[vec.size()];
  memcpy(data, vec.data(), sizeof(T) * vec.size());
  return data;
}

class PathRepairPlanner
{
public:
    PathRepairPlanner();
    PathRepairPlanner(const size_t num_of_knots, const double delta_s,
                  const std::array<double, 3>& x_init);
    bool plan(const vector<GPSInfo> &main_path, const vector<GPSInfo> &left_edge_points, const vector<GPSInfo> &right_edge_points, const int start_key, const int end_key,
              VehicleParameter m_vehicle, vector<GPSInfo> *refined_trajectory,vector<GPSInfo>*left_edge_trajectory,vector<GPSInfo>*right_edge_trajectory,
              std::vector<std::pair<double,double>> &edge_bound, double boundSetting, int &result_flag, 
              std::vector<GPSInfo> &error_point_se, std::vector<GPSInfo> &error_point_mid);
    bool generatePath(int matched_point_index,const std::vector<GPSInfo>& raw_path,
                      const std::vector<double> frenet_path, std::vector<Coffient> &coefficients_,
                      std::vector<GPSInfo> &optimal_path, const std::vector<float> raw_points_curvature);


      bool resultCheck(vector<GPSInfo> *refined_trajectory, int &res_flag);

      void CalculateKernel(std::vector<c_float>* P_data,
                           std::vector<c_int>* P_indices,
                           std::vector<c_int>* P_indptr);
      void CalculateAffineConstraint(std::vector<c_float>* A_data,
                                     std::vector<c_int>* A_indices,
                                     std::vector<c_int>* A_indptr,
                                     std::vector<c_float>* lower_bounds,
                                     std::vector<c_float>* upper_bounds,
                                     std::vector<float> curvature,
                                     bool &checkflag);
      void CalculateOffset(std::vector<c_float>* q);
      OSQPSettings* SolverDefaultSettings();
      OSQPData* FormulateProblem(vector<float> curvature, bool &boundErrorCheck);
      bool Optimize(const int max_iter, vector<float> curvature, int &res_flag);
      void FreeData(OSQPData* data);

      inline void set_weight_x(double w) { weight_x_ = w; }
      inline void set_weight_dx(double w) { weight_dx_ = w; }
      inline void set_weight_ddx(double w) { weight_ddx_ = w; }
      inline void set_weight_dddx(double w) { weight_dddx_ = w; }

      void set_end_state_ref(const std::array<double, 3>& weight_end_state,
                             const std::array<double, 3>& end_state_ref);

      void set_x_ref(const double weight_x_ref, std::vector<double> x_ref);

      inline void set_scale_factor(const std::array<double, 3>& factor) { scale_factor_ = factor; }

      inline void set_x_bounds(const std::vector<std::pair<double, double>>& x_bounds) { x_bounds_ = x_bounds; }
      inline void set_dx_bounds(const std::vector<std::pair<double, double>>& dx_bounds) { dx_bounds_ = dx_bounds; }

      void set_dx_bounds(double low_bound, double up_bound);
      void set_ddx_bounds(double low_bound, double up_bound);
      void set_dddx_bounds(double low, double up);
//      inline std::vector<double> x(){return x_;}
      void lateral_distance_min();
      double GetPointHeading(double p_x, double p_y, double ori_x, double ori_y);

      //叉乘分离左右边界点
      void SpiltLRPoints(const vector<GPSInfo> &main_path, const vector<GPSInfo> edge_points, const int start_key,
                         vector<GPSInfo> *left_edge_points, vector<GPSInfo> *right_edge_points);
      //验证修正路径后曲率
      void VerifyResultCurvature(GPSInfo start_gps, vector<GPSInfo> refined_trajectory);
      //验证修正路径后曲率-公司方式
      void VerifyResultCurvatureCompany(GPSInfo start_gps, vector<GPSInfo> refined_trajectory);
      //计算待修正路径曲率-后续曲率约束
      std::vector<float> CalRawPathCurvature(GPSInfo start_gps, vector<GPSInfo> raw_path_points);
      //修正后x，y坐标输出
      void CoutResultXY(GPSInfo start_gps, vector<GPSInfo> result);

private:
      size_t num_of_knots_;

      // output
      std::vector<double> x_;
      std::vector<double> dx_;
      std::vector<double> ddx_;


      std::array<double, 3> x_init_;
      std::array<double, 3> scale_factor_ = {{1.0, 1.0, 1.0}};

      std::vector<std::pair<double, double>> x_bounds_;
      std::vector<std::pair<double, double>> dx_bounds_;
      std::vector<std::pair<double, double>> ddx_bounds_;
      std::pair<double, double> dddx_bound_;

      double weight_x_ = 0.0;
      double weight_dx_ = 0.0;
      double weight_ddx_ = 0.0;
      double weight_dddx_ = 0.0;

      double delta_s_ = 0.5;

      bool has_x_ref_ = false;
      double weight_x_ref_ = 0.0;
      std::vector<double> x_ref_;

      bool has_end_state_ref_ = false;
      std::array<double, 3> weight_end_state_ = {{0.0, 0.0, 0.0}};
      std::array<double, 3> end_state_ref_;

};

#endif // PATHREPAIRPLANNER_H
