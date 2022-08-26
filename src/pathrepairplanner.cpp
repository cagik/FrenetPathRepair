#include "pathrepairplanner.h"
PathRepairPlanner::PathRepairPlanner()
{

}

namespace {
constexpr double kMaxVarRange = 1.0e10;
}

PathRepairPlanner::PathRepairPlanner(const size_t num_of_knots, const double delta_s,
              const std::array<double, 3>& x_init)
{
    if (num_of_knots <= 2) {
      std::cout << "Error when init problem." << std::endl;
      return;
    }
    num_of_knots_ = num_of_knots;
    x_init_ = x_init;
    delta_s_ = delta_s;
    x_bounds_.resize(num_of_knots_, std::make_pair(-kMaxVarRange, kMaxVarRange));
    dx_bounds_.resize(num_of_knots_, std::make_pair(-kMaxVarRange, kMaxVarRange));
    ddx_bounds_.resize(num_of_knots_,
                       std::make_pair(-kMaxVarRange, kMaxVarRange));
}

void generate_coff(GPSInfo start_pose,GPSInfo next_pose,double s0, double sf,std::vector<Coffient> &ref_coff)
{
    double x0 = start_pose.xEN;
    double y0 = start_pose.yEN;
    double theta0 = start_pose.theta;
    double xf = next_pose.xEN;
    double yf = next_pose.yEN;
    double thetaf = next_pose.theta;

    Eigen::Matrix4d A;
    A << 1.0, s0, s0*s0, s0*s0*s0,
            1.0, sf, sf*sf, sf*sf*sf,
            0.0, 1.0, 2.0 * s0, 3 * s0 * s0,
            0.0, 1.0, 2.0 * sf, 3 * sf * sf;
    Eigen::Vector4d B(x0, xf, std::cos(theta0), std::cos(thetaf));
    Eigen::Vector4d a_vec = A.colPivHouseholderQr().solve(B);

    Eigen::Vector4d C(y0, yf, std::sin(theta0), std::sin(thetaf));
    Eigen::Vector4d b_vec = A.colPivHouseholderQr().solve(C);

    Coffient p;
    p.s0 = s0;
    p.x_a0 = a_vec[0];
    p.x_a1 = a_vec[1];
    p.x_a2 = a_vec[2];
    p.x_a3 = a_vec[3];
    p.y_a0 = b_vec[0];
    p.y_a1 = b_vec[1];
    p.y_a2 = b_vec[2];
    p.y_a3 = b_vec[3];

    ref_coff.push_back(p);

}

GPSInfo frenet_to_cartesian(double s, double point_l, double delta_theta, std::vector<Coffient> referline_coff)
{

    GPSInfo result_point;

    for(std::size_t h = 0; h < referline_coff.size();h++)
    {

        if(s < referline_coff.at(h).s0)
        {
            Coffient param = referline_coff.at(h-1);
            double x = param.x_a0 + param.x_a1 * s + param.x_a2 * s * s + param.x_a3 * s * s * s;
            double d_x = param.x_a1 + 2 * param.x_a2 * s + 3 * param.x_a3 * s * s;
            double y = param.y_a0 + param.y_a1 * s + param.y_a2 * s * s + param.y_a3 * s * s * s;
            double d_y = param.y_a1 + 2 * param.y_a2 * s + 3 * param.y_a3 * s * s;
            double theta = std::atan2(d_y, d_x);
            result_point.xEN = x;
            result_point.yEN = y;
            result_point.theta = theta;

            result_point.xEN = result_point.xEN + point_l*std::cos(result_point.theta + M_PI*0.5);
            result_point.yEN = result_point.yEN + point_l*std::sin(result_point.theta + M_PI*0.5);
            break;
        }
    }

    result_point.theta = result_point.theta + delta_theta;

    double heading;
    double res_theta = result_point.theta;

    if(res_theta *180/M_PI < 90)
    {
        heading = 90 - res_theta *180/M_PI;
    }
    else
    {
        heading = 450 - res_theta*180/M_PI;
    }

    result_point.iHead = heading;

    return result_point;
}


GPSInfo frenet_to_cartesian(double point_l, double delta_theta, GPSInfo raw_point)
{

    GPSInfo result_point;

    result_point = raw_point;
    result_point.xEN = raw_point.xEN - point_l*std::sin(raw_point.theta);
    result_point.yEN = raw_point.yEN + point_l*std::cos(raw_point.theta);
    result_point.theta = raw_point.theta + delta_theta;

    double heading;
    double res_theta = result_point.theta;

    if(res_theta *180/M_PI < 90)
    {
        heading = 90 - res_theta * 180/M_PI;
    }
    else
    {
        heading = 450 - res_theta * 180/M_PI;
    }

    if(heading < 0.0)
        heading += 360.0;

    result_point.iHead = heading/ 180*M_PI;
    return result_point;
}


bool PathRepairPlanner::generatePath(int matched_point_index,const std::vector<GPSInfo>& raw_path,
                  const std::vector<double> frenet_path, std::vector<Coffient> &coefficients_, std::vector<GPSInfo> &optimal_path,
                  const std::vector<float> raw_points_curvature)
{

    GPSInfo cartesian_pose;
    double pre_x,pre_y,pre_head;

    for(std::size_t i = 0; i< frenet_path.size()-1 ; i++)
    {

        double heading = atan2(frenet_path[i+1] - frenet_path[i], 0.2); //这里 0.2代表相邻两个点的单位纵向长度
        cartesian_pose = frenet_to_cartesian(frenet_path[i], heading , raw_path[i]);//frenet_path heading

        if(i == 0)
        {
            pre_x = cartesian_pose.xEN;
            pre_y = cartesian_pose.yEN;
            pre_head = cartesian_pose.iHead;
        }else
        {
            double idistance = std::hypot(pre_x - cartesian_pose.xEN , cartesian_pose.yEN - pre_y);
            cartesian_pose.iDistance = idistance;
        }

        cartesian_pose.iCurvature = fabs(raw_points_curvature[i]/(1-raw_points_curvature[i]*frenet_path[i]));

        pre_x = cartesian_pose.xEN;
        pre_y = cartesian_pose.yEN;
        pre_head = cartesian_pose.iHead;
        optimal_path.push_back(cartesian_pose);
        if(i == frenet_path.size()-2)
        {
            cartesian_pose = frenet_to_cartesian(frenet_path[frenet_path.size()-1], heading , raw_path[frenet_path.size()-1]);//frenet_path heading
            cartesian_pose.iCurvature = raw_points_curvature[i]/(1-raw_points_curvature[i]*frenet_path[i]);
            optimal_path.push_back(cartesian_pose);
        }
    }
    return true;
}


void PathRepairPlanner::CalculateKernel(std::vector<c_float>* P_data,
                                  std::vector<c_int>* P_indices,
                                  std::vector<c_int>* P_indptr) {
      const int n = static_cast<int>(num_of_knots_);
      const int num_of_variables = 3 * n;
      const int num_of_nonzeros = num_of_variables + (n - 1);
      std::vector<std::vector<std::pair<c_int, c_float>>> columns(num_of_variables);
      int value_index = 0;

      //column列
      // x(i)^2 * (w_x + w_x_ref)
      for (int i = 0; i < n - 1; ++i) {
        columns[i].emplace_back(
            i, (weight_x_ + weight_x_ref_) / (scale_factor_[0] * scale_factor_[0]));
        ++value_index;
      }
      // x(n-1)^2 * (w_x + w_x_ref + w_end_x)
      columns[n - 1].emplace_back(
          n - 1, (weight_x_ + weight_x_ref_ + weight_end_state_[0]) /
                     (scale_factor_[0] * scale_factor_[0]));
      ++value_index;

      // x(i)'^2 * w_dx
      for (int i = 0; i < n - 1; ++i) {
        columns[n + i].emplace_back(
            n + i, weight_dx_ / (scale_factor_[1] * scale_factor_[1]));
        ++value_index;
      }
      // x(n-1)'^2 * (w_dx + w_end_dx)
      columns[2 * n - 1].emplace_back(2 * n - 1,
                                      (weight_dx_ + weight_end_state_[1]) /
                                          (scale_factor_[1] * scale_factor_[1]));
      ++value_index;

      auto delta_s_square = delta_s_ * delta_s_;
      // x(i)''^2 * (w_ddx + 2 * w_dddx / delta_s^2)
      columns[2 * n].emplace_back(2 * n,
                                  (weight_ddx_ + weight_dddx_ / delta_s_square) /
                                      (scale_factor_[2] * scale_factor_[2]));
      ++value_index;
      for (int i = 1; i < n - 1; ++i) {
        columns[2 * n + i].emplace_back(
            2 * n + i, (weight_ddx_ + 2.0 * weight_dddx_ / delta_s_square) /
                           (scale_factor_[2] * scale_factor_[2]));
        ++value_index;
      }


      columns[3 * n - 1].emplace_back(
          3 * n - 1,
          (weight_ddx_ + 1.0*weight_dddx_ / delta_s_square + weight_end_state_[2]) /
              (scale_factor_[2] * scale_factor_[2]));
      ++value_index;

      // -2 * w_dddx / delta_s^2 * x(i)'' * x(i + 1)''

      for (int i = 0; i < n-1; ++i) {
        columns[2 * n + i].emplace_back(2 * n + i + 1,
                                        (-2.0 * weight_dddx_ / delta_s_square) /
                                            (scale_factor_[2] * scale_factor_[2]));
        ++value_index;
      }

      if (value_index != num_of_nonzeros) {
        std::cout << "Error in calculate kernel!" << std::endl;
      }

      int ind_p = 0;
      for (int i = 0; i < num_of_variables; ++i) {
        P_indptr->push_back(ind_p);
        for (const auto& row_data_pair : columns[i]) {
          P_data->push_back(row_data_pair.second * 2.0);
          P_indices->push_back(row_data_pair.first);
          ++ind_p;
        }
      }
      P_indptr->push_back(ind_p);
}

void PathRepairPlanner::CalculateAffineConstraint(
    std::vector<c_float>* A_data, std::vector<c_int>* A_indices,
    std::vector<c_int>* A_indptr, std::vector<c_float>* lower_bounds,
    std::vector<c_float>* upper_bounds,std::vector<float> curvature, bool &checkflag) {

    const int n = static_cast<int>(num_of_knots_);
    const int num_of_variables = 3 * n;
    // ************************   compared with the initial, this add n kappa constrains *********************
    const int num_of_constraints = num_of_variables + 3 * (n - 1) + 3 + n;
    lower_bounds->resize(num_of_constraints);
    upper_bounds->resize(num_of_constraints);

    std::vector<std::vector<std::pair<c_int, c_float>>> variables(
      num_of_variables);

    int constraint_index = 0;
    //变量约束，通过赋值变量上下边界，完成约束设置
    // set x, x', x'' bounds
    for (int i = 0; i < num_of_variables; ++i) {
        if (i < n) {
            variables[i].emplace_back(constraint_index, 1.0);
            lower_bounds->at(constraint_index) = x_bounds_[i].first * scale_factor_[0];
            upper_bounds->at(constraint_index) = x_bounds_[i].second * scale_factor_[0];

            if(lower_bounds->at(constraint_index) > upper_bounds->at(constraint_index)){
                cout << "x error" << endl;
            }

        } 
        else if (i < 2 * n) {
            variables[i].emplace_back(constraint_index, 1.0);
            lower_bounds->at(constraint_index) = dx_bounds_[i - n].first * scale_factor_[1];
            upper_bounds->at(constraint_index) = dx_bounds_[i - n].second * scale_factor_[1];

            if(lower_bounds->at(constraint_index) > upper_bounds->at(constraint_index)){
                cout << "dx error" << endl;
            }
        } 
        else {
            variables[i].emplace_back(constraint_index, 1.0);
            lower_bounds->at(constraint_index) = ddx_bounds_[i - 2 * n].first * scale_factor_[2];
            upper_bounds->at(constraint_index) =    ddx_bounds_[i - 2 * n].second * scale_factor_[2];

            if(lower_bounds->at(constraint_index) > upper_bounds->at(constraint_index)){
                cout << "ddx error" << endl;
            }
        }
    ++constraint_index;
    }

    if (constraint_index != num_of_variables) {
    std::cout << "Error in Calculate Affine Constraint 1" << std::endl;
    }

    //运动学约束
    // x(i->i+1)''' = (x(i+1)'' - x(i)'') / delta_s
    for (int i = 0; i + 1 < n; ++i) {
        variables[2 * n + i].emplace_back(constraint_index, -1.0);
        variables[2 * n + i + 1].emplace_back(constraint_index, 1.0);
        lower_bounds->at(constraint_index) = dddx_bound_.first * delta_s_ * scale_factor_[2];
        upper_bounds->at(constraint_index) = dddx_bound_.second * delta_s_ * scale_factor_[2];

        if(lower_bounds->at(constraint_index) > upper_bounds->at(constraint_index)){
                cout << "dddx error" << endl;
        }

        ++constraint_index;
        if(lower_bounds->at(constraint_index) > 0){
            std::cout<<"error"<<std::endl;
        }
    }

    //连续性约束
    // x(i+1)' - x(i)' - 0.5 * delta_s * x(i)'' - 0.5 * delta_s * x(i+1)'' = 0
    for (int i = 0; i + 1 < n; ++i) {
        variables[n + i].emplace_back(constraint_index, -1.0 * scale_factor_[2]);
        variables[n + i + 1].emplace_back(constraint_index, 1.0 * scale_factor_[2]);
        variables[2 * n + i].emplace_back(constraint_index, -0.5 * delta_s_ * scale_factor_[1]);
        variables[2 * n + i + 1].emplace_back(constraint_index, -0.5 * delta_s_ * scale_factor_[1]);
        lower_bounds->at(constraint_index) = 0.0;
        upper_bounds->at(constraint_index) = 0.0;
        if(lower_bounds->at(constraint_index) > upper_bounds->at(constraint_index)){
                cout << "contiu error" << endl;
        }
        ++constraint_index;
    }

    // x(i+1) - x(i) - delta_s * x(i)'
    // - 1/3 * delta_s^2 * x(i)'' - 1/6 * delta_s^2 * x(i+1)''
    auto delta_s_sq_ = delta_s_ * delta_s_;
    for (int i = 0; i + 1 < n; ++i) {
        variables[i].emplace_back(constraint_index, -1.0 * scale_factor_[1] * scale_factor_[2]);
        variables[i + 1].emplace_back(constraint_index, 1.0 * scale_factor_[1] * scale_factor_[2]);
        variables[n + i].emplace_back(constraint_index, -delta_s_ * scale_factor_[0] * scale_factor_[2]);
        variables[2 * n + i].emplace_back(constraint_index, -delta_s_sq_ / 3.0 * scale_factor_[0] * scale_factor_[1]);
        variables[2 * n + i + 1].emplace_back(constraint_index, -delta_s_sq_ / 6.0 * scale_factor_[0] * scale_factor_[1]);

        lower_bounds->at(constraint_index) = 0.0;
        upper_bounds->at(constraint_index) = 0.0;
        if(lower_bounds->at(constraint_index) > upper_bounds->at(constraint_index)){
                cout << "ddt error" << endl;
        }
        ++constraint_index;
    }

    //曲率约束
    // ***********************  kappa constain *************************************************
    float k_max_limit = 0.075; // vehicle maximum kappa
//    float k_ref = 0.14; //this can be replaced by a ref_kappa_set according to the longitudinal offset

    for (int i= 0; i< n; ++i)
    {
//      variables[i].emplace_back(constraint_index,k_max_limit*k_ref);
        variables[i].emplace_back(constraint_index,k_max_limit*curvature[i]);
        lower_bounds->at(constraint_index) = -10000000*scale_factor_[0];
//      upper_bounds->at(constraint_index) = (k_max_limit-fabs(k_ref)) * scale_factor_[0];
        upper_bounds->at(constraint_index) = (k_max_limit-fabs(curvature[i])) * scale_factor_[0];
        ++constraint_index;
    }

    // constrain on x_init 初始状态约束，即最终轨迹的初始状态要与当前车辆状态保持一致
    variables[0].emplace_back(constraint_index, 1.0);
    lower_bounds->at(constraint_index) = x_init_[0] * scale_factor_[0];
    upper_bounds->at(constraint_index) = x_init_[0] * scale_factor_[0];
    ++constraint_index;

    variables[n].emplace_back(constraint_index, 1.0);
    lower_bounds->at(constraint_index) = x_init_[1] * scale_factor_[1];
    upper_bounds->at(constraint_index) = x_init_[1] * scale_factor_[1];
    ++constraint_index;

    variables[2 * n].emplace_back(constraint_index, 1.0);
    lower_bounds->at(constraint_index) = x_init_[2] * scale_factor_[2];
    upper_bounds->at(constraint_index) = x_init_[2] * scale_factor_[2];
    ++constraint_index;

    if (constraint_index != num_of_constraints) {
    std::cout << "Error in Calculate Affine Constraint 2" << std::endl;
    }

    for(int i = 0; i < lower_bounds->size(); i++){
        if(lower_bounds->at(i) > upper_bounds->at(i)){
            checkflag = false;
            return;
        }
    }

    checkflag = true;

    int ind_p = 0;
    for (int i = 0; i < num_of_variables; ++i) {
    A_indptr->push_back(ind_p);
    for (const auto& variable_nz : variables[i]) {
      // coefficient
      A_data->push_back(variable_nz.second);

      // constraint index
      A_indices->push_back(variable_nz.first);
      ++ind_p;
    }
    }

    A_indptr->push_back(ind_p);
}

void PathRepairPlanner::CalculateOffset(std::vector<c_float>* q) {
    if (q == nullptr) {
    std::cout << "q should not be nullptr!!" << std::endl;
    return;
    }
    const int n = static_cast<int>(num_of_knots_);
    const int kNumParam = 3 * n;
    q->resize(kNumParam, 0.0);

    if (has_x_ref_) {
    for (int i = 0; i < n; ++i) {
      q->at(i) += -2.0 * weight_x_ref_ * x_ref_[i] / scale_factor_[0];
    }
    }

    if (has_end_state_ref_) {
    q->at(n - 1) +=
        -2.0 * weight_end_state_[0] * end_state_ref_[0] / scale_factor_[0];
    q->at(2 * n - 1) +=
        -2.0 * weight_end_state_[1] * end_state_ref_[1] / scale_factor_[1];
    q->at(3 * n - 1) +=
        -2.0 * weight_end_state_[2] * end_state_ref_[2] / scale_factor_[2];
    }
}



OSQPData* PathRepairPlanner::FormulateProblem(vector<float> curvature, bool &boundErrorCheck) {
    // calculate kernel
    std::vector<c_float> P_data;
    std::vector<c_int> P_indices;
    std::vector<c_int> P_indptr;
    CalculateKernel(&P_data, &P_indices, &P_indptr);

    // calculate affine constraints
    std::vector<c_float> A_data;
    std::vector<c_int> A_indices;
    std::vector<c_int> A_indptr;
    std::vector<c_float> lower_bounds;
    std::vector<c_float> upper_bounds;
    CalculateAffineConstraint(&A_data, &A_indices, &A_indptr, &lower_bounds,
                            &upper_bounds, curvature, boundErrorCheck);

    OSQPData* data = reinterpret_cast<OSQPData*>(c_malloc(sizeof(OSQPData)));

    if(boundErrorCheck == false){
        cout << "bound error" << endl;
        return data;
    }

    // calculate offset
    //最后的参考线要考虑referline这一因素，即初始解
    //保证尽可能不要有太大偏差，这样有可能给车辆带来不稳定因素，这里主要是给目标函数进行补偿，目标函数的 ref一项
    std::vector<c_float> q;
    CalculateOffset(&q);

    
    if (lower_bounds.size() != upper_bounds.size()) {
    std::cout << "Formulate data failed!" << std::endl;
    }

    size_t kernel_dim = 3 * num_of_knots_;
    size_t num_affine_constraint = lower_bounds.size();

    data->n = kernel_dim;
    data->m = num_affine_constraint;
    data->P = csc_matrix(kernel_dim, kernel_dim, P_data.size(), CopyData(P_data),
                       CopyData(P_indices), CopyData(P_indptr));
    data->q = CopyData(q);
    data->A = csc_matrix(num_affine_constraint, kernel_dim, A_data.size(),
                 CopyData(A_data), CopyData(A_indices), CopyData(A_indptr));
    data->l = CopyData(lower_bounds);
    data->u = CopyData(upper_bounds);
    return data;
}

OSQPSettings* PathRepairPlanner::SolverDefaultSettings() {
    // Define Solver default settings
    OSQPSettings* settings =
      reinterpret_cast<OSQPSettings*>(c_malloc(sizeof(OSQPSettings)));
    osqp_set_default_settings(settings);
    settings->eps_abs = 5e-3;
    settings->eps_rel = 5e-3;
    settings->eps_prim_inf = 5e-4;
    settings->eps_dual_inf = 5e-4;
    settings->polish = true;
    settings->verbose = false;
    settings->scaled_termination = true;
    return settings;
}

bool PathRepairPlanner::Optimize(const int max_iter, vector<float> curvature, int &res_flag) {

    bool boundErrorCheck_ = true;

    OSQPData* data = FormulateProblem(curvature, boundErrorCheck_);//建造二次规划问题框架

    if(boundErrorCheck_ == false){
        res_flag = 2;
        return false;
    }

    OSQPSettings* settings = SolverDefaultSettings();
    settings->max_iter = max_iter;

    c_int status = 0;
    OSQPWorkspace* osqp_work = osqp_setup(data, settings);
    status = osqp_solve(osqp_work);

    cout << "status" << status << endl;

    // extract primal results
    x_.resize(num_of_knots_);
    dx_.resize(num_of_knots_);
    ddx_.resize(num_of_knots_);

    

    for (size_t i = 0; i < num_of_knots_; ++i) {
    x_.at(i) = osqp_work->solution->x[i] / scale_factor_[0];
    dx_.at(i) = osqp_work->solution->x[i + num_of_knots_] / scale_factor_[1];
    ddx_.at(i) = osqp_work->solution->x[i + 2 * num_of_knots_] / scale_factor_[2];
    }


    // Cleanup
    osqp_cleanup(osqp_work);
    FreeData(data);
    c_free(settings);

    res_flag = 1;

    return true;
}

void PathRepairPlanner::FreeData(OSQPData* data) {
    delete[] data->q;
    delete[] data->l;
    delete[] data->u;

    delete[] data->P->i;
    delete[] data->P->p;
    delete[] data->P->x;

    delete[] data->A->i;
    delete[] data->A->p;
    delete[] data->A->x;
}

void PathRepairPlanner::set_end_state_ref(
    const std::array<double, 3>& weight_end_state,
    const std::array<double, 3>& end_state_ref) {
    weight_end_state_ = weight_end_state;
    end_state_ref_ = end_state_ref;
    has_end_state_ref_ = true;
}

void PathRepairPlanner::set_x_ref(const double weight_x_ref,
                            std::vector<double> x_ref) {
    if (x_ref.size() != num_of_knots_) {
    std::cout << "x_ref size is wrong!" << std::endl;
    return;
    }
    weight_x_ref_ = weight_x_ref;
    x_ref_ = x_ref;
    has_x_ref_ = true;
}

void PathRepairPlanner::set_dx_bounds(double low_bound, double up_bound) {
    for (auto& x : dx_bounds_) {
    x.first = low_bound;
    x.second = up_bound;
    }
}

void PathRepairPlanner::set_ddx_bounds(double low_bound, double up_bound) {
    for (auto& x : ddx_bounds_) {
    x.first = low_bound;
    x.second = up_bound;
    }
}

void PathRepairPlanner::set_dddx_bounds(double low, double up) {
    dddx_bound_.first = low;
    dddx_bound_.second = up;
}

double PathRepairPlanner::GetPointHeading(double p_x, double p_y, double ori_x, double ori_y)//根据坐标x、y值计算其方位角
{
    double hAngle = 0;
    double dy = p_y - ori_y;
    double dx = p_x - ori_x;
    if (dx==0 && dy>0)
    {
        hAngle = 0;
    }
    else if(dx==0 && dy<0)
    {
        hAngle = 180;
    }
    else if(dy==0 && dx>0)
    {
        hAngle = 90;
    }
    else if(dy==0 && dx<0)
    {
        hAngle = 270;
    }
    else if(dx>0 && dy>0)//第一象限
    {
        hAngle = atan2(dx,dy)*180/M_PI;
    }
    else if(dx>0 && dy<0)//第二象限
    {
        hAngle = 180 - atan2(dx,-dy)*180/M_PI;
    }
    else if(dx<0 && dy<0)//第三象限
    {
        hAngle = 180 + atan2(-dx,-dy)*180/M_PI;
    }
    else if(dx<0 && dy>0)//第四象限
    {
        hAngle = 360 - atan2(-dx,dy)*180/M_PI;
    }

    return hAngle;
}

void PathRepairPlanner::SpiltLRPoints(const vector<GPSInfo> &main_path, const vector<GPSInfo> edge_points, const int start_key,
                                      vector<GPSInfo> *left_edge_points, vector<GPSInfo> *right_edge_points)
{
    //叉积划分左右边界点 上负下正 0306
    //若左右边界点按顺序排列 可以可视化观察出临界点 再手动进行区分
    pointVec main_path_Vec;
    point_t pt_main;
    GPSInfo start_gps = main_path.at(start_key);
    for (size_t i = 0; i < main_path.size(); i++) {
        float main_point_x, main_point_y;
        Tage::transGpsIntoInertia(start_gps, main_path[i],&main_point_x,&main_point_y);
        pt_main.push_back(main_point_x);
        pt_main.push_back(main_point_y);
        main_path_Vec.push_back(pt_main);
        pt_main.clear();
    }
    KDTree MainPathKDTree(main_path_Vec);
    std::pair<double,double> vector1;//向量1
    std::pair<double,double> vector2;//向量2
    for (size_t i = 0;i < 10; i++){
        (*left_edge_points).push_back(edge_points[i]);
    }
    for (size_t i = 10;i < edge_points.size(); i++) {
        float edge_point_x, edge_point_y, main_point_x_20, main_point_y_20, main_point_x, main_point_y;
        Tage::transGpsIntoInertia(start_gps, edge_points[i], &edge_point_x, &edge_point_y);
        pt_main.push_back(edge_point_x);
        pt_main.push_back(edge_point_y);
        auto project_index = MainPathKDTree.nearest_index(pt_main); //边界路径投影到最近主路径点的索引值
        pt_main.clear();
        Tage::transGpsIntoInertia(start_gps, main_path[project_index-20], &main_point_x_20, &main_point_y_20);
        Tage::transGpsIntoInertia(start_gps, main_path[project_index], &main_point_x, &main_point_y);

        vector1.first = edge_point_x - main_point_x_20;
        vector1.second = edge_point_y - main_point_y_20;
        vector2.first = main_point_x - main_point_x_20;
        vector2.second = main_point_y - main_point_y_20;
        //叉积 逆时针为正 顺时针为负 右手定则
        if((vector1.first * vector2.second - vector2.first * vector1.second) > 0)
        {
            (*left_edge_points).push_back(edge_points[i]);
        }
        else
        {
            (*right_edge_points).push_back(edge_points[i]);
        }

    }
}
void PathRepairPlanner::VerifyResultCurvatureCompany(GPSInfo start_gps, vector<GPSInfo> refined_trajectory)
{
    double iTemp = 0, length = 0, curvature;
    //计算曲率
    for(size_t i =0;i<refined_trajectory.size()-3;i++)
    {
        for(int n=i;n<i+4;n++)
        {
            float sub=refined_trajectory[n+1].iHead-refined_trajectory[n].iHead;
            iTemp+=sub;
            float nx,ny,n1x,n1y;
            Tage::transGpsIntoInertia(start_gps, refined_trajectory[n], &nx, &ny);
            Tage::transGpsIntoInertia(start_gps, refined_trajectory[n+1], &n1x, &n1y);
            length+=sqrt(pow(n1x-nx,2) + pow(n1y-ny, 2));
        }
        if (iTemp > 180) iTemp = iTemp - 360;
        if (iTemp < -180) iTemp = iTemp + 360;
        curvature=iTemp*M_PI/180/length;
    }
}

void PathRepairPlanner::VerifyResultCurvature(GPSInfo start_gps, vector<GPSInfo> refined_trajectory)
{
    cout<<"验证曲率约束中..."<<endl;
    GPSInfo cal_curvature;
    float next_x, next_y, now_x, now_y;
    double pre_x,pre_y,pre_head;
    std::vector<float> result_points_curvature;
    result_points_curvature.resize(refined_trajectory.size());
    for(std::size_t i = 0; i< refined_trajectory.size()-1 ; i++)
    {
        double iTemp = 0;
        double lengthTmp = 0;
        if(i < refined_trajectory.size() - 4){
            for(int n = i; n < i + 4; n++){
                float sub = refined_trajectory[n + 1].iHead - refined_trajectory[n].iHead;
                iTemp += sub;
                Tage::transGpsIntoInertia(start_gps, refined_trajectory[n], &now_x, &now_y);
                Tage::transGpsIntoInertia(start_gps, refined_trajectory[n+1], &next_x, &next_y);
                double dis_tmp = hypot(next_x - now_x, next_y - now_y);
                lengthTmp += dis_tmp;
            }
            result_points_curvature[i] = iTemp * M_PI / 180 / lengthTmp;
        }
        else{
            result_points_curvature[i] = 0;
        }

    }

    result_points_curvature[0] = result_points_curvature[1];
    result_points_curvature[refined_trajectory.size()-1] = result_points_curvature[refined_trajectory.size()-2];
    size_t j = 0;
    double max_curvature = 0;
    cout<<"曲率大于0.1输出:"<<endl;
    for (size_t i = 0;i < refined_trajectory.size();i++) {
        if(result_points_curvature[i] > max_curvature)
        {
            max_curvature = result_points_curvature[i];
        }
        if(result_points_curvature[i]>0.1)
        {
            cout<<"注意:第"<<i<<"个点曲率为"<<result_points_curvature[i]<<endl;
        }
        else
        {
            j++;
        }
    }
    cout<<"最大曲率为:"<<max_curvature<<endl;
    cout<<"曲率约束验证完成";
    if(j==refined_trajectory.size())
    {
        cout<<" 所有曲率均小于等于0.1"<<endl;cout<<endl;
    }
    else
    {
        cout<<endl;cout<<endl;
    }
}

std::vector<float> PathRepairPlanner::CalRawPathCurvature(GPSInfo start_gps, vector<GPSInfo> raw_path_points)
{
    cout<<"开始计算raw_path曲率"<<endl;
    GPSInfo cal_curvature;
    float next_x, next_y, now_x, now_y;
    std::vector<float> raw_points_curvature;
    raw_points_curvature.resize(raw_path_points.size());
    for(std::size_t i = 0; i< raw_path_points.size()-1 ; i++)
    {
        double iTemp = 0;
        double lengthTmp = 0;
        if(i < 20)
        {
            raw_points_curvature[i] = 0;
        }
        else if(i < raw_path_points.size() - 20)
        {
            for(int n = i - 10; n < i + 10; n++){
                float sub = raw_path_points[n + 1].iHead - raw_path_points[n].iHead;
                iTemp += sub;
                Tage::transGpsIntoInertia(start_gps, raw_path_points[n], &now_x, &now_y);
                Tage::transGpsIntoInertia(start_gps, raw_path_points[n+1], &next_x, &next_y);
                double dis_tmp = hypot(next_x - now_x, next_y - now_y);
                lengthTmp += dis_tmp;
            }
            if (iTemp > 180) iTemp = iTemp - 360;
            if (iTemp < -180) iTemp = iTemp + 360;
            raw_points_curvature[i] = iTemp * M_PI / 180 / lengthTmp;
        }
        else{
            raw_points_curvature[i] = 0;
        }

    }
    raw_points_curvature[0] = raw_points_curvature[1];
    raw_points_curvature[raw_path_points.size()-1] = raw_points_curvature[raw_path_points.size()-2];
        
    double max_curvature = 0;
    for (size_t i = 0;i < raw_points_curvature.size();i++) {
        if(raw_points_curvature[i] > max_curvature)
        {
            max_curvature = raw_points_curvature[i];
        }
    }

    cout<<"优化前路径最大曲率为"<<max_curvature<<endl;
    return raw_points_curvature;
}

bool PathRepairPlanner::resultCheck(vector<GPSInfo> *refined_trajectory, int &res_flag)
{
    for(size_t i = 0; i < refined_trajectory->size(); i++){
        if(refined_trajectory->at(i).iLatitude > 180 || refined_trajectory->at(i).iLatitude < 10 ||
           refined_trajectory->at(i).iLongitude > 180  || refined_trajectory->at(i).iLongitude < 10 )
        {
            res_flag = 3;
            return false;
        }
    }
    res_flag = 1;
    return true;
}


void PathRepairPlanner::CoutResultXY(GPSInfo start_gps, vector<GPSInfo> result)
{
    
}

bool PathRepairPlanner::plan(const vector<GPSInfo> &main_path, const vector<GPSInfo> &left_edge_points, const vector<GPSInfo> &right_edge_points,
                             const int start_key, const int end_key, VehicleParameter m_vehicle,
                             vector<GPSInfo> *refined_trajectory,vector<GPSInfo>*left_edge_trajectory,vector<GPSInfo>*right_edge_trajectory,
                             std::vector<std::pair<double,double>> &edge_bound, double boundSetting, int &result_flag,
                             std::vector<GPSInfo> &error_point_se, std::vector<GPSInfo> &error_point_mid)
{
    double default_boundary = m_vehicle.vehicle_width / 2 + boundSetting;


    // use the start gps as original point
    GPSInfo start_gps = main_path.at(start_key);
    vector<GPSInfo>::const_iterator start_iter = main_path.begin() + start_key;
    vector<GPSInfo>::const_iterator end_iter = main_path.begin() + end_key;
    vector<GPSInfo> raw_path_points(start_iter, end_iter);//左闭右开

    for(auto& path_point: raw_path_points)
    {
        float path_point_x, path_point_y;
        Tage::transGpsIntoInertia(start_gps, path_point, &path_point_x, &path_point_y);
        path_point.xEN = path_point_x;
        path_point.yEN = path_point_y;

        float temp_theta = 0.0;
        float bias_head = path_point.iHead - start_gps.iHead;
        if(bias_head > 90)
        {
            temp_theta = (450 - bias_head)/180.0*M_PI;
        }else
        {
            temp_theta = (90 - bias_head)/180.0*M_PI;
        }

        if(temp_theta > M_PI)
            temp_theta -=2*M_PI;

        path_point.theta = temp_theta;
    }

/*****************************************   frenet-l  ***********************************************/
    //主路径建KDTree
    pointVec main_path_Vec;
    point_t pt_main;

    for (size_t i = 0; i < main_path.size(); i++) {
        float main_point_x, main_point_y;
        Tage::transGpsIntoInertia(start_gps, main_path[i],&main_point_x,&main_point_y);
        pt_main.push_back(main_point_x);
        pt_main.push_back(main_point_y);
        main_path_Vec.push_back(pt_main);
        pt_main.clear();
    }
    KDTree MainPathKDTree(main_path_Vec);

    //左边界建KDTree
    pointVec left_path_Vec;
    point_t pt_left;

    for (size_t i = 0; i < left_edge_points.size(); i++) {
        float left_point_x, left_point_y;
        Tage::transGpsIntoInertia(start_gps, left_edge_points[i],&left_point_x,&left_point_y);
        pt_left.push_back(left_point_x);
        pt_left.push_back(left_point_y);
        left_path_Vec.push_back(pt_left);
        pt_left.clear();
    }
    KDTree LeftPathKDTree(left_path_Vec);

    //右边界建KDTree
    pointVec right_path_Vec;
    point_t pt_right;

    for (size_t i = 0; i < right_edge_points.size(); i++) {
        float right_point_x, right_point_y;
        Tage::transGpsIntoInertia(start_gps, right_edge_points[i],&right_point_x,&right_point_y);
        pt_right.push_back(right_point_x);
        pt_right.push_back(right_point_y);
        right_path_Vec.push_back(pt_right);
        pt_right.clear();
    }
    KDTree RightPathKDTree(right_path_Vec);

    //建立待优化raw路径KDTree
    pointVec raw_path_Vec;
    point_t pt_raw;

    for (size_t i = 0; i <raw_path_points.size(); i++) {
        float raw_point_x, raw_point_y;
        Tage::transGpsIntoInertia(start_gps, raw_path_points[i],&raw_point_x,&raw_point_y);
        pt_raw.push_back(raw_point_x);
        pt_raw.push_back(raw_point_y);
        raw_path_Vec.push_back(pt_raw);
        pt_raw.clear();
    }
    KDTree RawPathKDTree(raw_path_Vec);
    cout<<"完成建立KDTree"<<endl;

    //初始化边界范围 车宽一半加1m
    edge_bound.resize(raw_path_points.size());
    cout<<"开始设置左右边界默认值——距离边界最近点距离"<<endl;
    //利用KDTree求出距离各个边界点最近的主路径点 投影计算边界值 frenet-l —— hzx0222 ——
    /************将始终为初始值未设置横向距离的主路径点在边界上找最近点进行距离替换(左)****************/
    point_t pt_project_left;
    double left_dis;
    for (size_t i = 0; i < raw_path_points.size(); i++) {
            float path_point_x, path_point_y;
            Tage::transGpsIntoInertia(start_gps, raw_path_points[i], &path_point_x, &path_point_y);
            pt_project_left.push_back(path_point_x);
            pt_project_left.push_back(path_point_y);
            auto edge_nearest_point = LeftPathKDTree.nearest_point(pt_project_left);
            pt_project_left.clear();
            left_dis =(-1)*sqrt(pow(edge_nearest_point[0]-path_point_x,2) + pow(edge_nearest_point[1]-path_point_y, 2));
            edge_bound[i].first = left_dis;
        }
    cout<<"左--设置默认值完成"<<endl;

    /************将始终为初始值未设置横向距离的主路径点在边界上找最近点进行距离替换(右)****************/
    point_t pt_project_right;
    double right_dis;
    for (size_t i = 0; i < raw_path_points.size(); i++) {
            float path_point_x, path_point_y;
            Tage::transGpsIntoInertia(start_gps, raw_path_points[i], &path_point_x, &path_point_y);
            pt_project_right.push_back(path_point_x);
            pt_project_right.push_back(path_point_y);
            auto edge_nearest_point = RightPathKDTree.nearest_point(pt_project_right);
            pt_project_right.clear();
            right_dis =sqrt(pow(edge_nearest_point[0]-path_point_x,2) + pow(edge_nearest_point[1]-path_point_y, 2));
            edge_bound[i].second = right_dis;
    }
    cout<<"右--设置默认值完成"<<endl;cout<<endl;

    vector<double> tmp;
    tmp.resize(edge_bound.size());
    for(size_t i = 0; i < tmp.size(); ++i){
        tmp[i] = edge_bound[i].first;
    }
    DataToTxt("./data/ref_edge_l.txt", tmp);
    for(size_t i = 0; i < tmp.size(); ++i){
        tmp[i] = edge_bound[i].second;
    }
    DataToTxt("./data/ref_edge_r.txt", tmp);

    cout << "boundary" << default_boundary << endl;

    vector<int> error_point_index_sAnde;
    vector<int> error_point_index_mid;

    for(int i = 0; i < 50; i++){
        if(fabs(edge_bound[i].first) < default_boundary || fabs(edge_bound[i].second) < default_boundary){
            error_point_index_sAnde.push_back(i);
        }
    }

    for(size_t i = 0;i < raw_path_points.size();i++){
        double width = fabs(edge_bound[i].first) + fabs(edge_bound[i].second);
        if(width < 2 * default_boundary - 1){
            error_point_index_mid.push_back(i);
        }
    }

    for(int i = raw_path_points.size() - 50; i < raw_path_points.size(); i++){
        if(fabs(edge_bound[i].first) < default_boundary || fabs(edge_bound[i].second) < default_boundary){
            error_point_index_sAnde.push_back(i);
        }
    }

    cout << "error size se " << error_point_index_sAnde.size() << endl;
    cout << "error size mid " << error_point_index_mid.size() << endl;

    if(error_point_index_mid.size()!=0 || error_point_index_sAnde.size()!=0){

        for(size_t i = 0; i < error_point_index_mid.size(); i++){
            error_point_mid.push_back(raw_path_points[error_point_index_mid[i]]);
        }

        for(size_t i = 0; i < error_point_index_sAnde.size(); i++){
            error_point_se.push_back(raw_path_points[error_point_index_sAnde[i]]);
        }
        result_flag = 4;
        return false;
    }

    cout<<"开始碰撞检测"<<endl;
    //检测碰撞，设置边界
    for(size_t i = 0;i < raw_path_points.size();i++)
    {
        //左边界在左
        if(edge_bound[i].first < 0)
        {
            if( edge_bound[i].first > -1 * default_boundary )
            {
                edge_bound[i].first = (default_boundary) - fabs(edge_bound[i].first);
            }
        }
        //右边界在右
        if(edge_bound[i].second > 0)
        {
            if( edge_bound[i].second < default_boundary )
            {
                edge_bound[i].second = (-1)* ( (default_boundary) - edge_bound[i].second );
            }
        }
        //右边界在左
    }
    cout<<"完成碰撞检测，并设置碰撞点边界值"<<endl;cout<<endl;


    tmp.resize(edge_bound.size());
    for(size_t i = 0; i < tmp.size(); ++i){
        tmp[i] = edge_bound[i].first;
    }
    DataToTxt("./data/update_edge_l.txt", tmp);
    for(size_t i = 0; i < tmp.size(); ++i){
        tmp[i] = edge_bound[i].second;
    }
    DataToTxt("./data/update_edge_r.txt", tmp);


/*****************************************   frenet-l   **********************************************/

    std::vector<std::pair<double,double>> edge_v_bound; // raw velocity point set in the cartian coordinate system

    std::vector<double> init_x_solution;

    edge_v_bound.resize(raw_path_points.size());
    init_x_solution.resize(raw_path_points.size());

    for (size_t i =0; i <raw_path_points.size() ;++i) {
        edge_v_bound[i]=std::pair<double,double> (-1.0, 1.0);
        init_x_solution[i]=0.0;
    }

    double init_l = 0.0;
    double end_l = 0.0;
    std::array<double, 3> init_state = {init_l, 0.0, 0.0};   // 这个代表了车辆的初始横向状态，这里因为是主路径，所以都设为0
    PathRepairPlanner prob(raw_path_points.size(), 0.2, init_state);// 初始赋值
    std::array<double, 3> end_state = {end_l, 0.0, 0.0};   //这里代表车辆规划结束的终止状态，这里因为是主路径，所以都设为0

    prob.set_end_state_ref({1000.0, 0.0, 0.0}, end_state); //这个代表了不同目标的权重，这里期望终止点的横向位移为0

    prob.set_weight_x(2.0);
    prob.set_weight_dx(10.0);
    prob.set_weight_ddx(100.0);
    prob.set_weight_dddx(1000.0);

    prob.set_scale_factor({1.0, 10.0, 100.0});

    prob.set_x_bounds(edge_bound);
    prob.set_dx_bounds(edge_v_bound);
    prob.set_ddx_bounds(-2.0, 2.0);
    prob.set_dddx_bounds(-4.0, 4.0);

    prob.set_x_ref(10, init_x_solution);

    std::vector<Coffient> road_coffient;

    //计算待修正路径曲率-后续曲率约束
    std::vector<float> raw_points_curvature = CalRawPathCurvature(start_gps,raw_path_points);

    DataToTxt("./data/curvate.txt", raw_points_curvature);

    if(prob.Optimize(10000, raw_points_curvature, result_flag)){

        std::cout<<"优化成功！！！"<<std::endl;

        DataToTxt("./data/refine_s.txt", prob.x_);

        if (refined_trajectory->empty())
        {
            generatePath(0, raw_path_points, prob.x_, road_coffient, *refined_trajectory, raw_points_curvature);
        }
     
    }
    else
    {
      std::cout<<"Optimize failed!!"<<std::endl;cout<<endl;
      result_flag = 2;
      return false;
    }

    //验证优化后曲率
    //VerifyResultCurvature(start_gps, *refined_trajectory);

    for (size_t i =0; i < refined_trajectory->size() ;++i) {

        GPSInfo temp_point = Tage::transXYToGps(start_gps.iLatitude,start_gps.iLongitude, start_gps.iHead,
                                                refined_trajectory->at(i).xEN,refined_trajectory->at(i).yEN);
        refined_trajectory->at(i).iLongitude = temp_point.iLongitude;
        refined_trajectory->at(i).iLatitude = temp_point.iLatitude;
    }
    if(!resultCheck(refined_trajectory, result_flag)){
            return false;
    }

    std::cout<<"完成主路径修正工作"<<std::endl;
    return true;
}
