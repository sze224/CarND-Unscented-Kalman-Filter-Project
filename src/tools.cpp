#include <iostream>
#include "tools.h"

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;
using namespace std;
Tools::Tools() {}

Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) {
    
    VectorXd RMSE(4);
    RMSE << 0., 0., 0., 0.;
    
    /* check if estimation is empty or if ground_truth and estimations are different length */
    if (estimations.size() == 0 || estimations.size() != ground_truth.size()){
        cout << "Error: Estimation is a size 0 or Estimations and ground truth are different size" << endl;
    }else{
        /* if the inputs look good, compute the RMSE */
        for (int i = 0; i < estimations.size(); i++){
            VectorXd residual = estimations[i] - ground_truth[i];
            residual = residual.array() * residual.array();
            RMSE = RMSE + residual;
        }
        RMSE = RMSE / estimations.size();
        RMSE = RMSE.array().sqrt();
    }
    return RMSE;
}


VectorXd Tools::Polar2Cartesian(const Eigen::VectorXd& Polar_data){
    
    VectorXd c_state(4);
    
    double rho = Polar_data(0);
    double phi = Polar_data(1);
    double rho_dot = Polar_data(2);
    
    double px = rho * cos(phi);
    double py = rho * sin(phi);
    double vx = rho_dot * cos(phi);
    double vy = rho_dot * sin(phi);
    
    c_state << px, py, vx, vy;
    
    return c_state;
    
}
