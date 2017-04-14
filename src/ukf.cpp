#include "ukf.h"
#include "tools.h"
#include "Eigen/Dense"
#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

/**
 * Initializes Unscented Kalman filter
 */
UKF::UKF() {
    
  // if this is false, this will be the first measurement that is coming in
  is_initialized_ = false;
    
  // if this is false, laser measurements will be ignored (except during init)
  use_laser_ = true;

  // if this is false, radar measurements will be ignored (except during init)
  use_radar_ = true;

  // initial state vector
  x_ = VectorXd(5);

  // initial covariance matrix
  P_ = MatrixXd(5, 5);

  // Process noise standard deviation longitudinal acceleration in m/s^2
  std_a_ = 0.5;//0.085;

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = 0.5;//0.008;

  // Laser measurement noise standard deviation position1 in m
  std_laspx_ = 0.15;

  // Laser measurement noise standard deviation position2 in m
  std_laspy_ = 0.15;

  // Radar measurement noise standard deviation radius in m
  std_radr_ = 0.3;

  // Radar measurement noise standard deviation angle in rad
  std_radphi_ = 0.03;

  // Radar measurement noise standard deviation radius change in m/s
  std_radrd_ = 0.3;

  /**
  TODO:

  Complete the initialization. See ukf.h for other member properties.

  Hint: one or more values initialized above might be wildly off...
  */
    
  // dimension of the states vector (px, py, v, psi, psi_d)
  n_x_ = 5;

  // dimension of the states vector with augmented states (px, py, v, psi, psi_d, siga, sigpsid)
  n_aug_ = 7;
    
  // dimension of radar measurement dimension (rho, phi, rho_dot)
  n_r_ = 3;
    
  // dimension of lidar measurement dimension (px, py)
  n_l_ = 2;
    
  // tunable parameter for spreading sigma points
  lambda_ = 3 - n_aug_;
 
  // initialize weights
  weights_ = VectorXd(2 * n_aug_ + 1);
  
  // state sigma prediction
  Xsig_pred = MatrixXd(n_x_, 2 * n_aug_ + 1);
    
}

UKF::~UKF() {}

/**
 * @param {MeasurementPackage} meas_package The latest measurement data of
 * either radar or laser.
 */
void UKF::ProcessMeasurement(MeasurementPackage meas_package) {
  /**
  TODO:

  Complete this function! Make sure you switch between lidar and radar
  measurements.
  */
    Tools tools;
    
    if (!is_initialized_){
        x_ << 0,0,0,0,0;
        
        P_ << 1,0,0,0,0,
              0,1,0,0,0,
              0,0,1,0,0,
              0,0,0,1,0,
              0,0,0,0,1;
        
        if (meas_package.sensor_type_ == MeasurementPackage::LASER){
            float px_i = meas_package.raw_measurements_[0];
            float py_i = meas_package.raw_measurements_[1];
            if (fabs(px_i) < 0.001){
              px_i = 0.001;
            }
            if (fabs(py_i) < 0.001){
              py_i = 0.001;
            }
            x_ << px_i, py_i, 0, 0, 0;
        }
        
        if (meas_package.sensor_type_ == MeasurementPackage::RADAR){
            VectorXd c_state = tools.Polar2Cartesian(meas_package.raw_measurements_);
               x_ << c_state(0), c_state(1), 0, 0, 0;
            
        }
        
        time_us_ = meas_package.timestamp_;
        is_initialized_ = true;
        return;
    }
    
    double dt = (meas_package.timestamp_ - time_us_) / 1000000.0 ;
    time_us_ = meas_package.timestamp_;

    while (dt > 0.1){
      const double del_t = 0.05; 
      Prediction(del_t);
      dt = dt - del_t;
    }
    Prediction(dt);
    
    if (meas_package.sensor_type_ == MeasurementPackage::LASER){
      UpdateLidar(meas_package);
    }
    if (meas_package.sensor_type_ == MeasurementPackage::RADAR){
      UpdateRadar(meas_package);
    }
}

/**
 * Predicts sigma points, the state, and the state covariance matrix.
 * @param {double} delta_t the change in time (in seconds) between the last
 * measurement and this one.
 */
void UKF::Prediction(double delta_t) {
  /**
  TODO:

  Complete this function! Estimate the object's location. Modify the state
  vector, x_. Predict sigma points, the state, and the state covariance matrix.
  */

    // generate augmented sigma points
    VectorXd x_aug = VectorXd(n_aug_);
    MatrixXd P_aug = MatrixXd(n_aug_, n_aug_);
    MatrixXd Xsig_aug = MatrixXd(n_aug_, 2 * n_aug_ + 1);
    x_aug.fill(0.0);
    Xsig_aug.fill(0.0);
    P_aug.fill(0.0);

    x_aug.head(5) = x_;
    x_aug(5) = 0;
    x_aug(6) = 0;
    while (x_aug(3)> M_PI) x_aug(3) -= 2.*M_PI;
    while (x_aug(3)<-M_PI) x_aug(3) += 2.*M_PI;


    P_aug.topLeftCorner(5,5) = P_;
    P_aug(5,5) = std_a_ * std_a_;
    P_aug(6,6) = std_yawdd_ * std_yawdd_;
    
    MatrixXd L = P_aug.llt().matrixL();
         
    Xsig_aug.col(0) = x_aug;
    for(int i = 0; i < n_aug_; i++){
        Xsig_aug.col(i+1) = x_aug + sqrt(lambda_ + n_aug_) * L.col(i);
        Xsig_aug.col(i+1+n_aug_) = x_aug - sqrt(lambda_ + n_aug_) * L.col(i);
    }    

    // predict sigma points
    for (int i = 0; i < 2 * n_aug_ + 1; i++){
        double p_x = Xsig_aug(0, i);
        double p_y = Xsig_aug(1, i);
        double v  = Xsig_aug(2, i);
        double yaw = Xsig_aug(3, i);
        double yawd = Xsig_aug(4, i);
        double nu_a = Xsig_aug(5, i);
        double nu_yawdd = Xsig_aug(6, i);
        
        double px_p, py_p;
        
        if (fabs(yawd) > 0.001){
            px_p = p_x + v/yawd * (sin(yaw + yawd * delta_t) - sin(yaw));
            py_p = p_y + v/yawd * (cos(yaw) - cos(yaw + yawd * delta_t));
        }else{
            px_p = p_x + v * cos(yaw) * delta_t;
            py_p = p_y + v * sin(yaw) * delta_t;
        }

        double v_p = v; 
        double yaw_p = yaw + yawd * delta_t;
        double yawd_p = yawd;

        px_p = px_p + 0.5 * delta_t * delta_t * cos(yaw) * nu_a;
        py_p = py_p + 0.5 * delta_t * delta_t * sin(yaw) * nu_a;
        v_p = v_p + delta_t * nu_a;
        yaw_p = yaw_p + 0.5 * delta_t * delta_t * nu_yawdd;
        yawd_p = yawd_p + delta_t * nu_yawdd;
        
        Xsig_pred(0, i) = px_p;
        Xsig_pred(1, i) = py_p;
        Xsig_pred(2, i) = v_p;
        Xsig_pred(3, i) = yaw_p;
        Xsig_pred(4, i) = yawd_p;
    }
    
    //cout << Xsig_pred << endl;
    // Predict mean and convariance of state sigma points
    
    // set weights
    double weight_0 = lambda_/(lambda_ + n_aug_);
    weights_(0) = weight_0;
    for (int i = 1; i < 2 * n_aug_ + 1; i++){
        double weight = 0.5 / (lambda_ + n_aug_);
        weights_(i) =  weight;
    }
    
    // get predicted state mean
    x_.fill(0.0);
    for (int i = 0; i < 2 * n_aug_ + 1; i++){
        x_ = x_ + weights_(i) * Xsig_pred.col(i);
    }
    while (x_(3)> M_PI) x_(3) -= 2.*M_PI;
    while (x_(3)<-M_PI) x_(3) += 2.*M_PI;
    
    // get predicted state convariance
    P_.fill(0.0);
    for (int i = 0; i < 2 * n_aug_ + 1; i++){
        VectorXd x_diff = Xsig_pred.col(i) - x_;
        while (x_diff(3)> M_PI) x_diff(3) -= 2.*M_PI;
        while (x_diff(3)<-M_PI) x_diff(3) += 2.*M_PI;
        P_ = P_ + weights_(i) * x_diff * x_diff.transpose();
    }
}


/**
 * Updates the state and the state covariance matrix using a laser measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateLidar(MeasurementPackage meas_package) {
  /**
  TODO:

  Complete this function! Use lidar data to update the belief about the object's
  position. Modify the state vector, x_, and covariance, P_.

  You'll also need to calculate the lidar NIS.
  */
  // get the next Lidar measurement
  VectorXd z_Lidar = VectorXd(2);
  double px = meas_package.raw_measurements_[0];
  double py = meas_package.raw_measurements_[1];

  z_Lidar << px, py;
  
  MatrixXd Zsig = MatrixXd(n_l_, 2 * n_aug_ + 1);
  Zsig.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; i++){
    Zsig(0, i) = Xsig_pred(0, i);
    Zsig(1, i) = Xsig_pred(1, i);
  }

  VectorXd z_pred = VectorXd(n_l_);
  z_pred.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; i++){
    z_pred = z_pred + weights_(i) * Zsig.col(i);
  }

  MatrixXd S = MatrixXd(n_l_, n_l_);
  S.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; i++){
    VectorXd z_diff = Zsig.col(i) - z_pred;
    S = S + weights_(i) * z_diff * z_diff.transpose();
  }

  MatrixXd R = MatrixXd(n_l_,n_l_);
  R << std_laspx_*std_laspx_, 0,
       0, std_laspy_*std_laspy_;
  S = S + R;
  
  MatrixXd Tc = MatrixXd(n_x_, n_l_);
  Tc.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; i++){
    VectorXd z_diff = Zsig.col(i) - z_pred;
    
    VectorXd x_diff = Xsig_pred.col(i) - x_;
    while (x_diff(3)> M_PI) x_diff(3) -= 2.*M_PI;
    while (x_diff(3)<-M_PI) x_diff(3) += 2.*M_PI;
    
    Tc = Tc + weights_(i) * x_diff * z_diff.transpose();
  }

  VectorXd z_diff = z_Lidar -  z_pred;

  MatrixXd K = Tc * S.inverse();
  x_ = x_ + K * z_diff;
  while(x_(3)> M_PI) x_(3) -= 2.*M_PI;
  while(x_(3)<-M_PI) x_(3) += 2.*M_PI;
  P_ = P_ - K * S * K.transpose();


  NIS_laser_ = z_diff.transpose() * S.inverse() * z_diff;
}

/**
 * Updates the state and the state covariance matrix using a radar measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateRadar(MeasurementPackage meas_package) {
  /**
  TODO:

  Complete this function! Use radar data to update the belief about the object's
  position. Modify the state vector, x_, and covariance, P_.

  You'll also need to calculate the radar NIS.
  */
    // get current reading (rho, phi, rho_dot)
    VectorXd z_radar = VectorXd(3);
    z_radar << meas_package.raw_measurements_[0], meas_package.raw_measurements_[1], meas_package.raw_measurements_[2];
    
    // convert state sigma points into measurement pts
    MatrixXd Zsig = MatrixXd(n_r_, 2 * n_aug_ + 1);
    Zsig.fill(0.0);
    for (int i = 0; i < 2 * n_aug_ + 1; i++){
        double p_x = Xsig_pred(0, i);
        double p_y = Xsig_pred(1, i);
        double v = Xsig_pred(2, i);
        double yaw = Xsig_pred(3, i);
        
        double v1 = v * cos(yaw);
        double v2 = v * sin(yaw);

        if (fabs(p_x) < 0.001 && fabs(p_y) < 0.001){
          Zsig(0,i) = 0;
          Zsig(1,i) = 0;
          Zsig(2,i) = 0;
        }else{
          Zsig(0,i) = sqrt(p_x * p_x + p_y * p_y);
          Zsig(1,i) = atan2(p_y,p_x);
          Zsig(2,i) = (p_x * v1 + p_y * v2) / sqrt(p_x * p_x + p_y * p_y);
        } 
    }
    
    // get measurement means from measurement pts
    VectorXd z_pred = VectorXd(n_r_);
    z_pred.fill(0.0);
    for (int i = 0; i < 2 * n_aug_ + 1; i++){
        z_pred = z_pred + weights_(i) * Zsig.col(i);
    }
    
    // get measurement convariance from measurement pts
    MatrixXd S = MatrixXd(n_r_, n_r_);
    S.fill(0.0);
    for (int i = 0; i < 2 * n_aug_ + 1; i++){
        VectorXd z_diff = Zsig.col(i) - z_pred;
        while (z_diff(1)> M_PI) z_diff(1) -= 2.*M_PI;
        while (z_diff(1)<-M_PI) z_diff(1) += 2.*M_PI;
        S = S + weights_(i) * z_diff * z_diff.transpose();
    }

    MatrixXd R = MatrixXd(n_r_,n_r_);
    R << std_radr_*std_radr_,0,0,
         0,std_radphi_*std_radphi_,0,
         0,0,std_radrd_*std_radrd_;
    S = S + R;
    
    // cross correlation
    MatrixXd Tc = MatrixXd(n_x_, n_r_);
    Tc.fill(0.0);
    for (int i = 0; i < 2 * n_aug_ + 1; i++){
        VectorXd z_diff = Zsig.col(i) - z_pred;
        while (z_diff(1)> M_PI) z_diff(1) -= 2.*M_PI;
        while (z_diff(1)<-M_PI) z_diff(1) += 2.*M_PI;
        
        VectorXd x_diff = Xsig_pred.col(i) - x_;
        while (x_diff(3)> M_PI) x_diff(3) -= 2.*M_PI;
        while (x_diff(3)<-M_PI) x_diff(3) += 2.*M_PI;
        
        Tc = Tc + weights_(i) * x_diff * z_diff.transpose();
    }
    
    VectorXd z_diff = z_radar - z_pred;
    
    while (z_diff(1)> M_PI) z_diff(1) -= 2.*M_PI;
    while (z_diff(1)<-M_PI) z_diff(1) += 2.*M_PI;
    
    MatrixXd K = Tc * S.inverse();
    x_ = x_ + K * z_diff;
    while (x_(3)> M_PI) x_(3) -= 2.*M_PI;
    while (x_(3)<-M_PI) x_(3) += 2.*M_PI;
    P_ = P_ - K * S * K.transpose();

    NIS_radar_ = z_diff.transpose() * S.inverse() * z_diff;
}
