#include <math.h>
#include "kalman_filter.h"

using Eigen::MatrixXd;
using Eigen::VectorXd;
//using std::cout;
//using std::endl;
/*
 * Please note that the Eigen library does not initialize 
 *   VectorXd or MatrixXd objects with zeros upon creation.
 */
#define PI acos(-1)
const float PI2 = 2*PI;

KalmanFilter::KalmanFilter() {}

KalmanFilter::~KalmanFilter() {}

void KalmanFilter::Init(VectorXd &x_in, MatrixXd &P_in, MatrixXd &F_in,
                        MatrixXd &H_in, MatrixXd &R_in, MatrixXd &Q_in) {
  x_ = x_in;
  P_ = P_in;
  F_ = F_in;
  H_ = H_in;
  R_ = R_in;
  Q_ = Q_in;
}

void KalmanFilter::Predict() {
  /**
   * TODO: predict the state
   */
  x_ = F_ * x_;
  MatrixXd Ft = F_.transpose();
  P_ = F_ * P_ * Ft + Q_;
}

void KalmanFilter::Update(const VectorXd &z) {
  /**
   * TODO: update the state by using Kalman Filter equations
   */
  VectorXd z_pred = H_ * x_;
  VectorXd y = z - z_pred;
  MatrixXd Ht = H_.transpose();
  MatrixXd S = H_ * P_ * Ht + R_;
  MatrixXd Si = S.inverse();
  MatrixXd PHt = P_ * Ht;
  MatrixXd K = PHt * Si;
  
  x_ = x_ + (K * y);
  long x_size = x_.size();
  MatrixXd I = MatrixXd::Identity(x_size, x_size);
  P_ = (I - K * H_) * P_;
}

VectorXd CartesianToPolar(const VectorXd &x_state){
  //h(x) function to convert cartesian coordinates in x_ to polar coordinates
  VectorXd z_pred(3);
  float px, py, vx, vy;
  px = x_state[0];
  py = x_state[1];
  vx = x_state[2];
  vy = x_state[3];
  
  float rho, phi, rho_dot;
  rho = sqrt(px*px + py*py);
  phi = atan2(py, px);
  
  if (rho < 0.000001){
    //cout << "rho_dot - Error - Division by Zero" << endl;
    rho = 0.000001;
  }
  rho_dot = (px*vx + py*vy) / rho;
  z_pred << rho, phi, rho_dot;
  
  return z_pred;
}

void KalmanFilter::UpdateEKF(const VectorXd &z) {
  /**
   * TODO: update the state by using Extended Kalman Filter equations
   */
  // use h(x) function to calculate z_pred, which is CartesianToPoloar here
  VectorXd z_pred = CartesianToPolar(x_);
  VectorXd y = z - z_pred;
  
  // normalize the phi between -pi and pi
  while(y[1] > PI){
    y[1] -= PI2;
  }
  while(y[1] < -PI){
    y[1] += PI2;
  }
  // use the same process map function as KF
  MatrixXd Ht = H_.transpose();
  MatrixXd S = H_*P_*Ht + R_;
  MatrixXd Si = S.inverse();
  MatrixXd K = P_*Ht*Si;
  
  // update
  x_ = x_ + (K*y);
  long x_size = x_.size();
  MatrixXd I = MatrixXd::Identity(x_size, x_size);
  P_ = (I - K*H_) * P_;
  
}


