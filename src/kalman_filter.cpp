#include "kalman_filter.h"
#include <iostream>

using Eigen::MatrixXd;
using Eigen::VectorXd;

using namespace std;

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
      DONE:
      * predict the state
    */
    x_ = F_ * x_;
    MatrixXd Ft = F_.transpose();
    P_ = F_ * P_ * Ft + Q_;
}

void KalmanFilter::Update(const VectorXd &z) {
    cout << "KalmanFilter::Update ..." << endl;
    /**
    DONE:
    * update the state by using Kalman Filter equations
    */
    VectorXd z_pred = H_ * x_;
    VectorXd y = z - z_pred;
    // there is a requirement to in (-pi < y < pi), aproximating
    //while (y(1) < -M_PI) {y(1) += 2 * M_PI;}
    //while (y(1) > M_PI) {y(1) -= 2 * M_PI;}

    MatrixXd Ht = H_.transpose();
    MatrixXd S = H_ * P_ * Ht + R_;
    MatrixXd Si = S.inverse();
    MatrixXd PHt = P_ * Ht;
    MatrixXd K = PHt * Si;

    //new estimate
    x_ = x_ + (K * y);
    long x_size = x_.size();
    MatrixXd I = MatrixXd::Identity(x_size, x_size);
    P_ = (I - K * H_) * P_;
}

void KalmanFilter::UpdateEKF(const VectorXd &z) {
    cout << "KalmanFilter::UpdateEKF ..." << endl;

    cout << "P_: " << endl;
    cout << P_ << endl;
    cout << "H_: " << endl;
    cout << H_ << endl;
    cout << "R_: " << endl;
    cout << R_ << endl;

    /**
    TODO:
    * update the state by using Extended Kalman Filter equations
    * for RADAR
    *
    * x' = [px, py, vx, vy]
    * range:      ρ = sqrt(p​x^2 + py^2) 
    * bearing:    φ = atan(py/px)
    * range_rate: ρ˙= (px * vx + py * vy) / ρ
    */
    double range = sqrt(x_[0] * x_[0] + x_[1] * x_[1]);
    double bearing;
    double range_rate;
    // be sure we don't divide by zero
    if (fabs(range) > 0.001) {
      bearing = atan2(x_[1] , x_[0]);
      range_rate = ((x_[0] * x_[2] + x_[1] * x_[3]) / range);
    } else{
      bearing = 0;
      range_rate = 0;
    }
    //MatrixXd z_pred(3, 1);
    VectorXd h = VectorXd(3);
    h << range, bearing, range_rate;

    cout << "calculating y ..." << endl;
    VectorXd y = z - h;
    /**
      * One important point when calculating y with radar sensor data:
      * the second value in the polar coordinate vector is the angle ϕ.
      * You'll need to make sure to normalize ϕ in the y vector so that its angle is between −pi and pi;
      * in other words, add or subtract 2pi from ϕ until it is between −pi and pi.
      */
    //while (y(1) < -M_PI) {y(1) += 2 * M_PI;}
    //while (y(1) > M_PI) {y(1) -= 2 * M_PI;}

    MatrixXd Ht = H_.transpose();
    MatrixXd S = H_ * P_ * Ht + R_;
    cout << "S: " << endl;
    cout << S << endl;
    MatrixXd Si = S.inverse();
    MatrixXd PHt = P_ * Ht;
    cout << "PHt: " << endl;
    cout << PHt << endl;
    cout << "Si: " << endl;
    cout << Si << endl;
    MatrixXd K = PHt * Si;

    // new estimate
    cout << "New Estimate ..." << endl;
    cout << "x_: " << endl;
    cout << x_ << endl;
    cout << "K: " << endl;
    cout << K << endl;
    cout << "y: " << endl;
    cout << y << endl;

    x_ = x_ + (K * y);
    long x_size = x_.size();
    MatrixXd I = MatrixXd::Identity(x_size, x_size);
    P_ = (I - K * H_) * P_;

}
