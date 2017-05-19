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

    // if this is false, laser measurements will be ignored (except during init)
    use_laser_ = true;

    // if this is false, radar measurements will be ignored (except during init)
    use_radar_ = true;

    // initial state vector
    x_ = VectorXd(5);

    // initial covariance matrix
    P_ = MatrixXd(5, 5);

    // Process noise standard deviation longitudinal acceleration in m/s^2
    std_a_ = 3;

    // Process noise standard deviation yaw acceleration in rad/s^2
    std_yawdd_ = 0.9;

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

    // set initialization to false
    is_initialized_ = false;

    // set the length of state dimensions
    n_x_ = x_.rows();

    // set the length of the augmented state dimensions
    n_aug_ = 7;

    // initialize lambda
    lambda_ = 3 - n_aug_;

    // initialize the sigma points matrix
    Xsig_pred_ = MatrixXd(n_x_, 2 * n_aug_ + 1);
    Xsig_pred_.fill(0.0);

    // initialize the state covariance matrix
    P_.diagonal() << 1, 1, 10, 1, 1;

    // initialize the process noise matrix
    Q_ = MatrixXd(2, 2);
    Q_ << std_a_ * std_a_, 0,
            0, std_yawdd_ * std_yawdd_;

    // initialize the weights
    double w_others = 0.5 / (lambda_ + n_aug_);
    weights_ = VectorXd(2 * n_aug_ + 1);
    weights_(0) = lambda_ / (lambda_ + n_aug_);
    for (int j = 1; j < weights_.rows(); j++) {
        weights_(j) = w_others;
    }

    // initialize the NIS values
    NIS_radar_ = 0.0;
    NIS_laser_ = 0.0;

}

UKF::~UKF() {}

/**
 * @param {MeasurementPackage} meas_package The latest measurement data of
 * either radar or laser.
 */
void UKF::ProcessMeasurement(MeasurementPackage meas_package) {

    // initialize the state if this is the first measurement
    if (!is_initialized_) {

        x_ << 1, 1, 0, 0, 0;

        if (meas_package.sensor_type_ == MeasurementPackage::LASER) {
            x_ << meas_package.raw_measurements_[0], meas_package.raw_measurements_[1], 0, 0, 0;
        } // if measurement is from radar convert to state co-ordinates
        else if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {
            float rho = meas_package.raw_measurements_[0];
            float phi = meas_package.raw_measurements_[1];
            x_ << rho * cos(phi), rho * sin(phi), 0, 0, 0;
        }

        is_initialized_ = true;

        // save the current package's timestamp
        time_us_ = meas_package.timestamp_;

        return;
    }

    // skip if the sensor data should not be included
    if ((meas_package.sensor_type_ == MeasurementPackage::LASER && !use_laser_) ||
        (meas_package.sensor_type_ == MeasurementPackage::RADAR && !use_radar_)) {
        return;
    }

    // calculate the time difference and predict the new state
    double delta_t = (meas_package.timestamp_ - time_us_) / 1000000.0;
    time_us_ = meas_package.timestamp_;

    // predict the new state using the process model
    Prediction(delta_t);

    // update the state and covariance based on the type of sensor
    if (meas_package.sensor_type_ == MeasurementPackage::LASER) {
        UpdateLidar(meas_package);
    } else if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {
        UpdateRadar(meas_package);
    }

}

/**
 * Predicts sigma points, the state, and the state covariance matrix.
 * @param {double} delta_t the change in time (in seconds) between the last
 * measurement and this one.
 */
void UKF::Prediction(double delta_t) {

    // create the augmented state matrix
    VectorXd x_aug(n_aug_);
    x_aug.head(n_x_) = x_;
    x_aug(n_x_) = 0;
    x_aug(n_x_ + 1) = 0;

    // create the augmented covariance matrix
    MatrixXd P_aug(n_aug_, n_aug_);
    P_aug.fill(0.0);
    P_aug.topLeftCorner(n_x_, n_x_) = P_;
    P_aug.bottomRightCorner(n_aug_ - n_x_, n_aug_ - n_x_) = Q_;

    // get the sigma points
    MatrixXd Xsig_mat = UKF::GetSigmaPoints(x_aug, P_aug, lambda_);

    //predict sigma points
    long long delta_t_a = 0.5 * delta_t * delta_t;
    for (int i = 0; i < Xsig_mat.cols(); i++) {

        // extract the state and all parameters
        VectorXd aug_state = Xsig_mat.col(i);
        double p_x = aug_state(0);
        double p_y = aug_state(1);
        double v = aug_state(2);
        double yaw = aug_state(3);
        double yaw_rate = aug_state(4);
        double long_accel = aug_state(5);
        double yaw_accel = aug_state(6);

        Xsig_pred_.col(i) = aug_state.head(5);

        // save reusable calculations
        float sin_si = sin(yaw);
        float cos_si = cos(yaw);
        float angle_a = yaw + (yaw_rate * delta_t);
        VectorXd Xsig_state(5);
        VectorXd Xsig_noise(5);

        // calculate the predicted state
        if (fabs(yaw_rate) < 0.0001) {
            Xsig_state << v * cos_si * delta_t,
                    v * sin_si * delta_t,
                    0,
                    0,
                    0;
        } else {
            Xsig_state << (v / yaw_rate) * (sin(angle_a) - sin_si),
                    (v / yaw_rate) * (-cos(angle_a) + cos_si),
                    0,
                    yaw_rate * delta_t,
                    0;
        }

        // calculate the noise
        Xsig_noise << delta_t_a * cos_si * long_accel,
                delta_t_a * sin_si * long_accel,
                delta_t * long_accel,
                delta_t_a * yaw_accel,
                delta_t * yaw_accel;


        // calculate the predicted state points
        Xsig_pred_.col(i) = Xsig_pred_.col(i) + Xsig_state + Xsig_noise;
        Xsig_pred_(3, i) = NormalizeAngle(Xsig_pred_(3, i));
    }

    //predict state mean
    x_ = Xsig_pred_ * weights_;

    //predict state covariance matrix
    P_.fill(0.0);
    for (int i = 0; i < Xsig_pred_.cols(); i++) {

        VectorXd factor = Xsig_pred_.col(i) - x_;
        factor(3) = NormalizeAngle(factor(3));

        MatrixXd factor_p = weights_(i) * factor * factor.transpose();
        P_ = P_ + factor_p;
    }
}

/**
 * Updates the state and the state covariance matrix using a laser measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateLidar(MeasurementPackage meas_package) {

    // number of measurements in laser
    int n_z = 2;

    // initialize mean predicted measurement
    VectorXd z_pred = VectorXd(n_z);

    // initialize measurement covariance matrix S
    MatrixXd S = MatrixXd(n_z, n_z);

    // create the measurement sigma points matrix
    MatrixXd Zsig = MatrixXd(n_z, 2 * n_aug_ + 1);
    for (int i = 0; i < Xsig_pred_.cols(); i++) {
        Zsig.col(i) = Xsig_pred_.col(i).head(2);
    }

    //calculate mean predicted measurement
    z_pred = Zsig * weights_;

    //calculate measurement covariance matrix S
    MatrixXd R = MatrixXd(n_z, n_z);
    S.fill(0.0);

    //calculate sensor noise matrix R
    R.fill(0.0);
    R << std_laspx_*std_laspx_, 0,
            0, std_laspy_*std_laspy_;

    // calculate the innovation covariance
    for (int j = 0; j < 2 * n_aug_ + 1; j++) {
        VectorXd residual = Zsig.col(j) - z_pred;
        S = S + (weights_(j) * residual * residual.transpose());
    }
    S = S + R;

    //create matrix for cross correlation Tc
    MatrixXd Tc = MatrixXd(n_x_, n_z);
    Tc.fill(0.0);

    //calculate cross correlation matrix
    for (int i = 0; i < Xsig_pred_.cols(); i++) {
        VectorXd state_sig = Xsig_pred_.col(i) - x_;
        VectorXd meas_sig = Zsig.col(i) - z_pred;
        MatrixXd factor = weights_(i) * state_sig * meas_sig.transpose();
        Tc = Tc + factor;
    }

    // calculate Kalman gain K;
    MatrixXd K = Tc * S.inverse();

    // calculate residual
    VectorXd z_diff = meas_package.raw_measurements_ - z_pred;

    //update state mean and covariance matrix
    x_ = x_ + K * z_diff;
    P_ = P_ - (K * S * K.transpose());

    NIS_laser_ = z_diff.transpose() * S.inverse() * z_diff;
}

/**
 * Updates the state and the state covariance matrix using a radar measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateRadar(MeasurementPackage meas_package) {

    // number of measurements in laser
    int n_z = 3;

    //mean predicted measurement
    VectorXd z_pred = VectorXd(n_z);

    //measurement covariance matrix S
    MatrixXd S = MatrixXd(n_z, n_z);

    // create the measurement sigma points matrix
    MatrixXd Zsig = MatrixXd(n_z, 2 * n_aug_ + 1);
    for (int i = 0; i < Xsig_pred_.cols(); i++) {

        // extract the parameters
        VectorXd state = Xsig_pred_.col(i);
        double px = state(0);
        double py = state(1);
        double v = state(2);
        double yaw = state(3);
        double yaw_rate = state(4);

        Zsig.col(i) << sqrt((px * px) + (py * py)),
                atan2(py, px),
                ((px * cos(yaw) * v) + (py * sin(yaw) * v)) / sqrt((px * px) + (py * py));
    }

    //calculate mean predicted measurement
    z_pred = Zsig * weights_;

    //calculate measurement covariance matrix S
    MatrixXd R = MatrixXd(n_z, n_z);
    S.fill(0.0);

    //calculate sensor noise matrix R
    R.fill(0.0);
    R(0, 0) = std_radr_ * std_radr_;
    R(1, 1) = std_radphi_ * std_radphi_;
    R(2, 2) = std_radrd_ * std_radrd_;

    // calculate the innovation covariance matrix
    for (int j = 0; j < 2 * n_aug_ + 1; j++) {
        VectorXd meas = Zsig.col(j);
        VectorXd diff = meas - z_pred;

        //angle normalization
        while (diff(1) > M_PI) diff(1) -= 2. * M_PI;
        while (diff(1) < -M_PI) diff(1) += 2. * M_PI;

        MatrixXd state_factor = weights_(j) * diff * diff.transpose();
        S = S + state_factor;
    }
    S = S + R;

    //create matrix for cross correlation Tc
    MatrixXd Tc = MatrixXd(n_x_, n_z);
    Tc.fill(0.0);

    //calculate cross correlation matrix
    for (int i = 0; i < Xsig_pred_.cols(); i++) {
        VectorXd state_sig = Xsig_pred_.col(i) - x_;
        VectorXd meas_sig = Zsig.col(i) - z_pred;
        MatrixXd factor = weights_(i) * state_sig * meas_sig.transpose();
        Tc = Tc + factor;
    }

    //calculate Kalman gain K;
    MatrixXd K = Tc * S.inverse();

    // calculate residual
    VectorXd z_diff = meas_package.raw_measurements_ - z_pred;

    //update state mean and covariance matrix
    x_ = x_ + K * z_diff;
    P_ = P_ - (K * S * K.transpose());

    NIS_radar_ = z_diff.transpose() * S.inverse() * z_diff;

}

/**
 * Method to normalize the angle between -2pi and 2pi
 *
 * @param {VectorXd} x the state vector
 * @param {MatrixXd} P the covariance matrix
 * @param {double} lambda the lambda factor
 */
double UKF::NormalizeAngle(double value) {

    if ((value < -M_PI) || (value > M_PI)) {
        value -= (2.0 * M_PI) * floor(value / (2.0 * M_PI));
    }

    return value;
}

/**
 * Method to get the sigma points given the state, covariance and lambda factor
 *
 * @param {VectorXd} x the state vector
 * @param {MatrixXd} P the covariance matrix
 * @param {double} lambda the lambda factor
 */
MatrixXd UKF::GetSigmaPoints(VectorXd x, MatrixXd P, double lambda) {

    int n = x.rows();

    //create square root matrix
    MatrixXd A = P.llt().matrixL();
    MatrixXd sigma_a = sqrt(lambda + n) * A;

    // calculate the sigma points
    MatrixXd Xsig_mat(n, 2 * n + 1);
    Xsig_mat.col(0) = x;
    for (int i = 0; i < P.cols(); i++) {
        Xsig_mat.col(i + 1) = x + sigma_a.col(i);
        Xsig_mat.col(i + P.cols() + 1) = x - sigma_a.col(i);
    }

    return Xsig_mat;
}
