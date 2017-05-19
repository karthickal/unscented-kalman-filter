#include <iostream>
#include "tools.h"

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;

Tools::Tools() {}

Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) {

    // initialize the RMSE
    VectorXd rmse(4);
    rmse << 0.0, 0.0, 0.0, 0.0;

    // validate the estimation and ground truth values
    if (estimations.size() == 0 || estimations.size() != ground_truth.size()) {
        std::cout << "The estimation values are incorrect" << std::endl;
        return rmse;
    }

    // add the squared residual to the rmse
    for (int i=0; i<estimations.size(); i++ ) {
        VectorXd residual = estimations[i] - ground_truth[i];
        residual = residual.array() * residual.array();
        rmse = rmse + residual;
    }

    rmse = rmse / estimations.size();
    rmse = rmse.array().sqrt();

    // return the final result
    return rmse;


}
