//
// Created by Liang Zhuo on 5/9/17.
//

#ifndef SMOOTHFEATURELINE_MUTILS_H
#define SMOOTHFEATURELINE_MUTILS_H

#include <complex>
#include <Eigen/Dense>
#include <vector>

namespace GM {
    typedef std::complex <double> Cd;
    int sgn(double x);
    double dot(Cd a, Cd b);
    double det(Cd a, Cd b);
    Eigen::MatrixXd toEigenMatrix(const std::vector <Eigen::VectorXd> &a);
}


#endif //SMOOTHFEATURELINE_MUTILS_H
