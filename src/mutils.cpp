//
// Created by Liang Zhuo on 5/9/17.
//

#include <cmath>

#include "mutils.h"

namespace GM {
    int sgn(double x) {
        if (fabs(x) < 1e-6)
            return 0;
        return (x < 0) ? -1 : 1;
    }

    double dot(Cd a, Cd b) {
        return a.real() * b.real() + a.imag() * b.imag();
    }

    double det(Cd a, Cd b) {
        return a.real() * b.imag() - a.imag() * b.real();
    }

    Eigen::MatrixXd toEigenMatrix(const std::vector <Eigen::VectorXd> &a) {
        int m = (a.size() == 0) ? 3 : a[0].size();
        Eigen::MatrixXd result(a.size(), m);
        for (int i = 0; i < (int) a.size(); ++ i)
            result.row(i) = a[i];
        return result;
    }
}