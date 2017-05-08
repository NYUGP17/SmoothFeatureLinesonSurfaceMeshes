//
// Created by Liang Zhuo on 5/8/17.
//

#ifndef SMOOTHFEATURELINE_LAPLACE_SMOOTHING_H
#define SMOOTHFEATURELINE_LAPLACE_SMOOTHING_H


namespace GM {
    void compute_laplace_smoothing(const Eigen::MatrixXd &V, const Eigen::MatrixXi &F, const Eigen::MatrixXd &EV, const Eigen::VectorXd &E, Eigen::VectorXd &L);
}


#endif //SMOOTHFEATURELINE_LAPLACE_SMOOTHING_H
