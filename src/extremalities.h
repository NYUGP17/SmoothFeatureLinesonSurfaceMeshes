//
// Created by Liang Zhuo on 5/6/17.
//

#ifndef SMOOTHFEATURELINE_EXTREMALITIES_H
#define SMOOTHFEATURELINE_EXTREMALITIES_H


namespace GM {
    int sgn(double x);
    void compute_shape_operators(const Eigen::MatrixXd &V, const Eigen::MatrixXi &F, const Eigen::MatrixXd &VN, const Eigen::MatrixXd &FN, std::vector<Eigen::MatrixXd> &VS);
    void compute_eigens(const std::vector<Eigen::MatrixXd> &VS, std::vector <Eigen::VectorXd> &K, std::vector <Eigen::MatrixXd> &EV);
    void compute_area_star(const Eigen::MatrixXd &V, const Eigen::MatrixXi &F, Eigen::VectorXd &AS);
    void make_consistent(const Eigen::MatrixXi &F, Eigen::MatrixXd &EV, Eigen::VectorXi &is_regular);
    void compute_is_regular(const Eigen::MatrixXi &F, const Eigen::MatrixXd &EV, Eigen::VectorXi &is_regular);
    void compute_extremalities(const Eigen::MatrixXd &V, const Eigen::MatrixXi &F, const Eigen::MatrixXd &gradient, const Eigen::MatrixXd &EV, Eigen::VectorXd &E);
};


#endif //SMOOTHFEATURELINE_EXTREMALITIES_H
