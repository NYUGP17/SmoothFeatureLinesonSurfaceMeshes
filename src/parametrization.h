//
// Created by Liang Zhuo on 4/6/17.
//

#ifndef ASSIGNMENT4_PARAMETRIZATION_H
#define ASSIGNMENT4_PARAMETRIZATION_H

#include <igl/local_basis.h>

#include "mutils.h"

namespace GM {
    template <class MatrixType1, class MatrixType2> MatrixType1 *get_local_basis(const MatrixType1 &, const MatrixType2 &);
    Eigen::MatrixXd solve_vector_field_problem(const Eigen::MatrixXd &, const Eigen::MatrixXi &, const Eigen::VectorXi &, const Eigen::MatrixXd &);
    std::vector <Cd> solve_vector_field_problem(const Eigen::MatrixXd &, const Eigen::MatrixXi &, const Eigen::VectorXi &, const Eigen::MatrixXd &, const Eigen::MatrixXd *);
    std::vector <Cd> get_gradient(const Eigen::MatrixXd &, const Eigen::MatrixXi &, const Eigen::VectorXd &, const Eigen::MatrixXd *);
    template <class VectorType1, class VectorType2> Cd get_gradient(const VectorType1 &, const VectorType2 &);
    Eigen::MatrixXd to_vector_field(std::vector <Cd> &, const Eigen::MatrixXd *);
    Eigen::VectorXd get_signed_double_area(const Eigen::MatrixXd &, const Eigen::MatrixXi &);

    template <class MatrixType1, class MatrixType2> MatrixType1 *get_local_basis(const MatrixType1 &V, const MatrixType2 &F) {
        auto b = new MatrixType1[3];
        igl::local_basis(V, F, b[0], b[1], b[2]);
        return b;
    }

    template <class VectorType1, class VectorType2> Cd get_gradient(const VectorType1 &c, const VectorType2 &f) {
        return ((f[1] - f[0]) * (c[0] - c[2]) + (f[2] - f[0]) * (c[1] - c[0])) * Cd(0, 1) / det(c[1] - c[0], c[2] - c[0]);
    }
}

#endif //ASSIGNMENT4_PARAMETRIZATION_H
