//
// Created by Liang Zhuo on 5/8/17.
//

#include <igl/cotmatrix.h>
#include <igl/adjacency_list.h>

#include "laplace_smoothing.h"
#include "mutils.h"

namespace GM {
    void compute_laplace_smoothing(const Eigen::MatrixXd &V, const Eigen::MatrixXi &F, const Eigen::MatrixXd &EV, const Eigen::VectorXd &E, Eigen::VectorXd &L) {
        Eigen::SparseMatrix <double> cot_matrix;
        std::vector <std::vector <int>> adj;
        igl::cotmatrix(V, F, cot_matrix);
        igl::adjacency_list(F, adj);
        L.setZero(V.rows());
        for (int i = 0; i < V.rows(); ++ i)
            for (int j: adj[i])
                L(i) += cot_matrix.coeff(i, j) * (sgn(EV.row(i).dot(EV.row(j))) * E(j) - E(i));
    }
}