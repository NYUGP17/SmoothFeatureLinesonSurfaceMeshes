//
// Created by Liang Zhuo on 5/6/17.
//
#include <Eigen/Eigenvalues>
#include <igl/doublearea.h>
#include <igl/adjacency_list.h>
#include <igl/vertex_triangle_adjacency.h>
#include <igl/triangle_triangle_adjacency.h>
#include <set>

#include "extremalities.h"
#include "UFSet.h"

void _update_UFSet(const Eigen::MatrixXi &F, const Eigen::VectorXi &is_regular, const std::vector<std::vector<int>> &VF, UFSet &ufSet);

void _update_UFSet(const Eigen::MatrixXi &F, const Eigen::VectorXi &is_regular, const std::vector<std::vector<int>> &VF, UFSet &ufSet) {
    int n = VF.size();
    for (int i = 0; i < n; ++ i)
        for (int f: VF[i])
            if (is_regular(f))
                for (int j = 0; j < F.cols(); ++ j)
                    ufSet.merge(i, F(f, j));
}

int GM::sgn(double x) {
    if (fabs(x) < 1e-6)
        return 0;
    return (x < 0) ? -1 : 1;
}

void GM::compute_shape_operators(const Eigen::MatrixXd &V, const Eigen::MatrixXi &F, const Eigen::MatrixXd &VN, const Eigen::MatrixXd &FN, std::vector<Eigen::MatrixXd> &VS) {
    Eigen::MatrixXi TT;
    igl::triangle_triangle_adjacency(F, TT);
    /// initialization
    VS.clear();
    for (int i = 0; i < V.rows(); ++ i)
        VS.push_back(Eigen::MatrixXd::Zero(3, 3));
    for (int f1 = 0; f1 < F.rows(); ++ f1)
        for (int i = 0; i < F.cols(); ++ i) {
            int j = (i + 1 == F.cols()) ? 0 : i + 1;
            int vi1 = F(f1, i);
            int vi2 = F(f1, j);
            /// f2 is the opposite face of edge (vi1, vi2)
            int f2 = TT(f1, i);
            if (f2 < 0)
                continue;
            /// do once for each edge
            if (f1 < f2) {
                auto common_edge = V.row(vi2) - V.row(vi1);
                Eigen::Vector3d e_unit = common_edge.normalized();
                /// mean curvature = 2 |e| cos(x / 2)
                /// cos(x/2) = (+/-) sqrt((1 + cos(x)) / 2)
                /// first compute cos(x):
                auto fn1 = FN.row(f1);
                auto fn2 = FN.row(f2);
                auto cos_x = - fn1.dot(fn2) / (fn1.norm() * fn2.norm());
                auto cos_half = sqrt(std::max((1. + cos_x) / 2, 0.));
                if (cos_x < 0)
                    cos_half = - cos_half;
                auto mean_curvature = 2. * common_edge.norm() * cos_half;
                /// edge normal = (N1 + N2) / || N1 + N2 ||
                Eigen::Vector3d edge_normal = (VN.row(vi1) + VN.row(vi2)).normalized();
                auto eNe = e_unit.cross(edge_normal);
                auto edge_shape = mean_curvature * eNe * eNe.transpose();
                /// update vertex shape
                VS[vi1] += 0.5 * V.row(vi1).dot(edge_normal) * edge_shape;
                VS[vi2] += 0.5 * V.row(vi2).dot(edge_normal) * edge_shape;
            }
        }
}

void GM::compute_eigens(const std::vector<Eigen::MatrixXd> &VS, std::vector <Eigen::VectorXd> &K, std::vector <Eigen::MatrixXd> &EV) {
    /// initialization
    K.resize(2);
    EV.resize(2);
    for (int k = 0; k < 2; ++ k) {
        K[k].resize(VS.size());
        EV[k].resize(VS.size(), 3);
    }
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> solver;
    for (int i = 0; i < (int) VS.size(); ++ i) {
        solver.compute(VS[i]);
        auto eigenValues =  solver.eigenvalues();
        auto eigenVectors = solver.eigenvectors();
        int orders[3] = {0, 1, 2};
        std::sort(orders, orders + 3, [&eigenValues](int i, int j) {
            return fabs(eigenValues(i)) > fabs(eigenValues(j));
        });
        auto k_max = orders[0];
        auto k_min = orders[1];
        if (k_max < k_min)
            std::swap(k_max, k_min);
        K[0](i) = eigenValues(k_max);
        EV[0].row(i) = eigenVectors.col(k_max).normalized();
        K[1](i) = eigenValues(k_min);
        EV[1].row(i) = eigenVectors.col(k_min).normalized();
    }
}

void GM::compute_area_star(const Eigen::MatrixXd &V, const Eigen::MatrixXi &F, Eigen::VectorXd &AS) {
    Eigen::VectorXd A;
    igl::doublearea(V, F, A);
    AS.setZero(V.rows());
    for (int f = 0; f < F.rows(); ++ f)
        for (int j = 0; j < F.cols(); ++ j) {
            int i = F(f, j);
            AS(i) += 0.5 * fabs(A(f));
        }
}

void GM::make_consistent(const Eigen::MatrixXi &F, Eigen::MatrixXd &EV, Eigen::VectorXi &is_regular) {
    int n = EV.rows();
    std::vector<std::vector<int>> VF, VFi, A;
    igl::vertex_triangle_adjacency(n, F, VF, VFi);
    igl::adjacency_list(F, A);
    /// make consistency on edges
    for (int remain = n, next_time = 0; remain > 0; remain = next_time) {
        next_time = 0;
        for (int i = 0; i < n; ++ i) {
            int nImprovement = 0;
            for (int j: A[i]) {
                if (sgn(EV.row(i).dot(EV.row(j))) < 0)
                    ++ nImprovement;
                else if (sgn(EV.row(i).dot(EV.row(j))) > 0)
                    -- nImprovement;
            }
            if (nImprovement > 0) {
                next_time = remain - 1;
                EV.row(i) *= -1;
            }
        }
    }
    /// point-wise improvement
    for (bool admissible = true; admissible; ) {
        admissible = false;
        for (int i = 0; i < n; ++ i) {
            /// see if we can flip the sign?
            int nImprovement = 0;
            for (int j = 0; j < (int) VF[i].size(); ++ j) {
                int f = VF[i][j];
                int i1 = VFi[i][j];
                int i2 = (i1 + 1 == F.cols()) ? 0 : i1 + 1;
                int i3 = (i2 + 1 == F.cols()) ? 0 : i2 + 1;
                assert(F(f, i1) == i);
                if (GM::sgn(EV.row(F(f, i2)).dot(EV.row(F(f, i3)))) > 0) {
                    if (GM::sgn(EV.row(F(f, i1)).dot(EV.row(F(f, i2)))) < 0 && GM::sgn(EV.row(F(f, i1)).dot(EV.row(F(f, i3)))) < 0)
                        ++ nImprovement;
                    else if (sgn(EV.row(F(f, i1)).dot(EV.row(F(f, i2)))) > 0 && GM::sgn(EV.row(F(f, i1)).dot(EV.row(F(f, i3)))) > 0)
                        -- nImprovement;
                }
            }
            if (nImprovement > 0) {
                admissible = true;
                EV.row(i) *= -1;
            }
        }
    }
    /// connected component improvement
    UFSet ufSet(n);
    for (bool admissible = true; admissible; ) {
        admissible = false;
        GM::compute_is_regular(F, EV, is_regular);
        _update_UFSet(F, is_regular, VF, ufSet);
        /// compute connected components
        for (int f = 0; f < F.rows(); ++ f)
            if (!is_regular(f)) {
                std::set <int> owners;
                for (int j = 0; j < F.cols(); ++ j)
                    owners.insert(ufSet.find(F(f, j)));
                if (owners.size() == 1)
                    continue;
                if (owners.size() == 2) {
                    /// see if we could make it regular
                    bool admissible = true;
                    for (int fi1 = 0; fi1 < F.cols(); ++ fi1) {
                        int fi2 = (fi1 + 1 == F.cols()) ? 0 : fi1 + 1;
                        int i1 = F(f, fi1);
                        int i2 = F(f, fi2);
                        if (ufSet.find(i1) == ufSet.find(i2))
                            continue;
                        if (GM::sgn(EV.row(i1).dot(EV.row(i2))) > 0)
                            admissible = false;
                    }
                    if (!admissible)
                        continue;
                    /// great, flip one of the signs
                    for (int i = 0; i < n; ++ i)
                        if (ufSet.find(i) == *owners.begin())
                            EV.row(i) *= -1;
                    break;
                }
            }
    }
}

void GM::compute_is_regular(const Eigen::MatrixXi &F, const Eigen::MatrixXd &EV, Eigen::VectorXi &is_regular) {
    is_regular.setZero(F.rows());
    for (int f = 0; f < F.rows(); ++ f) {
        is_regular(f) = true;
        for (int fi1 = 0; fi1 < F.cols(); ++ fi1) {
            int fi2 = (fi1 + 1 == F.cols()) ? 0 : fi1 + 1;
            int i1 = F(f, fi1);
            int i2 = F(f, fi2);
            if (GM::sgn(EV.row(i1).dot(EV.row(i2))) < 0) {
                is_regular(f) = false;
                break;
            }
        }
    }
}

void GM::compute_extremalities(const Eigen::MatrixXd &V, const Eigen::MatrixXi &F, const Eigen::MatrixXd &gradient, const Eigen::MatrixXd &EV, Eigen::VectorXd &E) {
    Eigen::VectorXd A;
    Eigen::VectorXd AS;
    igl::doublearea(V, F, A);
    GM::compute_area_star(V, F, AS);
    E.setZero(V.rows());
    for (int f = 0; f < F.rows(); ++ f)
        for (int j = 0; j < F.cols(); ++ j) {
            int i = F(f, j);
            E(i) += 0.5 * fabs(A(f)) * gradient.row(f).dot(EV.row(i));
        }
    for (int i = 0; i < V.rows(); ++ i)
        E(i) /= AS(i);
}