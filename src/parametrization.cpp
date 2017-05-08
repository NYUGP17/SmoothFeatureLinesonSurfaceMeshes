//
// Created by Liang Zhuo on 4/6/17.
//

#include <igl/triangle_triangle_adjacency.h>
#include <Eigen/Sparse>
#include <Eigen/SparseCholesky>

#include "parametrization.h"

namespace GM {
    double _dot(Cd a, Cd b) {
        return a.real() * b.real() + a.imag() * b.imag();
    }

    double _det(Cd a, Cd b) {
        return a.real() * b.imag() - a.imag() * b.real();
    }

    Eigen::MatrixXd solve_vector_field_problem(const Eigen::MatrixXd &V, const Eigen::MatrixXi &F, const Eigen::VectorXi &constrainedF, const Eigen::MatrixXd &constrainedV) {
        auto my_basis = GM::get_local_basis(V, F);
        auto u = solve_vector_field_problem(V, F, constrainedF, constrainedV, my_basis);
        Eigen::MatrixXd field = to_vector_field(u, my_basis);
        delete[] my_basis;
        return field;
    }

    std::vector <Cd> solve_vector_field_problem(const Eigen::MatrixXd &V, const Eigen::MatrixXi &F, const Eigen::VectorXi &constrainedF, const Eigen::MatrixXd &constrainedV, const Eigen::MatrixXd *my_basis) {
        Eigen::MatrixXi TT;
        igl::triangle_triangle_adjacency(F, TT);
        auto mapping = std::vector <int> (F.rows());
        auto u = std::vector <Cd> (F.rows());
        for (int i = 0; i < constrainedF.size(); ++ i) {
            mapping[constrainedF(i)] = -1;
            u[constrainedF(i)] = Cd(constrainedV.row(i).dot(my_basis[0].row(constrainedF(i))), constrainedV.row(i).dot(my_basis[1].row(constrainedF(i))));
        }
        for (int i = 0, total_faces = 0; i < F.rows(); ++ i)
            if (mapping[i] != -1) {
                mapping[i] = total_faces;
                ++ total_faces;
            }
        std::vector <Eigen::Triplet <Cd>> coefficients_a, coefficients_b;
        int nr = 0;
        for (int f1 = 0; f1 < F.rows(); ++ f1)
            for (int i = 0; i < F.cols(); ++ i) {
                int j = (i + 1 == F.cols()) ? 0 : i + 1;
                int vi1 = F(f1, i);
                int vi2 = F(f1, j);
                /// f2 is the opposite face of edge (vi1, vi2)
                int f2 = TT(f1, i);
                if (f2 < 0)
                    continue;
                /// skip if both faces are constrained
                if (mapping[f1] == -1 && mapping[f2] == -1)
                    continue;
                if (f1 < f2) {
                    auto common_edge = V.row(vi2) - V.row(vi1);
                    auto ef1 = Cd(common_edge.dot(my_basis[0].row(f1)), common_edge.dot(my_basis[1].row(f1)));
                    auto ef2 = Cd(common_edge.dot(my_basis[0].row(f2)), common_edge.dot(my_basis[1].row(f2)));
                    ef1 = std::conj(ef1 / abs(ef1));
                    ef2 = std::conj(ef2 / abs(ef2));
                    if (mapping[f1] != -1)
                        coefficients_a.push_back(Eigen::Triplet <Cd>(nr, mapping[f1], ef1));
                    else
                        coefficients_b.push_back(Eigen::Triplet <Cd>(nr, 0, - u[f1] * ef1));
                    if (mapping[f2] != -1)
                        coefficients_a.push_back(Eigen::Triplet <Cd>(nr, mapping[f2], - ef2));
                    else
                        coefficients_b.push_back(Eigen::Triplet <Cd>(nr, 0, u[f2] * ef2));
                    ++ nr;
                }
            }
        Eigen::SparseMatrix <Cd> A(nr, F.rows() - constrainedF.size()), b(nr, 1);
        A.setFromTriplets(coefficients_a.begin(), coefficients_a.end());
        b.setFromTriplets(coefficients_b.begin(), coefficients_b.end());
        Eigen::SimplicialLDLT <Eigen::SparseMatrix <Cd>> solver;
        solver.compute(A.adjoint() * A);
        assert(solver.info() == Eigen::Success);
        Eigen::MatrixXcd z = solver.solve(A.adjoint() * Eigen::MatrixXcd(b)); // convert to dense matrix
        assert(solver.info() == Eigen::Success);
        for (int i = 0; i < F.rows(); ++ i)
            if (mapping[i] != -1)
                u[i] = z(mapping[i]);
        return u;
    }

    Eigen::MatrixXd to_vector_field(std::vector <Cd> &u, const Eigen::MatrixXd *my_basis) {
        int n = u.size();
        Eigen::MatrixXd field(n, my_basis[0].cols());
        for (int i = 0; i < n; ++ i)
            field.row(i) = u[i].real() * my_basis[0].row(i) + u[i].imag() * my_basis[1].row(i);
        return field;
    }

    std::vector <Cd> get_gradient(const Eigen::MatrixXd &V, const Eigen::MatrixXi &F, const Eigen::VectorXd &f_values, const Eigen::MatrixXd *my_basis) {
        std::vector <Cd> result;
        for (int i = 0; i < F.rows(); ++ i) {
            std::vector <double> f;
            std::vector <Cd> c;
            for (int j = 0; j < F.cols(); ++ j) {
                c.push_back(Cd(V.row(F(i, j)).dot(my_basis[0].row(i)), V.row(F(i, j)).dot(my_basis[1].row(i))));
                f.push_back(f_values(F(i, j)));
            }
            /// the following formula comes from the end of slide 8
            auto g = ((f[1] - f[0]) * (c[0] - c[2]) + (f[2] - f[0]) * (c[1] - c[0])) * Cd(0, 1) / _det(c[1] - c[0], c[2] - c[0]);
            result.push_back(g);
        }
        return result;
    }

    Eigen::VectorXd get_signed_double_area(const Eigen::MatrixXd &V, const Eigen::MatrixXi &F) {
        Eigen::VectorXd result(F.rows());
        assert(V.cols() == 2);
        for (int i = 0; i < F.rows(); ++ i) {
            result(i) = 0;
            for (int j = 0; j < F.cols(); ++ j) {
                int k = (j + 1 == F.cols()) ? 0 : j + 1;
                int v1 = F(i, j);
                int v2 = F(i, k);
                result(i) += V(v1, 0) * V(v2, 1) - V(v1, 1) * V(v2, 0);
            }
        }
        return result;
    }
}