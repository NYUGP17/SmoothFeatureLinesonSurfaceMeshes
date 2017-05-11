#include <igl/readOFF.h>
#include <igl/viewer/Viewer.h>
#include <igl/barycenter.h>
#include <igl/vertex_triangle_adjacency.h>
#include <igl/triangle_triangle_adjacency.h>
#include <igl/per_face_normals.h>
#include <igl/per_vertex_normals.h>
#include <igl/slice.h>
#include <igl/avg_edge_length.h>
#include <igl/file_dialog_open.h>
#include <igl/polyvector_field_poisson_reconstruction.h>
#include <igl/jet.h>
#include <igl/boundary_loop.h>
#include <igl/harmonic.h>
#include <igl/lscm.h>
#include <igl/map_vertices_to_circle.h>
#include <igl/serialize.h>
/*** insert any libigl headers here ***/
#include "mutils.h"
#include "extremalities.h"
#include "parametrization.h"
#include "laplace_smoothing.h"
#include <memory>
#include <iterator>

using namespace std;
using Viewer = igl::viewer::Viewer;

// Vertex array, #V x3
Eigen::MatrixXd V(0,3);
// Face array, #F x3
Eigen::MatrixXi F(0,3);
// Face barycenter array, #F x3
Eigen::MatrixXd MF(0,3);
// Face normal array, #F x3
Eigen::MatrixXd FN(0,3);
// Vertex normal array, #V x3
Eigen::MatrixXd VN(0,3);
// Vertex-to-face adjacency
std::vector<std::vector<int> > VF, VFi;
// Vertex based shape operators
std::vector<Eigen::MatrixXd> VS;
// curvature and eigenvectors
std::vector <Eigen::VectorXd> K;
std::vector <Eigen::MatrixXd> EV;
std::vector <Eigen::VectorXi> is_regular;
std::vector <Eigen::VectorXd> extremalities;
// feature line adj
std::vector <std::vector <std::vector <int>>> zero_adj;
std::vector <Eigen::MatrixXd> ZV;
std::vector <std::vector <double>> ZW;
// ...
Eigen::MatrixXd *my_basis;
// Scale for displaying vectors
double vScale = 0;
// Threshold
double threshold = 1.;

// Function declarations (see below for implementation)
bool callback_key_down  (Viewer &viewer, unsigned char key, int modifiers);

bool callback_key_down(Viewer& viewer, unsigned char key, int modifiers) {
    if (key == '1') {
        // Draw object only
        viewer.data.clear();
        viewer.data.set_mesh(V, F);

    }

    if (key == '2') {
        // Regular triangles
        viewer.data.clear();
        viewer.data.set_mesh(V, F);
        viewer.data.add_edges(V, V + vScale * EV[0], Eigen::RowVector3d(0.9, 0.1, 0.1));
        viewer.data.add_edges(V, V + vScale * EV[1], Eigen::RowVector3d(0.1, 0.9, 0.1));
        Eigen::MatrixXd face_colors = Eigen::MatrixXd::Constant(F.rows(), 3, 0.9);
        for (int f = 0; f < F.rows(); ++ f)
            if (!is_regular[0](f))
                face_colors.row(f) << 231. / 255, 99. / 255, 113. / 255.;
        viewer.data.set_colors(face_colors);
    }

    if (key == '3') {
        // Scalar field reconstruction
        viewer.data.clear();
        viewer.data.set_mesh(V, F);
        Eigen::MatrixXd colors;
        igl::jet(extremalities[0], true, colors);
        viewer.data.set_colors(colors);
    }

    if (key == '4') {
        // Laplacian Smoothing
        viewer.data.clear();
        viewer.data.set_mesh(V, F);

        std::vector <Eigen::VectorXd> smoothed;
        smoothed.resize(K.size());
        for (int k = 0; k < (int) K.size(); ++ k) {
            GM::compute_laplace_smoothing(V, F, EV[k], extremalities[k], smoothed[k]);
            extremalities[k] += 0.01 * smoothed[k];
        }

        Eigen::MatrixXd colors;
        igl::jet(extremalities[0], true, colors);
        viewer.data.set_colors(colors);
    }

    if (key == '5') {
        // Extract Feature Lines
        viewer.data.clear();
        viewer.data.set_mesh(V, F);

        zero_adj.clear();
        zero_adj.resize(K.size());
        ZV.resize(K.size());
        ZW.clear();
        ZW.resize(K.size());
        for (int k = 0; k < (int) K.size(); ++ k) {
            int n_points = 0;
            std::map <pair <int, int>, int> zero_id;
            std::vector <Eigen::VectorXd> z;
            for (int f = 0; f < F.rows(); ++ f)
                if (is_regular[k](f)) {
                    double sum_k_max = 0, sum_k_min = 0;
                    Eigen::Vector3d e_slice;
                    Eigen::Matrix3d ev_slice;
                    std::vector <GM::Cd> c;
                    for (int j = 0; j < F.cols(); ++ j) {
                        int i = F(f, j);
                        sum_k_max += K[0](i);
                        sum_k_min += K[1](i);
                        e_slice(j) = extremalities[k](i);
                        ev_slice.row(j) = EV[k].row(i);
                        c.push_back(GM::Cd(V.row(i).dot(my_basis[0].row(f)), V.row(i).dot(my_basis[1].row(f))));
                    }
                    /// make consistent
                    for (int j = 1; j < F.cols(); ++ j)
                        if (GM::sgn(ev_slice.row(0).dot(ev_slice.row(j))) < 0) {
                            e_slice(j) *= -1;
                            ev_slice.row(j) *= -1;
                        }
                    /// compute e_gradient
                    auto e_gradient_cd = GM::get_gradient(c, e_slice);
                    auto e_gradient = e_gradient_cd.real() * my_basis[0].row(f) + e_gradient_cd.imag() * my_basis[1].row(f);
                    /// check equation (6)
                    if (k == 0 && GM::sgn(fabs(sum_k_max) - fabs(sum_k_min)) <= 0)
                        continue;
                    if (k == 1 && GM::sgn(fabs(sum_k_max) - fabs(sum_k_min)) >= 0)
                        continue;
                    /// check equation (5)
                    auto desired_sign = (k == 0) ? -1 : 1;
                    double sum_inner_product = 0;
                    for (int j = 1; j < F.cols(); ++ j)
                        sum_inner_product += e_gradient.dot(ev_slice.row(j));
                    if (GM::sgn(sum_inner_product) != desired_sign)
                        continue;
                    /// pass two requirements, see if contains zero sets
                    auto last_found = -1;
                    for (int fi1 = 0; fi1 < F.cols(); ++ fi1) {
                        int fi2 = (fi1 + 1 == F.cols()) ? 0 : fi1 + 1;
                        if (GM::sgn(e_slice(fi1)) * GM::sgn(e_slice(fi2)) < 0) {
                            /// found: create a vertex
                            int point_id;
                            if (zero_id.count(make_pair(F(f, fi1), F(f, fi2))))
                                point_id = zero_id[make_pair(F(f, fi1), F(f, fi2))];
                            else {
                                /// new point
                                point_id = n_points;
                                zero_id[make_pair(F(f, fi1), F(f, fi2))] = zero_id[make_pair(F(f, fi2), F(f, fi1))] = point_id;
                                ++ n_points;
                                z.push_back((V.row(F(f, fi1)) * fabs(e_slice(fi2)) + V.row(F(f, fi2)) * fabs(e_slice(fi1))) / (fabs(e_slice(fi2)) + fabs(e_slice(fi1))));
                                ZW[k].push_back((K[k](F(f, fi1)) * fabs(e_slice(fi2)) + K[k](F(f, fi2)) * fabs(e_slice(fi1))) / (fabs(e_slice(fi2)) + fabs(e_slice(fi1))));
                                zero_adj[k].push_back(std::vector <int> ());
                            }
                            if (last_found != -1) {
                                /// connect edge
                                zero_adj[k][point_id].push_back(last_found);
                                zero_adj[k][last_found].push_back(point_id);
                            }
                            last_found = point_id;
                        }
                    }
                }
            for (int f = 0; f < F.rows(); ++ f)
                if (!is_regular[k](f)) {
                    auto total_zeros = 0;
                    for (int fi1 = 0; fi1 < F.cols(); ++ fi1) {
                        int fi2 = (fi1 + 1 == F.cols()) ? 0 : fi1 + 1;
                        if (zero_id.count(make_pair(F(f, fi1), F(f, fi2))))
                            ++ total_zeros;
                    }
                    if (total_zeros == 3) {
                        /// new point
                        int point_id = n_points;
                        ++ n_points;
                        z.push_back(MF.row(f));
                        zero_adj[k].push_back(std::vector <int> ());
                        double w_sum = 0;
                        for (int fi1 = 0; fi1 < F.cols(); ++ fi1) {
                            int fi2 = (fi1 + 1 == F.cols()) ? 0 : fi1 + 1;
                            int last_found = zero_id[make_pair(F(f, fi1), F(f, fi2))];
                            /// connect barycenter
                            zero_adj[k][point_id].push_back(last_found);
                            zero_adj[k][last_found].push_back(point_id);
                            w_sum += K[k](fi1);
                        }
                        ZW[k].push_back(w_sum / F.cols());
                    } else if (total_zeros == 2) {
                        int last_found = -1;
                        for (int fi1 = 0; fi1 < F.cols(); ++ fi1) {
                            int fi2 = (fi1 + 1 == F.cols()) ? 0 : fi1 + 1;
                            if (zero_id.count(make_pair(F(f, fi1), F(f, fi2)))) {
                                int point_id = zero_id[make_pair(F(f, fi1), F(f, fi2))];
                                if (last_found != - 1) {
                                    /// connect edge
                                    zero_adj[k][point_id].push_back(last_found);
                                    zero_adj[k][last_found].push_back(point_id);
                                }
                                last_found = point_id;
                            }
                        }
                    }
                }
            cerr << "#points = " << z.size() << endl;
            ZV[k] = GM::toEigenMatrix(z);
        }
        /// display
        for (int k = 0; k < (int) K.size(); ++ k) {
            int total_degree = 0;
            for (int i = 0; i < ZV[k].rows(); ++ i)
                total_degree += zero_adj[k][i].size();
            cerr << "#edges = " << total_degree << endl;
            Eigen::MatrixXd endpoints1(total_degree / 2, ZV[k].cols()), endpoints2(total_degree / 2, ZV[k].cols());
            for (int i = 0, count = 0; i < ZV[k].rows(); ++ i)
                for (int j: zero_adj[k][i])
                    if (i < j) {
                        endpoints1.row(count) = ZV[k].row(i);
                        endpoints2.row(count) = ZV[k].row(j);
                        ++ count;
                    }
            auto color = Eigen::RowVector3d(0.1, 0.1, 0.1);
            color(k) = 0.9;
            viewer.data.add_edges(endpoints1, endpoints2, color);
        }
    }

    if (key == '6') {
        // Remove small ridges by a threshold filter
        viewer.data.clear();
        viewer.data.set_mesh(V, F);
        for (int k = 0; k < (int) K.size(); ++ k) {
            double kScale = K[k].mean();
            std::vector <Eigen::VectorXd> endpoints1, endpoints2;
            for (int i = 0; i < ZV[k].rows(); ++ i)
                if (zero_adj[k][i].size() != 2) {
                    /// begin with one endpoint
                    for (int t: zero_adj[k][i]) {
                        std::vector <int> points_inside;
                        points_inside.push_back(i);
                        int v = t;
                        while (zero_adj[k][v].size() == 2) {
                            for (int u: zero_adj[k][v])
                                if (u != points_inside.back()) {
                                    points_inside.push_back(v);
                                    v = u;
                                    break;
                                }
                        }
                        points_inside.push_back(v);
                        /// find a line: check threshold
                        double total_weight = 0.;
                        for (int j = 0; j + 1 < (int) points_inside.size(); ++ j) {
                            total_weight += 0.5 * (ZW[k][points_inside[j]] + ZW[k][points_inside[j + 1]]) * (ZV[k].row(points_inside[j]) - ZV[k].row(points_inside[j + 1])).norm();
                        }
                        if ((k == 0 && total_weight > threshold * vScale * kScale) || (k == 1 && total_weight < threshold * vScale * kScale)) {
                            for (int j = 0; j + 1 < (int) points_inside.size(); ++ j)
                                if (points_inside[j] < points_inside[j + 1]) {
                                    endpoints1.push_back(ZV[k].row(points_inside[j]));
                                    endpoints2.push_back(ZV[k].row(points_inside[j + 1]));
                                }
                        }
                    }
                }
            auto color = Eigen::RowVector3d(0.1, 0.1, 0.1);
            color(k) = 0.9;
            viewer.data.add_edges(GM::toEigenMatrix(endpoints1), GM::toEigenMatrix(endpoints2), color);
        }
    }

    return true;
}

int main(int argc, char *argv[]) {
    if (argc != 2) {
        cout << "Usage smoothfeatureline mesh.obj" << endl;
        exit(0);
    }

    // Read mesh
    igl::readOFF(argv[1],V,F);

    // Plot the mesh
    Viewer viewer;
    viewer.callback_key_down = callback_key_down;
    callback_key_down(viewer, '1', 0);

    viewer.callback_init = [&](Viewer &v) {
        v.ngui->addVariable("Threshold", threshold);
        v.screen->performLayout();
        return false;
    };

    // Initialize scale for displaying vectors
    vScale = 0.5 * igl::avg_edge_length(V, F);

    // Compute face barycenters
    igl::barycenter(V, F, MF);

    // Compute face normals
    igl::per_face_normals(V, F, FN);

    // Compute vertex to face adjacency
    igl::vertex_triangle_adjacency(V, F, VF, VFi);

    // Compute vertex normals
    igl::per_vertex_normals(V, F, VN);

    GM::compute_shape_operators(V, F, VN, FN, VS);
    GM::compute_eigens(VS, K, EV);
    is_regular.resize(K.size());
    extremalities.resize(K.size());
    my_basis = GM::get_local_basis(V, F);
    for (int k = 0; k < (int) K.size(); ++ k) {
        Eigen::VectorXd area_star;
        GM::compute_is_regular(F, EV[k], is_regular[k]);
        GM::compute_area_star(V, F, is_regular[k], area_star);
        auto gradient_complex = GM::get_gradient(V, F, 3.0 * K[k].array() / area_star.array(), my_basis);
        auto gradient = GM::to_vector_field(gradient_complex, my_basis);
        GM::compute_extremalities(V, F, gradient, EV[k], is_regular[k], extremalities[k]);
//        K[k] = 3.0 * K[k].array() / area_star.array();
    }

    viewer.core.point_size = 10;

    viewer.launch();
}