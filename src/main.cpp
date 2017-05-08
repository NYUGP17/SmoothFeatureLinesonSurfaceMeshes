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
#include "extremalities.h"
#include "parametrization.h"
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
// Scale for displaying vectors
double vScale = 0;

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
        for (int k = 0; k < (int) K.size(); ++ k);

        Eigen::MatrixXd colors;
        igl::jet(smoothed[0], true, colors);
        viewer.data.set_colors(colors);
    }

    if (key == '5') {
        viewer.data.clear();
        viewer.data.set_mesh(V, F);
    }

    if (key == '6') {
        viewer.data.clear();
        viewer.data.set_mesh(V, F);
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
    auto my_basis = GM::get_local_basis(V, F);
    for (int k = 0; k < (int) K.size(); ++ k) {
        GM::make_consistent(F, EV[k], is_regular[k]);
        Eigen::VectorXd area_star;
        GM::compute_area_star(V, F, area_star);
        auto gradient_complex = GM::get_gradient(V, F, 3.0 * K[k].array() / area_star.array(), my_basis);
        auto gradient = GM::to_vector_field(gradient_complex, my_basis);
        GM::compute_extremalities(V, F, gradient, EV[k], extremalities[k]);
    }
    delete[] my_basis;

    viewer.core.point_size = 10;

    viewer.launch();
}