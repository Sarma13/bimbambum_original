#include <gsl/gsl_errno.h>
#include <gsl/gsl_integration.h>
#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <armadillo>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <functional>
#include <iomanip>
#include <iostream>
#include <string>


#include "cubic_spline.hpp"
#include "integrands_BIM.hpp"

namespace py = pybind11;
using namespace arma;

py::array_t<double> compute(py::array_t<double> grid_points_r,
                            py::array_t<double> grid_points_z,
                            py::array_t<double> r_nodes,
                            py::array_t<double> z_nodes,
                            py::array_t<double> phi_nodes,
                            py::array_t<double> un_nodes, int Nb, int Ns,
                            int N_points) {

    const int N = Nb + 1 + Ns + 1;

    // -- Converts the function parameters to Armadillo vectors --
    // -- for compatibility with the current BIM implementation --
    // -----------------------------------------------------------

    // Request buffer descriptors
    auto buf_grid_points_r = grid_points_r.request();
    auto buf_grid_points_z = grid_points_z.request();
    auto buf_r_nodes = r_nodes.request();
    auto buf_z_nodes = z_nodes.request();
    auto buf_phi_nodes = phi_nodes.request();
    auto buf_un_nodes = un_nodes.request();

    // Create solution array of type py::array_t<double>
    auto phi_domain = py::array(py::buffer_info(
            nullptr,        // Pointer to data
            sizeof(double), // Size of one item
            py::format_descriptor<double>::value, // Buffer format
            buf_grid_points_r.ndim,               // Dimensions
            {buf_grid_points_r.shape[0]}, // Elements on each dimension
            {sizeof(double)}              // Strides for each dimension
    ));

    auto buf_phi_domain = phi_domain.request();

    double *ptr_grid_points_r = (double *) buf_grid_points_r.ptr;
    double *ptr_grid_points_z = (double *) buf_grid_points_z.ptr;
    double *ptr_r_nodes = (double *) buf_r_nodes.ptr;
    double *ptr_z_nodes = (double *) buf_z_nodes.ptr;
    double *ptr_phi_nodes = (double *) buf_phi_nodes.ptr;
    double *ptr_un_nodes = (double *) buf_un_nodes.ptr;
    double *ptr_phi_domain = (double *) buf_phi_domain.ptr;

    vec grid_points_r_vec(N_points, fill::zeros);
    vec grid_points_z_vec(N_points, fill::zeros);
    vec r0_both(N, fill::zeros);
    vec z0_both(N, fill::zeros);
    vec phi_nodes_both(N, fill::zeros);
    vec un_nodes_both(N, fill::zeros);

    // Assign values to vectors
    for (int i = 0; i < N_points; i++) {
        grid_points_r_vec(i) = ptr_grid_points_r[i];
        grid_points_z_vec(i) = ptr_grid_points_z[i];
    }

    // Assign values to vectors
    for (int i = 0; i < N; i++) {
        r0_both(i) = ptr_r_nodes[i];
        z0_both(i) = ptr_z_nodes[i];
        phi_nodes_both(i) = ptr_phi_nodes[i];
        un_nodes_both(i) = ptr_un_nodes[i];
    }

    // Create Armadillo sub-vectors
    vec r0b = r0_both.subvec(0, Nb);
    vec r0s = r0_both.subvec(Nb + 1, N - 1);
    vec z0b = z0_both.subvec(0, Nb);
    vec z0s = z0_both.subvec(Nb + 1, N - 1);

    // Initialize coefficient matrix for BIM solution
    mat G_nodes(N_points, N, fill::zeros);
    mat H_nodes(N_points, N, fill::zeros);

    // ----- Set up cubic splines -----
    // --------------------------------

    double r0 = 0.0;
    double z0 = 0.0;
    double dr = 0.0;
    double dz = 0.0;

    vec drb = diff(r0b);
    vec dzb = diff(z0b);
    vec drs = diff(r0s);
    vec dzs = diff(z0s);

    vec pointsb = sqrt(drb % drb + dzb % dzb);
    vec pointss = sqrt(drs % drs + dzs % dzs);

    vec distance_b = cumsum(pointsb);
    vec distance_s = cumsum(pointss);
    vec begin(1, fill::zeros);

    distance_b = join_cols(begin, distance_b);
    distance_s = join_cols(begin, distance_s);

    // Spline requires c++ vectors
    std::vector<double> dist_b = conv_to < std::vector < double >> ::from(distance_b);
    std::vector<double> dist_s = conv_to < std::vector < double >> ::from(distance_s);

    std::vector<double> bubble_r_nodes = conv_to < std::vector < double >> ::from(r0b);
    std::vector<double> bubble_z_nodes = conv_to < std::vector < double >> ::from(z0b);
    std::vector<double> boundary_r_nodes =
            conv_to < std::vector < double >> ::from(r0s);
    std::vector<double> boundary_z_nodes =
            conv_to < std::vector < double >> ::from(z0s);

    // Cubic splines on bubble
    cubic_spline spline_rb;
    spline_rb.set_spline(dist_b, bubble_r_nodes, 1.0, -1.0);

    cubic_spline spline_zb;
    spline_zb.set_spline(dist_b, bubble_z_nodes, 0.0, 0.0);

    // Cubic splines on fluid-fluid interface
    cubic_spline spline_rs;
    spline_rs.set_spline(dist_s, boundary_r_nodes, 1, -1);

    cubic_spline spline_zs;
    spline_zs.set_spline(dist_s, boundary_z_nodes, 0, 0);

    gsl_set_error_handler_off();

    // ----- Loop over all node points -----
    // -------------------------------------

    for (int i = 0; i < N_points; ++i) { // nodes in the fluid domain
        r0 = grid_points_r_vec(i);
        z0 = grid_points_z_vec(i);
        for (int j = 0; j < N; ++j) { // node over the bubble and fluid-fluid interface


            // -- Extract integration intervals and cubic splines coefficients for integration --
            // ----------------------------------------------------------------------------------

            double a, a2, b, b2;  // two integration interval with subscript 1 and 2
            double m_ar, m_br, m_cr, m_dr, m_az, m_bz, m_cz,
                    m_dz;  // spline coefficients interval 1
            double m_ar2, m_br2, m_cr2, m_dr2, m_az2, m_bz2, m_cz2,
                    m_dz2;  // spline coefficients interval 2

            // The integration interval lies on the fluid-fluid interface
            if (j > Nb) {
                a = 0.0;
                a2 = 0.0;
                b = 0.0;
                b2 = 0.0;

                if (j < (N - 1)) {
                    b = distance_s[1 + j - (Nb + 1)] - distance_s[j - (Nb + 1)];
                    m_ar = spline_rs.m_a[j - (Nb + 1)];
                    m_br = spline_rs.m_b[j - (Nb + 1)];
                    m_cr = spline_rs.m_c[j - (Nb + 1)];
                    m_dr = r0_both(j);
                    m_az = spline_zs.m_a[j - (Nb + 1)];
                    m_bz = spline_zs.m_b[j - (Nb + 1)];
                    m_cz = spline_zs.m_c[j - (Nb + 1)];
                    m_dz = z0_both(j);
                }
                if (j > Nb + 1) {
                    b2 = distance_s[j - (Nb + 1)] - distance_s[j - 1 - (Nb + 1)];
                    m_ar2 = spline_rs.m_a[j - 1 - (Nb + 1)];
                    m_br2 = spline_rs.m_b[j - 1 - (Nb + 1)];
                    m_cr2 = spline_rs.m_c[j - 1 - (Nb + 1)];
                    m_dr2 = r0_both(j - 1);
                    m_az2 = spline_zs.m_a[j - 1 - (Nb + 1)];
                    m_bz2 = spline_zs.m_b[j - 1 - (Nb + 1)];
                    m_cz2 = spline_zs.m_c[j - 1 - (Nb + 1)];
                    m_dz2 = z0_both(j - 1);
                }
            }

                // The integration interval lies on the fluid-fluid interface
            else {
                a = 0.0;
                a2 = 0.0;
                b = 0.0;
                b2 = 0.0;

                if (j < Nb) {
                    b = distance_b[1 + j] - distance_b[j];
                    m_ar = spline_rb.m_a[j];
                    m_br = spline_rb.m_b[j];
                    m_cr = spline_rb.m_c[j];
                    m_dr = r0_both(j);
                    m_az = spline_zb.m_a[j];
                    m_bz = spline_zb.m_b[j];
                    m_cz = spline_zb.m_c[j];
                    m_dz = z0_both(j);
                }
                if (j > 0) {
                    b2 = distance_b[j] - distance_b[j - 1];
                    m_ar2 = spline_rb.m_a[j - 1];
                    m_br2 = spline_rb.m_b[j - 1];
                    m_cr2 = spline_rb.m_c[j - 1];
                    m_dr2 = r0_both(j - 1);
                    m_az2 = spline_zb.m_a[j - 1];
                    m_bz2 = spline_zb.m_b[j - 1];
                    m_cz2 = spline_zb.m_c[j - 1];
                    m_dz2 = z0_both(j - 1);
                }
            }

            // -- Compute the boundary integrals --
            // -----------------------------------

            // Special treatment for integration intervals at the bubble and
            // fluid-fluid interface endpoints
            if ((j == 0) || (j == (Nb + 1))) {
                gsl_integration_workspace *w_G1 = gsl_integration_workspace_alloc(2000);
                gsl_integration_workspace *w_H1 = gsl_integration_workspace_alloc(2000);

                double result_G1, error_G1;
                double result_H1, error_H1;

                G_H_params params;

                params.coeff_r0 = m_ar;
                params.coeff_r1 = m_br;
                params.coeff_r2 = m_cr;
                params.coeff_r3 = m_dr;
                params.coeff_z0 = m_az;
                params.coeff_z1 = m_bz;
                params.coeff_z2 = m_cz;
                params.coeff_z3 = m_dz;
                params.r0 = r0;
                params.z0 = z0;
                params.a = a;
                params.b = b;
                params.iter = j;
                params.N_b = Nb;
                params.N_s = Ns;

                gsl_function F1_1;
                F1_1.function = &G1;
                F1_1.params = &params;
                gsl_function F2_1;
                F2_1.function = &H1;
                F2_1.params = &params;

                gsl_integration_qags(&F1_1, a, b, 1e-12, 1e-12, 2000, w_G1, &result_G1,
                                     &error_G1);
                gsl_integration_qags(&F2_1, a, b, 1e-12, 1e-12, 2000, w_H1, &result_H1,
                                     &error_H1);

                G_nodes(i, j) = result_G1;
                H_nodes(i, j) = result_H1;

                gsl_integration_workspace_free(w_G1);
                gsl_integration_workspace_free(w_H1);

            }

                // Special treatment for integration intervals at the bubble and
                // fluid-fluid interface endpoints
            else if ((j == Nb) || (j == (Nb + Ns + 1))) {
                gsl_integration_workspace *w_G2 = gsl_integration_workspace_alloc(2000);
                gsl_integration_workspace *w_H2 = gsl_integration_workspace_alloc(2000);

                double result_G2, error_G2;
                double result_H2, error_H2;

                G_H_params params;

                params.coeff_r0 = m_ar2;
                params.coeff_r1 = m_br2;
                params.coeff_r2 = m_cr2;
                params.coeff_r3 = m_dr2;
                params.coeff_z0 = m_az2;
                params.coeff_z1 = m_bz2;
                params.coeff_z2 = m_cz2;
                params.coeff_z3 = m_dz2;
                params.r0 = r0;
                params.z0 = z0;
                params.a = a2;
                params.b = b2;
                params.iter = j;
                params.N_b = Nb;
                params.N_s = Ns;

                gsl_function F1_2;
                F1_2.function = &G2;
                F1_2.params = &params;
                gsl_function F2_2;
                F2_2.function = &H2;
                F2_2.params = &params;

                gsl_integration_qags(&F1_2, a2, b2, 1e-12, 1e-12, 2000, w_G2,
                                     &result_G2, &error_G2);
                gsl_integration_qags(&F2_2, a2, b2, 1e-12, 1e-12, 2000, w_H2,
                                     &result_H2, &error_H2);

                G_nodes(i, j) = result_G2;
                H_nodes(i, j) = result_H2;

                gsl_integration_workspace_free(w_G2);
                gsl_integration_workspace_free(w_H2);

            }

                // For all other integration intervals on the bubble and fluid-fluid
                // interface
            else {
                gsl_integration_workspace *w_G1 = gsl_integration_workspace_alloc(2000);
                gsl_integration_workspace *w_H1 = gsl_integration_workspace_alloc(2000);

                double result_G1, error_G1;
                double result_H1, error_H1;

                G_H_params params;

                params.coeff_r0 = m_ar;
                params.coeff_r1 = m_br;
                params.coeff_r2 = m_cr;
                params.coeff_r3 = m_dr;
                params.coeff_z0 = m_az;
                params.coeff_z1 = m_bz;
                params.coeff_z2 = m_cz;
                params.coeff_z3 = m_dz;
                params.r0 = r0;
                params.z0 = z0;
                params.a = a;
                params.b = b;
                params.iter = j;
                params.N_b = Nb;
                params.N_s = Ns;

                gsl_function F1_1;
                F1_1.function = &G1;
                F1_1.params = &params;
                gsl_function F2_1;
                F2_1.function = &H1;
                F2_1.params = &params;

                gsl_integration_qags(&F1_1, a, b, 1e-12, 1e-12, 2000, w_G1, &result_G1,
                                     &error_G1);
                gsl_integration_qags(&F2_1, a, b, 1e-12, 1e-12, 2000, w_H1, &result_H1,
                                     &error_H1);

                double G1_tmp = result_G1;
                double H1_tmp = result_H1;

                gsl_integration_workspace_free(w_G1);
                gsl_integration_workspace_free(w_H1);

                gsl_integration_workspace *w_G2 = gsl_integration_workspace_alloc(2000);
                gsl_integration_workspace *w_H2 = gsl_integration_workspace_alloc(2000);

                double result_G2, error_G2;
                double result_H2, error_H2;

                G_H_params params2;

                params2.coeff_r0 = m_ar2;
                params2.coeff_r1 = m_br2;
                params2.coeff_r2 = m_cr2;
                params2.coeff_r3 = m_dr2;
                params2.coeff_z0 = m_az2;
                params2.coeff_z1 = m_bz2;
                params2.coeff_z2 = m_cz2;
                params2.coeff_z3 = m_dz2;
                params2.r0 = r0;
                params2.z0 = z0;
                params2.a = a2;
                params2.b = b2;
                params2.iter = j;
                params2.N_b = Nb;
                params2.N_s = Ns;

                gsl_function F1_2;
                F1_2.function = &G2;
                F1_2.params = &params2;
                gsl_function F2_2;
                F2_2.function = &H2;
                F2_2.params = &params2;

                gsl_integration_qags(&F1_2, a2, b2, 1e-12, 1e-12, 2000, w_G2,
                                     &result_G2, &error_G2);
                gsl_integration_qags(&F2_2, a2, b2, 1e-12, 1e-12, 2000, w_H2,
                                     &result_H2, &error_H2);

                double G2_tmp = result_G2;
                double H2_tmp = result_H2;

                gsl_integration_workspace_free(w_G2);
                gsl_integration_workspace_free(w_H2);

                G_nodes(i, j) = G1_tmp + G2_tmp;
                H_nodes(i, j) = H1_tmp + H2_tmp;
            }
        }
    }

    vec V1(N_points, fill::zeros);
    vec V2(N_points, fill::zeros);
    vec phi_int(N_points, fill::zeros);

    V1 = H_nodes * phi_nodes_both;
    V2 = G_nodes * un_nodes_both;

    double pi = 3.14159265358979323846264338328;
    phi_int = -(V2 + V1) / (4 * pi);

    // Assign solution values to pybind array
    for (int i = 0; i < N_points; i++) {
        ptr_phi_domain[i] = phi_int(i);
    }
    return phi_domain;
}

PYBIND11_MODULE(flow_potential, module_handle
) {
module_handle.

doc() = "Compute the velocity potentials in the flow domain.";

module_handle.def("compute", &compute);
}
