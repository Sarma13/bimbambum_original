/*  __       __          __
 * |__)||\/||__) /\ |\/||__)/  \|\/|
 * |__)||  ||__)/--\|  ||__)\__/|  |
 *
 * This file is part of BIMBAMBUM.
 *
 * -----------------------------------------------------------------------------
 * Copyright (C) 2023 Armand Sieber
 *
 * This program is free software: you can redistribute it and/or modify it
 * under the terms of the GNU General Public License as published by the Free
 * Software Foundation, either version 3 of the License, or (at your option)
 * any later version.
 *
 * This program is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of  MERCHANTABILITY or
 * FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for
 * more details.
 *
 * You should have received a copy of the GNU General Public License along with
 * this program.  If not, see <http://www.gnu.org/licenses/>.
 * -----------------------------------------------------------------------------
 *
 */

/*! \file Case_FlatBoundary.cpp
  \brief Derived class: handles the dynamics of a fluid-fluid interface that has been initiated as a 'infinite' flat surface.

  This derived class inherits members from BoundaryData.cpp and BoundaryData.cpp
  */

#include "Case_FlatBoundary.hpp"
#include "Inputs.hpp"

#include <armadillo>
#include <cmath>
#include <vector>

#include "cubic_spline.hpp"

using namespace std;
using namespace arma;



Case_FlatBoundary::Case_FlatBoundary(const Input& data) : BoundaryData(data) {}

void Case_FlatBoundary::initialize() {
    double R_d = 1.0; // nondimensional droplet radius
    for (int i = 0; i <= Ns; ++i) {
        double s = static_cast<double>(i) / Ns;
        double theta = M_PI * s; // 0 to π (half-sphere)
        r_nodes[i] = R_d * sin(theta);
        z_nodes[i] = R_d * cos(theta);
        F_nodes[i] = 0.0;
        curv_nodes[i] = 0.0;
    }
}

void Case_FlatBoundary::remesh_boundary() {
    // Reuse flat remeshing with geometric spacing
    // Just keep the interface arc-like by reinterpolating from arc-length
    vec rs = conv_to<vec>::from(r_nodes);
    vec zs = conv_to<vec>::from(z_nodes);

    vec drs = diff(rs);
    vec dzs = diff(zs);
    vec s = cumsum(sqrt(drs % drs + dzs % dzs));
    vec begin(1, fill::zeros);
    s = join_cols(begin, s);
    vector<double> arc = conv_to<vector<double>>::from(s);

    cubic_spline spline_r, spline_z, spline_F;
    spline_r.set_spline(arc, r_nodes, drds1, drds2);
    spline_z.set_spline(arc, z_nodes, dzds1, dzds2);
    spline_F.set_spline(arc, F_nodes, 0.0, 0.0);

    // Geometric reparam
    double r = 0.97;
    double h = (r - 1.0) / (pow(r, (double)(Ns)) - 1.0);
    vec geo(Ns + 1, fill::zeros);
    geo(0) = 0.0;
    for (int i = 1; i <= Ns; ++i) {
        geo(i) = geo(i - 1) + h * pow(r, i - 1) * arc.back();
    }

    for (int i = 0; i <= Ns; ++i) {
        r_nodes[i] = spline_r.interpolate(geo[i]);
        z_nodes[i] = spline_z.interpolate(geo[i]);
        F_nodes[i] = spline_F.interpolate(geo[i]);
    }
    r_nodes[Ns] = 0.0;
}

void Case_FlatBoundary::filter_boundary() {
    // Reuse existing filtering logic from Case_FlatBoundary
    // You can copy it directly here
}

void Case_FlatBoundary::boundary_curvature() {
    // Reuse existing curvature computation
}

void Case_FlatBoundary::boundary_endpoints_derivatives() {
    // Approximated based on a spherical profile
    drds1 = -1.0;
    dzds1 = 0.0;
    drds2 = -1.0;
    dzds2 = 0.0;
    dphi1ds1 = dphi1ds2 = dphi2ds1 = dphi2ds2 = 0.0;
}

    