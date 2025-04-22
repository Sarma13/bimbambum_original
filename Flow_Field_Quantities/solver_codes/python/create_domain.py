""""
*  __       __          --
* |__)||\/||__) /\ |\/||__)/  \|\/|
* |__)||  ||__)/--\|  ||__)\__/|  |
*
* This file is part of BIMBAMBUM.
*
* create_domain.py
*
* Description:
* Inputs:   bubble coordinates (rb, zb), fluid-fluid interface coordinates
*           (rs, zs), coordinates of the rectangular fluid domain
*           (x_grid, y_grid) and an offseting constant shrink_value.
* Outputs:  in_domain, a vector indicating whether the point (x_grid, y_grid) is
*           bounded within the fluid domain.
*
* Based on (rb, zb), (rs, zs), and shrink_value, this script utilizes a winding
* number algorithm to solve a point-in-polygon problem. It determines whether the
* grid points (x_grid, y_grid) fall within the fluid domain bounded by the bubble and
* fluid-fluid interface boundaries that have been offseted by the constant
* quantity shrink_value.
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
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl

# ------------------
# --- Functions ---
# -----------------

# -- Find segments intersection to establish the vertices of the bounding polygon --
# ----------------------------------------------------------------------------------

# segment1 is defined by the points (x1, y1) and (x2, y2)
# segment2 is defined by the points (x3, y3) and (x4, y4)
def segments_intersection(segment1, segment2):

    x1 = segment1[0][0]
    x2 = segment1[1][0]
    x3 = segment2[0][0]
    x4 = segment2[1][0]
    y1 = segment1[0][1]
    y2 = segment1[1][1]
    y3 = segment2[0][1]
    y4 = segment2[1][1]

    dx1 = x1 - x2
    dy1 = y1 - y2
    dx2 = x3 - x4
    dy2 = y3 - y4

    det = dx1 * dy2 - dx2 * dy1

    if det == 0: # checks if the segment are parallel
       x_int = (x2 + x3) / 2.0
       y_int = (y2 + y3) / 2.0
    else:
       x_int = ((x1 * y2 - y1 * x2) * (x3 - x4) - (x1 - x2) * (x3 * y4 - y3 * x4)) / det
       y_int = ((x1 * y2 - y1 * x2) * (y3 - y4) - (y1 - y2) * (x3 * y4 - y3 * x4)) / det

    return x_int, y_int

# -- Checks if grid point is in the domain bounded by the polygon --
# ------------------------------------------------------------------

# Determine the signed shortest distance between a point and a line passing through vertex1 and vertex2
# Returns distance > 0 if the point lies on left of the line, distance < 0 if it's on the right
# and distance = 0 if the point lies on the line
def signed_distance(vertex1, vertex2, point):

    x_point = point[0]
    y_point = point[1]
    x1 = vertex1[0]
    y1 = vertex1[1]
    x2 = vertex2[0]
    y2 = vertex2[1]

    distance = ((x2 - x1) * (y_point - y1) - (y2 - y1) * (x_point - x1)) / ((x2 - x1)**2 + (y2 - y1)**2)**0.5

    return distance

# Winding number algorithm for point-in-polygon problem:
# Returns winding = 0 if the point is outside the polygon
def winding_number(point, polygon):

    y_point = point[1]
    polygon_tmp = polygon.copy()
    polygon_tmp.append(polygon[0]) #repeat first vertex of the polygon to close the contour

    n_vertices = len(polygon_tmp) - 1
    winding = 0   # winding number

    # loops over all segments of the polygon with vertices polygon[i] and polygon[i+1]
    for i in range(n_vertices):
        if polygon_tmp[i][1] <= y_point:
            if polygon_tmp[i+1][1]  > y_point:
                if signed_distance(polygon_tmp[i], polygon_tmp[i+1], point) > 0:
                    winding = winding + 1
        else:
            if polygon_tmp[i+1][1] <= y_point:
                if signed_distance(polygon_tmp[i], polygon_tmp[i+1], point) < 0:
                    winding = winding - 1
    return winding


# --- Create closed contour over ROI ---
# --------------------------------------

def setup_boundaries(rb, zb, rs, zs, x_grid, y_grid, shrink_value):

    # Computational domain limits (rectangular domain)
    minY = np.min(y_grid)

    # Additional points to close contour
    r_s_out = np.array([0, rs[0]])
    z_s_out = np.array([minY, minY])
    r_box = np.concatenate([rs, rb, r_s_out])
    z_box = np.concatenate([zs, zb, z_s_out])
    r_box = np.flipud(r_box)
    z_box = np.flipud(z_box)

    # build bounding polygon based on position of boundaries
    polygon = []
    for i in range(len(r_box)):
        vec = (r_box[i], z_box[i]) # (r, z) <--> (x, y) coordinates of polygon vertices
        polygon.append(vec)


    # extract polygon segments
    segments = []
    for i in range(len(polygon)):
        segments.append([polygon[i-1], polygon[i]])

    # Offset the vertices of the segments inwards
    new_segments = []
    for i in range(len(segments)):

        # If the segment lies on the axis of symmetry there is no need to perform the offsetting
        if (segments[i][1][0] == 0.0 and segments[i][0][0] == 0.0) or (segments[i][1][1] == minY and segments[i][0][1] == minY):

            new_segments.append([(segments[i][0][0], segments[i][0][1]), (segments[i][1][0], segments[i][1][1])])

        # Otherwise the vertices are offsetted inwards
        else:
            dx = segments[i][1][0] - segments[i][0][0]
            dy = segments[i][1][1] - segments[i][0][1]

            # this is to take into account the slopes of the segments
            offset_dx = dy * shrink_value / (dx*dx + dy*dy)**0.5
            offset_dy = dx * shrink_value / (dx*dx + dy*dy)**0.5

            new_segments.append([(segments[i][0][0] + offset_dx, segments[i][0][1] - offset_dy),
                              (segments[i][1][0] + offset_dx, segments[i][1][1] - offset_dy)])


    # Find vertices of new bounding polygon by checking intersection of new segments
    new_polygon = []
    for i in range(len(new_segments)):
        new_polygon.append((segments_intersection(new_segments[i-1], new_segments[i])))


    # Check whether grid points lay within the bounding polygon
    ###########################################################
    # in_domain = 1 if they are valid points (i.e. in polygon)
    # in_domain = 0 if they are outside of the polygon

    in_domain = []
    for i in range(len(y_grid)):
        for j in range(len(x_grid)):
            point = (x_grid[j], y_grid[i])
            point_in_polygon = winding_number(point, new_polygon)
            if point_in_polygon != 0:
                in_domain.append(1)
            else:
                in_domain.append(0)

    # Converts everything in numpy array
    in_domain = np.array(in_domain)

    return in_domain
