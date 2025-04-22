""""
*  __       __          --
* |__)||\/||__) /\ |\/||__)/  \|\/|
* |__)||  ||__)/--\|  ||__)\__/|  |
*
* This file is part of BIMBAMBUM.
*
* classify_domain.py
*
* Description:
* Inputs:   -   in_domain: vetcor indicating whether the fluid grid points are bounded
*               within the fluid domain or not.
*           -   x_domain, y_domain: 2d arrays containing the positions of the grid points.
*           -   Nx_domain, Ny_domain: number of grid points in each direction.
*           -   shrink_value: fluid domain offseting constant.
* Outputs:  -   domain_xy: 2d array indicating the nature of each grid point (see below)
*           -   add_coord: contains the indices and coordinates of the additional nodes
*               needed for the computation of the spatial derivatives
*
*
* The underlying meshing grid is rectangular and has evenly spaced nodes. However,
* the computation of velocity potential is limited to nodes located within the
* fluid domain bounded by the bubble and the fluid-fluid interface.
* This script assigns an integer value to each grid point based on the following
* classification:
*   -   Nodes outside of fluid domain bounded by the bubble and the fluid-fluid
*       interface (i.e. node discarded for the computation): 0
*   -   Interior nodes adjacent to the domain boundaries: 1
*   -   All other interior nodes: 2
*   -   Nodes on the axis of symmetry (except the ones adjacent to the bubble
*       or the fluid-fluid interface): 3
*   -   All other nodes on the axis of symmetry: 4
*
* To evaluate the spatial derivatives at nodes adjacent to the boundaries (i.e. nodes
* classified 1 or 4) additional "ghost" points a created at a distance shrink_value/4
* from the node under investigation.
* N.B. this will result in boundary integrals being evaluated at a distance shrink_value/4
* from the bubble surface or fluid-fluid interface, instead of shrink_value.
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

def classification(in_domain, x_domain, y_domain, Nx_domain, Ny_domain, shrink_value):


    domain_xy = in_domain.copy()
    domain_xy = domain_xy.reshape((Ny_domain, Nx_domain)) # node points classification
    add_coord = [] # additional nodes at domain boundaries

    #colors = ['red', 'green', 'orange', 'blue', 'yellow', 'purple']
    #cmap = mpl.colors.ListedColormap(colors)
    #plt.imshow(domain_xy, cmap = cmap)
    #plt.show()
    for i in range(Ny_domain):
        for j in range(Nx_domain):

            # Check if the grid point is bounded in the fluid domain
            if domain_xy[i][j] != 0:

                # Nodes away from the axis of symmetry
                if j > 0:
                    if i > 0 and i != (Ny_domain - 1) and j != (Nx_domain - 1):

                        # Checks if the adjacent nodes (left, right, up, down) are also in the fluid domain
                        if domain_xy[i-1][j] != 0 and domain_xy[i+1][j] != 0 and domain_xy[i][j-1] != 0 and domain_xy[i][j+1] != 0:
                            domain_xy[i][j] = 2

                        # If adjacent node missing --> create 4 new adjacent to compute the spatial derivatives
                        # New nodes created at a distance shrink_value/4.0 from node (i,j)
                        else:
                            domain_xy[i][j] = 1
                            x_ghost1 = x_domain[i][j] - shrink_value/4.0
                            x_ghost2 = x_domain[i][j]
                            x_ghost3 = x_domain[i][j] + shrink_value/4.0
                            x_ghost4 = x_domain[i][j]
                            y_ghost1 = y_domain[i][j]
                            y_ghost2 = y_domain[i][j] + shrink_value/4.0
                            y_ghost3 = y_domain[i][j]
                            y_ghost4 = y_domain[i][j] - shrink_value/4.0
                            points = [i, j, x_ghost1, x_ghost2, x_ghost3, x_ghost4, y_ghost1, y_ghost2, y_ghost3, y_ghost4]
                            add_coord.append(points)


                # Nodes on the axis of symmetry
                else:

                    # Checks if adjacent nodes (up, down) are in the fluid domain
                    if i > 0 and i != (Ny_domain - 1):
                        if domain_xy[i-1][j] != 0 and domain_xy[i+1][j] != 0:
                            domain_xy[i][j] = 3

                        # If adjacent node missing --> create 2 new adjacent to compute the spatial derivatives
                        # New nodes created at a distance shrink_value/4.0 from node (i,j)
                        else:
                            domain_xy[i][j] = 4
                            x_ghost1 = x_domain[i][j]
                            x_ghost2 = x_domain[i][j]
                            x_ghost3 = x_domain[i][j]
                            x_ghost4 = x_domain[i][j]
                            y_ghost1 = y_domain[i][j]
                            y_ghost2 = y_domain[i][j] + shrink_value/4.0
                            y_ghost3 = y_domain[i][j]
                            y_ghost4 = y_domain[i][j] - shrink_value/4.0
                            points = [i, j, x_ghost1, x_ghost2, x_ghost3, x_ghost4, y_ghost1, y_ghost2, y_ghost3, y_ghost4]
                            add_coord.append(points)

                    # If adjacent node missing --> create 2 new adjacent to compute the spatial derivatives
                    # New nodes created at a distance shrink_value/4.0 from node (i,j)
                    else:
                        domain_xy[i][j] = 4
                        x_ghost1 = x_domain[i][j]
                        x_ghost2 = x_domain[i][j]
                        x_ghost3 = x_domain[i][j]
                        x_ghost4 = x_domain[i][j]
                        y_ghost1 = y_domain[i][j]
                        y_ghost2 = y_domain[i][j] + shrink_value/4.0
                        y_ghost3 = y_domain[i][j]
                        y_ghost4 = y_domain[i][j] - shrink_value/4.0
                        points = [i, j, x_ghost1, x_ghost2, x_ghost3, x_ghost4, y_ghost1, y_ghost2, y_ghost3, y_ghost4]
                        add_coord.append(points)


    #colors = ['red', 'green', 'orange', 'blue', 'yellow', 'purple']
    #cmap = mpl.colors.ListedColormap(colors)
    #plt.imshow(domain_xy, cmap = "turbo")
    #plt.show()
    return domain_xy, add_coord
