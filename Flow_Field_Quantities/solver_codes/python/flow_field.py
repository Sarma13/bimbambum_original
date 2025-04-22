""""
*  __       __          --
* |__)||\/||__) /\ |\/||__)/  \|\/|
* |__)||  ||__)/--\|  ||__)\__/|  |
*
* This file is part of BIMBAMBUM.
*
* flow_field.py
*
* Description:
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
from mpl_toolkits.axes_grid1.inset_locator import inset_axes

import inputs
import flow_potential
import create_domain
import classify_domain


class flow_field:

    # initialize quantities
    def __init__(self, t_i):

        # --- Read input data form input structure ---
        # --------------------------------------------

        input_data = inputs.simulation_inputs # create object

        self.file_name = input_data.file_name
        self.dumper_name = input_data.dumper_name
        self.Nb = input_data.Nb
        self.Ns = input_data.Ns
        self.zeta = input_data.zeta
        self.time_step = t_i
        self.h_grid = input_data.h_grid
        self.s_grid = input_data.s_grid
        self.r_min = 0.0 # set to 0 since it is defined by the axis of symmetry
        self.r_max = input_data.r_max
        self.z_min = input_data.z_min
        self.z_max = input_data.z_max

        # --- Extract relevant data from simulation input file ---
        # --------------------------------------------------------

        # skip_header --> 10 first lines contain information about the simulation properties
        # skip_footer --> avoids errors if last result line is incompletely written
        coordinates = np.genfromtxt(self.file_name, skip_header = 10, skip_footer=1) # all data

        # Bubble
        self.r_b = coordinates[:, 2:self.Nb+3] # r-coordinate
        self.z_b = coordinates[:, self.Nb+3:2*self.Nb+4] # z-coordinate
        self.phi_b = coordinates[:, 2*self.Nb+4:3*self.Nb+5] # potential
        self.un_b = coordinates[:, 3*self.Nb+5:4*self.Nb+6] # normal velocity

        # Fluid-fluid interface
        self.r_s = coordinates[:, 4*self.Nb+6:4*self.Nb+6+self.Ns+1] # r-coordinate
        self.z_s = coordinates[:, 4*self.Nb+6+self.Ns+1:4*self.Nb+6+2*self.Ns+2] # z-coordinate
        self.phi_s = coordinates[:, 4*self.Nb+6+2*self.Ns+2:4*self.Nb+6+3*self.Ns+3] # potential
        self.un_s = coordinates[:, 4*self.Nb+6+3*self.Ns+3:4*self.Nb+6+4*self.Ns+4] # normal velocity

        # Time information
        self.time_sim = coordinates[:, 1] # all simulation times
        self.t_k = [self.time_step, self.time_step-1, self.time_step+1]  # time step and surrounding time step for temporal derivative

        # --- Computational grid ---
        # --------------------------

        # Grid nodes coordinates
        self.x_grid = np.arange(self.r_min, self.r_max, self.h_grid) # all nodes position along x
        self.y_grid = np.arange(self.z_min, self.z_max, self.h_grid) # all nodes position along y
        self.x_grid_coord = [] # list of nodes position along x that are bounded within the fluid domain
        self.y_grid_coord = [] # list of nodes position along y that are bounded within the fluid domain
        self.x_ghost_coord = [] # additional grid points needed to evaluate spatial derivatives at the boundaries of the domain
        self.y_ghost_coord = [] # additional grid points needed to evaluate spatial derivatives at the boundaries of the domain
        self.ghost_coord = []

        # Domain creation and additional information
        self.Nx_domain = len(self.x_grid) # number of elements along x
        self.Ny_domain = len(self.y_grid) # number of elements along y
        self.in_domain = np.zeros(self.Nx_domain*self.Ny_domain) # 1d array indicating whether a grid point is bounded within the fluid domain
        self.domain_xy = np.zeros((self.Ny_domain, self.Nx_domain)) # 2d array containing information about valid grid points (boundaries, axis of symmetry, ...)

        # --- Solution vectors -----
        # --------------------------
        self.phi_fluid = []
        self.ux_fluid = []
        self.uy_fluid = []
        self.U_fluid = []
        self.x_fluid = []
        self.y_fluid = []
        self.P_fluid = []

    # creates structured grid to evaluate flow-field quantities
    def domain_creation(self):

        # Loops over all three time_steps
        for t_step in range(3):

            # Time step data
            ti = self.t_k[t_step]
            rb = self.r_b[ti, :]
            zb = self.z_b[ti, :]
            rs = self.r_s[ti, :]
            zs = self.z_s[ti, :]

            # Create reference computational grid at "t = time_step"
            if t_step == 0:
                self.in_domain = create_domain.setup_boundaries(rb, zb, rs, zs, self.x_grid, self.y_grid, self.s_grid)

            # Checks if computational grid is valid at "t = time_step - 1" and "t = time_step + 1"
            else:

                in_domain_check = create_domain.setup_boundaries(rb, zb, rs, zs, self.x_grid, self.y_grid, self.s_grid)

                for i in range(len(in_domain_check)):
                    # Check if all points valid in the domain at time_step are also valid at time_step + 1 and time_step - 1
                    if self.in_domain[i] == 1 and in_domain_check[i] == 0:
                        self.in_domain[i] = 0 # point removed from valid domain

        # Classify domain node points and find necessary ghost nodes
        ############################################################

        x_domain , y_domain = np.meshgrid(self.x_grid, self.y_grid)
        self.domain_xy, self.ghost_coord = classify_domain.classification(self.in_domain, x_domain, y_domain, self.Nx_domain, self.Ny_domain, self.s_grid)


        for i in range(len(self.y_grid)):
            for j in range(len(self.x_grid)):
                if self.domain_xy[i][j] > 0: # if grid point lies within accepted fluid domain
                    self.x_grid_coord.append(self.x_grid[j])
                    self.y_grid_coord.append(self.y_grid[i])


        x_add_coord = np.zeros((len(self.ghost_coord), 4)) # extract solely the x-coordinate from ghost_coord
        y_add_coord = np.zeros((len(self.ghost_coord), 4)) # extract solely the y-coordinate from ghost_coord
        for i in range(len(self.ghost_coord)):
            coord = self.ghost_coord[i]
            x_add_coord[i][:] = coord[2:6]
            y_add_coord[i][:] = coord[6:10]

        self.x_ghost_coord = np.ravel(x_add_coord)
        self.y_ghost_coord = np.ravel(y_add_coord)

        in_domain_ghost_check = create_domain.setup_boundaries(rb, zb, rs, zs, self.x_ghost_coord, self.y_ghost_coord, self.s_grid)




    # evaluate flow field velocities: ux = dphi/dx and uy = dphi/dy
    def compute_velocities(self):

        x_domain , y_domain = np.meshgrid(self.x_grid, self.y_grid)

        # Compute potential over all three time steps (needed for pressure afterwards)
        ##############################################################################
        for t_step in range(3):

            ti = self.t_k[t_step]
            rb = self.r_b[ti, :]
            zb = self.z_b[ti, :]
            phib = self.phi_b[ti, :]
            unb = self.un_b[ti, :]
            rs = self.r_s[ti, :]
            zs = self.z_s[ti, :]
            phis = self.phi_s[ti, :]
            uns = self.un_s[ti, :]

            r_nodes = np.concatenate((rb, rs))
            z_nodes = np.concatenate((zb, zs))
            phi_nodes = np.concatenate((phib, phis))
            un_nodes = np.concatenate((unb, uns))

            #Compute the potentials at the grid nodes within the fluid domain
            N_points = int(np.sum(self.in_domain))
            phi_interior = flow_potential.compute(self.x_grid_coord, self.y_grid_coord, r_nodes, z_nodes, phi_nodes, un_nodes, self.Nb, self.Ns, N_points)

            #Compute the potentials at the ghost nodes within the fluid domain
            if t_step == 0:
                N_ghost = len(self.x_ghost_coord)
                phi_ghost = flow_potential.compute(self.x_ghost_coord, self.y_ghost_coord, r_nodes, z_nodes, phi_nodes, un_nodes, self.Nb, self.Ns, N_ghost)
                phi_ghost = np.reshape(phi_ghost, (len(self.ghost_coord), 4))




            #Reshape 1D vectors so that they match the 2D computational domain
            k = 0
            phi_domain = np.zeros(len(self.in_domain))
            for i in range(len(self.in_domain)):
                if self.in_domain[i] == 1:
                    phi_domain[i] = phi_interior[k]
                    k = k + 1

            phi_domain = np.reshape(phi_domain, (self.Ny_domain, self.Nx_domain))

            #plt.imshow(phi_domain)
            #plt.show()

            # Compute velocities
            ####################

            if t_step == 0:
                # in the fluid doamin
                ux = np.empty((self.Ny_domain, self.Nx_domain))
                ux[:] = np.NaN
                uy = np.empty((self.Ny_domain, self.Nx_domain))
                uy[:] = np.NaN
                for i in range(self.Ny_domain):
                    for j in range(self.Nx_domain):
                        if self.domain_xy[i][j] == 2: #all interior nodes
                            ux[i][j] = -(phi_domain[i][j-1] - phi_domain[i][j+1])/(2.0*self.h_grid)
                            uy[i][j] = -(phi_domain[i-1][j] - phi_domain[i+1][j])/(2.0*self.h_grid)
                        elif self.domain_xy[i][j] == 3: # on the axis of symmetry
                            ux[i][j] = 0.0
                            uy[i][j] = -(phi_domain[i-1][j] - phi_domain[i+1][j])/(2.0*self.h_grid)

                for k in range(len(phi_ghost)):
                    i = self.ghost_coord[k][0]
                    j = self.ghost_coord[k][1]
                    if self.domain_xy[i][j] == 1:
                        phi_ghost1 = phi_ghost[k][0]
                        phi_ghost2 = phi_ghost[k][1]
                        phi_ghost3 = phi_ghost[k][2]
                        phi_ghost4 = phi_ghost[k][3]
                        ux[i][j] = -(phi_ghost1 - phi_ghost3)/(2.0*self.s_grid/4.0)
                        uy[i][j] = -(phi_ghost4 - phi_ghost2)/(2.0*self.s_grid/4.0)
                    elif self.domain_xy[i][j] == 4:
                        phi_ghost2 = phi_ghost[k][1]
                        phi_ghost4 = phi_ghost[k][3]
                        ux[i][j] = 0.0
                        uy[i][j] = -(phi_ghost4 - phi_ghost2)/(2.0*self.s_grid/4.0)



                # Store intermediate results
                ############################
                self.phi_fluid.append(np.ravel(phi_domain))
                self.ux_fluid.append(np.ravel(ux))
                self.uy_fluid.append(np.ravel(uy))
                self.U_fluid.append(np.ravel(np.sqrt(ux*ux + uy*uy)))
                self.x_fluid.append(np.ravel(x_domain))
                self.y_fluid.append(np.ravel(y_domain))

            # stored for pressure computation
            else:
                self.phi_fluid.append(np.ravel(phi_domain))

    # evaluate flow field pressure
    def compute_pressure(self):

        # Compute pressure field
        ########################

        #Potential time derivative (second order polynomial fit)
        phi_fluid_array = np.array(self.phi_fluid)
        U_fluid_array = np.array(self.U_fluid)
        y_fluid_array = np.array(self.y_fluid)

        phi_fluid1 = phi_fluid_array[0,:] # t = i
        phi_fluid2 = phi_fluid_array[1,:] # t = i - 1
        phi_fluid3 = phi_fluid_array[2,:] # t = i + 1
        dt1 = self.time_sim[self.time_step] - self.time_sim[self.time_step-1]
        dt2 = self.time_sim[self.time_step+1] - self.time_sim[self.time_step]
        dphidt_fluid = phi_fluid3*dt1*dt1 - phi_fluid2*dt2*dt2 + phi_fluid1*(dt2*dt2 - dt1*dt1)
        dphidt_fluid = dphidt_fluid/(dt2*dt1*(dt1 + dt2))

        self.P_fluid = 1 - 0.5*U_fluid_array*U_fluid_array  - dphidt_fluid - self.zeta*y_fluid_array


    # write solution to user defined dumper files
    def write_solution(self):

        z_fluid = np.zeros_like(self.x_fluid) # included to allow post-processing in Paraview
        uz_fluid = np.zeros_like(self.ux_fluid) # included to allow post-processing in Paraview
        solution = np.vstack((self.x_fluid, self.y_fluid, z_fluid, self.ux_fluid, self.uy_fluid, uz_fluid, self.U_fluid, self.P_fluid))
        solution_name = self.dumper_name + "{0:0=4d}".format(self.time_step) + ".txt"
        np.savetxt(solution_name, solution.T, delimiter=',', newline = '\n', header="x, y, z, velocity_x, velocity_y, velocity_z, velocity_norm, pressure", comments='')

    # some basic visualizations of the flow field quantities
    def post_processing(self):

        x_domain, y_domain = np.meshgrid(self.x_grid, self.y_grid)
        ti = self.time_step


        # Velocity contour plot
        velocity_field = np.reshape(self.U_fluid, (self.Ny_domain, self.Nx_domain))

        fig1 = plt.figure(1)
        ax = fig1.add_subplot()
        plt.axis('equal')
        cont = ax.contourf(x_domain, y_domain, velocity_field, cmap ="turbo", levels = 100)
        ax.contourf(-x_domain, y_domain, velocity_field, cmap ="turbo", levels = 100)
        ax.plot(self.r_b[ti, :], self.z_b[ti, :], '-', color='black', linewidth = 1.5)
        ax.plot(-self.r_b[ti, :], self.z_b[ti, :], '-', color='black', linewidth = 1.5)
        ax.plot(self.r_s[ti, :], self.z_s[ti, :], '-', color='black', linewidth = 1.5)
        ax.plot(-self.r_s[ti, :], self.z_s[ti, :], '-', color='black', linewidth = 1.5)
        cbar = plt.colorbar(cont)
        cbar.set_label('U [-]')
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.set_xlim((-self.r_max, self.r_max))
        ax.set_xlabel('r [-]')
        ax.set_ylabel('z [-]')


        # Pressure contour plot
        pressure_field = np.reshape(self.P_fluid, (self.Ny_domain, self.Nx_domain))

        fig2 = plt.figure(2)
        ax = fig2.add_subplot()
        plt.axis('equal')
        cont = ax.contourf(x_domain, y_domain, pressure_field, cmap ="turbo", levels = 100)
        ax.contourf(-x_domain, y_domain, pressure_field, cmap ="turbo", levels = 100)
        ax.plot(self.r_b[ti, :], self.z_b[ti, :], '-', color='black', linewidth = 1.5)
        ax.plot(-self.r_b[ti, :], self.z_b[ti, :], '-', color='black', linewidth = 1.5)
        ax.plot(self.r_s[ti, :], self.z_s[ti, :], '-', color='black', linewidth = 1.5)
        ax.plot(-self.r_s[ti, :], self.z_s[ti, :], '-', color='black', linewidth = 1.5)
        cbar = plt.colorbar(cont)
        cbar.set_label('P [-]')
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.set_xlim((-self.r_max, self.r_max))
        ax.set_xlabel('r [-]')
        ax.set_ylabel('z [-]')

        # Velocity 3D scatter plot
        fig3 = plt.figure(3)
        ax = fig3.add_subplot(projection='3d')
        ax.set_title('U [-]')
        ax.scatter(self.x_fluid, self.y_fluid, self.U_fluid)
        ax.set_xlabel('r [-]')
        ax.set_ylabel('z [-]')
        ax.set_zlabel('U [-]')


        # Pressure 3D scatter plot
        fig4 = plt.figure(4)
        ax = fig4.add_subplot(projection='3d')
        ax.set_title('P [-]')
        ax.scatter(self.x_fluid, self.y_fluid, self.P_fluid)
        ax.set_xlabel('r [-]')
        ax.set_ylabel('z [-]')
        ax.set_zlabel('P [-]')


        plt.show()
