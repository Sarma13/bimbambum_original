""""
*  __       __          --
* |__)||\/||__) /\ |\/||__)/  \|\/|
* |__)||  ||__)/--\|  ||__)\__/|  |
*
* This file is part of BIMBAMBUM.
*
* inputs.py
*
* Description: class that contains all the user-defined inputs. Works as an
* argument parser.
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

class simulation_inputs:

    file_name = "rigid_boundary_Rayleigh_gamma095.txt" # Input file name (results of first processing phase)

    dumper_name = "rigid_boundary_gamma095_" # Output file name (no need to provide an extension --> the file will be named as: "dumper_name" + "time_step" + ".txt

    Nb = 80 # Number of elements on the bubble surface (must be that of the simulation results saved in "file_name")

    Ns = 60 # Number of elements on the fluid-fluid interface (must be that of the simulation results saved in "file_name")

    zeta = 0.0 # buoyancy parameter (must be that of the simulation results saved in "file_name")

    time_step = [850] # Time steps of interest (can be a list) for the computation of the flow field quantities
    # N.B. cannot be first nor last or second to last time step of computation owing to the temporal derivatives scheme

    h_grid = 0.02 # Distance between adjacent grid points

    s_grid = 0.015 # Shrinking value for the inwards ofsseting of the computational domain

    r_max = 2.1 # Maximum r-coordinate of the computational domain

    z_min = -2.1 # Minimium z-coordinate of the computational domain (fluid-fluid interface initially located at z = 0 and bubble initially located at z < 0)

    z_max = 2.0 # Maximum z-coordinate of the computational domain
