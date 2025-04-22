#File:		video_generator.py
#Date:		19.07.2023
#Author:	Armand Sieber
#Tag:		Short script to post-process key elements of the bubble dynamics

"""
This script is an example on how to read and use the the simulation results file
of the first processing stage of BIMBAMBUM v1.0.

The simulation results file is a '<file_name>.txt' file in which the 10 first lines contain key
information about the simulation parameters (spatial and temporal discretization, bubble properties, ...).
The remainder of the lines store quantities computed at the grid points on the bubble surface and on the
fluid-fluid interface. These may be accessed as follow:

    data = np.genfromtxt(file_name, skip_header = 10, skip_footer=1) # reads simulation file skipping the 10 first lines

    time_step = data[:, 0] # time steps id
    time = data[:, 1] # simulation time

    r_b = data[:, 2 : Nb + 3] # r-coordinate of the nodes on the bubble surface
    z_b = data[:, Nb + 3 : 2 * Nb + 4] # z-coordinate of the nodes on the bubble surface
    phi_b = data[:, 2 * Nb + 4 : 3 * Nb + 5] # value of the potentials on the bubble surface
    un_b = data[:, 3 * Nb + 5 : 4 * Nb + 6] # value of the normal velocities on the bubble surface

    r_s = data[:, 4 * Nb + 6 : 4 * Nb + 6 + Ns + 1] # r-coordinate of the nodes on the fluid-fluid interface
    z_s = data[:, 4 * Nb + 6 + Ns + 1 : 4 * Nb + 6 + 2 * Ns + 2] # z-coordinate of the nodes on the fluid-fluid interface
    phi_s = data[:, 4 * Nb + 6 + 2 * Ns + 2 : 4 * Nb + 6 + 3 * Ns + 3] # value of the potentials on the fluid-fluid interface (fluid 1)
    un_s = data[:, 4 * Nb + 6 + 3 * Ns + 3 : 4 * Nb + 6 + 4 * Ns + 4] # value of the normal velocities on the fluid-fluid interface (fluid 1)
    F_s = data[:, 4 * Nb + 6 + 4 * Ns + 4 : 4 * Nb + 6 + 5 * Ns + 5] # value of F on the fluid-fluid interface

    volume = data[:, -1] # bubble volume

where Nb is the number of elements used for the bubble discretization and number of elements for the fluid-fluid interface discretization.

The script below display the bubble and fluid-fluid interface shape at a time step of interest, as well as the temporal
evolution of the bubble volume and equivalent radius.
"""

import numpy as np
import matplotlib.pyplot as plt


# ---- User-defined inputs ----
# -----------------------------

# Must be provided by user
file_name = 'free_surface_Rayleigh_gamma075.txt' # simulation output filename
gamma = 0.75 # stand-off distance
Nb = 80 # number of elements on bubble surface
Ns = 60 # number of elements on fluid-fluid interface
time_step = -1 # time step of interest


# ---- Loads  solution vectors ----
# ---------------------------------
coordinates = np.genfromtxt(file_name, skip_header = 10, skip_footer=1)

time = coordinates[:, 1] # time vector

r_b = coordinates[:, 2:Nb+3] # r-coordinates bubble surface
z_b = coordinates[:, Nb+3:2*Nb+4] #z-coordinates bubble surface

r_s = coordinates[:, 4*Nb+6:4*Nb+6+Ns+1] # r-coordinates fluid-fluid interface
z_s = coordinates[:, 4*Nb+6+Ns+1:4*Nb+6+2*Ns+2] # z-coordinates fluid-fluid interface

volume = coordinates[:, -1] # bubble volume

# ---- Display solutions ----
# ---------------------------

# bubble volume and equivalent radius temporal evolution
# ------------------------------------------------------

pi = 3.1415926535897932
radius = (3.0 / 4.0 * volume / pi)**(1.0 / 3.0)

fig = plt.figure(1)
ax1 = fig.add_subplot(111)
plt.grid()
ax1.plot(time, volume, '-k', linewidth = 2, label = 'Volume')
ax1.set_ylabel('volume [-]',fontsize=12, color='black')
ax1.tick_params(axis='y', labelcolor='black')
ax1.set_xlabel('time [-]', fontsize=12)

ax2 = ax1.twinx()
ax2.plot(time, radius, '-r', linewidth = 2, label = 'Equ. radius')
ax2.set_ylabel('Equ. radius [-]',fontsize=12, color = 'red')
ax2.tick_params(axis='y', labelcolor='red')

# bubble and fluid-fluid interface shape
# --------------------------------------

# bubble
rb = r_b[time_step, :]
zb = z_b[time_step, :]
rb_flip = np.flipud(rb)
zb_flip = np.flipud(zb)
rb_tot = np.concatenate((rb, -rb_flip))
zb_tot = np.concatenate((zb, zb_flip))

# fluid-fluid interface
rs = r_s[-1, :]
zs = z_s[-1, :]
rs_flip = np.flipud(rs)
zs_flip = np.flipud(zs)
rs_tot = np.concatenate((rs, -rs_flip))
zs_tot = np.concatenate((zs, zs_flip))

# Region of interest
ymin = - gamma - 1.2
ymax = - gamma + 2.8
xmin = -2.0
xmax = 2.0

fig = plt.figure(2)
plt.axis('equal')
plt.ylim(ymin, ymax)
plt.xlim(xmin, xmax)
plt.grid()
plt.plot(rb_tot, zb_tot, '-r', linewidth = 2)
plt.plot(rs_tot, zs_tot, '-b', linewidth = 2)
plt.xlabel(r'$r$',fontsize=12)
plt.ylabel(r'$z$',fontsize=12)


plt.show()
