#File:		video_generator.py
#Date:		19.07.2023
#Author:	Armand Sieber
#Tag:		Create an animated view of the bubble and fluid-fluid interface dynamics


import numpy as np
import matplotlib.pyplot as plt

# --------- Function ---------
# ----------------------------
def find_nearest(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return idx

# ---- User-defined inputs ----
# -----------------------------

# Must be provided by user
file_name = 'free_surface_Rayleigh_gamma075.txt' # simulation output filename
gamma = 0.75 # stand-off distance
Nb = 80 # number of elements on bubble surface
Ns = 60 # number of elements on fluid-fluid interface


# ---- Loads  solution vectors ----
# ---------------------------------
coordinates = np.genfromtxt(file_name, skip_header = 10, skip_footer=1)

time = coordinates[:, 1] # time vector

r_b = coordinates[:, 2:Nb+3] # r-coordinates bubble surface
z_b = coordinates[:, Nb+3:2*Nb+4] #z-coordinates bubble surface

r_s = coordinates[:, 4*Nb+6:4*Nb+6+Ns+1] # r-coordinates fluid-fluid interface
z_s = coordinates[:, 4*Nb+6+Ns+1:4*Nb+6+2*Ns+2] # z-coordinates fluid-fluid interface

# ---- Display bubble dynamics ----
# ---------------------------------

# create time vector at constant time intervals
n_dt = len(coordinates[:,0])
indices =[]
time_max = time[-1]
iterations = int(n_dt/5)    # to increase or decrease the temporal resolution of
                            # the visualzation the value 5 can be changed
dt = time_max/iterations

for i in range(iterations):
    index = find_nearest(time, i*dt)
    indices.append(index)

# for loop to display bubble and fluid-fluid interface temporal evolution
for i in range(iterations):

    ti = int(indices[i]) # time step

    # bubble
    rb = r_b[ti, :]
    zb = z_b[ti, :]
    rb_flip = np.flipud(rb)
    zb_flip = np.flipud(zb)
    rb_tot = np.concatenate((rb, -rb_flip))
    zb_tot = np.concatenate((zb, zb_flip))

    # fluid-fluid interface
    rs = r_s[ti, :]
    zs = z_s[ti, :]
    rs_flip = np.flipud(rs)
    zs_flip = np.flipud(zs)
    rs_tot = np.concatenate((rs, -rs_flip))
    zs_tot = np.concatenate((zs, zs_flip))


    fig = plt.figure(1)

    # Region of interest
    ymin = - gamma - 1.2
    ymax = - gamma + 2.8
    xmin = -2.0
    xmax = 2.0

    plt.clf()
    plt.axis('equal')
    plt.ylim(ymin, ymax)
    plt.xlim(xmin, xmax)
    plt.grid()
    plt.plot(rb_tot, zb_tot, '-r', linewidth = 2)
    plt.plot(rs_tot, zs_tot, '-b', linewidth = 2)
    plt.xlabel(r'$r$',fontsize=12)
    plt.ylabel(r'$z$',fontsize=12)
    plt.pause(0.1)
